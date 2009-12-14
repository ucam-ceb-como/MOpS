/*
 * File:   stagflow.cpp
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the stagnation flow and twin flame model
 * Licence:
 *  This file is part of "Camflow".
 *
 *  Camflow is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * Contact:
 *  Dr Markus Kraft
 *  Dept of Chemical Engineering
 *  University of Cambridge
 *  New Museum Site
 *  Pembroke Street
 *  Cambridge
 *  CB2 3RA
 *  UK
 *
 *  Email   :   mk306@cam.ac.uk
 *  Website :   http://como.cheng.cam.ac.uk
 *
 * Created on 30 April 2009, 15:03
 */
#include "cam_reporter.h"
#include "cam_residual.h"
#include "cvode_wrapper.h"
#include "stagflow.h"
#include "array.h"
#include "cam_geometry.h"
#include "cam_setup.h"
#include "cam_math.h"
#include <cmath>
#include <algorithm>
#include "cam_params.h"
#include "cam_boundary.h"
#include <vector>

using namespace Camflow;
/*
 *solve the stagnation/twinflame model
 */
void StagFlow::solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg,CamProfile &cp,
             CamConfiguration &config, CamSoot &cs,  Mechanism &mech ){
    /*
     *function to set up the solver.
     *Initialisation of the solution vector is
     *done here
     */
    
    camMech = &mech;
    Thermo::Mixture mix(mech.Species());
    camMixture = &mix;
    spv = camMixture->Species();
    profile = &cp;
    camConfig = &config;

    
    admin = &ca;
    reacGeom = &cg;
    reacGeom->discretize();
    /*
     * 2 additional cells are padded to consider the
     * inlet and the exhaust
     */

    reacGeom->addZeroWidthCells();
    nCells = reacGeom->getnCells();
    configID = camConfig->getConfiguration();

    opPre = ca.getPressure();
    strainRate = ca.getStrainRate();
    profile->setGeometryObj(cg);
    reporter = new CamReporter();
    nSpc = camMech->SpeciesCount();
    /*
     *array offsets for ODE
     */
    ptrT = nSpc;
    //ptrG = ptrT+1;
    
    /*
     *number of finite volume cells and
     *loop initializations
     */
    
    iMesh_s = 1;
    iMesh_e = nCells-1;

    cellBegin = 0;
    cellEnd = nCells;
    

    /*
     *number of equations and the number of variables to
     *be slved by the ODE. (Species, Energy equations)
     */
    nVar = nSpc+1;
    nEqn = nVar*nCells;

    /*
     *init the solution vector
     */
    initSolutionVector(cc);

    reporter->header("StagFlow");    
    header();
    if(cc.getSolutionMode() == cc.COUPLED)
        csolve(cc);
    else{
        ssolve(cc);
        csolve(cc);
    }

}
/*
 *initialize the solution vector
 */
void StagFlow::initSolutionVector(CamControl& cc){

    /*
     *initialize the geometry
     */
    dz = reacGeom->getGeometry();    


    CamBoundary left, right;
    admin->getLeftBoundary(left);
    storeInlet(left,fuel);
    if(configID == camConfig->COUNTERFLOW){
        admin->getRightBoundary(right);
        storeInlet(right,oxid);
        oxid.Vel *= -1;
        oxid.FlowRate *= -1;
        right.setVelocity(-right.getVelocity());
        right.setFlowRate(-right.getFlowRate());
    }


    /*
     *initialize the ODE vector
     */
    solvect.resize(nEqn,0);
    std::vector<doublereal> vSpec, vT;
    initSpecies(left,right,cc,vSpec);        
    initTemperature(left,cc,vT);

    if(admin->getEnergyModel() == admin->ADIABATIC){
        if(configID == camConfig->STAGFLOW)
            vT[iMesh_e] = admin->getWallTemp();
    }

    initMomentum();
    initMassFlow();
    
    mergeEnergyVector(&vT[0]);
    mergeSpeciesVector(&vSpec[0]);    
    
   
}
/*
 *initialize the mass flow
 */
void StagFlow::initMassFlow(){
    //store the flow rates on the cell faces
    m_u.resize(cellEnd,fuel.Vel);    
    m_u[0]=fuel.Vel;
    m_u[cellEnd-1] = 0.0;
    m_flow.resize(cellEnd,0.001);
    m_flow[0] = fuel.FlowRate;
    m_flow[cellEnd-1] = 0.0;
    if(configID == camConfig->COUNTERFLOW) {
        m_u[cellEnd-1] = oxid.Vel;
        m_flow[cellEnd-1] = oxid.FlowRate;
    }
    
   
    /*
     *setting up TDMA coefficients
     */
    std::vector<doublereal> a, b,c;
    a.resize(nCells-2,0);
    b = c = a;

    doublereal oneBydelEP, oneBydelPW;
    //1st cell
    int i = 1;
    oneBydelEP = 1.0/(0.5*(dz[i]+dz[i+1]));
    oneBydelPW = 1.0/(0.5*dz[i]);
    b[0] = oneBydelEP+oneBydelPW;
    c[0] = -oneBydelEP;
    for(i=2; i<nCells-2; i++){
        oneBydelEP = 1.0/(0.5*(dz[i]+dz[i+1]));
        oneBydelPW = 1.0/(0.5*(dz[i]+dz[i-1]));
        b[i-1] = oneBydelEP+oneBydelPW;
        a[i-1] = -oneBydelPW;
        c[i-1] = -oneBydelEP;
    }
    //last cell
    i = nCells-2;
    oneBydelEP = 1.0/(0.5*dz[i]);
    oneBydelPW = 1.0/(0.5*(dz[i]+dz[i-1]));
    b[i-1] = oneBydelEP+oneBydelPW;
    a[i-1] = -oneBydelPW;

    //store the vectors to TDMA struct
    tdmaFlow.a = a;
    tdmaFlow.b = b;
    tdmaFlow.c = c;

}
/*
 *initialize momentum
 */
void StagFlow::initMomentum(){

    /*
     *Ref : G. Stahl and J. Warnatz
     *Combustion and flame 85:285-299
     */
    m_G.resize(nCells,0.0);

}

/*
 *coupled solver
 */
void StagFlow::csolve(CamControl& cc){
    
    CVodeWrapper cvw;
    eqn_slvd = EQN_ALL;
    int band = nVar*2;
    cvw.init(nEqn,solvect,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(),
                        cc.getMaxTime(),band,*this);
    cvw.solveDAE(CV_ONE_STEP,cc.getResTol());

    reportToFile(cc.getMaxTime(),&solvect[0]);
    cvw.destroy();

}

/*
 *segregated solver
 */
void StagFlow::ssolve(CamControl& cc){

    int seg_eqn, band;
    std::vector<doublereal> seg_soln_vec;
    CVodeWrapper cvw;

    //for(int i=0; i<cc.getNumIterations(); i++){

        /*
         *integrate species
         */
        //std::cout << "Solving species: " << i << std::endl;
        eqn_slvd = EQN_SPECIES;
        seg_eqn = nSpc*nCells;
        band = nSpc*2;
        extractSpeciesVector(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(),
                        cc.getMaxTime(),band,*this);
        //cvw.setMaxStep(1e-04);
        cvw.solveDAE(CV_ONE_STEP,1e-03);
        mergeSpeciesVector(&seg_soln_vec[0]);
        reportToFile(cc.getMaxTime(),&solvect[0]);
        cvw.destroy();
//        /*
//         *Integrate energy
//         */
//        if(admin->getEnergyModel()==admin->ADIABATIC){
//            std::cout << "Solving energy: " << i << std::endl;
//            eqn_slvd = EQN_ENERGY;
//            seg_eqn = nCells;
//            band = 1;
//            extractEnergyVector(seg_soln_vec);
//            cvw.init(seg_eqn,seg_soln_vec,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(),
//                    cc.getMaxTime(),band,*this);
//            cvw.solveDAE(CV_ONE_STEP,1e-03);
//            mergeEnergyVector(&seg_soln_vec[0]);
//            cvw.destroy();
//
//        }
//

    //}


    
}
/*
 *function called by the solver (DAEs and ODEs)
 */
int StagFlow::eval(doublereal t, doublereal* y, doublereal* ydot, bool jacEval){

    residual(t,y,ydot);
    return 0;
}
/*
 *residual function
 */
void StagFlow::residual(const doublereal& t, doublereal* y, doublereal* f){
    resSp.resize(cellEnd*nSpc,0);
    resT.resize(cellEnd,0);
    resFlow.resize(cellEnd,0);
    
    /*
     *residual evaluation function for Stagflow. Internal
     *cell residuals are evaluated by calling the base class
     *function
     */    
    if(eqn_slvd == EQN_ALL){
        if(admin->getEnergyModel() == admin->ADIABATIC){
            saveMixtureProp(t,y,true,true);
            updateThermo();
        }else{
            saveMixtureProp(t,y,false,true);
        }
        
        
        updateDiffusionFluxes();
        
        speciesBoundary(t,y,&resSp[0]);
        energyBoundary(t,y,&resT[0]);


        speciesResidual(t,y,&resSp[0]);
        energyResidual(t,y,&resT[0],true);
        
        
        
        for(int i=0; i<cellEnd; i++){
            for(int l=0; l<nSpc; l++){
                f[i*nVar+l] = resSp[i*nSpc+l];
                
            }                        
            f[i*nVar+ptrT] = resT[i];            
            
        }        
    }else{
        if(eqn_slvd == EQN_SPECIES){
            mergeSpeciesVector(y);
            saveMixtureProp(t,&solvect[0],false,true);
            updateDiffusionFluxes();            
            speciesBoundary(t,y,f);
            speciesResidual(t,y,f);
        }else if(eqn_slvd==EQN_ENERGY){
            mergeEnergyVector(y);
            saveMixtureProp(t,&solvect[0],true,false);
            updateDiffusionFluxes();
            updateThermo();         
            energyResidual(t,y,f,true);
            energyBoundary(t,y,f);

        }
    }


}


/*
 *species boundary condition implementation
 */
void StagFlow::speciesBoundary(const doublereal& t, doublereal* y, doublereal* f){

    //-------------------------------------------
    //
    //  Left Boundary Settings
    //
    //------------------------------------------

    doublereal convection, diffusion;

    for(int l=0; l<nSpc ; l++){
        convection = m_u[0]*dydx(fuel.Species[l],y[l],dz[0]);
        diffusion = dydx(0,s_jk(iMesh_s,l),dz[0])/m_rho[0];
        f[l] = convection + diffusion;

    }

     //---------------------------------------------------
     //
     // Right Boundary Settings
     //
     //---------------------------------------------------
    if(configID == camConfig->STAGFLOW){
        for(int l=0; l<nSpc; l++){
            f[iMesh_e*nSpc+l] = -m_u[iMesh_e-1]*dydx(s_mf(iMesh_e,l),s_mf(iMesh_e-1,l),dz[iMesh_e]);
        }
    }else{
        for(int l=0; l<nSpc; l++){

            convection = m_u[iMesh_e]*dydx(s_mf(iMesh_e,l),oxid.Species[l],dz[iMesh_e]);
            diffusion = 0.0;//dydx(s_jk(iMesh_e,l),0,dz[iMesh_e])/m_rho[iMesh_e];
            f[iMesh_e*nSpc+l] = convection+diffusion;
        }
    }

}
/*
 *energy boundary
 */
void StagFlow::energyBoundary(const doublereal& t, doublereal* y, doublereal* f){

    /*--------------------------------------------------------------------------
     * Left Boundary:
     *--------------------------------------------------------------------------
     * The temperature is fixed at the fuel inlet temperature
     */
    f[0] = 0.0;

    /*-------------------------------------------------------------------------
     * Right Boundary:
     *---------------------------------------------------------------------
     * The temperature at the right boundary is  fixed at the surface
     * temperature in case of stagnation flow, and oxidizer inlet temperature
     * in case of counter flow flames
     */
     
     f[iMesh_e] = 0;
    
    
}

/*
 *calculate the axial velocity
 */
void StagFlow::calcFlowField(const doublereal& time, doublereal* y){

        
    /*
     *prepare
     */
    if(eqn_slvd == EQN_ALL){
        saveMixtureProp(time, y,true,true);
    }else{
        if(eqn_slvd == EQN_SPECIES){
            mergeSpeciesVector(y);
            saveMixtureProp(time, &solvect[0],false,true);
        }else if(eqn_slvd == EQN_ENERGY){
            mergeEnergyVector(y);
            saveMixtureProp(time,&solvect[0],true,true);
        }
    }

    for(int i=0; i<10; i++){
        //------------------------------------------------------------
        //           Inlet face
        //----------------------------------------------------------
        //m_flow[0] = fuel.FlowRate; // this is the vanishing cell
        m_flow[0] = fuel.Vel*m_rho[0];
        fuel.FlowRate = m_flow[0];
        //----------------------------------------------------------
        //           Exit face
        //------------------------------------------------------------
           // this is the vanishing cell
        if(configID == camConfig->STAGFLOW){
            m_flow[iMesh_e] = 0.0;
        }else{
            m_flow[iMesh_e] = oxid.Vel*m_rho[iMesh_e];
            oxid.FlowRate = m_flow[iMesh_e];
        }
        std::vector<doublereal> flow;
        calcVelocity(flow);
        for(int i=1; i<iMesh_e; i++){
            m_flow[i] = flow[i-1];            
        }
        

        m_G[0] = sqrt(m_rho[iMesh_e]/m_rho[0]);
        if(configID==camConfig->COUNTERFLOW){
            m_G[iMesh_e] = 1.0;
        }
        std::vector<doublereal> mom;
        calcMomentum(mom);
        for(int i=1; i<iMesh_e;i++){
            m_G[i] = mom[i-1];
        }
        if(configID==camConfig->STAGFLOW){
            m_G[iMesh_e] = m_G[iMesh_e-1];
        }else{
            m_G[iMesh_e] = 1.0;
        }
    }

    for(int i=0; i<iMesh_e+1; i++){
        m_u[i]=m_flow[i]/m_rho[i];        
    }
    
   
}

void StagFlow::calcVelocity(std::vector<doublereal>& u){
    
    std::vector<doublereal> r;
    u.resize(tdmaFlow.b.size(),0);
    r = u;
    /*
     *right hand side for tridiagonal matrix
     */
    //first cell
    int i = 1;
    doublereal rho_e = 0.5*(m_rho[i+1]+m_rho[i]);
    doublereal rho_w = m_rho[i-1];
    doublereal Ge = 0.5*(m_G[i]+m_G[i+1]);
    doublereal Gw = m_G[i-1];
    doublereal rhoGe = rho_e*Ge;
    doublereal rhoGw = rho_w*Gw;
    r[i-1] = 2*strainRate*(rhoGe - rhoGw) +
            (2*fuel.FlowRate/dz[i]);
    //last cell
    i = iMesh_e-1;
    rho_e = m_rho[i+1];
    rho_w = 0.5*(m_rho[i]+m_rho[i-1]);
    Ge = m_G[i+1];
    Gw = 0.5*(m_G[i]+m_G[i-1]);
    rhoGe = rho_e*Ge;
    rhoGw = rho_w*Gw;
    r[i-1] = 2*strainRate*(rhoGe-rhoGw);
    if(configID == camConfig->COUNTERFLOW){
        r[i-1] += 2*oxid.FlowRate/dz[i];
    }
    //other cells
    for( i = 2; i<iMesh_e-1; i++){
        rho_e = 0.5*(m_rho[i+1]+m_rho[i]);
        rho_w = 0.5*(m_rho[i]+m_rho[i-1]);
        Ge = 0.5*(m_G[i]+m_G[i+1]);
        Gw = 0.5*(m_G[i]+m_G[i-1]);
        rhoGe = rho_e*Ge;
        rhoGw = rho_w*Gw;
        r[i-1] = 2*strainRate*(rhoGe-rhoGw);
    }
    
    CamMath cm;//    
    cm.TDMA(tdmaFlow.a,tdmaFlow.b,tdmaFlow.c,r,u);
    
   
}
/*
 *calculate the momentum
 */
void StagFlow::calcMomentum(std::vector<doublereal>& mom){
    /*
     *this evaluates the g equation using TDMA.
     *At the stagnation plane the boundary zero gradient boundary
     *condition is implemented. At the hot edge (fuel inlet) the boundary
     * condition is g=sqrt(rho_ub/rho_b). At the cold edge the boundary
     * condition is g=1
     * (Ref: JG. Stahl and J. Warnatz Cobust. Flame 85 (1991) 285-299
     */
    /*
     *preparation
     */
    int i;
    doublereal De, fe, Dw, fw;
    doublereal delPW, delPE;
    doublereal mu_e, mu_w;
    std::vector<doublereal> a,b,c,r;
    c.resize(nCells-2,0.0);
    mom = a = b = r = c;
    /*---------------------------------------------------
     *        First cell
     *--------------------------------------------------/
     */
    
    i = 1;
    delPE = 0.5*(dz[i]+dz[i+1]);
    delPW = 0.5*dz[i];
    mu_e = 0.5*(m_mu[i]+m_mu[i+1]);
    De = mu_e/delPE;
    fe = 0.5*(m_flow[i]+m_flow[i+1]);
    b[i-1] = De + 0.5*fuel.FlowRate + m_mu[i-1]/delPW;
    c[i-1] = De - 0.5*fe ;
    r[i-1] = (m_mu[i-1]/delPW + fuel.FlowRate) *m_G[i-1] +
            strainRate*m_rho[iMesh_e]*dz[i];


    /*------------------------------------------------
     * interior cells
     *------------------------------------------------/
     */
    for(i=2; i<iMesh_e-1;i++){
        delPE = 0.5*(dz[i]+dz[i+1]);
        delPW = 0.5*(dz[i]+dz[i-1]);
        mu_e = 0.5*(m_mu[i]+m_mu[i+1]);
        mu_w = 0.5*(m_mu[i]+m_mu[i-1]);
        fe = 0.5*(m_flow[i]+m_flow[i+1]);
        fw = 0.5*(m_flow[i]+m_flow[i-1]);
        De = mu_e/delPE;
        Dw = mu_w/delPW;
        a[i-1] = -(0.5*fw+Dw);
        b[i-1] = De+Dw;
        c[i-1] = -(De-0.5*fe);
        r[i-1] = strainRate*m_rho[iMesh_e]*dz[i];

    }

    /*------------------------------------------------
     * last cell
     *------------------------------------------------/
     */
    i = iMesh_e-1;
    mu_w = 0.5*(m_mu[i]+m_mu[i-1]);
    delPW = 0.5*(dz[i]+dz[i-1]);
    Dw = mu_w/delPW;
    fw = 0.5*(m_flow[i]+m_flow[i-1]);

    b[i-1] = Dw;
    a[i-1] = -(Dw+0.5*fw);
    r[i-1] = strainRate*m_rho[iMesh_e]*dz[i];

    if(configID == camConfig->COUNTERFLOW){
        mu_e = m_mu[iMesh_e];
        delPE = 0.5*dz[i];
        fe = 0.5*oxid.FlowRate;
        De = mu_e/delPE;
        b[i-1] = Dw+De-fe;
        r[i-1] += (De-fe)*m_G[iMesh_e];

    }
    

    CamMath cm;
    cm.TDMA(a,b,c,r,mom);        

}

/*
 *update diffusion fluxes
 */
void StagFlow::updateDiffusionFluxes(){

    CamResidual::updateDiffusionFluxes();
    for(int l=0; l<nSpc; l++)
        s_jk(iMesh_e+1,l)=0.0;
    /*
     *flux at the oxidizer inlet
     */
//    doublereal delta = 0.5*(dz[iMesh_e]+dz[iMesh_e-1]);
//    std::vector<doublereal> flx;
//    flx.resize(nSpc,0);
//    //preperation for flux correction
//    doublereal jCorr = 0;
//    for(int l=0; l<nSpc; l++){
//        doublereal grad = dydx(s_mf(iMesh_e,l),s_mf(iMesh_e-1,l),delta);
//        doublereal avgRho = 0.5*(m_rho[iMesh_e]+m_rho[iMesh_e-1]);
//        doublereal avgD = (s_Diff(iMesh_e,l)+s_Diff(iMesh_e-1,l))/2.0;
//        flx[l] = -avgD*avgRho*grad;
//
//        jCorr += flx[l];
//
//    }
//   //correction
//    for(int l=0; l<nSpc; l++){
//        flx[l] -= s_mf(iMesh_e,l)*jCorr;
//        s_jk(iMesh_e,l) = flx[l];
//    }

}
/*
 *report functions
 */
void StagFlow::report(doublereal t, doublereal* solution){
    
}

void StagFlow::report(doublereal t, doublereal* solutio, doublereal& res){
    static int nStep=0;
    std::cout.width(5);
    std::cout.setf(std::ios::scientific);
    if(nStep%10==0) reporter->consoleHead("time(s) \t residual");
    std::cout << t <<"\t" << res << std::endl;
    nStep++;
    
    
}

void StagFlow::header(){
    headerData.clear();
    headerData.push_back("int_time");
    headerData.push_back("x");
    headerData.push_back("rho");
    headerData.push_back("u");
    headerData.push_back("G");
    headerData.push_back("T");
    for (int l = 0; l < nSpc; l++) {
        headerData.push_back( (*spv)[l]->Name() );
    }
    headerData.push_back("sumfracs");

}

void StagFlow::reportToFile(doublereal t, doublereal* soln){

    saveMixtureProp(t,soln,false,true);
    doublereal sum;
    reporter->openFiles();
    reporter->writeHeader(headerData);
    std::vector<doublereal> data, axpos;
    std::vector<doublereal> molfrac, massfrac;
    axpos = reacGeom->getAxpos();
    int len = axpos.size();
    for(int i=0; i<len; i++){
        data.clear();
        data.push_back(t);
        data.push_back(axpos[i]);
        data.push_back(m_rho[i]);
        data.push_back(m_u[i]);
        data.push_back(m_G[i]);
        data.push_back(soln[i*nVar+ptrT]);

        massfrac.clear();
        molfrac.clear();
        for(int l=0; l<nSpc; l++){
            massfrac.push_back(soln[i*nVar+l]);
        }
        if(admin->getSpeciesOut() == admin->MASS){
            sum = 0;
            for(int l=0; l<nSpc; l++){
                data.push_back(fabs(massfrac[l]));
                sum += massfrac[l];
            }
        }else{
            CamConverter cc;
            cc.mass2mole(massfrac,molfrac,*camMech);
            sum = 0;
            for(int l=0; l<nSpc; l++){
                data.push_back(fabs(molfrac[l]));
                sum += molfrac[l];
            }
        }
        data.push_back(sum);
        reporter->writeStdFileOut(data);

    }

    reporter->closeFiles();
}
