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
#include <stdlib.h>
#include <vector>

using namespace Camflow;
/*
 *solve the stagnation/twinflame model
 */
void StagFlow::solve(CamControl& cc, CamAdmin& ca, CamGeometry& cg,
                            CamProfile& cp, Mechanism& mech){
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

    
    admin = &ca;
    reacGeom = &cg;
    reacGeom->discretize();
    /*
     * 2 additional cells are padded to consider the
     * inlet and the exhaust
     */

    reacGeom->addZeroWidthCells();
    nCells = reacGeom->getnCells();

    opPre = ca.getPressure();
    profile->setGeometryObj(cg);
    reporter = new CamReporter();
    nSpc = camMech->SpeciesCount();

    /*
     *array offsets for newton solver
     */
    ptrF = 0;
    /*
     *array offsets for ODE
     */
    ptrT = nSpc;
    ptrG = ptrT+1;
    
    /*
     *number of finite volume cells and
     *loop initializations
     */
    
    iMesh_s = 1;
    iMesh_e = nCells-1;

    cellBegin = 0;
    cellEnd = nCells;
    


    /*
     *number of variables and the number of equations
     *to be solved by the newton solver
     */
    alg_nVar = 1; // F(x)
    alg_nEqn =  alg_nVar*nCells;
    alg_band = alg_nVar*2;
    

    /*
     *number of equations and the number of variables to
     *be slved by the ODE
     */
    nVar = nSpc+2;
    nEqn = nVar*nCells;

    /*
     *init the solution vector
     */
    initSolutionVector(cc);

    /*
     *newton solver pointer
     */
    newton = new KinsolWrapper();
    tol_res = cc.getResTol();

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
    admin->getRightBoundary(right);

    if(left.getVelocity()==0)
        throw CamError("Stagflow error: inlet velocity not specified\n");

    /*
     *store the inlets and change the
     *sign of flow for the oxidizer
     */
    storeInlet(left,fuel);
    storeInlet(right,oxid);

    oxid.FlowRate *= -1;
    oxid.Vel *= -1;
    right.setVelocity(-right.getVelocity());
    right.setFlowRate(-right.getFlowRate());


    /*
     *initialize the ODE vector
     */
    solvect.resize(nEqn,0);
    vector<doublereal> vSpec, vT, vMom, vU;
    initSpecies(left,right,cc,vSpec);
    initMomentum(fuel.FlowRate,oxid.FlowRate,vMom);
    initTemperature(left,cc,vT);
    //initMassFlow(fuel.Vel, oxid.Vel,vU);

    mergeEnergyVector(&vT[0]);
    mergeSpeciesVector(&vSpec[0]);
    mergeMomentum(&vMom[0]);
    //mergeContinuity(&vU[0]);
    //mergeContinuity(&m_u[0]);

    /*
     *initialize the algebraic equations vector
     */
    alg_solvect.resize(alg_nEqn,0.0);
    initMassFlow(fuel.Vel, oxid.Vel,m_u);
    /*
     *merge the dependent vaiables into
     *algebraic solution vector
     */
    for(int i=0; i<cellEnd; i++){
        alg_solvect[i*alg_nVar+ptrF] = m_u[i];

    }

}
/*
 *initialize the mass flow
 */
void StagFlow::initMassFlow(const doublereal fLeft, const doublereal fRight,
                            vector<doublereal>& soln){
    soln.resize(cellEnd,0);
    doublereal flowGrad = (fRight-fLeft)/reacGeom->getLenth();
    doublereal xPos=0;
    for(int i=0; i<cellEnd; i++){
        xPos += dz[i];
        soln[i] = flowGrad*(xPos-0.0) + fLeft;
    }

    //readVelocity();

}
/*
 *initialize momentum
 */
void StagFlow::initMomentum(const doublereal fLeft, const doublereal fRight,
                            vector<doublereal>& soln){

    
    doublereal flowGrad = 0.01;//(fRight-fLeft)/reacGeom->getLenth();
    soln.resize(cellEnd,flowGrad);
    soln[0] = fuel.rVelGrad;
    soln[iMesh_e] = oxid.rVelGrad;

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
    cvw.solve(CV_ONE_STEP,cc.getResTol());
    cvw.destroy();

}

/*
 *segregated solver
 */
void StagFlow::ssolve(CamControl& cc){

    int seg_eqn, band;
    vector<doublereal> seg_soln_vec;
    CVodeWrapper cvw;

    for(int i=0; i<cc.getNumIterations(); i++){
        /*
         *solve momentum
         */
        cout << "solving momentum\n";
        eqn_slvd = nCells;
        band = 1;
        extractMomentum(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,cc.getFlowAbsTol(),cc.getFlowRelTol(),
                        cc.getMaxTime(),band,*this);
        cvw.solve(CV_ONE_STEP,1e-03);
        cvw.destroy();

        /*
         *integrate species
         */
        cout << "Solving species: " << i << endl;
        eqn_slvd = EQN_SPECIES;
        seg_eqn = nSpc*nCells;
        band = nSpc*2;
        extractSpeciesVector(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(),
                        cc.getMaxTime(),band,*this);
        //cvw.setMaxStep(1e-04);
        cvw.solve(CV_ONE_STEP,1e-03);
        mergeSpeciesVector(&seg_soln_vec[0]);
        reportToFile(cc.getMaxTime(),&solvect[0]);
        cvw.destroy();
        /*
         *Integrate energy
         */
        if(admin->getEnergyModel()==admin->ADIABATIC){
            cout << "Solving energy: " << i << endl;
            eqn_slvd = EQN_ENERGY;
            seg_eqn = nCells;
            band = 1;
            extractEnergyVector(seg_soln_vec);
            cvw.init(seg_eqn,seg_soln_vec,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(),
                    cc.getMaxTime(),band,*this);
            cvw.solve(CV_ONE_STEP,1e-03);
            mergeEnergyVector(&seg_soln_vec[0]);
            cvw.destroy();

        }


    }


    
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
    resAxVel.resize(cellEnd,0);
    /*
     *residual evaluation function for Stagflow. Internal
     *cell residuals are evaluated by calling the base class
     *function
     */
    if(eqn_slvd == EQN_ALL){
        if(admin->getEnergyModel() == admin->ADIABATIC){
            saveMixtureProp(y,true,true);
            updateThermo();
        }else{
            saveMixtureProp(y,false,true);
        }
        updateMomentum(y);        
        updateDiffusionFluxes();


        solveAES();
        speciesBoundary(t,y,&resSp[0]);
        energyBoundary(t,y,&resT[0]);


        speciesResidual(t,y,&resSp[0]);
        energyResidual(t,y,&resT[0]);
        
        momentumResidual(t,y,&resFlow[0]);
        
        for(int i=0; i<cellEnd; i++){
            for(int l=0; l<nSpc; l++){
                f[i*nVar+l] = resSp[i*nSpc+l];
                
            }                        
            f[i*nVar+ptrT] = resT[i];
            f[i*nVar+ptrG] = resFlow[i];
            
        }        
    }else{
        if(eqn_slvd == EQN_SPECIES){
            mergeSpeciesVector(y);
            saveMixtureProp(&solvect[0],false,true);
            updateFlow(&solvect[0]);
            updateMomentum(&solvect[0]);
            updateDiffusionFluxes();
            //solveAES();
            speciesBoundary(t,y,f);
            speciesResidual(t,y,f);
        }else if(eqn_slvd == EQN_CONTINUITY) {
            mergeContinuity(y);
            saveMixtureProp(&solvect[0],false,true);
            updateFlow(&solvect[0]);
            continuity(t,y,f);
        }else if(eqn_slvd == EQN_MOMENTUM){
            mergeMomentum(y);
            saveMixtureProp(&solvect[0],false,true);
            updateMomentum(&solvect[0]);
            //solveAES();
            momentumResidual(t,y,f);
        }else if(eqn_slvd==EQN_ENERGY){
            mergeEnergyVector(y);
            saveMixtureProp(&solvect[0],true,false);
            updateMomentum(&solvect[0]);
            updateDiffusionFluxes();
            updateThermo();
            //solveAES();
            energyResidual(t,y,f);
            energyBoundary(t,y,f);

        }
    }


}
/*
 *function called by newton solver
 */
int StagFlow::eval(doublereal* y, doublereal* ydot){
    /*
     *save the dependent variables
     */
    saveAxVel(y);
    /*
     *mass flow residual
     */    
    continuity(0.0,y,ydot);
    return 0;
}


void StagFlow::solveAES(){

//    for(int i=0; i<cellEnd; i++)
//        alg_solvect[i] = m_u[i];
    newton->init(alg_nEqn,alg_solvect,tol_res,alg_band,*this);
    newton->solve();    
    saveNewton();
    newton->destroy();
//    for(int i=0; i<alg_nEqn; i++)
//        cout << alg_solvect[i] << endl;
    
}
/*
 *save newton results
 */
void StagFlow::saveNewton(){
    for(int i=0; i<cellEnd; i++){
        m_u[i] = alg_solvect[i];        
    }
}
/*
 *mass flow residual: node storage
 */
void StagFlow::continuity(const doublereal& time, doublereal* y,
                                                doublereal* f){
/*
 *This is for Newton Solver
 */
    //---------------------------------------------------------------
    //  Fuel inlet
    //---------------------------------------------------------------
    f[0] = fuel.Vel - m_u[0];
    //---------------------------------------------------------------
    //  Interior mesh
    //---------------------------------------------------------------
    doublereal grad, phiP,phiE,phiW;
    for(int i=iMesh_s; i<iMesh_e; i++){
        phiP = m_rho[i]*m_u[i];
        phiE = m_rho[i+1]*m_u[i+1];
        phiW = m_rho[i-1]*m_u[i-1];

        grad = (m_u[i] > 0) ? (phiP-phiW)/dz[i] : (phiE-phiP)/dz[i];
        //grad = (phiE-phiP)/dz[i];

        f[i] = grad + 2*m_rho[i]*m_G[i];

    }

    //---------------------------------------------------------------
    //   Last cell
    //---------------------------------------------------------------
    f[iMesh_e] = oxid.Vel - m_u[iMesh_e];

/*
 *ODE implementation
 */
//    //---------------------------------------------------------------
//    //  Fuel inlet
//    //---------------------------------------------------------------
//    f[0] = 0;
//    //---------------------------------------------------------------
//    //  Interior mesh
//    //---------------------------------------------------------------
//    doublereal grad, sh_e, sh_w;
//    doublereal muAvg_e,muAvg_w, delta;
//    doublereal rhs1, rhs2, rhs3, rhs4;
//    doublereal fourthird = 4.0/3.0;
//    for(int i=iMesh_s; i<iMesh_e; i++){
//        //east
//        muAvg_e = 0.5*(m_mu[i]+m_mu[i+1]);
//        delta = 0.5*(dz[i]+dz[i+1]);
//        sh_e = muAvg_e*(m_u[i+1]-m_u[i])/delta;
//
//        muAvg_w = 0.5*(m_mu[i]+m_mu[i-1]);
//        delta = 0.5*(dz[i]+dz[i-1]);
//        sh_w = muAvg_w*(m_u[i]-m_u[i-1])/delta;
//
//        doublereal Ge = 0.5*(m_G[i]+m_G[i+1]);
//        doublereal Gw = 0.5*(m_G[i]+m_G[i-1]);
//
//        doublereal muGe = muAvg_e*Ge;
//        doublereal muGw = muAvg_w*Gw;
//
//        rhs1 = (fourthird/m_rho[i])*(muGe-muGw)/dz[i];
//
//        rhs2 = (2*m_mu[i]/m_rho[i])*(Ge-Gw)/dz[i];
//
//        rhs3 = (fourthird/m_rho[i])*(sh_e-sh_w)/dz[i];
//
//        rhs4 = (m_u[i] > 0) ? (m_u[i] -m_u[i-1]) : (m_u[i+1] - m_u[i]);
//        rhs4 *= (m_u[i]/dz[i]);
//
//
//        f[i] = -rhs1 + rhs2 + rhs3 - rhs4;
//
//    }
//
//    //---------------------------------------------------------------
//    //   Last cell
//    //---------------------------------------------------------------
//    f[iMesh_e] = 0;
//


}
/*
 *momentum residual
 */
void StagFlow::momentumResidual(const doublereal& time, doublereal* y,
                                                     doublereal* f){

    eigen = -m_rho[iMesh_e]*100*100;
    //---------------------------------------------------------------
    //  Inlet cell
    //---------------------------------------------------------------
    f[0] = 0;

    //---------------------------------------------------------------
    //  Interior mesh
    //---------------------------------------------------------------
    doublereal shear, velGrad,pgrad;
    for(int i=iMesh_s; i<iMesh_e; i++){

        shear = (m_shear[i+1]-m_shear[i])/(m_rho[i]*dz[i]);
        //velGrad = m_u[i]*(m_G[i+1]-m_G[i-1])/(2*dz[i]);
        velGrad = ( m_u[i] > 0 )? (m_G[i]-m_G[i-1]) : (m_G[i+1]-m_G[i]) ;
        velGrad *= (m_u[i]/dz[i]);
        pgrad = eigen/m_rho[i];
        f[i] = shear - velGrad - (m_G[i]*m_G[i]) - pgrad;
    }
    //---------------------------------------------------------------
    //  Last cell
    //---------------------------------------------------------------
    f[iMesh_e] = 0;
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

    for(int l=0; l<nSpc; l++){

        convection = m_u[iMesh_e]*dydx(s_mf(iMesh_e,l),oxid.Species[l],dz[iMesh_e]);
        diffusion = dydx(s_jk(iMesh_e,l),0,dz[iMesh_e])/m_rho[iMesh_e];
        f[iMesh_e*nSpc+l] = convection+diffusion;
    }

}
/*
 *energy boundary
 */
void StagFlow::energyBoundary(const doublereal& t, doublereal* y, doublereal* f){
    //-------------------------------------------
    //
    //  Left Boundary Settings
    //
    //------------------------------------------

    f[0] = 0.0;

     //---------------------------------------------------
     //
     // Right Boundary Settings
     //
     //---------------------------------------------------
    f[iMesh_e] = 0;
    
}
/*
 *calculate the eigen value pressure gradient
 */
void StagFlow::calcEigenValuePGad(){
    eigen = -(m_shear[iMesh_e]/dz[iMesh_e]) -
            (m_G[iMesh_e]*m_G[iMesh_e]*m_rho[iMesh_e]) -
        (m_u[iMesh_e]*m_rho[iMesh_e]*(oxid.rVelGrad-m_G[iMesh_e])/dz[iMesh_e]);
   
}
/*
 *save axial velocity
 */
void StagFlow::saveAxVel(doublereal* y){
    for(int i=0; i<cellEnd; i++){
        m_u[i] = y[i*alg_nVar+ptrF];
    }
}
/*
 *update flow
 */
void StagFlow::updateFlow(doublereal* y){
    m_u.resize(cellEnd,0);
    m_G.resize(cellEnd,0);
    for(int i=0; i<cellEnd; i++){
        m_u[i] = y[i*nVar+ptrF];
        m_G[i] = y[i*nVar+ptrG];
    }
}
/*
 *update momentum
 */
void StagFlow::updateMomentum(doublereal* y){
    //save the radial velocity gradient (Vr/r)
    m_G.resize(cellEnd,0);
    for(int i=0; i<cellEnd; i++){
        m_G[i] = y[i*nVar+ptrG];
    }
    //shear calculation
    m_shear.resize(cellEnd+1,0);
    doublereal delta,muAvg, grad;
    for(int i= iMesh_s; i<iMesh_e; i++){
        
        delta = (dz[i]+dz[i-1])/2.0;
        muAvg = (m_mu[i]+m_mu[i-1])/2.0;
        if((i-1) == 0)muAvg = m_mu[i];
        
        grad = dydx(m_G[i],m_G[i-1],delta);
        m_shear[i] = muAvg*grad;
        
    }

    /*
     *oxidizer inlet
     */
    delta = 0.5*(dz[iMesh_e]+dz[iMesh_e-1]);
    muAvg = (m_mu[iMesh_e]+m_mu[iMesh_e-1])/2.0;
    grad = dydx(m_G[iMesh_e],m_G[iMesh_e-1],delta);
    m_shear[iMesh_e] = muAvg*grad;
    
}
/*
 *update diffusion fluxes
 */
void StagFlow::updateDiffusionFluxes(){

    CamResidual::updateDiffusionFluxes();
    /*
     *flux at the oxidizer inlet
     */
    doublereal delta = 0.5*(dz[iMesh_e]+dz[iMesh_e-1]);
    vector<doublereal> flx;
    flx.resize(nSpc,0);
    //preperation for flux correction
    doublereal jCorr = 0;
    for(int l=0; l<nSpc; l++){
        doublereal grad = dydx(s_mf(iMesh_e,l),s_mf(iMesh_e-1,l),delta);
        doublereal avgRho = 0.5*(m_rho[iMesh_e]+m_rho[iMesh_e-1]);
        doublereal avgD = (s_Diff(iMesh_e,l)+s_Diff(iMesh_e-1,l))/2.0;
        flx[l] = -avgD*avgRho*grad;

        jCorr += flx[l];

    }
   //correction
    for(int l=0; l<nSpc; l++){
        flx[l] -= s_mf(iMesh_e,l)*jCorr;
        s_jk(iMesh_e,l) = flx[l];
    }

}
/*
 *report functions
 */
void StagFlow::report(doublereal t, doublereal* solution){
    
}

void StagFlow::report(doublereal t, doublereal* solutio, doublereal& res){
    static int nStep=0;
    cout.width(5);
    cout.setf(ios::scientific);
    if(nStep%10==0) reporter->consoleHead("time(s) \t residual");
    cout << t <<"\t" << res << endl;
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

}

void StagFlow::reportToFile(doublereal t, doublereal* soln){

    saveMixtureProp(soln,false,true);
    doublereal sum;
    reporter->openFiles();
    reporter->writeHeader(headerData);
    vector<doublereal> data, axpos;
    vector<doublereal> molfrac, massfrac;
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
        //data.push_back(sum);
        reporter->writeStdFileOut(data);

    }

    reporter->closeFiles();
}

void StagFlow::readVelocity(){
    vector<doublereal> grid;
    ifstream inf;
    string velFile = "vel.inp";
    inf.open(velFile.c_str(),ios::in);
    if(inf.good()){
        std::string position;

        while(!inf.eof()){
            getline(inf,position);
            if(! isEmpty(position)){
                grid.push_back(cdble(position));
            }
        }
    }

    inf.close();
    int len = grid.size();
    cout << len << endl;
    m_u.resize(cellEnd,0);
    m_u[0] = grid[0]/100.0;
    for(int i=1; i<len; i++){
        //cout << i << "  " << grid[i-1] << "  " << grid[i] << endl;
        m_u[i] = 0.5*(grid[i-1]+grid[i])/100;
    }
    m_u[len] = grid[len-1]/100;

//    for(int i=0; i<= iMesh_e; i++)
//        cout << i << "  " << m_u[i] << endl;
}