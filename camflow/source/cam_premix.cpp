/*
 * File:   cam_premix.cpp
 * Author: vinod
  * File purpose:
 *  This class implements the premix model
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
 *  Website
 *
 * Created on January 18, 2009, 7:57 PM
 */


#include <vector>
#include <map>
#include <math.h>
#include <cstring>

#include "cam_soot.h"
#include <iostream>
#include <sstream>
#include "cam_params.h"
#include "cam_admin.h"
#include "cam_profile.h"
#include "cam_residual.h"
#include "cam_setup.h"
#include "cam_premix.h"
#include "cam_reporter.h"
#include "kinsol_wrapper.h"
#include "cvode_wrapper.h"
//#include "ida_wrapper.h"
using namespace Camflow;
using namespace std;

int CamPremix::eval(doublereal t, doublereal* y, doublereal* ydot, bool jacEval){

    /*
     *this is called by the DAE wrapper object. Given y, ydot is returned
     */
    residual(t,y,ydot);
    return 0;
}


void CamPremix::residual(const doublereal& t, doublereal* y, doublereal* f){

    resSp.resize(cellEnd*nSpc,0);
    resT.resize(cellEnd,0);
    resFlow.resize(cellEnd,0);
    /*
     *residual evaluation function for premix. Internal
     *cell residuals are evaluated by calling the base class
     *function
     */
    if(eqn_slvd == EQN_ALL){
        if(admin->getEnergyModel() == admin->ADIABATIC){
            saveMixtureProp(y,true,false);
            updateThermo();
        }else{
            saveMixtureProp(y,false,false);
        }
        saveFlowVariables(y);
        updateDiffusionFluxes();


        //massFlowBoundary(t,y,&resM[0]);
        speciesBoundary(t,y,&resSp[0]);
        energyBoundary(t,y,&resT[0]);

        speciesResidual(t,y,&resSp[0]);
        energyResidual(t,y,&resT[0]);
        //massFlowResidual(t,y,&resM[0]);

        for(int i=0; i<cellEnd; i++){
            for(int l=0; l<nSpc; l++){
                f[i*nVar+l] = resSp[i*nSpc+l];
            }
            f[i*nVar+ptrF] = resFlow[i];
            f[i*nVar+ptrT] = resT[i];
        }

        //soot moments
        if(sootMom->active()){
            resMoment.resize(cellEnd*nMoments,0);
            sootMom->momentResidual(t,iMesh_s,iMesh_e,nVar,nSpc,dz,m_u,m_rho,y,&resMoment[0]);
            momentBoundary(t,y,&resMoment[0]);

            for(int i=0; i<cellEnd; i++){
//                for(int l=nSpc; l<(nSpc+nMomemts); l++){
//                    f[i*nVar+l] = resMoment[l-nSpc];
//
//                }
                for(int l=0; l<nMoments; l++){
                    f[i*nVar+l+nSpc] = resMoment[i*nMoments+l];
                }
            }
        }


    }else{
        if(eqn_slvd == EQN_SPECIES){
            mergeSpeciesVector(y);
            saveMixtureProp(&solvect[0],false,false);
            saveFlowVariables(&solvect[0]);
            updateDiffusionFluxes();
            speciesBoundary(t,y,f);
            speciesResidual(t,y,f);
        }else if(eqn_slvd == EQN_CONTINUITY){
            mergeContinuity(y);
            saveMixtureProp(&solvect[0],false,false);
            saveFlowVariables(&solvect[0]);
            massFlowBoundary(t,y,f);
            massFlowResidual(t,y,f);
        }else if(eqn_slvd==EQN_ENERGY){
            mergeEnergyVector(y);
            saveMixtureProp(&solvect[0],true,false);
            saveFlowVariables(&solvect[0]);
            updateDiffusionFluxes();
            updateThermo();
            energyResidual(t,y,f);
            energyBoundary(t,y,f);
        }else if(eqn_slvd==EQN_MOMENTS){
            mergeSootMoments(y);
            saveMixtureProp(&solvect[0],false,false);
            saveFlowVariables(&solvect[0]);
            momentBoundary(t,y,f);
            sootMom->momentResidual(t,iMesh_s,iMesh_e,dz,m_u,m_rho,y,f);
        }
    }

}


void CamPremix::massFlowBoundary(const doublereal& t, doublereal* y, doublereal* f){
    //f[0] = ud_inlet.FlowRate - y[0]
    //f[iMesh_e] = y[iMesh_e-1] - y[iMesh_e];

    f[0] = ud_inlet.Vel*(ud_inlet.FlowRate -m_flow[0])/dz[0];

    /*
     *Mass flow
     */

    f[iMesh_e] = -m_u[iMesh_e]*(m_flow[iMesh_e]-m_flow[iMesh_e-1])/dz[iMesh_e];

}

void CamPremix::speciesBoundary(const doublereal& t, doublereal* y, doublereal* f){

    //-------------------------------------------
    //
    //  Left Boundary Settings
    //
    //------------------------------------------

    doublereal convection, diffusion;

    for(int l=0; l<nSpc ; l++){
        convection = m_u[0]*dydx(ud_inlet.Species[l],y[l],dz[0]);
        diffusion = dydx(0,s_jk(iMesh_s,l),dz[0])/m_rho[0];
        f[l] = convection + diffusion;

    }

     //---------------------------------------------------
     //
     // Right Boundary Settings
     //
     //---------------------------------------------------

    for(int l=0; l<nSpc; l++){
        f[iMesh_e*nSpc+l] = -m_u[iMesh_e]*dydx(s_mf(iMesh_e,l),s_mf(iMesh_e-1,l),dz[iMesh_e]);
    }

}

/*
 *soot moment boundary conditions
 */
void CamPremix::momentBoundary(const doublereal& t, doublereal* y, doublereal* f){
    //-------------------------------------------
    //
    //  Left Boundary Settings
    //
    //------------------------------------------
    for(int l=0; l<nMoments; l++){
        f[l] = 0;
    }
    //-------------------------------------------
    //
    //  Right Boundary Settings
    //
    //------------------------------------------
    doublereal convection;
    for(int l=0; l<nMoments; l++){
        doublereal phi_e = y[iMesh_e*nVar+nSpc+l];
        doublereal phi_w = y[(iMesh_e-1)*nVar+nSpc+l];
        convection = -m_u[iMesh_e]*(phi_e-phi_w)/dz[iMesh_e];
        f[iMesh_e*nMoments+l] = convection;
    }

}

void CamPremix::energyBoundary(const doublereal& t, doublereal* y, doublereal* f){
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
    if(admin->getEnergyModel() == admin->ADIABATIC){

        f[iMesh_e] = -m_u[iMesh_e]*(dydx(m_T[iMesh_e],m_T[iMesh_e-1],dz[iMesh_e]));

    }else{

        f[iMesh_e] = 0;
    }

}



void CamPremix::solve(CamControl& cc, CamAdmin& ca, CamGeometry& cg,
        CamProfile& cp, CamSoot &cs, Mechanism& mech){
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


    CamBoundary cb;
    admin = &ca;
    reacGeom = &cg;
    admin->getLeftBoundary(cb);
    sootMom = &cs;
    reacGeom->discretize();
    /*
     * 2 additional cells are padded to consider the
     * inlet and the exhaust
     */

    reacGeom->addZeroWidthCells();

    opPre = ca.getPressure();
    profile->setGeometryObj(cg);
    reporter = new CamReporter();
    nSpc = camMech->SpeciesCount();

    nVar = nSpc + 2;
    /*
     *check whether to solve for the soot moments
     */
    if(cs.active()){
        nMoments = cs.getNumMoments();
        nVar += nMoments;
        ptrF = nSpc + nMoments;
        ptrT = ptrF + 1;
    }else{
        ptrF = nSpc;
        ptrT = ptrF + 1;
    }

    nEqn = nVar *(reacGeom->getnCells());

    /*
     *this loop is for residual evaluation by the base class
     */
    iMesh_s = 1;
    iMesh_e = reacGeom->getnCells()-1;
    cellBegin = 0;
    cellEnd = reacGeom->getnCells();

    initSolutionVector(cb,cc);

    reporter->header("PREMIX");
    reporter->problemDescription(cb,*this);
    header();

    /*
     *solution call
     */
    if(cc.getSolutionMode() == cc.COUPLED)
        csolve(cc);
    else{
        ssolve(cc);
        csolve(cc);
    }


}

void CamPremix::initSolutionVector(CamBoundary &cb, CamControl &cc){

    

    dz = reacGeom->getGeometry();

    //admin->getLeftBoundary(cb);
    opPre = admin->getPressure();

//    /*
//     *assert the inlet flow rate specification
//     */
//    if(cb.getFlowRate() == 0)
//        throw CamError("premix error: inlet flow rate not specified\n");

    storeInlet(cb,ud_inlet);
    /*
     *store the inlet species mass fraction, and flow rate for
     *use in the boundary residual definitions
     */

    //createSolnVector(cb,cc,solvect);
    solvect.resize(nEqn,0);
    vector<doublereal> vSpec, vMass, vT, vSoot;
    initMassFlow(cb,cc,vMass);
    initSpecies(cb,cc,vSpec);
    initTemperature(cb,cc,vT);

    mergeEnergyVector(&vT[0]);
    mergeContinuity(&vMass[0]);
    mergeSpeciesVector(&vSpec[0]);

    if(sootMom->active()){
        //nMoments = sootMom->getNumMoments();
        vector<doublereal> temp;
        sootMom->initMoments(*camMech,temp,cellEnd);
        vSoot.resize(nMoments*cellEnd,0.0);
        for(int i=0; i<cellEnd; i++){
            for(int l=0; l<nMoments; l++){
                vSoot[i*nMoments+l] = temp[l];
            }
        }
        mergeSootMoments(&vSoot[0]);
    }

   

}
/*
 *couples solver
 */
void CamPremix::csolve(CamControl& cc){

    cout << "\nPremix coupled solver\n";
    eqn_slvd = EQN_ALL;
    int solver = cc.getSolver();
    int band = nVar*2;

    if(solver == cc.CVODE){
        CVodeWrapper cvw;
        cvw.init(nEqn,solvect,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(), cc.getMaxTime(),band,*this);
        if(cc.getResidualMonitor()){

            cvw.solve(CV_ONE_STEP,cc.getResTol());
        }else{
            resNorm = cvw.solve(CV_ONE_STEP);
            cout << "ressidual at the end of maxTime " << resNorm << endl;

        }
        reportToFile(cc.getMaxTime(),&solvect[0]);
        cvw.destroy();
    }
//    if(solver == cc.IDA){
//        IDAWrapper ida;
//        ida.init(nEqn,solvect,cc.getResTol(),cc.getMaxTime(),band,*this);
//        ida.solve();
//        reportToFile(cc.getMaxTime(),&solvect[0]);
//
//    }

}

/*
 *segregated solver
 */
void CamPremix::ssolve(CamControl& cc){
    /*
     *preperations
     */
    int nCell = reacGeom->getnCells();
    int seg_eqn, band;
    vector<doublereal> seg_soln_vec;
    CVodeWrapper cvw;
    for(int i=0; i<cc.getNumIterations(); i++){
        /*
         *integrate mass flow
         */
        cout << "solving mass flow: " << i << endl;
        eqn_slvd = EQN_CONTINUITY;
        seg_eqn = nCell;
        band = 1;
        extractContinuity(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(), cc.getMaxTime(),band,*this);
        cvw.solve(CV_ONE_STEP,1e-03);
        mergeContinuity(&seg_soln_vec[0]);
        cvw.destroy();
        /*
         *integrate soot moemnts
         */
        if(sootMom->active()){
            /*
             *integrate moments
             */
            cout << "Solving moments: " << i << endl;
            eqn_slvd = EQN_MOMENTS;
            seg_eqn = nCell*nMoments;
            band = nMoments*2;
            extractSootMoments(seg_soln_vec);
            cvw.init(seg_eqn,seg_soln_vec,1e-03,1e-03,cc.getMaxTime(),band,*this);
            cvw.solve(CV_ONE_STEP,1e-03);
            mergeSootMoments(&seg_soln_vec[0]);
            cvw.destroy();
        }

        /*
         *integrate species
         */
        cout << "Solving species: " << i << endl;
        eqn_slvd = EQN_SPECIES;
        seg_eqn = nSpc*nCell;
        band = nSpc*2;
        extractSpeciesVector(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(), cc.getMaxTime(),band,*this);
        cvw.solve(CV_ONE_STEP,1e-03);
        mergeSpeciesVector(&seg_soln_vec[0]);
        cvw.destroy();

        if(admin->getEnergyModel()==admin->ADIABATIC){
            /*
             *integrate energy
             */
            cout << "Solving energy: " << i << endl;
            eqn_slvd = EQN_ENERGY;
            seg_eqn = nCell;
            band = 1;
            extractEnergyVector(seg_soln_vec);
            cvw.init(seg_eqn,seg_soln_vec,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(), cc.getMaxTime(),band,*this);
            cvw.solve(CV_ONE_STEP,1e-03);
            mergeEnergyVector(&seg_soln_vec[0]);
            cvw.destroy();
        }



    }
}
/*
 *console putput
 */
void CamPremix::report(doublereal t, doublereal* soln){
    static int nSteps=0;
    cout.width(5);
    cout.setf(ios::scientific);
    cout << t << endl;
    if(nSteps%10==0)reporter->consoleHead("time(s)");
    nSteps++;
}
/*
 *console output with residual monitoring
 */
void CamPremix::report(doublereal t, doublereal* soln, doublereal& res){
    static int nStep=0;
    cout.width(5);
    cout.setf(ios::scientific);
    if(nStep%10==0) reporter->consoleHead("time(s) \t residual");
    cout << t <<"\t" << res << endl;
    nStep++;
}


void CamPremix::reportToFile(doublereal t, doublereal* soln){
    saveMixtureProp(soln,false,false);
    doublereal sum =0;

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
        data.push_back(axpos[i]/m_u[i]);
        data.push_back(m_rho[i]);
        data.push_back(m_u[i]);
        data.push_back(soln[i*nVar+ptrF]);
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

        //soot moments
        if(sootMom->active()){
            for(int l=nSpc; l<(nMoments+nSpc); l++){
                data.push_back(soln[i*nVar+l]);
            }
        }
        reporter->writeStdFileOut(data);

    }

    reporter->closeFiles();

    /*
     *moment output
     */
    //sootMom->report(len);

}

//return the initial solution vector
void CamPremix::getInitial(vector<doublereal>& initial){

    initial = solvect;
}


/*
 *update the diffusion fluxes. The east side of
 *last cell flux need to be set as zero
 */

void CamPremix::updateDiffusionFluxes(){


    /*
     *call the base class function and then
     *pad the species outlet boundary condition
     *at the exit
     */
    CamResidual::updateDiffusionFluxes();
    for(int l=0; l<nSpc; l++)
        s_jk(iMesh_e+1,l)=0.0;
}

/*
 *update the thermal fluxes by calling
 *the base class function
 */
void CamPremix::updateThermo(){
    CamResidual::updateThermo();

}
/*
 *save flow variables
 */
void CamPremix::saveFlowVariables(doublereal* y){
    m_u.clear();
    m_flow.clear();
    for(int i=cellBegin; i<cellEnd; i++){
        m_u.push_back(y[i*nVar+ptrF]/m_rho[i]);
        m_flow.push_back(y[i*nVar+ptrF]);
    }

}

void CamPremix::header(){
    headerData.clear();
    headerData.push_back("int_time");
    headerData.push_back("x");
    headerData.push_back("tau");
    headerData.push_back("rho");
    headerData.push_back("u");
    headerData.push_back("mdot");
    headerData.push_back("T");
    for (int l = 0; l < nSpc; l++) {
        headerData.push_back( (*spv)[l]->Name() );
    }
    headerData.push_back("sumfracs");

    //moment headeres
    if(sootMom->active()){

        for(int l=0; l<nMoments; l++){
            stringstream int2str;
            int2str << l;
            string moment = "M$"+int2str.str();
            headerData.push_back(moment);
        }
    }

}



////mass matrix evaluation
//
//void CamPremix::massMatrix(doublereal** M){
//
//    /*
//     *inlet boundary: all governing equations are algebraic except energy
//     */
//
//    /*
//     * on all interior mesh points except the mass flow
//     * all other governing equations are ODE. Two additional
//     * imaginary cells are padded to handle the inlet and exhaust
//     */
//    int nCells = reacGeom->getnCells()+2;
//
//    for(int i=0; i< nEqn; i++)
//        M[0][i] = 1.0;
//    /*
//     *at the inlet species conservation eqauations are algebraic
//     */
//    for(int l=0; l<nSpc; l++)
//        M[0][l] = 0.0;
//    /*
//     *mass flow
//     */
//    for(int i=0; i< nCells; i++){
//        int ptr = nVar*i + ptrF;
//        M[0][ptr] = 0.0;
//    }
//    /*
//     *at the exit both species and energy eqn are algebraic
//     */
//    int ptr;
//    for(int l=0; l<nSpc; l++){
//        ptr = nVar*iMesh_e+l;
//        M[0][ptr] = 0.0;
//    }
//    if(admin->getEnergyModel() == admin->ADIABATIC){
//
//        M[0][ptrT] = 0.0;
//
//        ptr = nVar*iMesh_e + ptrT;
//        M[0][ptr] = 0.0;
//    }
//
//}
//

