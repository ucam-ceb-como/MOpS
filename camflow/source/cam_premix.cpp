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

#include "cam_params.h"
#include "cam_admin.h"
#include "cam_profile.h"
#include "cam_residual.h"
#include "cam_setup.h"
#include "radau_wrapper.h"
#include "cam_premix.h"
#include "cam_reporter.h"
#include "kinsol_wrapper.h"
#include "cvode_wrapper.h"
#include "ida_wrapper.h"
using namespace Camflow;



int CamPremix::eval(doublereal t, doublereal* y, doublereal* ydot, bool jacEval){

    /*
     *this is called by the DAE wrapper object. Given y, ydot is returned
     */   
    residual(t,y,ydot);
    return 0;
}
//mass matrix evaluation

void CamPremix::massMatrix(doublereal** M){

    /*
     *inlet boundary: all governing equations are algebraic except energy
     */

    /*
     * on all interior mesh points except the mass flow
     * all other governing equations are ODE. Two additional
     * imaginary cells are padded to handle the inlet and exhaust
     */
    int nCells = reacGeom->getnCells()+2;

    for(int i=0; i< nEqn; i++)
        M[0][i] = 1.0;
    /*
     *at the inlet species conservation eqauations are algebraic
     */
    for(int l=0; l<nSpc; l++)
        M[0][l] = 0.0;
    /*
     *mass flow
     */
    for(int i=0; i< nCells; i++){
        int ptr = nVar*i + ptrC;
        M[0][ptr] = 0.0;
    }
    /*
     *at the exit both species and energy eqn are algebraic
     */
    int ptr;
    for(int l=0; l<nSpc; l++){
        ptr = nVar*loopEnd+l;
        M[0][ptr] = 0.0;
    }
    if(admin->getEnergyModel() == admin->ADIABATIC){
        
        M[0][ptrT] = 0.0;

        ptr = nVar*loopEnd + ptrT;
        M[0][ptr] = 0.0;
    }

}




void CamPremix::residual(const doublereal& t, doublereal* y, doublereal* f){
    /*
     *residual evaluation function for premix. Internal
     *cell residuals are evaluated by calling the base class
     *function
     */
    saveMixtureProp(y,true);
    updateDiffusionFluxes();

    if(admin->getEnergyModel() == admin->ADIABATIC) updateThermo();
    
    boundaryCondition(t, y,f);
   
    speciesResidual(t,y,f);
   
    energyResidual(t,y,f);
    massFlowResidual(t,y,f);
        
}


void CamPremix::boundaryCondition(const doublereal& t,  doublereal* y, doublereal* f){

    //-------------------------------------------
    //
    //  Left Boundary Settings
    //
    //------------------------------------------

    /*
     *species boundary conditions
     */        
    doublereal convection, diffusion;
    
    for(int l=0; l<nSpc ; l++){
        convection = y[ptrC]*dydx(ud_inlet.Species[l],y[l],dz[0])/m_rho[0];
        diffusion = dydx(0,s_jk(loopBegin,l),dz[0])/m_rho[0];
        f[l] = convection + diffusion;

    }
    /*
     *mass flow
     */
    f[ptrC] = ud_inlet.Vel*(ud_inlet.FlowRate -y[ptrC])/dz[0];
    /*
     *energy: irrespective of the condition, inlet is fixed at the
     *user defined temperature
     */
     // once newton is implmented uncommnet this and comment the one below
     //f[ptrT] = ud_inlet.T -y[ptrT];
     f[ptrT] = 0.0;


     //---------------------------------------------------
     //
     // Right Boundary Settings
     //
     //---------------------------------------------------

    /**
     * Species boundary condition
     */
    
    int ptrSP, ptrSW, ptrCP, ptrCW;
    ptrCP = loopEnd*nVar + ptrC;
    for(int l=0; l<nSpc; l++){
        ptrSP = loopEnd*nVar + l;
        ptrSW = (loopEnd-1)*nVar + l;

        f[ptrSP] = (-y[ptrCP]/m_rho[loopEnd])*dydx(y[ptrSP],y[ptrSW],dz[loopEnd]);
    }

    /*
     *Mass flow
     */
    ptrCW = (loopEnd-1)*nVar + ptrC;
    f[ptrCP] = -(y[ptrCP]/m_rho[loopEnd])*(y[ptrCP]-y[ptrCW])/dz[loopEnd];

    /*
     * Energy
     */
    int ptrTP, ptrTW;
    ptrTP = loopEnd*nVar + ptrT;
    if(admin->getEnergyModel() == admin->ADIABATIC){

        ptrTW = (loopEnd-1)*nVar + ptrT;
        f[ptrTP] = (-y[ptrCP]/m_rho[loopEnd])*(dydx(y[ptrTP],y[ptrTW],dz[loopEnd]));

    }else{

        f[ptrTP] = 0;
    }

    
    
}



void CamPremix::solve(CamControl& cc, CamAdmin& ca, CamGeometry& cg, CamProfile& cp, Mechanism& mech){
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
    reacGeom->descretize();
    reacGeom->addZeroWidthCells();

    opPre = ca.getPressure();
    profile->setGeometryObj(cg);
    reporter = new CamReporter();
    nSpc = camMech->SpeciesCount();
    nVar = nSpc + 2;
    /*
     * 2 additional cells are padded to consider the
     * inlet and the exhaust
     */
    nEqn = nVar *(reacGeom->getnCells());
    ptrC = nSpc;
    ptrT = ptrC + 1;

    /*
     *this loop is for residual evaluation by the base class
     */
    loopBegin = 1;
    loopEnd = reacGeom->getnCells()-1;

    cellBegin = 0;
    cellEnd = reacGeom->getnCells();

    setupSolutionVector(cb,cc);

}

void CamPremix::setupSolutionVector(CamBoundary &cb, CamControl &cc){

    dz = reacGeom->getGeometry();
   
    admin->getLeftBoundary(cb);
    opPre = admin->getPressure();

    storeInlet(cb);
    /*
     *store the inlet species mass fraction, and flow rate for
     *use in the boundary residual definitions
     */

    createSolnVector(cb,cc,solvect);
     
    reporter->header("PREMIX");
    reporter->problemDescription(cb,*this);
    reporter->consoleHead("time (s) ");
    header();

    int solver = cc.getSolver();
    int band = nVar*2;

    if(solver == cc.RADAU){
        RadauWrapper rw;
        rw.setBandWidth(band);
        rw.setControl(cc);
        rw.initSolver(nEqn, 0.0, cc.getMaxTime() ,solvect,rTol,aTol,*this);
        rw.Integrate();
        reportToFile(cc.getMaxTime(),&solvect[0]);
    }
    if(solver == cc.CVODE){
        CVodeWrapper cvw;
        cvw.init(nEqn,solvect,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(), cc.getMaxTime(),band,*this);
        //cvw.setIniStep(1e-06);
        cvw.solve(CV_ONE_STEP);
        reportToFile(cc.getMaxTime(),&solvect[0]);
        cvw.destroy();
    }
    if(solver == cc.IDA){
        IDAWrapper ida;
        ida.init(nEqn,solvect,cc.getResTol(),cc.getMaxTime(),band,*this);
        ida.solve();
        reportToFile(cc.getMaxTime(),&solvect[0]);

    }

    

    
}


void CamPremix::report(doublereal t, doublereal* soln){
    cout.width(5);
    cout.setf(ios::scientific);    
    cout << t << endl;   
	

}

void CamPremix::reportToFile(doublereal t, doublereal* soln){
    saveMixtureProp(soln,false);
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
        data.push_back(m_rho[i]);
        data.push_back(m_u[i]);
        data.push_back(soln[i*nVar+ptrC]);
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
        s_jk(loopEnd+1,l)=0.0;
}

/*
 *update the thermal fluxes by calling
 *the base class function
 */
void CamPremix::updateThermo(){
    CamResidual::updateThermo();

}

void CamPremix::storeInlet(CamBoundary& cb){
    /*
     *stote the inlet properties in the
     *structure to use with the inlet boundary
     */

    vector<doublereal> temp;
    getInletMassFrac(cb,temp);
    camMixture->SetMassFracs(temp);
    doublereal T = getInletTemperature(cb);
    camMixture->SetTemperature(T);
    ud_inlet.T = T;
    ud_inlet.FlowRate = getInletFlowRate(cb);
    ud_inlet.Vel = getInletVelocity(cb);
    ud_inlet.Dens = ud_inlet.FlowRate/ud_inlet.Vel;
    ud_inlet.Species = temp;
    ud_inlet.Dk = camMixture->getMixtureDiffusionCoeff(opPre);

}

void CamPremix::header(){
    headerData.clear();
    headerData.push_back("time");
    headerData.push_back("x");
    headerData.push_back("rho");
    headerData.push_back("u");
    headerData.push_back("mdot");
    headerData.push_back("T");
    for (int l = 0; l < nSpc; l++) {
        headerData.push_back( (*spv)[l]->Name() );
    }
    headerData.push_back("sumfracs");
}
