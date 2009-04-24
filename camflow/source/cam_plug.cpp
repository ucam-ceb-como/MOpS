/*
 * File:   cam_plug.cpp
 * Author: vinod (vj231@cam.ac.uk)
 * File purpose:
 *  This class implements the plug flow model
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
 * Created on January 24, 2009, 12:08 PM
 */

#include "cam_profile.h"


#include "cam_admin.h"

#include <vector>
//#include "RadauWrapper.h"
#include "cam_residual.h"
#include "cam_setup.h"
#include "cam_plug.h"
//#include "ida_wrapper.h"
#include "cam_math.h"
#include "cvode_wrapper.h"
using namespace Camflow;
using namespace Sprog;


//return the initial solution vector
void CamPlug::getInitial(vector<doublereal>& initial){
    initial = solvect;
}

//plug solve main call
void CamPlug::solve(CamControl& cc, CamAdmin& ca, CamGeometry &cg, CamProfile& cp,
                                Mechanism &mech){

    camMech = &mech;
    Thermo::Mixture mix(mech.Species());
    camMixture = &mix;
    spv = camMixture->Species();
    profile = &cp;
    CamBoundary cb;
    admin = &ca;

    reporter = new CamReporter();
    //cg.descretize();
    reacGeom = &cg;
    Ac = cg.getArea();
    As = cg.getSurfAres_l();

    

    nSpc = mech.SpeciesCount();
    nVar = nSpc+3; //species + temperature + massflow + residence time
    nEqn = nVar;
    ptrT = nSpc;
    ptrC = ptrT +1;
    ptrR = ptrC +1;
    
    /*get the fuel inlet conditions and the fill the
     *solution vector with species mass fractions
     */
    real avgMolWt=0;
    ca.getLeftBoundary(cb);
    getInletMassFrac(cb,solvect);
    //set the relative tolerance for species    
    rTol.resize(solvect.size(),cc.getSpeciesRelTol());
    //set the absolute tolerance for species
    aTol.resize(solvect.size(),cc.getSpeciesAbsTol());
    
    camMixture->SetMassFracs(solvect);

    //fill the temparature
    doublereal temp = cb.getTemperature();
    camMixture->SetTemperature(temp);
    solvect.push_back(cb.getTemperature());
    rTol.push_back(cc.getTempRelTol());
    aTol.push_back(cc.getTempAbsTol());


    //fill the mass flow rate
    opPre = ca.getPressure();
    avgMolWt = camMixture->getAvgMolWt();
    doublereal rho = opPre*avgMolWt/(R*temp);
    if(cb.getFlowRate() == 0.0){
        cb.setFlowRate(rho*cb.getVelocity());
    }
    solvect.push_back(cb.getFlowRate());
    rTol.push_back(cc.getFlowRelTol());
    aTol.push_back(cc.getFlowAbsTol());

    //residence time
    solvect.push_back(0.0);
    rTol.push_back(cc.getFlowRelTol());
    aTol.push_back(cc.getFlowAbsTol());

    reporter->header("PLUG  ");
    reporter->problemDescription(cb,*this);
    reporter->consoleHead("axial position (m)");
    reporter->openFiles();
    header();
    reporter->writeHeader(headerData);

    int solver = cc.getSolver();

    if(solver == cc.RADAU){

        RadauWrapper rw;
        rw.setBandWidth(nEqn);
        rw.setControl(cc);
        rw.initSolver(nEqn,0.0,cg.getLenth(),solvect,rTol,aTol, *this);
        rw.Integrate();
    }else if(solver == cc.CVODE ){

        CVodeWrapper cvw;
        cvw.init(nEqn,solvect,cc.getSpeciesAbsTol(), cc.getSpeciesRelTol(), cg.getLenth(),nEqn,*this);
        cvw.setIniStep(cc.getIniStep());
        cvw.solve(CV_ONE_STEP);
    }

    reporter->closeFiles();
}


int CamPlug::eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval){
    /*
     *this is called by the DAE wrapper object. Given y, ydot is returned
     */
    residual(x,y,ydot);
    return 0;
}

/*
 *mass matrix evaluation
 */

void CamPlug::massMatrix(doublereal** M){

    for(int i=0; i< nEqn; i++)
        M[0][i] = 1.0;

}


//residual definitions
void CamPlug::residual(const doublereal& x, doublereal* y, doublereal* f){
    /*
     *@y - solution vector
     *@f - dy/dx
     */
    updateMixture(x,y);
    speciesResidual(x,y,f);
    energyResidual(x,y,f);
    massFlowResidual(x,y,f);
    residenceTime(x,y,f);
}

void CamPlug::updateMixture(const doublereal& x, doublereal* y){
    /*
     *update the mixture with the current mass fraction
     *temperature and density
     */
    doublereal tmptr;
    vector<doublereal> massfracs;
    massfracs.resize(nSpc,0.0);
    for (int l = 0; l < nSpc; l++) {
        massfracs[l] = y[l];
    }
    camMixture->SetMassFracs(massfracs);
    if(admin->getEnergyModel() == admin->USERDEFINED)
        tmptr = profile->getUserDefTemp(x);
    else
        tmptr = y[ptrT];
    camMixture->SetTemperature(tmptr);

    /*
     *evaluate density and hence the velocity
     */
    real avgMolWt = camMixture->getAvgMolWt();
    rho = avgMolWt*opPre/(R*tmptr);
    camMixture->SetMassDensity(rho);
    vel = y[ptrC]/rho;

}
//species residual definition
void CamPlug::speciesResidual(const doublereal& x, doublereal* y, doublereal* f){

    camMech->Reactions().GetMolarProdRates(*camMixture,wdot);
    for (int l = 0; l < nSpc; l++) {
        f[l]= wdot[l]*(*spv)[l]->MolWt()/(y[ptrC]);
    }

}

//temperature residual
void CamPlug::energyResidual(const doublereal& x, doublereal* y, doublereal* f){
    int engModel = admin->getEnergyModel();
    if(engModel == admin->ISOTHERMAL || engModel == admin->USERDEFINED){
        f[ptrT]= 0.0;
    }else{
        //get the molar enthalpy
        CamMath cm;
        vector<doublereal> eth = camMixture->getMolarEnthalpy();
        doublereal cp = camMixture->getSpecificHeatCapacity();
        //heat release due to chemical reactions
        doublereal heat = cm.sumVector(eth,wdot);
        doublereal extSource = 0.0;
       
        /*
         * Calculation of heat flux from the reactor wall
         * The heat transfer coefficient is calculated based
         * on Nusselt number correlation based on Graets number
         * This is only an approximation. For exact results the user
         * need to specify the overall heat transfer coefficient
         */
        if(engModel == admin->NONISOTHERMAL){

            doublereal lambda = camMixture->getThermalConductivity(opPre);
            doublereal eta = camMixture->getViscosity();
            doublereal dia = reacGeom->getDia();
            doublereal ht = admin->getHeatTransferCoeff(x,vel,dia,rho,eta,lambda,cp);
            extSource = ht*reacGeom->getSurfAres_l()*(admin->getWallTemp() - y[ptrT]);

        }

        f[ptrT] = (-heat*Ac + extSource)/(y[ptrC]*cp*Ac);

    }
}

//mass flow residual
void CamPlug::massFlowResidual(const doublereal& x, doublereal* y, doublereal* f){
    f[ptrC] = 0;
}

//residence time
void CamPlug::residenceTime(const doublereal& x, doublereal* y, doublereal* f){
    f[ptrR] = 1.0/vel;
}

//generate the header data
void CamPlug::header(){
    headerData.clear();
    headerData.push_back("x");
    headerData.push_back("rho");
    headerData.push_back("u");
    headerData.push_back("mdot");
    headerData.push_back("tau");
    headerData.push_back("T");
    for (int l = 0; l < nSpc; l++) {
        headerData.push_back( (*spv)[l]->Name() );
    }
    headerData.push_back("sumfracs");

}
//report the solution
void CamPlug::report(doublereal x, doublereal* soln){
    //std::setw(5); std::setprecision(4);
    cout.width(5);
    cout.setf(ios::scientific);
    updateMixture(x,soln);
    cout << x << endl;

    //prepare to report
    doublereal sum =0;
    vector<doublereal> data;
    data.clear();
    data.push_back(x);
    data.push_back(rho);
    data.push_back(vel);
    data.push_back(soln[ptrC]);
    data.push_back(soln[ptrR]);
    data.push_back(soln[ptrT]);

    if(admin->getSpeciesOut()==admin->MOLE){
        vector<doublereal> molefracs;
        molefracs = camMixture->MoleFractions();
        sum =0.0;
        for (int l = 0; l < nSpc; l++) {
            data.push_back(fabs(molefracs[l]));
            sum += molefracs[l];
        }

    }else{
        sum =0;
        for (int l = 0; l < nSpc; l++) {
            data.push_back(fabs(soln[l]));
            sum += soln[l];
        }
    }
    data.push_back(sum);
    reporter->writeStdFileOut(data);

}

/*
 *this is a dummy function
 */
void CamPlug::report(doublereal x, doublereal* soln, doublereal& res){
    cout << "nothing implemented\n";
}