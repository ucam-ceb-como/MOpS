/*
 * File:   batch.h
 * Author: vj231
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the batch reactor model
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
 * Created on 04 June 2009, 11:43
 */
#include "batch.h"
#include "cam_math.h"
#include "cvode_wrapper.h"
#include <vector>
#include <cstring>
#include <iostream>
#include <sstream>
using namespace Camflow;

/*
 *set the reactor type const volume or const pressure
 */
void Batch::setType(int n){
    batchType = n;
}
/*
 *return the reactor type
 */
int Batch::getType(){
    return batchType;
}

/*
 *Solve the batch reactor. This is called by CamModel.
 *Entering call to batch reactor model
 */

/*
 * batch solve main call. This function stores the mechanism object and the
 * mixture object to the corresponding pointers, sets the array offsets, the
 * number of equations to be solved, the total number of variables in the
 * system etc.
 */
void Batch::solve(CamControl& cc, CamAdmin& ca, CamGeometry &cg, CamProfile& cp,
                      CamSoot &cs, Mechanism &mech){

    camMech = &mech;
    Thermo::Mixture mix(mech.Species());
    camMixture = &mix;
    spv = camMixture->Species();
    profile = &cp;
    CamBoundary cb;
    admin = &ca;
    sootMom = &cs;

    reporter = new CamReporter();

    nSpc = mech.SpeciesCount();
    nVar = nSpc+1; //species + temperature
    nEqn = nVar;
    ptrT = nVar-1;
    

    /*get the fuel inlet conditions and the fill the
     *solution vector with species mass fractions
     */
    
    ca.getLeftBoundary(cb);
    getInletMassFrac(cb,solvect);
    //set the initial pressure
    opPre = ca.getPressure();

    
    /*
     *initialize the solution vector with species
     *mass fractions
     */
    camMixture->SetMassFracs(solvect);
    //fill the temparature
    doublereal temp = cb.getTemperature();
    camMixture->SetTemperature(temp);

    real avgMolWt = camMixture->getAvgMolWt();
    rho = avgMolWt*opPre/(R*temp);
    camMixture->SetMassDensity(rho);
    camMixture->GetConcs(solvect);
    solvect.push_back(temp);
   


    reporter->header("BATCH ");
    reporter->problemDescription(cb,*this);
    reporter->openFiles();
    header();
    reporter->writeHeader(headerData);
    /*
     *report the values at the inlet
     */
    report(0.0,&solvect[0]);



    CVodeWrapper cvw;
    cvw.init(nEqn,solvect,cc.getSpeciesAbsTol(), cc.getSpeciesRelTol(),
                        cc.getMaxTime(),nEqn,*this);
    
    cvw.solve(CV_ONE_STEP,cc.getResTol());
    reporter->closeFiles();
}

/*
 * function called by the ODE solver
 */
int Batch::eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval){
    /*
     *this is called by the DAE wrapper object. Given y, ydot is returned
     */
    residual(x,y,ydot);
    return 0;
}


//residual definitions
void Batch::residual(const doublereal& time, doublereal* y, doublereal* f){

    updateMixture(time,y);          //saves the dependent variables
    speciesResidual(time,y,f);      // solve species residual
    energyResidual(time,y,f);       // solve energy residual
    /*
     *moment residuals
     */
    if(sootMom->active()){
        vector<doublereal> mom;
        for(int m=0; m<nMoments; m++)
            mom.push_back(y[nSpc+m]);

        resMoment.resize(nMoments,0.0);
        /*
         *evaluate soot residual
         */
        sootMom->residual(time,momRates,&mom[0],&resMoment[0]);
        
        for(int m=0; m<nMoments; m++){
            f[nSpc+m] = resMoment[m];
            //cout << f[nSpc+m] << endl;int dd; cin >> dd;
        }
        
    }
}

void Batch::updateMixture(const doublereal& x, doublereal* y){
    /*
     *update the mixture with the current mass fraction
     *temperature and density
     */
    doublereal tmptr;
    vector<doublereal> molefracs;
    molefracs.resize(nSpc,0.0);
    doublereal cb = 0;
    for (int l = 0; l < nSpc; l++) {
        cb += y[l];
    }
    for(int l=0; l< nSpc; l++){
        molefracs[l] = y[l]/cb;
    }
    //camMixture->SetMassFracs(massfracs);
    camMixture->SetFracs(molefracs);
    if(admin->getEnergyModel() == admin->USERDEFINED)
        tmptr = profile->getUserDefTemp(x);
    else
        tmptr = y[ptrT];
    camMixture->SetTemperature(tmptr);

    opPre = cb*R*tmptr;

}
//species residual definition
void Batch::speciesResidual(const doublereal& x, doublereal* y, doublereal* f){

    camMech->Reactions().GetMolarProdRates(*camMixture,wdot);
    for (int l = 0; l < nSpc; l++) {
        f[l]= wdot[l];
        
    }
    

}

//temperature residual
void Batch::energyResidual(const doublereal& x, doublereal* y, doublereal* f){
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

        f[ptrT] = (-heat*Ac + extSource)/(y[ptrF]*cp*Ac);

    }
}

//generate the header data
void Batch::header(){
    headerData.clear();
    headerData.push_back("time");
    headerData.push_back("pre");
    headerData.push_back("rho");
    headerData.push_back("T");
    for (int l = 0; l < nSpc; l++) {
        headerData.push_back( (*spv)[l]->Name() );
    }
    headerData.push_back("sumfracs");
    /*
     *additional output if soot moment are active
     */

}
//report the solution
void Batch::report(doublereal x, doublereal* soln){
    //std::setw(5); std::setprecision(4);
    static int nStep =0;
    cout.width(5);
    cout.setf(ios::scientific);
    updateMixture(x,soln);
    if(nStep%20==0) reporter->consoleHead("time (s)");
    cout << x << endl;
    reportToFile(x,soln);
    nStep++;

}

/*
 *this is a dummy function
 */
void Batch::report(doublereal x, doublereal* soln, doublereal& res){
    static int nStep = 0;
    cout.width(5);
    cout.setf(ios::scientific);
    updateMixture(x,soln);
    if(nStep%20==0) reporter->consoleHead("time(s) \t residual");
    cout << x << "     " << res << endl;
    reportToFile(x,soln);
    nStep++;

}

void Batch::reportToFile(doublereal time, doublereal* soln){
    //prepare to report
    doublereal sum =0;
    vector<doublereal> data;
    data.clear();
    data.push_back(time);
    data.push_back(opPre);
    data.push_back(rho);
    data.push_back(soln[ptrT]);

    if(admin->getSpeciesOut()==admin->MOLE){
        vector<doublereal> molefracs;
        molefracs = camMixture->MoleFractions();
        sum =0.0;
        for (int l = 0; l < nSpc; l++) {
            if(fabs(molefracs[l])> 1e-99)
                data.push_back(fabs(molefracs[l]));
            else
                data.push_back(0.0);
            sum += molefracs[l];
        }

    }else{
        vector<doublereal> massfracs;
        camMixture->GetMassFractions(massfracs);
        sum =0.0;
        for (int l = 0; l < nSpc; l++) {
            if(fabs(massfracs[l])> 1e-99)
                data.push_back(fabs(massfracs[l]));
            else
                data.push_back(0.0);
            sum += massfracs[l];
        }
    }
    data.push_back(sum);

  

    reporter->writeStdFileOut(data);

}