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

using namespace Camflow;
using namespace Gadgets;

Batch::Batch
(
    CamAdmin& ca,
    CamConfiguration& config,
    CamControl& cc,
    CamGeometry& cg,
    CamProfile& cp,
    CamSoot& cs,
    Mechanism& mech
)
:
  CamSetup(ca, config, cc, cg, cp, cs, mech)
{}

Batch::~Batch()
{}

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

void Batch::checkSetup()
{

	if(mCord != 1) throw std::invalid_argument("Grid size does not equal one\n");
	if(nVar != nSpc + 1) throw std::invalid_argument("nVar != nSpc + 1\n");
	if(nEqn != nVar) throw std::invalid_argument("nEqn != nVar\n");

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
void Batch::solve()
{

    CamBoundary cb;

    // Check that the problem has been setup properly.
    checkSetup();

    /*get the fuel inlet conditions and the fill the
     *solution vector with species mass fractions
     */
    admin_.getLeftBoundary(cb);
    getInletMassFrac(cb,solvect);
    solvect.push_back(cb.getTemperature());

    // Update the mixture properties.
    updateMixture(&solvect[0]);

    reporter_->header("BATCH");
    reporter_->problemDescription(cb,*this);
    reporter_->openFiles();
    reporter_->writeHeader(header());

    /*
     *report the values at the inlet
     */
    report(0.0,&solvect[0]);

    if (solverID == control_.CVODE) {

    	CVodeWrapper cvw;
        cvw.init(nEqn,solvect,control_.getSpeciesAbsTol(), control_.getSpeciesRelTol(),
            control_.getMaxTime(),nEqn,*this);

        cvw.solve(CV_ONE_STEP,control_.getResTol());
        reporter_->closeFiles();

    } else if (solverID == control_.RADAU) {

      RadauWrapper radauWrapper;

        std::vector<doublereal> relTolVector;
        std::vector<doublereal> absTolVector;

        relTolVector.push_back(control_.getSpeciesRelTol());
        absTolVector.push_back(control_.getSpeciesAbsTol());

        radauWrapper.setBandWidth(nEqn);

        radauWrapper.initSolver(nEqn,
                                0.0,
                                control_.getMaxTime(),
                                solvect,
                                relTolVector,
                                absTolVector,
                                *this);

        radauWrapper.Integrate();

        /*
       *write the output to file only if the call is not
       *from the interface
       */
      //if(!interface) {
      //    reportToFile(cc.getMaxTime(),&solvect[0]);
      //}
      // radauWrapper.destroy();

         reporter_->closeFiles();

    } else if (solverID == control_.LIMEX) {
          throw std::logic_error("Error -- Limex is not yet supported");
      }
}

/*
 * function called by the ODE solver
 */
int Batch::eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval)
{
    residual(x,y,ydot);
    return 0;
}


//residual definitions
void Batch::residual(const doublereal& time, doublereal* y, doublereal* f)
{
    updateMixture(y);          // saves the dependent variables
    speciesResidual(time,y,f);      // solve species residual
    energyResidual(time,y,f);       // solve energy residual
    sootResidual(time,y,f);         // solve soot residual
}

void Batch::updateMixture(doublereal* y)
{
    /*
     *update the mixture with the current mass fraction,
     *temperature and density
     */
    doublereal temperature = y[ptrT];
    std::vector<doublereal> massfracs(nSpc,0.0);
    for (int l=0; l< nSpc; ++l)
    {
        massfracs[l] = y[l];
    }

    camMixture_->SetMassFracs(massfracs);
    camMixture_->SetTemperature(temperature);

    // This line MUST come after SetMassFracs/SetTemperature!!!!!
    // Otherwise getAvgMolWt is wrong.
    rho = camMixture_->getAvgMolWt() * opPre / (R * temperature);
    camMixture_->SetMassDensity(rho);

    // Put the molar production rate (mol per m3 per sec) in wdot.
    camMech_->Reactions().GetMolarProdRates(*camMixture_,wdot);

}

//species residual definition
void Batch::speciesResidual(const doublereal& x, doublereal* y, doublereal* f)
{

    for (int l = 0; l < nSpc; ++l)
    {
        f[l] = (1.0/camMixture_->MassDensity()) * (*spv_)[l]->MolWt() * wdot[l];
    }

}

//temperature residual
void Batch::energyResidual(const doublereal& x, doublereal* y, doublereal* f)
{

    int engModel = admin_.getEnergyModel();

    if (engModel == admin_.ISOTHERMAL)
    {
        f[ptrT]= 0.0;
    }
    else
    {
        //get the molar enthalpy
        CamMath cm;
        std::vector<doublereal> eth = camMixture_->getMolarEnthalpy();
        doublereal cp = camMixture_->getSpecificHeatCapacity();
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
        if(engModel == admin_.NONISOTHERMAL)
        {
            doublereal lambda = camMixture_->getThermalConductivity(opPre);
            doublereal eta = camMixture_->getViscosity();
            doublereal dia = reacGeom_.getDia();
            doublereal ht = admin_.getHeatTransferCoeff(x,vel,dia,rho,eta,lambda,cp);
            extSource = ht*reacGeom_.getSurfAres_l()*(admin_.getWallTemp() - y[ptrT]);
        }

        //f[ptrT] = (-heat*Ac + extSource)/(y[ptrF]*cp*Ac);
        f[ptrT] = -heat/(camMixture_->MassDensity()*cp);

    }

}


// soot residual definition
void Batch::sootResidual(const doublereal& x, doublereal* y, doublereal* f)
{
    /*
     *moment residuals
     */
    if (sootMom_.active())
    {
        std::vector<doublereal> mom;
        for (int m=0; m<nMoments; ++m)
            mom.push_back(y[nSpc+m]);

        resMoment.resize(nMoments,0.0);
        sootMom_.residual(x,momRates,&mom[0],&resMoment[0]);

        for (int m=0; m<nMoments; ++m)
        {
            f[nSpc+m] = resMoment[m];
            std::cout << f[nSpc+m] << std::endl;//int dd; cin >> dd;
        }

    }

}

//generate the header data
std::vector<std::string> Batch::header(){


    std::vector<std::string> headerData;

    headerData.clear();
    headerData.push_back("time");
    headerData.push_back("pre");
    headerData.push_back("rho");
    headerData.push_back("T");
    for (int l = 0; l < nSpc; l++) {
        headerData.push_back( (*spv_)[l]->Name() );
    }
    headerData.push_back("sumfracs");
    /*
     *additional output if soot moment are active
     */
    return headerData;

}
//report the solution
void Batch::report(doublereal x, doublereal* soln){
    //std::setw(5); std::setprecision(4);
    static int nStep =0;
    std::cout.width(5);
    std::cout.setf(std::ios::scientific);
    //updateMixture(soln);
    if(nStep%20==0) reporter_->consoleHead("time (s)");
    std::cout << x << std::endl;
    reportToFile(x,soln);
    nStep++;

}

/*
 *this is a dummy function
 */
void Batch::report(doublereal x, doublereal* soln, doublereal& res){
    static int nStep = 0;
    std::cout.width(5);
    std::cout.setf(std::ios::scientific);
    //updateMixture(soln);
    if(nStep%20==0) reporter_->consoleHead("time(s) \t residual");
    std::cout << x << "     " << res << std::endl;
    reportToFile(x,soln);
    nStep++;

}

void Batch::reportToFile(doublereal time, doublereal* soln)
{

    //prepare to report
    doublereal sum =0;
    std::vector<doublereal> data;
    data.clear();
    data.push_back(time);
    data.push_back(opPre);
    data.push_back(rho);
    data.push_back(soln[ptrT]);

    if(admin_.getSpeciesOut()==admin_.MOLE){
        std::vector<doublereal> molefracs;
        molefracs = camMixture_->MoleFractions();
        sum =0.0;
        for (int l = 0; l < nSpc; l++) {
            if(fabs(molefracs[l])> 1e-99)
                data.push_back(fabs(molefracs[l]));
            else
                data.push_back(0.0);
            sum += molefracs[l];
        }

    }else{
        std::vector<doublereal> massfracs;
        camMixture_->GetMassFractions(massfracs);
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

    reporter_->writeStdFileOut(data);

}

doublereal Batch::getResidual()
const
{

    doublereal resNorm=0;

    for (int i=0; i<nSpc*mCord; ++i)
    {
        resNorm += resSp[i]*resSp[i];
    }
    for (int i=0; i<mCord; ++i)
    {
        resNorm += resT[i]*resT[i];
    }

    return std::sqrt(resNorm);

}

/*
 *mass matrix evaluation
 */
void Batch::massMatrix(doublereal** M)
{}
