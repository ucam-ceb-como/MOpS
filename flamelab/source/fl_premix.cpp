/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This is the premix class. This class manages all functions related
	to premix problem solutions. This class contains implemention of 
	solver call for both steady state and transient problems. The residual
	functions are also defined in this class. 
  Licence:
    This file is part of "flameLab".

    flameLab is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#include "fl_premix.h"
#include "gpc_idealgas.h"
#include <iostream>
#include <cassert>
#include<cvode/cvode_spgmr.h>
#include<cvode/cvode_bandpre.h>
using namespace std;
using namespace FlameLab;
using namespace Sprog;
using namespace Strings;

Premix* Premix::p;

//premix constructor 
Premix::Premix(const Sprog::Mechanism &mech):MFLX(mech.SpeciesCount()),
											TEMP(mech.SpeciesCount()+1),
											DENS(mech.SpeciesCount()+2),
											VEL(mech.SpeciesCount()+3){
	
	
}

//initializes premix variables and updates the mixture with mass fraction temperature and mass density
void Premix::initVariables(Reactor &reac)
{
	
	int index; // species index
	real avgMolarMass = 0.0; // avg molecular mass
	real rho=0.0 ; // mass density
	InitialConditions ic = reac.getFuelInletConditions();
	const Mechanism mech = *reac.getMechanism();
	
	SpeciesPtrVector spv;
	// initialization for the first cell for marching problem
	// resize the vector to hold the species+temp+mass and initialize with zero
	variables.resize(mech.SpeciesCount()+4,0.0); 
	// get the inlet mass/mole fractions from the fuel object
	map<std::string,real> fuelSpecies = ic.getFraction();
	// define iterator for the inlet species map
	map<std::string, real>::iterator mi;
	mi = fuelSpecies.begin();
	//iterate through the map
	while ( mi != fuelSpecies.end()){
		index = mech.FindSpecies(convertToCaps(trim(mi->first)));
		if(index < 0 )
			throw ErrorHandler("Species "+mi->first+" not found in species list\n",200);
		else
			variables[index] = mi->second;		
		mi++;
	}
	 //convert to mass fractions if the inputs are in mole fractions.		
	if(ic.Molefraction == ic.getMassOrMole()){
		unsigned int i;
		spv = mech.Species();		
		//calculate everage molecular mass
		for(i=0; i<mech.SpeciesCount(); i++)
			avgMolarMass += spv[i]->MolWt()*variables[i];
		// convert to mass fractions
		for(i=0;i<mech.SpeciesCount(); i++)
			variables[i] = variables[i]*spv[i]->MolWt()/avgMolarMass;

		//mix.SetMassFracs(variables);
		
		reac.getMixture().SetMassFracs(variables);
		//avgMolarMass = mix.getAvgMolWt(variables);
		avgMolarMass = reac.getMixture().getAvgMolWt(variables);

	}else{ // mass fraction
		//mix.SetMassFracs(variables);
		reac.getMixture().SetMassFracs(variables);
		//avgMolarMass = mix.getAvgMolWt(variables);
		avgMolarMass = reac.getMixture().getAvgMolWt(variables);
	}

	// set the initial mass fractions---------------------
	vector<real> massfracs;
	reac.getMixture().GetMassFractions(massfracs);
	reac.getFuelInletConditions().setFraction(massfracs);
	
	// end of initial mass fraction setting---------------
	
	reac.setSpeciesCount(mech.SpeciesCount());
	if(reac.getEnergyModel() == reac.Isothermal) {
		rho = reac.getPressure()*avgMolarMass/(R*reac.getTemperature());
		reac.getFuelInletConditions().setTemperature(reac.getTemperature());
		variables[TEMP] = reac.getTemperature();
	}
	else if(reac.getEnergyModel() == reac.Adiabatic){
		rho = reac.getPressure()*avgMolarMass/(R*ic.getTemperature());	
		variables[TEMP] = ic.getTemperature();
	}else if(reac.getEnergyModel() == reac.UserDefined){
		vector<real> dz = reac.getGeometry();
		variables[TEMP] = reac.getUserTemperature(0.5*dz[0]);
		//rho = reac.getPressure()*avgMolarMass/(R*variables[TEMP]);	
		real inletTemp = ic.getTemperature();		
		rho = reac.getPressure()*avgMolarMass/(R*inletTemp);	
	}
	//mix.SetMassDensity(rho);	
	reac.getMixture().SetMassDensity(rho);
	//set the initial density
	reac.getFuelInletConditions().setDensity(rho);

	variables[DENS] = rho;
	if(ic.getFlowRate() != 0.0){
		variables[MFLX] = ic.getFlowRate();
		variables[VEL] = variables[MFLX]/rho;
		reac.getFuelInletConditions().setVelocity(variables[VEL]);
	}else{
		variables[MFLX] = rho*ic.getVelocity();
		reac.getFuelInletConditions().setFlowRate(rho*ic.getVelocity());
		variables[VEL] = ic.getVelocity();
	}
	
	
	
	//mix.SetTemperature(variables[TEMP]);
	reac.getMixture().SetTemperature(variables[TEMP]);

	//set the fuel nozzle inlet mixture
	reac.getFuelInletConditions().setFuelMixture(reac.getMixture());
	
	
}
// Initialize transient solution vector
void Premix::initTransVector(FlameLab::Reactor &reac){
	vector<real> massFrac;	
	variables.clear();
	if(reac.getSpaceToTime() == reac.ON){
		for(int i = 0; i<reac.getnCells(); i++){
			cells[i].getMixture().GetMassFractions(massFrac);
			// store the mass fraction of cell "i" into solution vector
			for(int l=0;l<nSpecies; l++){
				variables.push_back(massFrac[l]);
			}
			variables.push_back(cells[i].getMassFlux());
			variables.push_back(cells[i].getMixture().Temperature());

		}					
	}else{

		initVariables(reac);		
		//initialize the cells
		
		vector<real> dz = reac.getGeometry();
		real axPos = 0.0;		
		for(int i=0; i<reac.getnCells(); i++){			
			SingleCell sc(i);
			//set the cell centroid
			if(i==0){
				sc.setCentroid(0.5*dz[i]);
			}else{
				axPos += dz[i-1];
				sc.setCentroid(axPos+0.5*dz[i]);
			}
			//get the user defined temperature profile
			if(reac.getEnergyModel() == reac.UserDefined){
				if( i==0){
					real temp = reac.getUserTemperature(0.5*dz[i]);
					reac.getMixture().SetTemperature(temp);
					
				}else{
					real temp = reac.getUserTemperature(axPos+0.5*dz[i]);
					reac.getMixture().SetTemperature(temp);					
					
				}
			}

			sc.setMixture(reac.getMixture());
			sc.setVelocity(variables[VEL]);
			sc.setPressure(reac.getPressure());
			sc.setMassFlux(variables[MFLX]);
			

			cells.push_back(sc);
		}
		//reset the variables based on cells
		variables.clear();
		for(int i = 0; i<reac.getnCells(); i++){
			cells[i].getMixture().GetMassFractions(massFrac);
			// store the mass fraction of cell "i" into solution vector
			for(int l=0;l<nSpecies; l++){
				variables.push_back(massFrac[l]);
			}
			variables.push_back(cells[i].getMassFlux());
			variables.push_back(cells[i].getMixture().Temperature());
		}


	}

}
			

//call premix to solve the problem
void Premix::solve(Sprog::Mechanism &mech, FlameLab::SolverControl &sc, Reactor &reac, FlameLabIO &io){

	ptrToReactor = &reac;
	// set the reactor run model
	if(sc.getSolMode() ==sc.steadyState){
		if(reac.getReactorModel() == reac.Plug) reac.setReactorRunModel(reac.NDS);
		if(reac.getReactorModel() == reac.PremixFlame ||
			reac.getReactorModel() == reac.PremixFree ){
			if(reac.getDiffusion())
				reac.setReactorRunModel(reac.WDS);
			else
				reac.setReactorRunModel(reac.NDS);
		}
	}else if(sc.getSolMode() ==sc.transient){
		if(reac.getReactorModel() ==reac.Plug) reac.setReactorRunModel(reac.NDT);
		if(reac.getReactorModel() == reac.PremixFlame ||
			reac.getReactorModel() == reac.PremixFree ){
			if(reac.getDiffusion())
				reac.setReactorRunModel(reac.WDT);
			else
				reac.setReactorRunModel(reac.NDT);
		}
	}

	
	int nCells;
	nCells = reac.getnCells();	
	Thermo::Mixture mix(mech.Species());
	reac.setMixture(mix);
	reac.setMechanism(mech);
	reac.setSpaceToTime(reac.OFF);


	// this call is for a marching problem
	if(sc.getMaxTime()<=0)
		tMax = 100000.0;
	else
		tMax = sc.getMaxTime();



	

	/////////////////////////////////////////////////////////
	///////////////    Solver Calls  ////////////////////////
	//////////////// Single cell sweep //////////////////////
	if( reac.getReactorRunModel() == reac.NDS ||
		reac.getReactorRunModel() == reac.WDS){
		
		reac.setSpaceToTime(reac.OFF);

		//loop over the number of sweeps
		for(currentLoop = 0; currentLoop <sc.getNSweep(); currentLoop++){
			cout << "\n";
			cout << " Premix problem statistics\n";
			cout << " Number of species :" << mech.SpeciesCount() << endl;
			cout << " Number of cells :"<< nCells << endl;
			cout << " Sweep number :" << (currentLoop+1) << endl;
			cout << "\n";
			//initialize the solver
			initSolver(sc,reac);
			// prepare for console out
			if(io.getMonitorSwitch() == io.ON) io.prepareConsole(mech,*this);
					
			for(int i=0; i< nCells; i++){
				reac.setAxialPosition(i);
				while (currentTime < tMax ){

					int CVodeError = CVode(cvode_mem,tMax,solVect,&currentTime,CV_ONE_STEP);
					if(CVodeError < 0) {
						string error = CVodeGetReturnFlagName(CVodeError);
						error += "\n";
						throw ErrorHandler(error, CVodeError);
					}
				}
				// save the solution for time marching
				saveSolution(NV_DATA_S(solVect),reac);
				//set the current time to zero
				currentTime = 0.0;
				// output to console
				if(io.getMonitorSwitch()==io.ON) io.writeToConsole(reac);			
				//reinitialize the solver for next cell
				reInitSolver(sc,reac);
				
			}
			//calculate the fluxes 
			calcFluxes();
			//if(reac.getSolveEnergy())
			//	reac.setEnergModel(reac.Adiabatic);
			
		}
		
		if(reac.getReactorRunModel()==reac.WDS) 
			reac.setSpaceToTime(reac.ON);
		io.writeToFile(0,cells,reac);

		if(reac.getSolveEnergy()) 
			reac.setEnergModel(reac.Adiabatic);
		//if(reac.getReactorModel() == reac.PremixFree){
		//	cout << "Enabling energy equation \n";
		//	reac.setEnergModel(reac.Adiabatic);
		//}

	}
	// migrate from space marching to time marching in case of steady state with diffusion
	if (reac.getReactorRunModel() == reac.WDS || // with diffusion steady state
		reac.getReactorRunModel() == reac.WDT || // with diffusion transient
		reac.getReactorRunModel() == reac.NDT )  // no diffusion transient
	
	{

		cout << "Premix problem statistics\n";
		cout << "Number of species :" << mech.SpeciesCount() << endl;
		cout << "Number of cells :"<< nCells << endl;


		initSolver(sc,reac);	

		cout.setf(ios::scientific);
		cout.width(5);
		integrate(sc,reac,io);

		//if(reac.getReactorModel() == reac.PremixFree){ // include energy equation
		//	cout << "Enabling energy equation \n";
		//	reac.setEnergModel(reac.Adiabatic);
		//	reInitSolver(sc,reac);
		//	integrate(sc,reac,io);

		//}


		
	}

	io.writeToFile(currentTime, cells,reac);
	CVodeFree(&cvode_mem);
	if(sc.getMethod() == sc.KRYLOV) CVBandPrecFree(&bp_data);
	
}

void Premix::integrate(SolverControl &sc, Reactor &reac, FlameLabIO &io){

		///*****************************************************

		map<map<real,real>,real> outPut = sc.getOutputInterval();
		map<map<real,real>,real>::iterator outPutIter;
		outPutIter = outPut.begin();

		map<real,real> fromTo;
		map<real,real>::iterator ftIter;
		real nextOut = 0.0;
		//real prevTime = 0.0;

		cout << "Calculating premix flame\n";
		while(outPutIter != outPut.end()){
			
			fromTo = outPutIter->first;
			ftIter = fromTo.end();
			ftIter--;
			real tOut = ftIter->second;
			if(tOut == 0) tOut=tMax;
			real delta_t = outPutIter->second;
			nextOut += delta_t;
			
			int count = 0;
			do{
				if((count%20)==0){
					cout << "\n";
					cout << " Time(s)       Residual" <<endl;
					cout <<"----------     ----------\n";
					
				}
				count++;
				CVodeSetStopTime(cvode_mem,nextOut);
				int CVodeError = CVode(cvode_mem,nextOut,solVect,&currentTime,CV_NORMAL_TSTOP);
				if(CVodeError < 0) {
					string error = CVodeGetReturnFlagName(CVodeError);
					error += "\n";
					throw ErrorHandler(error, CVodeError);
				}else{
					CVodeGetDky(cvode_mem,currentTime,1,derivative);
					
					real resid = getResidual(NV_DATA_S(derivative));
					
					if(io.getMonitorSwitch() == io.ON) 
						cout << " "<<currentTime << "  " << resid << endl;
						//cout << " "<<currentTime << endl;
					if(resid < sc.getResTol() ) break;
					if(sc.getReportMode() == sc.Intermediate)io.writeToFile(currentTime, cells,reac);
					io.writeToFile(currentTime,cells,reac,true);
					nextOut += delta_t;
				}
			}while(currentTime <= tOut);
			
			outPutIter++;
		}			

}

// premix solver initialization. This is called from solve
void Premix::initSolver(SolverControl &sc, Reactor &reac){

	// calculate the number of equations: Always number of species + 2
	// one for mass continuity and the other one for temperature. In case of 
	// isothermal calculation set the right side of temperature solution to zero.
	const Mechanism mech = *reac.getMechanism();
	// nuber of variables
	nVar = mech.SpeciesCount() +2;
	// number of species
	nSpecies = mech.SpeciesCount();
	// current integration time
	currentTime = 0.0;	
	// set the reactor pointer
	//ptrToReactor = &reac;

	
	if ((reac.getReactorRunModel() == reac.WDS ||
		reac.getReactorRunModel() == reac.NDS) &&
		reac.getSpaceToTime() == reac.OFF )
		// this is a marching problem
	{
		if(currentLoop > 0 ){
			// just reinitialize the problem
			variables.clear();
			variables.resize(mech.SpeciesCount()+4,0.0);
			vector<SingleCell>::iterator scp = cells.begin();
			vector<real> massfracs;
			scp->getMixture().GetMassFractions(massfracs);
			for(unsigned int l=0; l<mech.SpeciesCount(); l++)
				variables[l] = massfracs[l];

			variables[MFLX] = scp->getMassFlux();
			variables[VEL] = scp->getVelocity();
			variables[TEMP] = scp->getMixture().Temperature();
			variables[DENS] = scp->getMixture().Density();

			solVect = NULL;
			solVect = N_VMake_Serial(nEq,&variables[0]);
			aTol = sc.getATol();
			CVodeReInit(cvode_mem,&residual,currentTime,solVect,CV_SS,sc.getRTol(),(void*)&aTol);


		}else{
			nEq = mech.SpeciesCount() + 2;
			//initialize the solution vector and update mixture mixture
			initVariables(reac);	
			// create memory for integration
			cvode_mem = NULL;
			cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);	

			/*vriable size is larger than the number of variables in the solution vector
			However since the initialization of the solution vector is done according to the
			Number of equation, only the required variables will be mapped into the solution vector
			*/
			solVect = NULL;
			solVect = N_VMake_Serial(nEq,&variables[0]);
		
			aTol = sc.getATol();	
			//initialize integrator memory
			CVodeMalloc(cvode_mem,&residual,currentTime,solVect,CV_SS,sc.getRTol(),(void*)&aTol);
			
			CVDense(cvode_mem,nEq);
		}


	}else if(reac.getSpaceToTime() == reac.ON){
		// this is a continuation to transient from marching
		if( ! (reac.getReactorRunModel() == reac.WDT ||	reac.getReactorRunModel() == reac.NDT)){
			CVodeFree(&cvode_mem);
			cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);
		}
		nEq = (mech.SpeciesCount()+2)*reac.getnCells();		
		initTransVector(reac);
		//ptrToSlnVector = &variables[0];
		solVect = NULL;
		solVect = N_VMake_Serial(nEq,&variables[0]);

		//NV_DATA_S(solVect);
		aTol = sc.getATol();	
		CVodeMalloc(cvode_mem,&trResidual,currentTime,solVect,CV_SS,sc.getRTol(),(void*)&aTol);
		// Jacobian set up
		int nUp = nVar*2;
		int nLo = nVar*2;

		derivative = NULL;
		derivative = N_VNew_Serial(nEq);

		if(sc.getMethod() == sc.DIRECT){
			CVBand(cvode_mem,nEq,nUp,nLo);
			CVBandSetJacFn(cvode_mem,Jacob,(void*)this);
		}else{
			bp_data = CVBandPrecAlloc(cvode_mem,nEq,nVar,nVar);
			CVBPSpgmr(cvode_mem,PREC_RIGHT,0,bp_data);
		}

	}else if( (reac.getSpaceToTime() == reac.OFF && 
		reac.getReactorRunModel() == reac.WDT) || (reac.getReactorRunModel() == reac.NDT) ){
		// this is a marching problem
		cvode_mem = NULL;
		cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);	
		nEq = (mech.SpeciesCount()+2)*reac.getnCells();		
		initTransVector(reac);
		//ptrToSlnVector = &variables[0];
		solVect = NULL;
		solVect = N_VMake_Serial(nEq,&variables[0]);

		//NV_DATA_S(solVect);
		aTol = sc.getATol();	
		CVodeMalloc(cvode_mem,&trResidual,currentTime,solVect,CV_SS,sc.getRTol(),(void*)&aTol);

		// Jacobian set up
		int nUp = nVar*2;
		int nLo = nVar*2;

		if(sc.getMethod() == sc.DIRECT){
			CVBand(cvode_mem,nEq,nUp,nLo);
			CVBandSetJacFn(cvode_mem,Jacob,(void*)this);
		}else{
			bp_data = CVBandPrecAlloc(cvode_mem,nEq,nVar,nVar);
			CVBPSpgmr(cvode_mem,PREC_RIGHT,0,bp_data);
		}

		derivative = NULL;
		derivative = N_VNew_Serial(nEq);
		

	}

	// set the memory pointer to user defined data
	CVodeSetFdata(cvode_mem,(void*)this);
	//set initial step size
	if(sc.getIniStep() > 0 )
		CVodeSetInitStep(cvode_mem,sc.getIniStep());
	// set the max step size
	if (sc.getMaxStep() != 0)
		CVodeSetMaxStep(cvode_mem,sc.getMaxStep());

	if(sc.getMinStep () != 0)
		CVodeSetMinStep(cvode_mem,sc.getMinStep());

	CVodeSetMaxNumSteps(cvode_mem,5000);
	//CVodeSetMaxNonlinIters(cvode_mem,5);
	
//============================= NEWTON SOLVER IMPLEMENTATION =================//
	//kin_mem = KINCreate();
	//KINSetFdata(kin_mem,(void*)this);
	//KINSetNumMaxIters(kin_mem,5000);
	//templ = NULL;
	//templ = N_VNew_Serial(nEq);
	//KINMalloc(kin_mem,&newtonRes,templ);
	//KINSpgmr(kin_mem,0);


}
// Reinitialization of the solver
void Premix::reInitSolver(SolverControl &sc, Reactor &reac){
	
	aTol = sc.getATol();	
	if(reac.getEnergyModel() == reac.UserDefined){
		vector<real> dz = reac.getGeometry();
		int cellNr = reac.getAxialPosition();
		real axPos = reac.getAxialPosition(cellNr);
		if(cellNr+1 < reac.getnCells())
			axPos += dz[cellNr+1]/2.0;

		real temp = reac.getUserTemperature(axPos);
		variables[TEMP] = temp;
		solVect = NULL;
		solVect = N_VMake_Serial(nEq,&variables[0]);
		NV_DATA_S(solVect);
	}

	if(sc.getSolMode() == sc.transient){
		currentTime = 0.0;
		CVodeReInit(cvode_mem,&trResidual,currentTime,solVect,CV_SS,sc.getRTol(),(void*)&aTol);
	}else{
		CVodeReInit(cvode_mem,&residual,currentTime,solVect,CV_SS,sc.getRTol(),(void*)&aTol);
	}
	
}
// routine for updating the variables in user memory
void Premix::updateVariables(FlameLab::real *y){
	

	real avgMolarMass;
	vector<real> massFracs;
	//p = static_cast<Premix*>(object);



	if( p->ptrToReactor->getSpaceToTime() == p->ptrToReactor->OFF &&
		(p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->WDS ||
		p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->NDS ))
	{

		/*-------------------------------------
		 This is a single cell marching problem
		 --------------------------------------*/		
		massFracs.resize(p->nSpecies);
		//real sumFracs = 0.0;
		for(int i=0; i<p->nSpecies; i++){			
			massFracs[i] = y[i];
			//if(massFracs[i] < 0 ) 
				//cout << "-ve mass fraction found " << massFracs[i] << endl;
		}


		//massFracs[p->nSpecies-1] = 1-sumFracs;
		// store the current solutiuon to user memory
		for(int i=0; i< p->nEq; i++)
			p->variables[i] = y[i];	

		// update the mixture with latest mass fractions
		p->ptrToReactor->getMixture().SetMassFracs(massFracs);
		// set the mixture temperature
		real Tmax  = p->ptrToReactor->getMaxTemperature();
		p->ptrToReactor->getMixture().SetTemperature(min(y[p->TEMP],Tmax));
		// get the avg mol wt
		avgMolarMass = p->ptrToReactor->getMixture().getAvgMolWt(variables);
		// calculate the density		
		variables[DENS] = p->ptrToReactor->getPressure()*avgMolarMass/(R*variables[TEMP]);
		// set the mass density to mixture
		p->ptrToReactor->getMixture().SetMassDensity(variables[DENS]);
		variables[MFLX] = y[p->MFLX];
		// calculate the velocity
		variables[VEL] = variables[MFLX]/variables[DENS];
	}
	if( ((p->ptrToReactor->getSpaceToTime() == p->ptrToReactor->ON) && 
		(p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->WDS) )||
		 p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->WDT  ||
		 p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->NDT   ) {
		
		
		{

			//update the cells
			massFracs.resize(p->nSpecies);
			for(int i=0; i< p->ptrToReactor->getnCells(); i++){
				// extract the mass fraction for each cell
				for(int l=0; l<p->nSpecies; l++){
					int spIndx = (i*p->nVar) + l;
					massFracs[l] = y[spIndx];
				}
				// set the mass fraction to the mixture
				p->cells[i].getMixture().SetMassFracs(massFracs);
				// set the temperature to the mixture
				real Tmax = p->ptrToReactor->getMaxTemperature();
				

				p->cells[i].getMixture().SetTemperature(min(y[(i*p->nVar)+p->TEMP],Tmax));
				// set the mass flux to the cell
				p->cells[i].setMassFlux(y[(i*p->nVar)+p->MFLX]);
				// get the avergae molar mass
				avgMolarMass = p->cells[i].getMixture().getAvgMolWt(massFracs);
				// calculate the mass density
				real rho = p->cells[i].getPressure()*avgMolarMass/(R*y[(i*p->nVar)+p->TEMP]);
				// set the mixture density
				p->cells[i].getMixture().SetMassDensity(rho);
				// calculate the velocity
				real vel = y[(i*p->nVar)+p->MFLX]/rho;
				// set the velocity to the cell
				p->cells[i].setVelocity(vel);
				
			}

			p->calcFluxes();

		}
	}
}

void Premix::calcFluxes(){
	//p = static_cast<Premix*>(object);
	vector<real> dz = p->ptrToReactor->getGeometry();
	real mfW, mfP;
	for(int i=0; i<p->ptrToReactor->getnCells(); i++){
		real pre = p->ptrToReactor->getPressure();
		if(i==0){
			mfW = p->ptrToReactor->getFuelInletConditions().getFlowRate();
			mfP = p->cells[i].getMassFlux();
			p->cells[i].evaluateFluxes(pre,mfW,mfP,dz,p->ptrToReactor->getFuelInletConditions());
		}else{
			mfW = p->cells[i-1].getMassFlux();
			mfP = p->cells[i].getMassFlux();
			p->cells[i].evaluateFluxes(pre,mfW,mfP,dz);
		}
	}
}
// premix residual
int Premix::residual(double time,N_Vector y, N_Vector ydot, void *object){

	p = static_cast<Premix*>(object);		
	p->updateVariables(NV_DATA_S(y));	

	scMassFlux(time, NV_DATA_S(y), NV_DATA_S(ydot));
	scSpeciesResidual(time, NV_DATA_S(y), NV_DATA_S(ydot));
	scEnergyResidual(time, NV_DATA_S(y), NV_DATA_S(ydot));

	return 0;
	
}

// newton residual
int Premix::newtonRes(N_Vector y, N_Vector fval, void *object){
	//p= static_cast<Premix*>(object);
	//p->updateVariables(NV_DATA_S(y));
	//real time = 0.0;
	//mcResidual(time,NV_DATA_S(y), NV_DATA_S(ydot), object);
	//boundary(time,NV_DATA_S(y), NV_DATA_S(ydot), object);
	return 0;
}

int Premix::trResidual(double time,N_Vector y, N_Vector ydot, void *object){

	p = static_cast<Premix*>(object);

	p->updateVariables(NV_DATA_S(y));
	
	mcResidual(time,NV_DATA_S(y), NV_DATA_S(ydot));
	boundary(time,NV_DATA_S(y), NV_DATA_S(ydot));
	return 0;
}

int Premix::jacRes(double time,N_Vector y, N_Vector ydot, void *object, bool firstCall){

	p = static_cast<Premix*>(object);
	p->updateVariables(NV_DATA_S(y));
	mcResidual(time,NV_DATA_S(y), NV_DATA_S(ydot));
	boundary(time,NV_DATA_S(y), NV_DATA_S(ydot));
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// SINGLE CELLS RESIDUAL DEFINTIONS FOLLOW/////////////////////////////////

// single cell mass flux residual definition. Only convective terms needs to be handled for space
// marching problem
void Premix::scMassFlux(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot){
	int axPos;
	real convect;
	vector<real> dz;

	real dens, vel;
	
	axPos = p->ptrToReactor->getAxialPosition();
	dz = p->ptrToReactor->getGeometry();

	InitialConditions ic = p->ptrToReactor->getFuelInletConditions();

	dens = ic.getDensity();
	vel = ic.getVelocity();

	if(axPos == 0){
		convect = (ic.getDensity()*ic.getVelocity() - y[p->MFLX])/(0.5*dz[axPos]);
	}else{
		// previsous cell index for mass flux
		//int pcIndxFlx= p->MFLX + ((axPos-1)*p->nVar);
		real delta = 0.5*(dz[axPos-1]+dz[axPos]);			


		//convect = (p->mcVariables[pcIndxFlx] - y[p->MFLX])/delta;	
		convect = (p->cells[axPos-1].getMassFlux() - y[p->MFLX])/delta;	
			
	
	}

	ydot[p->MFLX] = p->variables[p->VEL] * convect;
	
}

// single cell residual definition for species transport equations. Only convecting term needs to be handled 
// for space marching problem
void Premix::scSpeciesResidual(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot){

	real convect, source, diffusion;
	vector<real> dz;
	static vector<real> wdot, jW,jE;
	int axPos, nSpec;	
	vector<SingleCell>::iterator scp;
	int nCells = p->ptrToReactor->getnCells();



	InitialConditions ic = p->ptrToReactor->getFuelInletConditions();	
	
	
	axPos = p->ptrToReactor->getAxialPosition();
	nSpec = p->ptrToReactor->getSpeciesCount();
	dz = p->ptrToReactor->getGeometry();


	
	// get the molar production rate
	const Mechanism *mech = p->ptrToReactor->getMechanism();
	mech->Reactions().GetMolarProdRates(p->ptrToReactor->getMixture(),wdot);
	const Sprog::SpeciesPtrVector *spv = p->ptrToReactor->getMixture().Species();
		
	vector<real> massFracs;
	if(axPos == 0) {

		//get the diffusion flux terms if sweep is more than one
		if(p->currentLoop > 0){
			scp = p->cells.begin();
			jW = scp->getFaceSpFluxes();
			jE = (scp+1)->getFaceSpFluxes();
		}


		massFracs = ic.getMassFractions();
		real density = p->variables[p->DENS];
		for(int l=0; l<nSpec; l++){
			convect = (ic.getDensity()*ic.getVelocity()*massFracs[l]- y[p->MFLX]*y[l])/
				(dz[axPos]*density);	

			diffusion = 0.0;
			if(p->currentLoop > 0){
				diffusion = (jW[l]-jE[l]) / (dz[axPos]*density);	
			}

			source = (*spv)[l]->MolWt()*wdot[l]/density;			
			ydot[l] = convect + diffusion + source;
		}

	}else{

		//calculate diffusion terms if sweep is more than one
		if(p->currentLoop > 0){
			scp = p->cells.begin() + axPos;
			jW = scp->getFaceSpFluxes();
			jE = (scp+1)->getFaceSpFluxes();
		}

		//int pcIndxSp, pcIndxFlx;
		//real delta = 0.5*(dz[axPos-1]+dz[axPos]);
		

		//Get mass fractions for the previous cell
		p->cells[axPos-1].getMixture().GetMassFractions(massFracs);
		real sumSource = 0.0;
		real density = p->variables[p->DENS];
		for(int l=0; l<nSpec; l++){
		
			convect = (p->cells[axPos-1].getMassFlux()* massFracs[l] - y[p->MFLX]*y[l])/
				(dz[axPos]*density);

			diffusion = 0.0;
			if(p->currentLoop >0){
				if(axPos == nCells-1) 
					diffusion = jW[l]/ (dz[axPos]*density);
				else
					diffusion = (jW[l]-jE[l]) / (dz[axPos]*density);
			}

			source = (*spv)[l]->MolWt()*wdot[l]/density;

			sumSource += (*spv)[l]->MolWt()*wdot[l];

			ydot[l] = convect + diffusion+ source;
		}

	}

	
}
// single cell residual definition for temperature solution
void Premix::scEnergyResidual(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot){
	
	if( p->ptrToReactor->getEnergyModel() == p->ptrToReactor->Isothermal ||
		p->ptrToReactor->getEnergyModel() == p->ptrToReactor->UserDefined) {
		
		ydot[p->TEMP] = 0.0;
	}else if(p->ptrToReactor->getEnergyModel() == p->ptrToReactor->Adiabatic ){
		
		//get the axial position
		int axPos = p->ptrToReactor->getAxialPosition();
		
		//get the number of cells
		int nCells = p->ptrToReactor->getnCells();
		

		real convect, cpP, cpW, source;

		vector<real> dz,hs;		
		static vector<real> wdot;

		// storage for species fluxes inorder to calculate enthalpy flux
		vector<real> jE,jW;
		vector<real> hW,hE;
		// single cell iterator
		vector<SingleCell>::iterator scp;
		// species pointer for calculation of enthaly fluxes
		const SpeciesPtrVector *spv = p->ptrToReactor->getMixture().Species();




		InitialConditions ic = p->ptrToReactor->getFuelInletConditions();

		// get the molar production rate
		const Mechanism *mech = p->ptrToReactor->getMechanism();
		mech->Reactions().GetMolarProdRates(p->ptrToReactor->getMixture(),wdot);
		

		dz = p->ptrToReactor->getGeometry();

		real conduction = 0.0; //conduction term
		real enthalpyFlx = 0.0; // enthalpy flux term

		if(axPos == 0){

			//get the diffusion flux terms if sweep is more than one
			if(p->currentLoop > 0){
				scp = p->cells.begin();
				jW = scp->getFaceSpFluxes();
				jE = (scp+1)->getFaceSpFluxes();
				hW = scp->getFaceMolarEnthalpy();
				hE = (scp+1)->getFaceMolarEnthalpy();
			}
			//define local density for efficiency reason
			real density = p->variables[p->DENS];

			// west cell center specific heat
			cpW = ic.getFuelMixture().getSpecificHeatCapacity(ic.getTemperature());

			// present cell center specific heat
			real temp = min(p->ptrToReactor->getMaxTemperature(),y[p->TEMP]);
			cpP = p->ptrToReactor->getMixture().getSpecificHeatCapacity(temp);
			// get the molar enthalpy of all species
			hs = p->ptrToReactor->getMixture().getMolarEnthalpy(temp);
			

			convect = (ic.getDensity()*ic.getVelocity()*cpW*ic.getTemperature() -
				y[p->MFLX]*cpP*y[p->TEMP])/(dz[axPos]*density*cpP);

			//calculate the heat source terms
			source = 0;
			real sumEnthFlxE = 0.0;
			real sumEnthFlxW = 0.0;			
			for(unsigned int l=0; l!=mech->SpeciesCount(); l++){
				source += wdot[l]*hs[l];

				//calculate the enthalpy flux if sweep is more than ONE
				if(p->currentLoop > 0){
					sumEnthFlxW += jW[l]*hW[l]* (*spv)[l]->MolWt();
					sumEnthFlxE += jE[l]*hE[l]* (*spv)[l]->MolWt();
				}
			
			}

			// calculate the conduction and enthalpy flux terms in case sweep is more than 1
			if(p->currentLoop > 0){
				real qW = scp->getFaceThermalCondFluxes();
				real qE = (scp+1)->getFaceThermalCondFluxes();
				conduction= (qE-qW)/(density*cpP*dz[axPos]);
				enthalpyFlx = (sumEnthFlxW - sumEnthFlxE)/(density*cpP*dz[axPos]);

			}


			source /= (p->variables[p->DENS]*cpP);

		}else{


			//get the diffusion flux terms if sweep is more than one
			if(p->currentLoop > 0){
				scp = p->cells.begin() + axPos;				
				jW = scp->getFaceSpFluxes();
				jE = (scp+1)->getFaceSpFluxes();
				hW = scp->getFaceMolarEnthalpy();
				hE = (scp+1)->getFaceMolarEnthalpy();
			}

			//real delta = 0.5*(dz[axPos-1]+dz[axPos]);

			//define local copy of density
			real density = p->variables[p->DENS];

			// present cell center specific heat
			real temp = min(y[p->TEMP],p->ptrToReactor->getMaxTemperature());
			cpP = p->ptrToReactor->getMixture().getSpecificHeatCapacity(temp);
			// get the molar enthalpy of all species
			hs = p->ptrToReactor->getMixture().getMolarEnthalpy(temp);
			

			// west cell center specific heat
			//int cellId = axPos-1;			
			real temperature = p->cells[axPos-1].getMixture().Temperature();
			cpW = p->cells[axPos-1].getMixture().getSpecificHeatCapacity(temperature);

			
			real flxW = p->cells[axPos-1].getMassFlux();
			convect = (flxW*cpW*temperature - y[p->MFLX]*cpP*y[p->TEMP])/(cpP*density*dz[axPos]);

			source = 0;
			real sumEnthFlxE = 0.0;
			real sumEnthFlxW = 0.0;				
			for(unsigned int l=0; l!=mech->SpeciesCount(); l++){
				source += wdot[l]*hs[l];

				//calculate the enthalpy flux if sweep is more than ONE
				if(p->currentLoop > 0){
					sumEnthFlxW += jW[l]*hW[l]* (*spv)[l]->MolWt();
					if( axPos == nCells-1)
						sumEnthFlxE += 0.0;
					else
						sumEnthFlxE += jE[l]*hE[l]* (*spv)[l]->MolWt();
				}
			
			}

			source /= (p->variables[p->DENS]*cpP);


			if(p->currentLoop > 0){
				real qE =0.0;
				real qW = scp->getFaceThermalCondFluxes();
				if( axPos != nCells-1) qE = (scp+1)->getFaceThermalCondFluxes();
				conduction = (qE-qW)/(density*cpP*dz[axPos]);
				enthalpyFlx = (sumEnthFlxW - sumEnthFlxE)/(density*cpP*dz[axPos]);
			}


		}
		
		p->heatConvection = convect;
		p->heatSource = source;

		ydot[p->TEMP] = convect - source;

	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// MULTI CELL RESIDUAL DEFINITIONS FOLLOW////////////////////////////////////////

void Premix::mcResidual(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot){


	
	const Mechanism *mech = p->ptrToReactor->getMechanism();
	const Sprog::SpeciesPtrVector *spv = p->ptrToReactor->getMixture().Species();

	

	// index variables for species,mass and temperature
	int  mIndex, mIndex_;
	int  spIndex, spIndex_;
	int  TIndex, TIndex_;


	vector<real> jW,jE; // vector for species flux
	vector<real> hW,hE,hP; // vecotor of enthalpies
	vector<real> wdot; // vector for rate of production

	vector<real> dz = p->ptrToReactor->getGeometry();// Geometry information
	int nCells = p->ptrToReactor->getnCells();	// number of finite volume cells

	real diffusion, convection, source, heatSource; // source terms
	real sumEnthFlxW, sumEnthFlxE; // sum of energy flux terms
	real density; // present cell density
	
	vector<SingleCell>::iterator scp;
	scp = p->cells.begin()+1;	

	for(int axPos =1; axPos < nCells-1; axPos++){
		//species fluxes
		jW = scp->getFaceSpFluxes();	// west cell face flux for the present cell
		jE = (scp+1)->getFaceSpFluxes(); // east cell face flux for the present cell
		mech->Reactions().GetMolarProdRates(scp->getMixture(),wdot);// molar production rates
		// get the cell density
		density = scp->getMixture().MassDensity();

		// get the face mass fluxes
		real mfW = scp->getFaceMassFlux();
		real mfE = (scp+1)->getFaceMassFlux();


		mIndex_ = ((axPos-1)*p->nVar)  + p->MFLX;
		mIndex  = (axPos*p->nVar) + p->MFLX;
		
		// ========================Species residual definitions ====================
		for(int l=0; l<p->nSpecies; l++){
			spIndex =  (axPos*p->nVar) +l;
			spIndex_ = ((axPos-1)*p->nVar) + l;

			diffusion = (jW[l]-jE[l])/(density*dz[axPos]);
			
			convection = (mfW*y[spIndex_] - mfE*y[spIndex])/(density*dz[axPos]);
			

			source = (*spv)[l]->MolWt()*wdot[l]/density;

			ydot[spIndex] =  diffusion + convection + source;	

		}//========== End of Species residuals=======================================
		
		//===============      mass flux  residual ==================================

		ydot[mIndex] = scp->getVelocity()*( y[mIndex_] - y[mIndex])/dz[axPos];
		
		//===================  temperature residual ================================
		TIndex  = (axPos*p->nVar) + p->TEMP;			
		TIndex_ = ((axPos-1)*p->nVar) + p->TEMP;
		if (p->ptrToReactor->getEnergyModel() == p->ptrToReactor->Isothermal ||
			p->ptrToReactor->getEnergyModel() == p->ptrToReactor->UserDefined) {
			ydot[TIndex] = 0.0; //isothermal
		}else if(p->ptrToReactor->getEnergyModel() == p->ptrToReactor->Adiabatic) {
			sumEnthFlxW = 0.0;
			sumEnthFlxE = 0.0;
			heatSource = 0.0;
			// enthalpy flux		
			hW = scp->getFaceMolarEnthalpy();	// west cell face enthalpy for the present cell
			hE = (scp+1)->getFaceMolarEnthalpy(); // east cell face enthalpy for the present cell
			// get the molar enthalpy for heat source calculation
			real temp = min(y[TIndex],p->ptrToReactor->getMaxTemperature());
			hP = scp->getMixture().getMolarEnthalpy(temp);

			// cell center specific heat for the west cell
			real cpW =  (scp-1)->getMixture().getSpecificHeatCapacity(min(y[TIndex_],p->ptrToReactor->getMaxTemperature()));
			// cell center specific heat for the present cell
			real cpP =	scp->getMixture().getSpecificHeatCapacity(temp);

			sumEnthFlxW = 0.0;
			sumEnthFlxE = 0.0;
			heatSource = 0.0;

			for(int l=0; l< p->nSpecies; l++){
				sumEnthFlxW += jW[l]*hW[l]* (*spv)[l]->MolWt();
				sumEnthFlxE += jE[l]*hE[l]* (*spv)[l]->MolWt();
				heatSource += wdot[l]*hP[l];
			}

			heatSource /= (density*cpP);

			convection = (mfW*cpW*y[TIndex_] - mfE*cpP*y[TIndex])/
					(density*cpP*dz[axPos]);


			real qW = scp->getFaceThermalCondFluxes(); // west cell face flux for the present cell
			real qE = (scp+1)->getFaceThermalCondFluxes(); // east cell face flux for the present cell
			real conduction = (qE- qW)/(density*cpP*dz[axPos]);
			real enthalpyFlx = (sumEnthFlxW - sumEnthFlxE)/(density*cpP*dz[axPos]);

			ydot[TIndex] = convection - conduction + enthalpyFlx -heatSource; //adiabatic
			
		}

		scp++;

	}

}
	
void Premix::boundary(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot){

	
	
	const Mechanism *mech = p->ptrToReactor->getMechanism();
	const SpeciesPtrVector *spv = p->ptrToReactor->getMixture().Species();

	int nCells = p->ptrToReactor->getnCells();

	int  mIndex, mIndex_;
	int  spIndex, spIndex_;
	int  TIndex, TIndex_;

	vector<real> jW,jE; // vector for species flux
	vector<real> hW,hE,hP; // vecotor of enthalpies
	vector<real> wdot; // vector for rate of production

	vector<real> dz = p->ptrToReactor->getGeometry();// Geometry information
	

	real diffusion, convection, source, heatSource; // source terms
	real sumEnthFlxW, sumEnthFlxE; // sum of energy flux terms
	real density; // present cell density
	vector<SingleCell>::iterator scp;	

	//******************************first cell********************************************//
	{
		scp = p->cells.begin();//initialize the cell iterator

		int axPos = 0;
		density = scp->getMixture().MassDensity();
		vector<real> inletMassFracs;
		InitialConditions ic = p->ptrToReactor->getFuelInletConditions();
		ic.getFuelMixture().GetMassFractions(inletMassFracs);
		jE = (scp+1)->getFaceSpFluxes();
		jW = scp->getFaceSpFluxes();
		mech->Reactions().GetMolarProdRates(scp->getMixture(),wdot); // get the molar production rate

		real mfE = (scp+1)->getFaceMassFlux();

		for(int l=0; l<p->nSpecies; l++){
			spIndex = (axPos*p->nVar) +l;	
			diffusion = (jW[l]-jE[l]) / (density*dz[axPos]);
			convection = (ic.getDensity()*ic.getVelocity()*inletMassFracs[l] - mfE*y[l])/(density*dz[axPos]);			

			source = (*spv)[l]->MolWt()*wdot[l]/density;

			ydot[spIndex] = diffusion + convection + source;										
		}
		mIndex = (axPos*p->nVar) + p->MFLX;		

		ydot[mIndex] = scp->getVelocity()* (ic.getVelocity()*ic.getDensity() - y[p->MFLX])/(dz[axPos]);
					
		TIndex= (axPos*p->nVar) + p->TEMP;

		if (p->ptrToReactor->getEnergyModel() == p->ptrToReactor->Isothermal ||
			p->ptrToReactor->getEnergyModel() == p->ptrToReactor->UserDefined) {

			ydot[TIndex] = 0.0;

		}else if(p->ptrToReactor->getEnergyModel() == p->ptrToReactor->Adiabatic){

			real cpW = ic.getFuelMixture().getSpecificHeatCapacity(ic.getTemperature());
			real cpP = scp->getMixture().getSpecificHeatCapacity(y[TIndex]);
			hW = scp->getFaceMolarEnthalpy();
			hE = (scp+1)->getFaceMolarEnthalpy();
			hP = scp->getMixture().getMolarEnthalpy(y[TIndex]);// molar enthalpy for heat source calc
			sumEnthFlxE = 0.0;
			sumEnthFlxW = 0.0;
			heatSource = 0.0;
			for(int l=0; l<p->nSpecies; l++){
				sumEnthFlxW += jW[l]*hW[l]* (*spv)[l]->MolWt();
				sumEnthFlxE += jE[l]*hE[l]* (*spv)[l]->MolWt();
				heatSource += wdot[l]*hP[l];
			}

			heatSource /= (density*cpP);
			convection = (ic.getVelocity()*ic.getDensity()*cpW*ic.getTemperature()-
					mfE*cpP*y[p->TEMP])/(dz[axPos]*cpP*density);

			real qE = (scp+1)->getFaceThermalCondFluxes();
			real qW = (scp)->getFaceThermalCondFluxes();
			real conduction = (qE-qW)/(scp->getMixture().MassDensity()*cpP*dz[axPos]);

			real enthalpyFlx = (sumEnthFlxW - sumEnthFlxE)/(density*cpP*dz[axPos]);


			ydot[TIndex] = convection - conduction + enthalpyFlx-heatSource;
			
		}

	}
	//**********************************last cell****************************************//
	{
		int axPos = nCells-1;
		scp = p->cells.begin();
		scp+= nCells-1;
		mIndex = (axPos*p->nVar) + p->MFLX;
		mIndex_ = (axPos-1)*p->nVar + p->MFLX;
		jW = scp->getFaceSpFluxes();

		real mfW = scp->getFaceMassFlux();

		mech->Reactions().GetMolarProdRates(scp->getMixture(),wdot); //get the molar production rate
		density = scp->getMixture().MassDensity();
		for(int l=0; l<p->nSpecies; l++){
			spIndex = (axPos*p->nVar) +l;
			spIndex_ = ((axPos-1)*p->nVar) + l;

			diffusion = jW[l]/(density*dz[axPos]);
			convection = (mfW*y[spIndex_] - y[mIndex]*y[spIndex])/(density*dz[axPos]);
			
			source = (*spv)[l]->MolWt()*wdot[l]/density;
			ydot[spIndex] = diffusion + convection + source;										
		}
		
		ydot[mIndex] = scp->getVelocity()*( y[mIndex_]-y[mIndex] )/dz[axPos];
					
		TIndex= (axPos*p->nVar) + p->TEMP;	
		TIndex_ = ((axPos-1)*p->nVar) + p->TEMP;

		if (p->ptrToReactor->getEnergyModel() == p->ptrToReactor->Isothermal ||
			p->ptrToReactor->getEnergyModel() == p->ptrToReactor->UserDefined) {

			ydot[TIndex] = 0.0;

		}else if(p->ptrToReactor->getEnergyModel() == p->ptrToReactor->Adiabatic){

			real cpW = (scp-1)->getMixture().getSpecificHeatCapacity(y[TIndex_]);
			real cpP = scp->getMixture().getSpecificHeatCapacity(y[TIndex]);

			hP = scp->getMixture().getMolarEnthalpy(y[TIndex]);// molar enthalpy for heat source calc
			hW = scp->getFaceMolarEnthalpy();

			sumEnthFlxW = 0.0;
			heatSource = 0.0;
			for(int l=0; l< p->nSpecies; l++){
				sumEnthFlxW += jW[l]*hW[l];
				heatSource += wdot[l]*hP[l];
			}

			heatSource /= (density*cpP);
			convection = (mfW*cpW*y[TIndex_] - y[mIndex]*cpP*y[TIndex])/(density*cpP*dz[axPos]);

			real qW = scp->getFaceThermalCondFluxes();
			real conduction = -qW/(density*cpP*dz[axPos]);

			real enthalpyFlx = sumEnthFlxW/(density*cpP*dz[axPos]);

			ydot[TIndex] = convection - conduction + enthalpyFlx-heatSource;
			

		}

	}
		
}


	
				





////////////////////////////////////////////////////////////////////////////////////////////////////////////

// save the solution to use with time marching problem. Solution vector contains
// species mass fractions, mass flux, and temperature
void Premix::saveSolution(FlameLab::real *y, Reactor &reac){

	//int oneDIndx, densIndx, velIndx, axPos;
	int  axPos;
	axPos = reac.getAxialPosition();
	vector<real> dz = reac.getGeometry();
	real centr = 0.0;
	for(int i=0; i<axPos; i++){
		centr += dz[i];
	}

	//save mixture
	Sprog::Thermo::Mixture newMix = reac.getMixture();
	//sc.setCellId(axPos);
	SingleCell sc(axPos);
	//set the centroid
	if(axPos == 0)
		sc.setCentroid(0.5*dz[axPos]);
	else
		sc.setCentroid(centr + 0.5*dz[axPos]);	
	//cstrMix.push_back(newMix);
	sc.setMixture(newMix);	
	//set the velocity
	sc.setVelocity(variables[VEL]);
	//set the massflux
	sc.setMassFlux(variables[MFLX]);
	sc.setPressure(reac.getPressure());
	cells.push_back(sc);
	
}

// returns the number of equations
int Premix::getNeq(){
	return nEq;
}
// returns true if converged
FlameLab::real Premix::getResidual(FlameLab::real *dy){
	real sum = 0.0;
	for(int i=0; i<nEq; i++)
		sum += fabs(dy[i]);
	return sum;

}


int Premix::Jacob(long int N, long int mu,  long int ml,
				  BandMat J, 
				  FlameLab::real time, 
				  N_Vector y, 
				  N_Vector fy, 
				  void *jac_data, 
				  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{

	p = static_cast<Premix*>(jac_data);	
	//update the local copy to evaluate fluzes
	//p->updateVariables(NV_DATA_S(y));

	N_Vector ftemp, ytemp, ewt;
	/* Rename work vectors for use as temporary values of y and f */
	ftemp = tmp1;
	ytemp = tmp2;
	ewt = tmp3;

	real *ytemp_data, *ftemp_data, *ewt_data, *y_data, *fy_data, *col_j;

	CVodeGetErrWeights(p->cvode_mem,ewt);

	ytemp_data = N_VGetArrayPointer(ytemp);
	ftemp_data = N_VGetArrayPointer(ftemp);
	y_data = N_VGetArrayPointer(y);
	fy_data = N_VGetArrayPointer(fy);
	ewt_data = N_VGetArrayPointer(ewt);

	 /* Load ytemp with y = predicted y vector */
	N_VScale(1.0,y,ytemp);


		//calcuate the machine unit round off

	real uround = UNIT_ROUNDOFF;
	real srur = sqrt(uround);
	real fnorm = N_VWrmsNorm(fy,ewt);
	real *hcur, curh;
	curh = 1.e-05;
	hcur = &curh;
	CVodeGetCurrentStep(p->cvode_mem,hcur);
	real minInc = (fnorm != 0.0) ? ( 1000.0*fabs(*hcur)* uround * fnorm *N) : 1.0;


	  /* Set bandwidth and number of column groups for band differencing */
	long int width = ml + mu + 1;
	long int ngroups = min(width, N);
	int mupper,mlower;
	mupper = mu;
	mlower = ml;
	bool firstCall = true;
	for (int group=1; group <= ngroups; group++) {
		/* Increment all y_j in group */
		for(int j=group-1;j < N; j+=width) {
			real inc = max(srur*fabs(y_data[j]), minInc/ewt_data[j]);
			ytemp_data[j] += inc;
		}
		/* Evaluate residual with incremented y */
		jacRes(time,ytemp, ftemp, jac_data, firstCall);
		firstCall = false;
		
		/* Restore ytemp, then form and load difference quotients */
		for (int j=group-1; j < N; j+=width) {
		  ytemp_data[j] = y_data[j];
		  col_j = BAND_COL(J,j);
		  real inc = max(srur*fabs(y_data[j]), minInc/ewt_data[j]);
		  real inc_inv = 1.0/inc;
		  int i1 = max(0, j-mupper);
		  int i2 = min(j+mlower, (int(N)-1) );
		  for (int i=i1; i <= i2; i++)
				BAND_COL_ELEM(col_j,i,j) =
				inc_inv * (ftemp_data[i] - fy_data[i]);
		}

	}

	return 0;
}



	
