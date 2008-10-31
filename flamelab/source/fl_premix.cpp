#include "fl_premix.h"
#include "gpc_idealgas.h"
#include <iostream>
using namespace std;
using namespace FlameLab;
using namespace Sprog;
using namespace Strings;


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
	real rho ; // mass density
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
	 //conver to mass fractions if the inputs are in mole fractions.		
	if(ic==ic.getMassOrMole()){ // true if mole fraction
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
	rho = reac.getPressure()*avgMolarMass/(R*ic.getTemperature());	
	//mix.SetMassDensity(rho);	
	reac.getMixture().SetMassDensity(rho);
	//set the initial density
	reac.getFuelInletConditions().setDensity(rho);

	variables[MFLX] = rho*ic.getVelocity();
	variables[TEMP] = ic.getTemperature();
	variables[VEL] = ic.getVelocity();
	variables[DENS] = rho;
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
		for(int i=0; i<reac.getnCells(); i++){
			SingleCell sc(i);
			sc.setMixture(reac.getMixture());
			sc.setVelocity(variables[VEL]);
			sc.setPressure(reac.getPressure());
			sc.setMassFlux(variables[MFLX]);
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
	
	int nCells;
	nCells = reac.getnCells();	
	Thermo::Mixture mix(mech.Species());
	reac.setMixture(mix);
	reac.setMechanism(mech);
	reac.setSpaceToTime(reac.OFF);
	initSolver(sc,reac);
	// prepare for console out
	io.prepareConsole(mech,*this);
	// this call is for a marching problem
	if(sc.getMaxTime()<=0)
		tMax = 100000.0;
	else
		tMax = sc.getMaxTime();

	// set the reactor run model
	if(sc==sc.steadyState){
		if(reac == reac.Plug) reac.setReactorRunModel(reac.NDS);
		if(reac == reac.PremixFlame){
			if(reac.getDiffusion())
				reac.setReactorRunModel(reac.WDS);
			else
				reac.setReactorRunModel(reac.NDS);
		}
	}else if(sc==sc.transient){
		if(reac ==reac.Plug) reac.setReactorRunModel(reac.NDT);
		if(reac == reac.PremixFlame){
			if(reac.getDiffusion())
				reac.setReactorRunModel(reac.WDT);
			else
				reac.setReactorRunModel(reac.NDT);
		}
	}

	/////////////////////////////////////////////////////////
	///////////////    Solver Calls  ////////////////////////
	//////////////// Single cell sweep //////////////////////
	if( reac.getReactorRunModel() == reac.NDS ||
		reac.getReactorRunModel() == reac.WDS){
		
		reac.setSpaceToTime(reac.OFF);
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
			//save the current mixture
			Sprog::Thermo::Mixture newMix = mix;
			currentTime = 0.0;
			// output to console
			io.writeToConsole(reac);
			//reinitialize the solver for next cell
			reInitSolver(sc);
			
		}

		if(reac.getReactorRunModel()==reac.WDS) 
			reac.setSpaceToTime(reac.ON);

	}
	// migrate from space marching to time marching in case of steady state with diffusion
	if (reac.getReactorRunModel() == reac.WDS ||
		reac.getReactorRunModel() == reac.WDT ||
		reac.getReactorRunModel() == reac.NDT )
	{

		initSolver(sc,reac);	
		currentTime = 0.0;
		while (currentTime < tMax ){

			int CVodeError = CVode(cvode_mem,tMax,solVect,&currentTime,CV_ONE_STEP);
			if(CVodeError < 0) {
				string error = CVodeGetReturnFlagName(CVodeError);
				error += "\n";
				throw ErrorHandler(error, CVodeError);
			}
		}	
	}

	CVodeFree(&cvode_mem);
	
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
	ptrToReactor = &reac;

	
	if (reac.getSpaceToTime()== reac.OFF)
	{
		nEq = mech.SpeciesCount() + 2;
		//initialize the solution vector and update mixture mixture
		initVariables(reac);	
		// create memory for integration
		cvode_mem = NULL;
		cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);	
		ptrToSlnVector = &variables[0];
		//solVect = N_VMake_Serial(nEq,ptrToSlnVector);

		/*vriable size is larger than the number of variables in the solution vector
		However since the initialization of the solution vector is done according to the
		Number of equation, only the required variables will be mapped into the solution vector
		*/
		solVect = NULL;
		solVect = N_VMake_Serial(nEq,&variables[0]);

	}else if(reac.getSpaceToTime() == reac.ON){
		if( ! (reac.getReactorRunModel() == reac.WDT ||	reac.getReactorRunModel() == reac.NDT)){
			CVodeFree(&cvode_mem);
			cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);
		}
		nEq = (mech.SpeciesCount()+2)*reac.getnCells();		
		initTransVector(reac);
		ptrToSlnVector = &variables[0];
		solVect = NULL;
		solVect = N_VMake_Serial(nEq,&variables[0]);

	}
			
	NV_DATA_S(solVect);
	// nVector for derivatives
	//derivative = N_VNew_Serial(nEq);	
	aTol = sc.getATol();	
	//initialize integrator memory
	CVodeMalloc(cvode_mem,&residual,currentTime,solVect,CV_SS,sc.getRTol(),(void*)&aTol);

	// set the memory pointer to user defined data
	CVodeSetFdata(cvode_mem,(void*)this);

	if(reac.getSpaceToTime() == reac.OFF)
		// dense need to be set for space marching		
		CVDense(cvode_mem,nEq);
	else{
		int nUp = (nVar*2) + 7;
		int nLo = (nVar*2) + 7;
		CVBand(cvode_mem,nEq,nUp,nLo);
	}


	
	//set the solver controllers
	//set initial step size
	CVodeSetInitStep(cvode_mem,sc.getIniStep());
	// set the max step size
	if (sc.getMaxStep() != 0)
		CVodeSetMaxStep(cvode_mem,sc.getMaxStep());


}
// Reinitialization of the solver
void Premix::reInitSolver(SolverControl &sc){
	aTol = sc.getATol();
	CVodeReInit(cvode_mem,&residual,currentTime,solVect,CV_SS,sc.getRTol(),(void*)&aTol);
}
// routine for updating the variables in user memory
void Premix::updateVariables(FlameLab::real *y, void *object){
	

	real avgMolarMass;
	vector<real> massFracs;
	Premix *p = static_cast<Premix*>(object);

	//if (p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->WDS ||
	//	p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->NDS){

	if( p->ptrToReactor->getSpaceToTime() == p->ptrToReactor->OFF){

		/*-------------------------------------
		 This is a single cell marching problem
		 --------------------------------------*/		
		massFracs.resize(p->nSpecies);
		for(int i=0; i<p->nSpecies; i++)
			massFracs[i] = y[i];
		// store the current solutiuon to user memory
		for(int i=0; i< p->nEq; i++)
			p->variables[i] = y[i];	

		// update the mixture with latest mass fractions
		p->ptrToReactor->getMixture().SetMassFracs(massFracs);
		// set the mixture temperature
		p->ptrToReactor->getMixture().SetTemperature(y[p->TEMP]);
		// get the avg mol wt
		avgMolarMass = p->ptrToReactor->getMixture().getAvgMolWt(variables);
		// calculate the density
		variables[DENS] = p->ptrToReactor->getPressure()*avgMolarMass/(R*variables[TEMP]);
		// set the mass density to mixture
		p->ptrToReactor->getMixture().SetMassDensity(variables[DENS]);
		// calculate the velocity
		variables[VEL] = variables[MFLX]/variables[DENS];
	}
	if( (p->ptrToReactor->getSpaceToTime() == p->ptrToReactor->ON) && 
		(p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->WDS) ||
		 p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->WDT  ||
		 p->ptrToReactor->getReactorRunModel() == p->ptrToReactor->NDT   ) {
		
		
		{

			//update the cells
			massFracs.resize(p->nSpecies);
			for(int i=0; i< p->ptrToReactor->getnCells(); i++){
				for(int l=0; l<p->nSpecies; l++){
					int spIndx = (i*p->nVar) + l;
					massFracs[l] = y[spIndx];
				}
				p->cells[i].getMixture().SetMassFracs(massFracs);
				p->cells[i].getMixture().SetTemperature(y[(i*p->nVar)+p->TEMP]);
				p->cells[i].setMassFlux(y[(i*p->nVar)+p->MFLX]);

				avgMolarMass = p->cells[i].getMixture().getAvgMolWt(massFracs);
				real rho = p->cells[i].getPressure()*avgMolarMass/(R*y[(i*p->nVar)+p->TEMP]);
				p->cells[i].getMixture().SetMassDensity(rho);
				real vel = y[(i*p->nVar)+p->MFLX]/rho;
				p->cells[i].setVelocity(vel);
				
			}

			//loop over cells and evaluate fluxes
			// No inlet diffusion supported
			vector<real> dz = p->ptrToReactor->getGeometry();
			for(int i=1; i<p->ptrToReactor->getnCells(); i++){
				real pre = p->ptrToReactor->getPressure();
				p->cells[i].evaluateFluxes(pre,dz);
			}
		}
	}
}

// premix residual
int Premix::residual(double time,N_Vector y, N_Vector ydot, void *object){

	Premix *p = static_cast<Premix*>(object);		
	

	p->updateVariables(NV_DATA_S(y), object);
	
	switch (p->ptrToReactor->getSpaceToTime()){

		case p->ptrToReactor->OFF :

			scMassFlux(time, NV_DATA_S(y), NV_DATA_S(ydot), object);
			scSpeciesResidual(time, NV_DATA_S(y), NV_DATA_S(ydot), object);
			scEnergyResidual(time, NV_DATA_S(y), NV_DATA_S(ydot), object);
			break;

		default :
			mcResidual(time,NV_DATA_S(y), NV_DATA_S(ydot), object);
			boundary(time,NV_DATA_S(y), NV_DATA_S(ydot), object);
	}	
	return 0;
	
}


////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// SINGLE CELLS RESIDUAL DEFINTIONS FOLLOW/////////////////////////////////

// single cell mass flux residual definition. Only convective terms needs to be handled for space
// marching problem
void Premix::scMassFlux(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot, void *obj){
	int axPos;
	real convect;
	vector<real> dz;

	real dens, vel;
	
	
	Premix *p = static_cast<Premix*>(obj);

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

		real massflx = p->cells[axPos-1].getMassFlux();

		//convect = (p->mcVariables[pcIndxFlx] - y[p->MFLX])/delta;	
		convect = (p->cells[axPos-1].getMassFlux() - y[p->MFLX])/delta;	
			
	
	}

	ydot[p->MFLX] = p->variables[p->VEL] * convect;
	
}

// single cell residual definition for species transport equations. Only convecting term needs to be handled 
// for space marching problem
void Premix::scSpeciesResidual(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot, void *object){

	real convect, source;
	vector<real> dz;
	static vector<real> wdot;
	int axPos, nSpec;	
	Premix *p = static_cast<Premix*>(object);
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
		massFracs = ic.getMassFractions();
		for(int l=0; l<nSpec; l++){
			convect = (ic.getDensity()*ic.getVelocity()*massFracs[l]- y[p->MFLX]*y[l])/
				(0.5*dz[axPos]*p->variables[p->DENS]);	
			source = (*spv)[l]->MolWt()*wdot[l]/p->variables[p->DENS];			
			ydot[l] = convect + source;
		}

	}else{
		//int pcIndxSp, pcIndxFlx;
		real delta = 0.5*(dz[axPos-1]+dz[axPos]);
		//pcIndxFlx = p->MFLX + ((axPos-1)*p->nVar);

		//Get mass fractions for the previous cell
		p->cells[axPos-1].getMixture().GetMassFractions(massFracs);

		for(int l=0; l<nSpec; l++){
			//pcIndxSp = l + ((axPos-1)*p->nVar);
			/*convect = (p->mcVariables[pcIndxFlx] * p->mcVariables[pcIndxSp] - y[p->MFLX]*y[l])/
				(delta*p->variables[p->DENS]);*/

			convect = (p->cells[axPos-1].getMassFlux()* massFracs[l] - y[p->MFLX]*y[l])/
				(delta*p->variables[p->DENS]);


			source = (*spv)[l]->MolWt()*wdot[l]/p->variables[p->DENS];

			ydot[l] = convect + source;
		}

	}

	
}
// single cell residual definition for temperature solution
void Premix::scEnergyResidual(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot, void *object){
	Premix *p = static_cast<Premix*>(object);
	if( !(*p->ptrToReactor != p->ptrToReactor->Isothermal)) {
		
		ydot[p->TEMP] = 0.0;
	}else if(!(*p->ptrToReactor != p->ptrToReactor->Adiabatic)){
		int axPos = p->ptrToReactor->getAxialPosition();
		

		real convect, cpP, cpW, source;

		vector<real> dz,hs;		
		static vector<real> wdot;

		InitialConditions ic = p->ptrToReactor->getFuelInletConditions();

		// get the molar production rate
		const Mechanism *mech = p->ptrToReactor->getMechanism();
		mech->Reactions().GetMolarProdRates(p->ptrToReactor->getMixture(),wdot);
		

		dz = p->ptrToReactor->getGeometry();

		if(axPos == 0){
			// west cell center specific heat
			cpW = ic.getFuelMixture().getSpecificHeatCapacity(ic.getTemperature());

			// present cell center specific heat

			cpP = p->ptrToReactor->getMixture().getSpecificHeatCapacity(y[p->TEMP]);
			// get the molar enthalpy of all species
			hs = p->ptrToReactor->getMixture().getMolarEnthalpy(y[p->TEMP]);
			

			convect = (ic.getDensity()*ic.getVelocity()*cpW*ic.getTemperature() -
				y[p->MFLX]*cpP*y[p->TEMP])/(0.5*dz[axPos]*p->variables[p->DENS]*cpP);
			//calculate the heat source terms
			source = 0;
			for(int l=0; l!=mech->SpeciesCount(); l++){
				source += wdot[l]*hs[l];
			}

			source /= p->variables[p->DENS]*cpP;

		}else{
			real delta = 0.5*(dz[axPos-1]+dz[axPos]);
			// present cell center specific heat
			cpP = p->ptrToReactor->getMixture().getSpecificHeatCapacity(y[p->TEMP]);
			// get the molar enthalpy of all species
			hs = p->ptrToReactor->getMixture().getMolarEnthalpy(y[p->TEMP]);
			

			// west cell center specific heat
			//int cellId = axPos-1;			
			real temperature = p->cells[axPos-1].getMixture().Temperature();
			cpW = p->cells[axPos-1].getMixture().getSpecificHeatCapacity(temperature);

			
			real flxW = p->cells[axPos-1].getMassFlux();
			convect = (flxW*cpW*temperature - y[p->MFLX]*cpP*y[p->TEMP])/(cpP*p->variables[p->DENS]*delta);
			//calculate the heat source terms
			source = 0;
			for(int l=0; l!=mech->SpeciesCount(); l++){
				source += wdot[l]*hs[l];
			}

			source /= p->variables[p->DENS]*cpP;


		}
		ydot[p->TEMP] = convect - source;

	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// MULTI CELL RESIDUAL DEFINITIONS FOLLOW////////////////////////////////////////

void Premix::mcResidual(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot, void *object){

	Premix *p = static_cast<Premix*>(object);
	vector<real> jW,jE;
	vector<real> hW,hE,hP;
	vector<real> wdot;

	const Mechanism *mech = p->ptrToReactor->getMechanism();
	const Sprog::SpeciesPtrVector *spv = p->ptrToReactor->getMixture().Species();
	// get the number of cells
	int nCells = p->ptrToReactor->getnCells();	
	vector<real> dz = p->ptrToReactor->getGeometry();
	
	real diffusion, convection, source, heatSource;
	real sumEnthFlxW, sumEnthFlxE;
	int mIndex_, mIndex;
	int spIndex_, spIndex;
	int TIndex_, TIndex;

	real density; // present cell density
	vector<SingleCell>::iterator scp;
	scp = p->cells.begin();	
	
	// loop over the cells
	for(int axPos =0; axPos <nCells; axPos++){
		// for the present cell only the west face properties are defined.
		// inorder to get the east cell properties, one need to get the 
		// face property of the next cell

		// check for boundary cells
		if(axPos > 0 && axPos <(nCells-1)){
			// mass flux storage indexes
			mIndex_ = ( (axPos-1)*p->nVar ) + p->MFLX;
			mIndex = (axPos*p->nVar) + p->MFLX;
			// get the species diffusion fluxes
			
			jW = scp->getFaceSpFluxes();	// west cell face flux for the present cell
			jE = (scp+1)->getFaceSpFluxes(); // east cell face flux for the present cell

			// species residual definitions
			// get the molar production rates for the present cell
			mech->Reactions().GetMolarProdRates(scp->getMixture(),wdot);


			// enthalpy flux		
			hW = scp->getFaceMolarEnthalpy();	// west cell face enthalpy for the present cell
			hE = (scp+1)->getFaceMolarEnthalpy(); // east cell face enthalpy for the present cell

			// get the molar enthalpy for heat source calculation
			hP = scp->getMixture().getMolarEnthalpy(y[p->TEMP]);
			// get the cell density
			density = scp->getMixture().MassDensity();
			sumEnthFlxW = 0.0;
			sumEnthFlxE = 0.0;
			heatSource = 0.0;
			for(int l=0; l<p->nSpecies; l++){

				diffusion = (jW[l]-jE[l])/(density*dz[axPos]);

				spIndex_ = ((axPos-1)*p->nVar) + l;
				spIndex = (axPos*p->nVar) +l;

				convection = (y[mIndex_]*y[spIndex_] - y[mIndex]*y[mIndex])/
					(density*dz[axPos]);
				source = wdot[l]*(*spv)[l]->MolWt()/density;
				ydot[spIndex] = diffusion + convection + source;

				// enthalpy transport
				sumEnthFlxW += jW[l]*hW[l]* (*spv)[l]->MolWt();
				sumEnthFlxE += jE[l]*hE[l]* (*spv)[l]->MolWt();

				heatSource += wdot[l]*hP[l];
				
				
			}
			// mass flx
			ydot[mIndex] = scp->getVelocity()*( y[mIndex_] - y[mIndex])/dz[axPos];


			if ( !(*p->ptrToReactor != p->ptrToReactor->Isothermal) ) {
				ydot[(axPos*p->nVar) + p->TEMP] = 0.0;
			}else {
				// energy balance
				TIndex = (axPos*p->nVar) + p->TEMP;
				TIndex_ = ((axPos-1)*p->nVar) + p->TEMP;

				// cell center specific heat for the west cell
				real cpW =  (scp-1)->getMixture().getSpecificHeatCapacity(y[TIndex_]);
				// cell center specific heat for the present cell
				real cpP =	scp->getMixture().getSpecificHeatCapacity(y[TIndex]);
				

				convection = (y[mIndex_]*cpW*y[TIndex_] - y[mIndex]*cpP*y[TIndex])/
					(density*cpP*dz[axPos]);

				real qW = scp->getFaceThermalCondFluxes(); // west cell face flux for the present cell
				real qE = (scp+1)->getFaceThermalCondFluxes(); // east cell face flux for the present cell

				real conduction = (qE- qW)/(density*cpP*dz[axPos]);

				real enthalpyFlx = (sumEnthFlxW - sumEnthFlxE)/(density*cpP*dz[axPos]);

				heatSource /= (density*cpP);
				
				ydot[TIndex] = convection + conduction + enthalpyFlx - heatSource;
			}				

		}
		scp++;
	}

}
	
void Premix::boundary(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot, void *object){

	Premix *p = static_cast<Premix*>(object);
	const Mechanism *mech = p->ptrToReactor->getMechanism();
	const SpeciesPtrVector *spv = p->ptrToReactor->getMixture().Species();

	vector<real> jW,jE;
	vector<real> hW,hE,hP;
	vector<real> wdot,dz;
	real diffusion, convection, source, heatSource, density;
	real sumEnthFlxE, sumEnthFlxW;
	vector<SingleCell>::iterator scp;
	scp = p->cells.begin();

	dz = p->ptrToReactor->getGeometry();

	int nCells = p->ptrToReactor->getnCells();	
	for(int axPos = 0; axPos<nCells; axPos+=(nCells-1)){

		if(axPos==0){
			// first cell
			vector<real> inletMassFracs;
			// get the east cell face fluxes
			InitialConditions ic = p->ptrToReactor->getFuelInletConditions();
			ic.getFuelMixture().GetMassFractions(inletMassFracs);
			jE = (scp+1)->getFaceSpFluxes();
			// get face molar enthalpy
			hE = (scp+1)->getFaceMolarEnthalpy();
			// get the molar production rate
			mech->Reactions().GetMolarProdRates(scp->getMixture(),wdot);
			density = scp->getMixture().MassDensity(); // get the present cell density
			hP = scp->getMixture().getMolarEnthalpy(y[p->TEMP]);// molar enthalpy for heat source calc
			sumEnthFlxE = 0.0;
			heatSource = 0.0;
			for(int l=0; l<p->nSpecies; l++){
				diffusion = -jE[l] / (density*dz[axPos]);

				convection = (ic.getDensity()*ic.getVelocity()*inletMassFracs[l]  - 
					y[p->MFLX]*y[l])/(density*dz[axPos]);

				source = wdot[l]*(*spv)[l]->MolWt()/density;

				ydot[l] = diffusion + convection + source;

				sumEnthFlxE += jE[l]*hE[l]* (*spv)[l]->MolWt();
				heatSource += wdot[l]*hP[l];
			}
			// mass flux
			ydot[p->MFLX] = scp->getVelocity()* 
				(ic.getVelocity()*ic.getDensity() - y[p->MFLX])/(dz[axPos]);
			//energy balance
			if( !(*p->ptrToReactor != p->ptrToReactor->Isothermal)) {

				ydot[p->TEMP] = 0.0;
			}else{
				real cpW = ic.getFuelMixture().getSpecificHeatCapacity(ic.getTemperature());
				real cpP = scp->getMixture().getSpecificHeatCapacity(y[p->TEMP]);

				convection = (ic.getVelocity()*ic.getDensity()*cpW*ic.getTemperature()-
					y[p->MFLX]*cpP*y[p->TEMP])/(dz[axPos]*cpP*density);

				real qE = (scp+1)->getFaceThermalCondFluxes();

				real conduction = qE/(scp->getMixture().MassDensity()*cpP*dz[axPos]);

				real enthalpyFlx = - sumEnthFlxE/(density*cpP*dz[axPos]);

				heatSource /= (density*cpP);
				ydot[p->TEMP] = convection + conduction + enthalpyFlx - heatSource;
			}

		}else if(axPos == nCells-1){
			scp += nCells-1;
			//get the molar production rates
			mech->Reactions().GetMolarProdRates(scp->getMixture(),wdot);
			jW = scp->getFaceSpFluxes();
			density = scp->getMixture().MassDensity();
			
			int mIndex_, mIndex;
			int spIndex_, spIndex;
			int TIndex_, TIndex;
			
			mIndex_ = (axPos-1)*p->nVar + p->MFLX;
			mIndex = axPos*p->nVar + p->MFLX;

			hP = scp->getMixture().getMolarEnthalpy(y[p->TEMP]);// molar enthalpy for heat source calc
			hW = scp->getFaceMolarEnthalpy();
			sumEnthFlxW = 0.0;
			heatSource = 0.0;
	
			for(int l=0; l<p->nSpecies; l++){
				spIndex = axPos * p->nVar + l;
				spIndex_ = (axPos-1)* p->nVar + l;

				diffusion = jW[l]/(density*dz[axPos]);
				convection = (y[mIndex_]*y[spIndex_] - y[mIndex]*y[spIndex])/(density*dz[axPos]);
				source = wdot[l]*(*spv)[l]->MolWt()/density;
				ydot[spIndex] = diffusion + convection + source;

				sumEnthFlxW += jW[l]*hW[l];
				heatSource += wdot[l]*hP[l];
			}
				// mass balance
			ydot[mIndex] = scp->getVelocity()*( y[mIndex_]-y[mIndex] )/dz[axPos];
			//energy balance
			if( !(*p->ptrToReactor != p->ptrToReactor->Isothermal)) {
				ydot[(axPos-1)*p->nVar + p->TEMP] = 0.0;
			}else{
				TIndex_ = (axPos-1)*p->nVar + p->TEMP;
				TIndex = axPos*p->nVar + p->TEMP;

				real cpW = (scp-1)->getMixture().getSpecificHeatCapacity(y[TIndex_]);
				real cpP = scp->getMixture().getSpecificHeatCapacity(y[TIndex]);

				convection = (y[mIndex_]*cpW*y[TIndex_] - y[mIndex]*cpP*y[TIndex])/(density*cpP*dz[axPos]);

				real qW = scp->getFaceThermalCondFluxes();

				real conduction = -qW/(density*cpP*dz[axPos]);

				real enthalpyFlx = sumEnthFlxW/(density*cpP*dz[axPos]);

				heatSource /= (density*cpP);

				ydot[TIndex] = convection + conduction + enthalpyFlx + heatSource;
			}
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
	// store species, flx, and temperature
	//for(int i=0; i<nEq; i++){
	//	oneDIndx = i + (axPos*nVar); // convert 2D index to 1D
	//	mcVariables[oneDIndx] = y[i];
	//}

	// save density
	//densIndx = DENS + (axPos*nVar); // convert 2D index to 1D
	//mcVariables[densIndx] = variables[DENS];

	// save velocity
	//velIndx = VEL + (axPos*nVar);// convert 2D index to 1D
	//mcVariables[velIndx] = variables[VEL];
	
	//save mixture
	Sprog::Thermo::Mixture newMix = reac.getMixture();
	//cstrMix.push_back(newMix);
	SingleCell sc(axPos);
	//sc.setCellId(axPos);
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
bool Premix::checkConvergance(FlameLab::real *dy){
	real sum = 0.0;
	for(int i=0; i<nEq; i++)
		sum += fabs(dy[i]);

	if(sum < 1.0e-04)
		return 1;
	else
		return 0;
}

	
	
