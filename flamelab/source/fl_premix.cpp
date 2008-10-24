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
	
	
}

//call premix to solve the problem
void Premix::solve(Sprog::Mechanism &mech, FlameLab::SolverControl &sc, Reactor &reac, FlameLabIO &io){
	
	int nCells;
	nCells = reac.getnCells();	
	Thermo::Mixture mix(mech.Species());
	reac.setMixture(mix);
	reac.setMechanism(mech);
	initSolver(sc,reac);
	// prepare for console out
	io.prepareConsole(mech,*this);
	// this call is for a marching problem
	if(sc.getMaxTime()<=0)
		tMax = 100000.0;
	else
		tMax = sc.getMaxTime();

	/////////////////////////////////////////////////////////
	//////////////  Main call to solver /////////////////////
	//////////////// Single cell sweep //////////////////////
	if(sc == sc.steadyState || reac==reac.Plug){

		mcVariables.resize(nVar*nCells,0.0);

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
			currentTime = 0.0;
			// output to console
			io.writeToConsole(reac);
			reInitSolver(sc);
			
		}
	}else if(sc==sc.transient){
		while (currentTime < tMax ){
			CVode(cvode_mem,tMax,solVect,&currentTime,CV_ONE_STEP);
		}	
	}		
	
}

// premix solver initialization. This is called from solve
void Premix::initSolver(SolverControl &sc, Reactor &reac){

	// calculate the number of equations: Always number of species + 2
	// one for mass continuity and the other one for temperature. In case of 
	// isothermal calculation set the right side of temperature solution to zero.
	const Mechanism mech = *reac.getMechanism();
	nEq = mech.SpeciesCount() + 2;
	// nuber of variables
	nVar = mech.SpeciesCount() +4;
	// number of species
	nSpecies = mech.SpeciesCount();
	// current integration time
	currentTime = 0.0;
	
	// set the reactor pointer
	ptrToReactor = &reac;

	//initialize the solution vector and update mixture mixture
	initVariables(reac);
	
	// create memory for integration
	cvode_mem = NULL;
	cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);	

	ptrToSlnVector = &variables[0];
	//solVect = N_VMake_Serial(nEq,ptrToSlnVector);
	solVect = NULL;
	solVect = N_VMake_Serial(nEq,&variables[0]);

	NV_DATA_S(solVect);

	// nVector for derivatives
	derivative = N_VNew_Serial(nEq);
	
	aTol = sc.getATol();	
	//initialize integrator memory
	CVodeMalloc(cvode_mem,&residual,currentTime,solVect,CV_SS,sc.getRTol(),(void*)&aTol);

	// set the memory pointer to user defined data
	CVodeSetFdata(cvode_mem,(void*)this);

	// dense need to be set only for the initialization call or in case of 
	// neglecting diffusion.
	CVDense(cvode_mem,nEq);
	
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

	
	massFracs.resize(p->nSpecies);
	for(int i=0; i<p->nSpecies; i++)
		massFracs[i] = y[i];
	// store the current solutiuon to user memory
	for(int i=0; i< p->nEq; i++)
		p->variables[i] = y[i];	

	// update the mixture with latest mass fractions
	p->ptrToReactor->getMixture().SetMassFracs(massFracs);
	// get the avg mol wt
	avgMolarMass = p->ptrToReactor->getMixture().getAvgMolWt(variables);
	// calculate the density
	variables[DENS] = p->ptrToReactor->getPressure()*avgMolarMass/(R*variables[TEMP]);
	// calculate the velocity
	variables[VEL] = variables[MFLX]/variables[DENS];
	
	
	
}

// premix residual
int Premix::residual(double time,N_Vector y, N_Vector ydot, void *object){
	//real *data;
	//data = NV_DATA_S(y);
	//SolverManager *sm = static_cast<SolverManager*>(object);	
	//Premix *p = static_cast<Premix*>(object);		
	//vector<real> dz;
	
	//dz = p->ptrToReactor->getGeometry();
	
	
	Premix *p = static_cast<Premix*>(object);		
	

	p->updateVariables(NV_DATA_S(y), object);

	scMassFlux(time, NV_DATA_S(y), NV_DATA_S(ydot), object);
	scSpeciesResidual(time, NV_DATA_S(y), NV_DATA_S(ydot), object);
	scEnergyResidual(time, NV_DATA_S(y), NV_DATA_S(ydot), object);

	
	//p->ptrToReactor->SpeciesResidual(time, NV_DATA_S(y),NV_DATA_S(ydot), object);
	//if(p->ptrToReactor->operator !=(p->ptrToReactor->Isothermal))
	//	p->ptrToReactor->TemperatureResidual();
	//exit(1);

	return 0;
}

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
		int pcIndxFlx= p->MFLX + ((axPos-1)*p->nVar);
		real delta = 0.5*(dz[axPos-1]+dz[axPos]);	

		convect = (p->mcVariables[pcIndxFlx] - y[p->MFLX])/delta;
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
		
	vector<real> massFracs = ic.getMassFractions();
	if(axPos == 0) {
		
		for(int l=0; l<nSpec; l++){
			convect = (ic.getDensity()*ic.getVelocity()*massFracs[l]- y[p->MFLX]*y[l])/
				(0.5*dz[axPos]*p->variables[p->DENS]);	
			source = (*spv)[l]->MolWt()*wdot[l]/p->variables[p->DENS];			
			ydot[l] = convect + source;
		}

	}else{
		int pcIndxSp, pcIndxFlx;
		real delta = 0.5*(dz[axPos-1]+dz[axPos]);
		pcIndxFlx = p->MFLX + ((axPos-1)*p->nVar);

		for(int l=0; l<nSpec; l++){
			pcIndxSp = l + ((axPos-1)*p->nVar);

			convect = (p->mcVariables[pcIndxFlx] * p->mcVariables[pcIndxSp] - y[p->MFLX]*y[l])/
				(delta*p->variables[p->DENS]);
			source = (*spv)[l]->MolWt()*wdot[l]/p->variables[p->DENS];
			ydot[l] = convect + source;
		}

	}

	
}
// single cell residual definition for temperature solution
void Premix::scEnergyResidual(FlameLab::real &time, FlameLab::real *y, FlameLab::real *ydot, void *object){
	Premix *p = static_cast<Premix*>(object);
	if( *p->ptrToReactor != p->ptrToReactor->Adiabatic ||
		*p->ptrToReactor != p->ptrToReactor->UserDefined){
		ydot[p->TEMP] = 0.0;
	}else{
		int axPos = p->ptrToReactor->getAxialPosition();
		real convect, cpP, cpW;
		vector<real> massFracs,dz;		
		InitialConditions ic = p->ptrToReactor->getFuelInletConditions();
		dz = p->ptrToReactor->getGeometry();
		if(axPos == 0){
			// west cell center specific heat
			massFracs = ic.getMassFractions();
			cpW = p->ptrToReactor->getMixture().getSpecificHeatCapacity(massFracs,ic.getTemperature());
			// present cell center specific heat
			p->ptrToReactor->getMixture().GetMassFractions(massFracs);
			cpP = p->ptrToReactor->getMixture().getSpecificHeatCapacity(massFracs,y[p->TEMP]);
			convect = (ic.getDensity()*ic.getVelocity()*cpW*ic.getTemperature() -
				y[p->MFLX]*cpP*y[p->TEMP])/(0.5*dz[axPos]*p->variables[p->DENS]*cpP);
		}else{
			// present cell center specific heat
			p->ptrToReactor->getMixture().GetMassFractions(massFracs);
			cpP = p->ptrToReactor->getMixture().getSpecificHeatCapacity(massFracs,y[p->TEMP]);

			// west cell center specific heat
			int cellId = axPos-1;
			p->getMassFracs(cellId,massFracs);
			cpW = p->ptrToReactor->getMixture().getSpecificHeatCapacity(massFracs,p->getCellTemp(cellId));

			real delta = 0.5*(dz[axPos-1]+dz[axPos]);
			real flxW = p->getFlux(cellId);
			convect = (flxW*cpW*p->getCellTemp(cellId) - y[p->MFLX]*cpP*y[p->TEMP])/(cpP*p->variables[p->DENS]*delta);

		}
		ydot[p->TEMP] = convect;

	}

}


// save the solution to use with time marching problem. Solution vector contains
// species mass fractions, mass flux, and temperature
void Premix::saveSolution(FlameLab::real *y, Reactor &reac){

	int oneDIndx, densIndx, velIndx, axPos;
	axPos = reac.getAxialPosition();
	// store species, flx, and temperature
	for(int i=0; i<nEq; i++){
		oneDIndx = i + (axPos*nVar); // convert 2D index to 1D
		mcVariables[oneDIndx] = y[i];
	}

	// save density
	densIndx = DENS + (axPos*nVar); // convert 2D index to 1D
	mcVariables[densIndx] = variables[DENS];

	// save velocity
	velIndx = VEL + (axPos*nVar);// convert 2D index to 1D
	mcVariables[velIndx] = variables[VEL];
	
}

// returns the mass fractions for the given cellid
void Premix::getMassFracs(int &cellId, std::vector<FlameLab::real> &fracs){
	fracs.resize(nSpecies);
	for(int i=0; i<nSpecies; i++){
		fracs[i]= mcVariables[i+ (cellId*nVar)];
	}
}
// returns the temperature given the cell id
FlameLab::real Premix::getCellTemp(int &cellId){
	return mcVariables[TEMP +(cellId*nVar)];
}
// returns the density for a given cell id
FlameLab::real Premix::getDensity(int &cellId){
	return mcVariables[DENS + (cellId*nVar) ];
}

// returns the velocity for a given cell id
FlameLab::real Premix::getVelocity(int &cellId){
	return mcVariables[VEL + (cellId*nVar) ];
}
// returns the mass flx given the cell id
FlameLab::real Premix::getFlux(int &cellId){
	return mcVariables[MFLX + (cellId*nVar)];
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



	
	
