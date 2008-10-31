#include "fl_reactor.h"
#include "fl_premix.h"
#include <iostream>
using namespace FlameLab;
using namespace Strings;
using namespace Sprog;
using namespace std;

int Reactor::nSpecies;
const Sprog::Mechanism* Reactor::mech;

Reactor::Reactor(){
	fuel = NULL;
	oxidizer = NULL;
}

Reactor::~Reactor(){
	if(fuel != NULL)
		delete fuel;
	if(oxidizer != NULL)
		delete oxidizer;
}

Reactor::Reactor(const Reactor &copy){
	*this = copy;
}
		
// sets the reactor with premix or cdflame model
void Reactor::setReactorModel(ReactorModel modelID){
	rModel = modelID;
}

// returns the reactor model PremixFlame or CDflame
Reactor::ReactorModel Reactor::getReactorModel(){
	return rModel;
}

// returns the fuel initial conditions object
FlameLab::InitialConditions& Reactor::getFuelInletConditions(){
	if(fuel == NULL) 
		fuel = new InitialConditions();
	return *fuel;
}
//returns the oxidizer initial conditions object
FlameLab::InitialConditions& Reactor::getOxidizerInletConditions(){
	if(oxidizer == NULL) 
		oxidizer = new InitialConditions();
	return *oxidizer;
}
// operator overloaded function for checking the reactor model (premix or cdflame)
bool Reactor::operator ==(Reactor::ReactorModel rm){
	return this->rModel == rm;
}
// fill the reactor with the mixture
void Reactor::setMixture(Sprog::Thermo::Mixture &mix){
	this->mix = &mix;
}
// return the mixture in the reactor
Sprog::Thermo::Mixture& Reactor::getMixture() const{
	return *mix;
}

// set the current mechsnism
void Reactor::setMechanism(const Sprog::Mechanism &mech){
	this->mech = &mech;
}
// returns the current mechanism
const Sprog::Mechanism *const Reactor::getMechanism() const{
	return mech;
}


// set the number of species in the mechanism
void Reactor::setSpeciesCount(int n){
	this->nSpecies = n;
}
// get the number of species
int Reactor::getSpeciesCount(){
	return this->nSpecies;
}

// sets the guess values for species mass/mole fraction
void Reactor::setFraction(const std::string &name, FlameLab::real fraction){
	guessSpecies.insert(pair<std::string, FlameLab::real>(name,fraction));
}