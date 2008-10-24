#include "fl_initial.h"

using namespace FlameLab;

// sets the velocity in (m/s)
void InitialConditions::setVelocity(real vel){
	this->velocity = vel;
}
// returns the velocity in (m/s)
real InitialConditions::getVelocity() const{
	return this->velocity;
}
// sets the temperature in (K)
void InitialConditions::setTemperature(real temp){
	this->temperature = temp;
}
// returns the temperature in (K)
real InitialConditions::getTemperature() const{
	return this->temperature;
}
// sets the inlet species conditions such as mass or mole fracs
void InitialConditions::setMassOrMole(FlameLab::InitialConditions::MassOrMole sp){
	this->mom = sp;
}
// returns the nature of inlet species values (mass/mole) fracs
InitialConditions::MassOrMole InitialConditions::getMassOrMole(){
	return this->mom;
}

// sets the species mass/mole fraction into species map for the calling nozzle
void InitialConditions::setFraction(const std::string &name, real fraction){
	species.insert(pair<std::string, real>(name,fraction));
}
// returns the species map for the calling inlet nozzle
map<std::string,real> InitialConditions::getFraction() const{
	return species;
}


void InitialConditions::setFraction(std::vector<real> fracs){
	massFracs.resize(fracs.size());
	for(int i = 0; i!= fracs.size(); i++)
		massFracs[i] = fracs[i];
}

vector<real> InitialConditions::getMassFractions() const{
	return massFracs;
}

//overloaded operator to check if the inlet inputs are in mass or mole
bool InitialConditions::operator ==(FlameLab::InitialConditions::MassOrMole sp){
	return this->Molefraction == sp;
}
// set the density (Kg/m3)
void InitialConditions::setDensity(real dens){
	this->density = dens;
}
// returns the density in Kg/m3

real InitialConditions::getDensity() const{
	return this->density;
}
