/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This class contains information about the initial conditions.
	Inlet conditions at the fuel and oxidizer nozzle are manged in
	this class. The reactor object contains objects of intial conditions
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
#include "fl_initial.h"
#include<string>
using namespace FlameLab;
using namespace Sprog::Thermo;
// sets the velocity in (m/s)
void InitialConditions::setVelocity(real vel){
	this->velocity = vel;
}
// returns the velocity in (m/s)
real InitialConditions::getVelocity() const{
	return this->velocity;
}
//set flow rate
void InitialConditions::setFlowRate(real flr){
	this->flowRate = flr;
}
// return the flow rate
real InitialConditions::getFlowRate() const{
	return this->flowRate;
}
//set radial velocity gradient
void InitialConditions::setRadialVelocityGrad(FlameLab::real vel){
	this->radVelGrad = vel;
}
//return the radial velocity gradient
real& InitialConditions::getRadialVelocityGrad() {
	return this->radVelGrad;
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
void InitialConditions::setMassOrMole(int sp){
	this->mom = sp;
}
// returns the nature of inlet species values (mass/mole) fracs
int InitialConditions::getMassOrMole() const{
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
	for(unsigned int i = 0; i!= fracs.size(); i++)
		massFracs[i] = fracs[i];
}

vector<real>& InitialConditions::getMassFractions(){
	return massFracs;
}

//overloaded operator to check if the inlet inputs are in mass or mole
//bool InitialConditions::operator ==(FlameLab::InitialConditions::MassOrMole sp){
//	return this->Molefraction == sp;
//}
// set the density (Kg/m3)
void InitialConditions::setDensity(real dens){
	this->density = dens;
}
// returns the density in Kg/m3

real InitialConditions::getDensity() const{
	return this->density;
}


//sets the inlet fuel mixture
void InitialConditions::setFuelMixture(Sprog::Thermo::Mixture mix){
	nozzleMixture.insert(pair<string,Sprog::Thermo::Mixture>("fuel",mix));
}
// sets the inlet oxidizer mixture
void InitialConditions::setOxidizerMixture(Sprog::Thermo::Mixture mix){
	nozzleMixture.insert(pair<string,Sprog::Thermo::Mixture>("oxidizer",mix));
}
// returns the initial fuel inlet mixture
Sprog::Thermo::Mixture& InitialConditions::getFuelMixture() {
	ni = nozzleMixture.find("fuel");
	return ni->second;
}
// retuns the initial oxidizer inlet mixture
Sprog::Thermo::Mixture& InitialConditions::getOxidizerMixture()  {
	ni = nozzleMixture.find("oxidizer");
	return ni->second;
}




