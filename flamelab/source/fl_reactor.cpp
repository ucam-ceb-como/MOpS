/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This is the reactor class. It contains information
	regarding the reactor under simulation. It holds information
	on the inlet condition, geometry of the reactor, and
	operating conditions defined for the reactor.
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
void Reactor::setReactorModel(int modelID){
	rModel = modelID;
}

// returns the reactor model PremixFlame or CDflame
int Reactor::getReactorModel(){
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
//bool Reactor::operator ==(Reactor::ReactorModel rm){
//	return this->rModel == rm;
//}

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