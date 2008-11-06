/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This class represents each finite volume in the descretized reactor.
	The flame object will contain a static instance of a vector of
	single cell object. Single cells contain the mixture, velocity, and pressure
	and madd flux. All properties are defined at the cell center.
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
#include "fl_single_cell.h"

using namespace Sprog;
using namespace Sprog::Thermo;
using namespace FlameLab;

//	faceId	0       1       2       3       4       5       6
//			|-------|-------|-------|-------|-------|-------|
//	cellid	|  0    |   1   |   2   |   3   |  4    |   5   |
//			|-------|-------|-------|-------|-------|-------|
//

vector<Sprog::Thermo::Mixture> SingleCell::cellMixture;


void SingleCell::setMixture(Sprog::Thermo::Mixture mix){
	cellMixture.push_back(mix);
}

void SingleCell::setCellId(int n){
	this->cellId = n;
}

int SingleCell::getCellId(){
	return this->cellId;
}

Sprog::Thermo::Mixture& SingleCell::getMixture(){
	return cellMixture[this->cellId];
}
// Evaluates the species diffusion fluxes and the thermal conduction fluxes
void SingleCell::evaluateFluxes(FlameLab::real &pre,
								FlameLab::real &mfW,
								FlameLab::real &mfP,
								vector<FlameLab::real> &dz){	
	wFace.calcFluxes(this->cellId,
		pre,
		mfW,
		mfP,
		cellMixture[this->cellId -1],
		cellMixture[this->cellId],
		dz);

}
// returns the species diffusion fluxes in kg/m2s
const vector<FlameLab::real>& SingleCell::getFaceSpFluxes() const{
	return this->wFace.getFaceSpeciesFlx();
}
// returns the thermal conduction fluxes in J/m2s
const FlameLab::real& SingleCell::getFaceThermalCondFluxes() const{
	return this->wFace.getFaceCondFlx();
}
//set the velocity
void SingleCell::setVelocity(FlameLab::real vel){
	this->velocity = vel;
}
//returns the velocity in m/s
const FlameLab::real& SingleCell::getVelocity() const{
	return this->velocity;
}
//set the mass flux
void SingleCell::setMassFlux(FlameLab::real mf){
	this->massFlux = mf;
}
// return the mass flux in Kg/m2s
const FlameLab::real& SingleCell::getMassFlux() const{
	return this->massFlux;
}
	
//return the molar enthapy for all species in J/mol
const vector<FlameLab::real>& SingleCell::getFaceMolarEnthalpy() const{
	return this->wFace.getFaceMolarEnthalpy();
}

// return the face mass flux in Kg/m2s
const FlameLab::real& SingleCell::getFaceMassFlux() const{
	return this->wFace.getFaceMassFlux();
}

//set the pressure
void SingleCell::setPressure(FlameLab::real pre){
	this->pressure =pre;
}

// get the pressure
const FlameLab::real& SingleCell::getPressure() const{
	return this->pressure;
}