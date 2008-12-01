/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	Cell Interface class defines the interfaces of finite volume 
	descretized cells. It is a container for face properties 
	such as mass flux, species flux, and thermal flux. Each single 
	cell object will contain a cell interface object. For the purpose
	of memory management the single cell objects contain only the west
	side face interface. To get the information on the east side face
	the next single cells interface object need to be accesssed.

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
#include "fl_cell_interface.h"
#include <cmath>

//	faceId	0       1       2       3       4       5       6
//			|-------|-------|-------|-------|-------|-------|
//	cellid	|  0    |   1   |   2   |   3   |  4    |   5   |
//			|-------|-------|-------|-------|-------|-------|
//
using namespace FlameLab;
using namespace Sprog::Thermo;
void CellInterface::calcFluxes(int cellId, 
							   real &pre, 
							   real &mfW,
							   real &mfP,
							   Sprog::Thermo::Mixture &wMix,
							   Sprog::Thermo::Mixture &pMix,
							   vector<real> &dz){

	vector<real> dCoeffW,dCoeffP;
	vector<real> mF_W, mF_P;
	vector<real> hW, hP;
	real delta;
	dCoeffW = wMix.getMixtureDiffusionCoeff(pre);
	dCoeffP = pMix.getMixtureDiffusionCoeff(pre);
	

	pMix.GetMassFractions(mF_P);
	wMix.GetMassFractions(mF_W);


	// molar enthalpies
	hW = wMix.getMolarEnthalpy(wMix.Temperature());
	hP = pMix.getMolarEnthalpy(pMix.Temperature());


	if(cellId == 0){
		delta = 0.5*dz[cellId];
		intfcRho = wMix.MassDensity();
		mFlux = mfW;
		intfVisc = wMix.getViscosity();
	}else if(cellId == dz.size()){		
		delta = 0.5*dz[cellId-1];
		intfcRho = pMix.MassDensity();
		mFlux = mfP;
		intfVisc = pMix.getViscosity();
	}else{
		delta = 0.5*(dz[cellId]+dz[cellId-1]);
		real rhoW = wMix.MassDensity();
		real rhoP = pMix.MassDensity();
		//arithmatic mean for density
		//real intfcRho = rhoW + (rhoW-rhoP)*(0.5*dz[cellId]/delta);
		//Harmonic mean for density
		intfcRho = rhoP*rhoW/(rhoP- ((rhoP-rhoW)*(0.5*dz[cellId]/delta)));
		// Harmonic average for mass flux
		mFlux = mfW*mfP/(mfP- ((mfP-mfW)*(0.5*dz[cellId]/delta)));
		real viscW = wMix.getViscosity();
		real viscP = pMix.getViscosity();
		intfVisc = viscP*viscW/(viscP - ( (viscP-viscW)*(0.5*dz[cellId]/delta)));

	}

	
	real jCorr = 0.0;
	real jAbsSum = 0.0;
	speciesFlx.resize(dCoeffP.size(),0.0);
	molarEnthalpy.resize(dCoeffP.size(),0.0);
	real Dmix_f;
	for(int l=0; l != dCoeffW.size(); l++){
		if(cellId == 0){
			Dmix_f = dCoeffW[l];
			molarEnthalpy[l] = hW[l];
		}else if(cellId = dz.size()){
			Dmix_f = dCoeffP[l];
			molarEnthalpy[l] = hP[l];
		}else{
		//calculate the inteface diffusion coefficient based on Harmonic average
			real temp = (1.0/dCoeffW[l]) - ( ((1/dCoeffW[l]) - (1/dCoeffP[l]) )*(0.5*dz[cellId]/delta));
			Dmix_f = 1/temp;
		// harmonic mean molar enthalpy
			molarEnthalpy[l] = hP[l]*hW[l]/(hP[l]- ((hP[l]-hW[l])*(0.5*dz[cellId]/delta)));

		}

		real grad = (mF_P[l] - mF_W[l])/delta;
		speciesFlx[l] = -Dmix_f*intfcRho*grad;
		
		// data for flux correction
		jCorr += speciesFlx[l];
		jAbsSum += fabs(speciesFlx[l]);

		
	}
	// Apply flux correction
	int last = speciesFlx.size()-1;
	
	if(fabs(jCorr) <= fabs(speciesFlx[last])*1e-2){
		speciesFlx[last]=speciesFlx[last]-jCorr;
	}else{
		for(int l=0;l!=speciesFlx.size(); l++){
			speciesFlx[l] -= jCorr*fabs(speciesFlx[l])/jAbsSum;
			
		}
	}
	//set the species fluxes to cell interface
	

	// Calcualte the thermal conduction flux
	real kW = wMix.getThermalConductivity(pre);
	real kP = pMix.getThermalConductivity(pre);
	real kF;
	if(cellId == 0)
		kF = kW;
	else if (cellId == dz.size())
		kF = kP;
	else
	//Harmonic average for thermal conductuvity
		kF = kW*kP/(kP- ((kP-kW)*(0.5*dz[cellId]/delta)));

	real TGrad = (pMix.Temperature()-wMix.Temperature())/delta;
	q = kF*TGrad;


	

	
}
// set the face species flux
//void CellInterface::setFaceSpecDiffFlx(int cellId, std::vector<real> spFlx){
//	spFaceFlx[cellId+1] = spFlx;
//}
// returns the face species flux in m2/s
const vector<real>& CellInterface::getFaceSpeciesFlx() const{
	//return spFaceFlx[cellId+1];
	return this->speciesFlx;
}
// set the thermal conduction flux for the given cell face
//void CellInterface::setFaceCondFlx(int cellId, real &condFlx){
//	qFlx[cellId+1] = condFlx;
//}
// returns the interface thermal conduction flux in J/m2 s
const real& CellInterface::getFaceCondFlx() const{
	//return qFlx[cellId+1];
	return this->q;
}

const vector<real>& CellInterface::getFaceMolarEnthalpy() const{
	return this->molarEnthalpy;
}

// returns the mass flux for the west face
const real& CellInterface::getFaceMassFlux() const{
	return this->mFlux;
}

// return the face density
const real& CellInterface::getFaceDensity() const{
	return this->intfcRho;
}
//return the viscosity
const real& CellInterface::getFaceViscosity() const{
	return this->intfVisc;
}