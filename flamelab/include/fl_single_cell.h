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

#ifndef FL_SINGLE_CELL
#define FL_SINGLE_CELL
#include "gpc.h"
#include "fl_cell_interface.h"
#include "fl_geometry.h"
#include "fl_initial.h"
#include <vector>

//	faceId	0       1       2       3       4       5       6
//			|-------|-------|-------|-------|-------|-------|
//	cellid	|  0    |   1   |   2   |   3   |  4    |   5   |
//			|-------|-------|-------|-------|-------|-------|
//
namespace FlameLab{
	class SingleCell {
		int cellId;
		real velocity;
		real radVelGrad;
		real massFlux;
		real pressure;
		real density;
		static std::vector<Sprog::Thermo::Mixture> cellMixture;
		//only one face for each cell is stored to save memory
		CellInterface wFace, eFace;
		
	public:
		SingleCell(int n){cellId=n;};

		~SingleCell(){};
		//set the cell id
		void setCellId(int n);
		//set the cell mixture
		void setMixture(Sprog::Thermo::Mixture mix);
		// get the cell mixture
		Sprog::Thermo::Mixture& getMixture();
		// get the cell id
		int getCellId();

		// calculate fluxes use for interior cells
		void evaluateFluxes(real &pre, // pressure
			real &mfW,					// west cell mass flux
			real &mfP,					// present cell mass flux
			vector<real>& dz);			// geometry information

		// calculate fluxes for boundary cells
		void evaluateFluxes(real &pre, // pressure
			real &mfW,					// west cell mass flux
			real &mfP,					// present cell mass flux
			vector<real>& dz,			// geometry information
			InitialConditions &ic);		// nozzle conditions


		//always rturns the east side face fluxes
		const vector<real>& getFaceSpFluxes(int eastFace=0) const;
		const real& getFaceThermalCondFluxes() const;
		const real& getFaceMassFlux() const;
		// set the cell center velocity
		void setVelocity(FlameLab::real vel);
		// return the cell center velocity m/s
		const FlameLab::real& getVelocity() const;
		//set the cell center mass flx in kg/m2s
		void setMassFlux(FlameLab::real mf);
		// return the cell center mass flux in kg/m2s
		const FlameLab::real& getMassFlux() const;
		//return the vector of molar enthalpies all species J/mol
		const vector<real>& getFaceMolarEnthalpy()const;
		//set the pressure
		void setPressure(FlameLab::real pre);
		//return the pressure in Pa
		const FlameLab::real& getPressure() const;
		//set the radial velocity gradient (1/s)
		void setRadialVelocityGrad(real vel);
		//return the radial velocity (1/s)
		const FlameLab::real& getRadVelGrad() const;
		// set the mass density
		void setDensity(FlameLab::real dens);
		// returns the mass density in Kg/m3
		const FlameLab::real& getDensity() const;
		// returns the interface mass density in Kg/m3
		const FlameLab::real& getFaceDensity() const;
		//rturn the west face viscosity in Kg/ms
		const FlameLab::real& getFaceViscosity() const;


	protected:
		SingleCell(){};
	};
}

#endif