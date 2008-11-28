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

#ifndef FL_CELL_INTERFACE_H
#define FL_CELL_INTERFACE_H
#include "fl_params.h"
#include "gpc_mixture.h"
#include <vector>

//	faceId	0       1       2       3       4       5       6
//			|-------|-------|-------|-------|-------|-------|
//	cellid	|  0    |   1   |   2   |   3   |  4    |   5   |
//			|-------|-------|-------|-------|-------|-------|
//
namespace FlameLab{
	class CellInterface{
	public:
		CellInterface(){ q = 0;};
		~CellInterface(){};
		void calcFluxes(int cellId, // cell number
				real &pre,			// cell pressure
				real &mfW,			// west cell mass flux
				real &mfP,			// present cell mass flux
				Sprog::Thermo::Mixture &lMix, // left cell mixture object
				Sprog::Thermo::Mixture &rMix, // right cell mixture object
				vector<real> &dz);		// descretization 
		// returns the cell inteface flux in m2/s
		const vector<real>& getFaceSpeciesFlx() const;

		// returns the face thermal conduction flux in J/m2-s
		const real& getFaceCondFlx() const;

		//returns the vector of molar enthalpy of all species the face in J/mol
		const vector<real>& getFaceMolarEnthalpy() const;
		// return the mass flux in Kg/m2s
		const real& getFaceMassFlux() const;
		//returns the interface density in Kg/m3
		const real& getFaceDensity() const;
		//returns the face viscosity in Kg/m-s
		const real& getFaceViscosity() const;
	
	protected:		
		real q;
		real mFlux,intfcRho, intfVisc;
		vector<real> speciesFlx;
		vector<real> molarEnthalpy;
	};
}

#endif