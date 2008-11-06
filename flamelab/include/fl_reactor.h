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

#ifndef FL_REACTOR_H
#define FL_REACTOR_H
#include "fl_params.h"
#include "fl_geometry.h"
#include "fl_process_conditions.h"
#include "fl_reactor_controls.h"
#include "fl_initial.h"
#include "fl_error_handler.h"
#include "gpc.h"
#include "string_functions.h"
#include <iostream>
#include <vector>
#include <map>
namespace FlameLab{
	class Reactor : public Geometry, public ProcessConditions, public ReactorControls{
	
		

	public:		

		enum ReactorModel{
			PremixFlame,
			CDflmae,
			Plug
		};
		

		Reactor();
		Reactor(const Reactor &copy); //copy constructor
		~Reactor();
		//sets the reactor model
		void setReactorModel(int modelID);
		// returns the reactor model
		int getReactorModel();
		// returns the fuel inlet conditions
		FlameLab::InitialConditions& getFuelInletConditions();
		//returns the oxidizer inlet conditions
		FlameLab::InitialConditions& getOxidizerInletConditions();	

		// sets the mixture object to the reactor
		void setMixture(Sprog::Thermo::Mixture &mix);
		// sets the mechanism object to the reactor
		void setMechanism(const Sprog::Mechanism &mech);
		// returns the mixture object
		Sprog::Thermo::Mixture& getMixture() const;
		// returns constant pointer to mechanism object
		const Sprog::Mechanism *const getMechanism() const;
		//sets the number of species
		void setSpeciesCount(int n);
		//returns the number of species
		int getSpeciesCount();
		//set the guess values for species
		void setFraction(const std::string &name, real fraction);
		
	protected:		
			
		int rModel;
		
		
		InitialConditions *fuel, *oxidizer;
		vector<real> variables;
		map<string,real> guessSpecies;
		Sprog::Thermo::Mixture *mix;		
		static const Sprog::Mechanism *mech;
		static int nSpecies;		
				

	};
}

#endif