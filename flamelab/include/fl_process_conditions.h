/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This class encapsulates the operating conditions. It stores
	information such as whether the reactor is isothermal/adiabatic.
	Operating pressure and temperature in case of isothermal operation
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

#ifndef FL_PROCESS_CONDITIONS_H
#define FL_PROCESS_CONDITIONS_H
#include "fl_params.h"
#include <vector>
namespace FlameLab{
	class ProcessConditions{

		//std::vector<real> species;
		//real vecocity;
		real temperature;
		real pressure;
		

	public:
		ProcessConditions(){}
		~ProcessConditions(){}

		enum EnergyEquation{
			Isothermal,
			Adiabatic,
			UserDefined
		};
		
		//set the process temperature
		void setTemperature(real temp);

		//set the reactor pressure. A constant pressure system is assumed in all case
		void setPressure(real pre);

		//set the energy model(isothermal, userdefined, or adiabatic)
		void setEnergModel(int n);

		//return the energy model
		int getEnergyModel();
		
		//return the process temperature in K
		real getTemperature() const;

		//return the reactor pressure in Pa
		real getPressure() const;

		//set the maximum allowed temperature
		void setMaxTemperature(real T_max);

		//reutrn the maximum allowed temperature. In certain cases if there is 
		//a sudden jump in temerature above certain limit may lead to failure
		//of thermodynamic property calculations
		real getMaxTemperature();


	protected:
		EnergyEquation temptrSolution;
		int energyModel; //adiabatic,isothermal,userdefined
		real Tmax;//maximum allowed temperature
	};


}
	
#endif
