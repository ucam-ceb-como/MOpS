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

#include "fl_process_conditions.h"

using namespace FlameLab;
// set the operating pressure
void ProcessConditions::setPressure(real pre){
	pressure = pre;
}

// sets the operating temperature
void ProcessConditions::setTemperature(real temp){
	temperature = temp;
}
// sets the temperature equation solution method
//void ProcessConditions::setTemptrSolution(ProcessConditions::EnergyEquation ee){
//	temptrSolution = ee;
//}
//set the energy model
void ProcessConditions::setEnergModel(int n){
	energyModel = n;
}
//return the energy model
int ProcessConditions::getEnergyModel(){
	return energyModel;
}
// returns the pressure in (Pa)
real ProcessConditions::getPressure() const{
	return pressure;
}
// returns the temperature in (K)
real ProcessConditions::getTemperature() const{
	return temperature;
}
//set the maximum allowed temperature
void ProcessConditions::setMaxTemperature(FlameLab::real T_max){
	Tmax = T_max;
}
//return the maximum allowed temperature
real ProcessConditions::getMaxTemperature(){
	return Tmax;
}
