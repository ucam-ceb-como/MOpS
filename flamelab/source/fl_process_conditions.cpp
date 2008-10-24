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
void ProcessConditions::setTemptrSolution(ProcessConditions::EnergyEquation ee){
	temptrSolution = ee;
}
// returns the pressure in (Pa)
real ProcessConditions::getPressure() const{
	return pressure;
}
// returns the temperature in (K)
real ProcessConditions::getTemperature() const{
	return temperature;
}
// returns the nature of energy soltuion methods
ProcessConditions::EnergyEquation ProcessConditions::getTemptrSolution(){
	return temptrSolution;
}

bool ProcessConditions::operator !=(FlameLab::ProcessConditions::EnergyEquation ee){
	return this->temptrSolution != ee;
}

