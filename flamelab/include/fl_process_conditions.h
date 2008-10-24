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
		
		void setTemperature(real temp);
		void setPressure(real pre);
		void setTemptrSolution(EnergyEquation ee);
		
		real getTemperature() const;
		real getPressure() const;
		EnergyEquation getTemptrSolution();
		bool operator != (EnergyEquation ee);

	protected:
		EnergyEquation temptrSolution;
	};


}
	
#endif