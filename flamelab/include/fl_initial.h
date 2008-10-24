#ifndef FL_INITIAL_H
#define FL_INITIAL_H
#include "fl_params.h"
#include <map>
#include <vector>
namespace FlameLab{
	class InitialConditions{

	public:
		enum MassOrMole{
			Massfraction,
			Molefraction
		};

		InitialConditions(){}
		~InitialConditions(){}


		void setVelocity(real vel);
		real getVelocity() const;

		void setTemperature(real temp);
		real getTemperature() const;

		void setDensity(real dens);
		real getDensity() const;

		void setFraction(const std::string &name, real fraction);
		map<std::string,real> getFraction() const;
		
		void setFraction(vector<real> fracs);
		vector<real> getMassFractions() const;

		void setMassOrMole(MassOrMole sp);
		MassOrMole getMassOrMole();
		bool operator==(MassOrMole sp);

	private:

		real velocity;
		real temperature;
		real density;

		MassOrMole mom;
		std::map<std::string,real> species;
		vector<real> massFracs;
		
	};
}

#endif