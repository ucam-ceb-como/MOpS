#ifndef FL_INITIAL_H
#define FL_INITIAL_H
#include "fl_params.h"
#include "gpc.h"
#include <map>
#include <vector>
namespace FlameLab{
	class InitialConditions{

	public:
		enum MassOrMole{
			Massfraction,
			Molefraction
		};

		enum{
			FUEL=0,
			OXIDIZER
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

		void setFuelMixture(Sprog::Thermo::Mixture mix);
		void setOxidizerMixture(Sprog::Thermo::Mixture mix);

		Sprog::Thermo::Mixture& getFuelMixture() ;
		Sprog::Thermo::Mixture& getOxidizerMixture() ;


	private:

		real velocity;
		real temperature;
		real density;

		MassOrMole mom;
		std::map<std::string,real> species;
		vector<real> massFracs;	
		map<string,Sprog::Thermo::Mixture> nozzleMixture;
		map<string,Sprog::Thermo::Mixture>::iterator ni;
		
	};
}

#endif