#ifndef FL_REACTOR_H
#define FL_REACTOR_H
#include "fl_params.h"
#include "fl_geometry.h"
#include "fl_process_conditions.h"
#include "fl_initial.h"
#include "fl_error_handler.h"
#include "gpc.h"
#include "string_functions.h"
#include <iostream>
#include <vector>
#include <map>
namespace FlameLab{
	class Reactor : public Geometry, public ProcessConditions{

	public:

		enum ReactorModel{
			PremixFlame,
			CDflmae,
			Plug
		};

		//Reactor(const Sprog::Mechanism &mech);
		Reactor();
		~Reactor();
		//sets the reactor model
		void setReactorModel(ReactorModel modelID);
		// returns the reactor model
		ReactorModel getReactorModel();
		// returns the fuel inlet conditions
		FlameLab::InitialConditions& getFuelInletConditions();
		//returns the oxidizer inlet conditions
		FlameLab::InitialConditions& getOxidizerInletConditions();	

		bool operator ==(ReactorModel rm);
		// sets the mixture object to the reactor
		void setMixture(Sprog::Thermo::Mixture &mix);
		// sets the mechanism object to the reactor
		void setMechanism(const Sprog::Mechanism &mech);
		// returns the mixture object
		Sprog::Thermo::Mixture& getMixture() const;
		// returns const pointer to mixture object
		Sprog::Thermo::Mixture *const getMixturePtr() const;
		// returns constant pointer to mechanism object
		const Sprog::Mechanism *const getMechanism() const;
		//sets the number of species
		void setSpeciesCount(int n);
		//returns the number of species
		int getSpeciesCount();


	protected:		
		ReactorModel rModel;
		InitialConditions *fuel, *oxidizer;
		vector<real> variables;
		Sprog::Thermo::Mixture *mix;
		const Sprog::Mechanism *mech;
		int nSpecies;


	};
}

#endif