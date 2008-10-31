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
		// returns constant pointer to mechanism object
		const Sprog::Mechanism *const getMechanism() const;
		//sets the number of species
		void setSpeciesCount(int n);
		//returns the number of species
		int getSpeciesCount();
		//set the guess values for species
		void setFraction(const std::string &name, real fraction);
		
	protected:		
			
		ReactorModel rModel;
		
		InitialConditions *fuel, *oxidizer;
		vector<real> variables;
		map<string,real> guessSpecies;
		Sprog::Thermo::Mixture *mix;		
		static const Sprog::Mechanism *mech;
		static int nSpecies;		
				

	};
}

#endif