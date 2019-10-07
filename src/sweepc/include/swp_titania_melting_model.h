/*
  Author(s):      Casper Lindberg (csl37)
  Project:        sweep (population balance solver)
  
  Copyright (C) 2019 Casper Lindberg.

  File purpose:
    Implement size dependent phase transformation from non-crystalline (liquid) titania 
	to solid crystal phases. Three models are implemented:
	1. "melting" point dependent transformation, T = [Tbulk]*(1-[A]/d)
	2. "Gibbs" model with fixed crossover size and transformation temperature
	3. O/Ti "composition" dependent crossover size and fixed transformation temperature

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Prof. Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    Philippa Fawcett Drive
    Cambridge
    CB3 0AS
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef SWEEP_MELTING_H
#define SWEEP_MELTING_H

#include "swp_params.h"
#include <iostream>

namespace Sweep
{
// Forward declare required classes.
class Cell;
class Particle;

namespace AggModels {
class Primary;
class PriPartPrimary;
}

namespace Processes
{
class MeltingModel
{
public:

	//! Model types.
	enum TransformationType {
		Melting,	// Melting temperature dependent transformation
		Gibbs,      // Gibbs energy of transformation: fixed crossover size and temperature
		Composition // (Oxygen/Titanium) composition dependent transformation
	};

    //! Constructors.
	MeltingModel(); // Default constructor.
	MeltingModel(const MeltingModel &copy); // Copy-constructor.
	MeltingModel(      // Stream-reading constructor.
        std::istream &in //  - Input stream.
        );        

    //! Destructor.
	~MeltingModel(void);

    //! Operators.
	MeltingModel &operator=(const MeltingModel &rhs);

    // MODEL ENABLE/DISABLE.

    //! Returns true if the model is enabled, otherwise false.
    bool IsEnabled(void) const;

    //! Enables the transformation model.
    void Enable(void);

    //! Disables the transformation model.
    void Disable(void);

    // READ/WRITE/COPY.

    //! Creates a copy of the model.
	MeltingModel *const Clone(void) const;

    //! Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    void Deserialize(
        std::istream &in // Input stream.
        );

	// MODEL TYPE

	//! Returns the model type.
	TransformationType Type(void) const;

	//! Sets the model type.
	void SetType(TransformationType t);

	// DATA ACCESS

	//! Set oxygen component indices
	void SetOxygenComp(std::vector<unsigned int> OIndex);

	//! Return oxygen component indices
	std::vector<unsigned int> GetOxygenComp() const;

	//! Set titanium component indices
	void SetTitaniumComp(std::vector<unsigned int> TiIndex);

	//! Return titanium component indices
	std::vector<unsigned int> GetTitaniumComp() const;

	//! Return number of phase transformations
	unsigned int PhaseChangeCount(void) const;

	//! Returns true if the phase transformation exists
	bool PhaseChangeValid(std::string id) const;

	// INITIALISATION

	//! Add a phase transformation
	void AddPhase(std::string name, double A, double T, std::vector<fvector> ddcomp, unsigned int iid);

	//! Add phase crossover diameter
	void AddCrossover(std::string name, std::string id_l, std::string id_r, double aa, double ab, double ac, double ad);

	// PERFORM TRANSFORMATION

	//! Returns composition change for melting (crystal -> non-crystalline)
	void MeltingCompositionChange(std::vector<fvector> &dcomp) const;

	//! Return composition change for non-crystalline -> crystal transformation
	void CompositionChange(unsigned int iid, std::vector<fvector> &dcomp);

	//! Return composition change for given non-crystalline -> crystal phase transformation
	void CompositionChange(const Cell &sys, const AggModels::Primary &p, std::vector<fvector> &dcomp);

	//! Return true if particle is above the transformation temperature
	bool IsLiquid(const Cell &sys, const AggModels::Primary &p) const;

private:

    //! Model on/off flag.
    bool m_enable;

	//! Transformation type
	TransformationType m_type;

	//! Oxygen component indices
	std::vector<unsigned int> m_OIndex;

	//! Titanium component indices
	std::vector<unsigned int> m_TiIndex;

	//! Melting composition changes (reverse of phase changes)
	// A separate change in composition must be specified for each type of element of each type of solid phase
	std::vector<fvector> m_dcompmelt;

	//! Return highest melting temperature for particle
	double MeltingTemp(const AggModels::Primary &p) const;

	//! Return O:Ti ratio
	double OTiRatio(const AggModels::Primary &p) const;

	// STRUCTURES

	//! Phase transformation
	struct PHASE
	{
		std::string name;	//! name (crystal phase name/symbol)
		// Melting temperature parameters T = [Tbulk]*(1-[A]/d)
		double A;		//! parameter A
		double T_bulk;	//! bulk melting temperature
		std::vector<fvector> dcomp; //! change in composition non-crystalline (liquid) -> solid crystal
		// The composition change must be separated by element so that adjust is applied 
		// the correct number of times for each element during a transformation
		unsigned int id;	//! Phase id

		// CONSTRUCTORS

		//! Default constructor.
		PHASE(void);
		
		//! Initialising constructor.
		PHASE(std::string aname, double aA, double aT, std::vector<fvector> adcomp, unsigned int iid);

		//! Copy constructor.
		PHASE(const PHASE &copy);

		// FUNCTIONS

		//! Clone
		PHASE* const Clone(void) const;
		
		//! Returns the melting point for given particle size
		double MeltingPoint(double d);
	};

	//! Crossover diameter between solid phases for the "Composition" and "Gibbs" models
	struct CROSSOVER
	{
		std::string name;	// name
		//pointers to phase transformations
		PHASE* phase_l;	//! pointer to phase transformation for dp < dcross
		PHASE* phase_r;	//! pointer to phase transformation for dp > dcross
		//! crossover size parameters c4e-238 Eq. (30) lambda_2...5
		double a, b, c, d;

		// CONSTRUCTORS
		
		//! Default constructor.
		CROSSOVER(void);

		//! Initialising constructor (Composition model).
		CROSSOVER(std::string aname, PHASE* aphase_l, PHASE* aphase_r,
			double aa, double ab, double ac, double ad);

		//! Initialising constructor (Gibbs model).
		CROSSOVER(std::string aname, PHASE* aphase_l, PHASE* aphase_r, double aa);

		//! Copy constructor.
		CROSSOVER(const CROSSOVER &copy);

		// FUNCTIONS

		//! Clone
		CROSSOVER* const Clone(void) const;

		//! Crossover diameter
		double CrossoverDiameter(double sigma_o) const; //composition dependent diameter

		double CrossoverDiameter() const;	//fixed crossover diameter

		//! Return phase transformation
		const PHASE* const GetTransformation(double dp, double sigma_o) const; //composition dependent
		
		const PHASE* const GetTransformation(double dp) const; //fixed crossover diameter
	};

	//! Define a vector of phases
	typedef std::vector<PHASE*> PhasePtrVector;
	PhasePtrVector m_phases;

	//The model currently supports only a single crossover
	CROSSOVER* m_crossover; //! Crossover diameter
	//	// Define a vector of crossovers
	//	typedef std::vector<CROSSOVER*> CrossoverPtrVector;
	//	CrossoverPtrVector m_crossovers;

};
};
};

#endif
