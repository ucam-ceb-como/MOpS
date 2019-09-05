/*
  Author(s):      Casper Lindberg (csl37)
  Project:        sweep (population balance solver)
  
  Copyright (C) 2019 Casper Lindberg.

  File purpose:
    Implement a melting point dependent phase change

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
/*! Model for rate at which primary particles merge with their
 *  neighbours.
 */
class MeltingModel
{
public:

    // Constructors.
	MeltingModel(); // Default constructor.
	MeltingModel(const MeltingModel &copy); // Copy-constructor.
	MeltingModel(      // Stream-reading constructor.
        std::istream &in //  - Input stream.
        );        

    // Destructor.
	~MeltingModel(void);

    // Operators.
	MeltingModel &operator=(const MeltingModel &rhs);


    // MODEL ENABLE/DISABLE.

    // Returns true is the model is enabled, otherwise false.
    bool IsEnabled(void) const;

    // Enables this sintering model.
    void Enable(void);

    // Disables this sintering model.
    void Disable(void);


    // READ/WRITE/COPY.

    // Creates a copy of the sintering model.
	MeltingModel *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in // Input stream.
        );

	// DATA ACCESS METHODS

	void EnableCrossover();

	bool IsEnableCrossover() const;

	void AddPhase(std::string name, double A, double T, double dmin,
		double dmax, std::vector<fvector> ddcomp, unsigned int iid);

	void MeltingCompositionChange(std::vector<fvector> &dcomp) const;

	void CompositionChange(unsigned int iid, std::vector<fvector> &dcomp);

	void CompositionChange(const Cell &sys, const AggModels::Primary &p, std::vector<fvector> &dcomp);

	bool IsLiquid(const Cell &sys, const AggModels::Primary &p) const;

private:

    // Model on/off flag.
    bool m_enable;

	// fixed crossover flag
	bool m_crossover;

	//Get highed melting temperature
	double MeltingTemp(const AggModels::Primary &p) const;

	//Melting composition changes
	//A separate change in composition must be specified for each type of element of each type of solid phase
	std::vector<fvector> m_dcompmelt;

	// phase changes
	struct PHASE
	{
		std::string name;
		double A; // parameter
		double T_bulk; // bulk melting temperature
		std::vector<fvector> dcomp; // change in composition
		double dmin;   // Diameter range for phase change
		double dmax;
		unsigned int id;	// Phase id

		// CONSTRUCTORS

		// Default constructor.
		PHASE(void);
		
		// Initialising constructor.
		PHASE(std::string aname, double aA, double aT, double admin, double admax, std::vector<fvector> adcomp, unsigned int iid);

		// Copy constructor.
		PHASE(const PHASE &copy);

		// FUNCTIONS

		// Clone
		PHASE* const Clone(void) const;
		
		// Returns melting point for given particle size
		double MeltingPoint(double d);
	};

	// define a vector of phases
	typedef std::vector<PHASE*> PhasePtrVector;
	PhasePtrVector m_phases;

};
};
};

#endif
