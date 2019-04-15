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

	void SetLiquidIndex(unsigned int index);

	unsigned int GetLiquidIndex(void) const;

	void AddPhase(std::string name, double A, double T, fvector dcomp);

	void CompositionChange(const Cell &sys, const AggModels::Primary &p, fvector &dcomp);

	bool IsLiquid(const Cell &sys, const AggModels::Primary &p) const;

private:

    // Model on/off flag.
    bool m_enable;

	// index of liquid phase
	unsigned int m_liquidindex; 

	//Get highed melting temperature
	double MeltingTemp(const AggModels::Primary &p) const;

	// phase changes
	struct PHASE
	{
		std::string name;
		double A; // parameter
		double T_bulk; // bulk melting temperature
		fvector dcomp; // change in composition

		// CONSTRUCTORS

		// Default constructor.
		PHASE(void);
		
		// Initialising constructor.
		PHASE(std::string aname, double aA, double aT, fvector adcomp);

		// Copy constructor.
		PHASE(const PHASE &copy);

		//FUNCTIONS

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
