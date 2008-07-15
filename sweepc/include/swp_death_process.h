/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of a uniform death process.  Particles are uniformly
    deleted at the calculated rate.  This provides a method for
    simulating Cell outflow events.

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

#ifndef SWEEP_DEATH_PROCESS_H
#define SWEEP_DEATH_PROCESS_H

#include "swp_params.h"
#include "swp_process.h"
#include "swp_process_type.h"
#include "swp_cell.h"
#include "sprog.h"
#include <iostream>

namespace Sweep
{
// Forward declare the Mechanism class.
class Mechanism;

namespace Processes
{
class DeathProcess : public Process
{
public: 
    // Constructors.
    DeathProcess(const Sweep::Mechanism &mech); // Initialising constructor.
    DeathProcess(const DeathProcess &copy);     // Copy constructor.
    DeathProcess(                               // Stream-reading constructor.
        std::istream &in,                       //  - Input stream.
        const Sweep::Mechanism &mech            //  - Parent mechanism.
        );

    // Destructors.
    ~DeathProcess(void);

    // Operators.
    DeathProcess &operator=(const DeathProcess &rhs);


    // RATE CONSTANT.

    // Returns the rate constant.
    real A(void) const;

    // Sets the rate constant.
    void SetA(real a);


	// TOTAL RATE CALCULATIONS.

    // Returns rate of the process for the given system.
    real Rate(
        real t,         // Time.
        const Cell &sys // System for which to calculate rate.
        ) const;

	// Calculates the process rate using the given 
    // chemical conditions, rather than those conditions in the
    // given system.
    real Rate(
        real t,                   // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys           // System for which to calculate rate.
        ) const;


	// RATE TERM CALCULATIONS.

    // Returns the number of rate terms for this process.
    unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all terms.
    real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  The given chemical conditions are used instead of those
    // in the given system object.
    real RateTerms(
        real t,                   // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys,          // System for which to calculate rate terms.
        fvector::iterator &iterm  // Iterator to the first term.
        ) const;


	// PERFORMING THE PROCESS.

    // Performs the process on the given system.  The responsible rate term is given
    // by index.  Returns 0 on success, otherwise negative.
    int Perform(
        real t,                // Time.
        Cell &sys,             // System to update.
        unsigned int iterm = 0 // The process term responsible for this event.
        ) const;


    // READ/WRITE/COPY.

    // Creates a copy of the inception.
    DeathProcess *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    ProcessType ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Rate parameters.
    real m_a; // Rate constant.

    // Default constructor is protected to prevent a process being
    // defined without knowledge of the parent mechanism.
    DeathProcess(void);
};
};
};

#endif
