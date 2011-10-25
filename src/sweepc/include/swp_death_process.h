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
#include "swp_process_type.h"
#include "swp_process.h"

#include <iostream>
#include <vector>

namespace Sweep
{
// Forward declare the Mechanism class.
class Mechanism;
// Forward declare the Cell class.
class Cell;

namespace Transport
{
    // Forward declare structure to hold details of particle for onward transport
    struct TransportOutflow;
}

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
        real t,          // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d &local_geom
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
        const Geometry::LocalGeometry1d &local_geom,                  // position information
        fvector::iterator &iterm // Iterator to the first term.
        ) const;


	// PERFORMING THE PROCESS.

    //! Kill one particle
    virtual int Perform(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng
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
typedef std::vector<DeathProcess*> DeathPtrVector;
};
};

#endif
