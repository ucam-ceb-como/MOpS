/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ParticleProcess class defines a process which occurs on single
    particles.

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

#ifndef SWEEP_PARTICLEPROCESS_H
#define SWEEP_PARTICLEPROCESS_H

#include "swp_params.h"
#include "swp_process.h"
#include "swp_cell.h"
#include "swp_particle.h"
#include "swp_process_type.h"
#include "sprog.h"
#include <iostream>

namespace Sweep
{
// Forward declare the Mechanism class.
class Mechanism;

namespace Processes
{
class ParticleProcess;
typedef std::vector<ParticleProcess*> PartProcPtrVector;

class ParticleProcess : public Process
{
public:
	/// Constructors.
    ParticleProcess(const Sweep::Mechanism &mech); // Default Constructor.
	ParticleProcess(const ParticleProcess &copy);  // Copy-constructor.
    ParticleProcess(                 // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    // Destructor.
    virtual ~ParticleProcess(void);

    // Operators.
    ParticleProcess &operator=(const ParticleProcess &rhs);


    // PROCESS DEFERMENT.

    // Returns TRUE if process should be deferred, otherwise false.
    bool IsDeferred(void) const;

    // Sets the process to be deferred or not.  Interface version
    // sets deferred off.
    virtual void SetDeferred(bool defer);


    // CHANGES TO PARTICLE ON PROCESS OCCURANCE.

    // Returns the composition vector of the new particle.
    const fvector &CompChange(void) const;

    // Returns the amount of the ith component of the new particle.
    real CompChange(unsigned int i) const;

    // Sets the particle composition vector.
    void SetCompChange(const fvector &comp);

    // Sets the amount of the ith component in the new particle.
    void SetCompChange(unsigned int i, real comp);

    // Returns the tracker variable vector of the new particle.
    const fvector &TrackChange(void) const;

    // Returns the value of the ith tracker variable of the
    // new particle.
    real TrackChange(unsigned int i) const;

    // Sets the new particle tracker variable vector.
    void SetTrackChange(const fvector &track);

    // Sets the value of the ith tracker variable in the
    // new particle.
    void SetTrackChange(unsigned int i, real track);

    
    // TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

    // Returns rate of the process for the given system.
    virtual real Rate(
        real t,         // Time.
        const Cell &sys // System for which to calculate rate.
        ) const = 0;

/*
	// Calculates the process rate using the given 
    // chemical conditions, rather than those conditions in the
    // given system.
    virtual real Rate(
        real t,                   // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys           // System for which to calculate rate.
        ) const = 0;
*/

    static real CalcRates(
        real t, // Time.
        const Cell &sys, // System for which to calculate rates.
        const PartProcPtrVector &proc, // Vector of processes to calculate.
        fvector &rates, // Return vector for rates.
        unsigned int start = 0 // Start index in output vector.
        );

	// SINGLE PARTICLE RATE CALCULATIONS.

    // Returns the rate of the process for the given particle in
    // the system. Process must be linear in particle number.
    virtual real Rate(
        real t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const = 0;

/*
	// Returns rate of the process for the given particle using the
    // given chemical conditions rather than those conditions in the
    // the given system.
    virtual real Rate(
        real t,                   // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys,          // System to which the particle belongs.
        const Particle &sp        // Particle for which to calculate rate.
        ) const = 0;
*/

	// RATE TERM CALCULATIONS.
    //   These routines return the individual rate terms for a 
    //   process, which may have multiple terms (e.g. condensation).

    // Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const = 0;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        fvector::iterator &iterm // Iterator to the first term.
        ) const = 0;

/*
    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  The given chemical conditions are used instead of those
    // in the given system object.
    virtual real RateTerms(
        real t,                   // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys,          // System for which to calculate rate terms.
        fvector::iterator &iterm  // Iterator to the first term.
        ) const = 0;
*/

	// PERFORMING THE PROCESS.

    // Performs the process on the given system.  The responsible rate term is given
    // by index.  Returns 0 on success, otherwise negative.
    virtual int Perform(
        real t,                // Time.
        Cell &sys,             // System to update.
        unsigned int iterm = 0 // The process term responsible for this event.
        ) const = 0;

    // Performs the process on a given particle in the system.  Particle
    // is given by index.  The process is performed n times.
    virtual int Perform(
        real t,        // Current time (s).
        Cell &sys,     // System to which the particle belongs.
        Particle &sp,  // Particle for which to perform process.
        unsigned int n // Number of times to perform the process.
        ) const = 0;


    // READ/WRITE/COPY.
    
    // Creates a copy of the particle process.
    virtual ParticleProcess *const Clone(void) const = 0;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    virtual ProcessType ID(void) const = 0;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    bool m_defer; // Is the process solved by LPDA?

    // Default constructor is protected to prevent processes being
    // defined without knowledge of the parent mechanism.
    ParticleProcess(void);

    // Changes to particle composition and tracker values
    // when this process occurs.
    fvector m_dcomp; // Composition change.
    fvector m_dvals; // Tracker variable value changes.
};
};
};

#endif
