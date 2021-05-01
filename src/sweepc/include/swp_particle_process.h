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

#include <iostream>
#include <vector>

namespace Sweep
{
// Forward declare the Mechanism class.
class Mechanism;

namespace Processes
{
class ParticleProcess;
typedef std::vector<ParticleProcess*> PartProcPtrVector;

/*!
 * \brief Processes that change the internal state of particles without transport or interaction
 */
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
    double CompChange(unsigned int i) const;

    // Sets the particle composition vector.
    void SetCompChange(const fvector &comp);

    // Sets the amount of the ith component in the new particle.
    void SetCompChange(unsigned int i, double comp);

    // Returns the tracker variable vector of the new particle.
    const fvector &TrackChange(void) const;

    // Returns the value of the ith tracker variable of the
    // new particle.
    double TrackChange(unsigned int i) const;

    // Sets the new particle tracker variable vector.
    void SetTrackChange(const fvector &track);

    // Sets the value of the ith tracker variable in the
    // new particle.
    void SetTrackChange(unsigned int i, double track);

    
    // TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

    // Returns rate of the process for the given system.
    virtual double Rate(
        double t,          // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d &local_geom
        ) const = 0;

    // Return rate constant and chemistry part for hybrid method
    virtual double Rate(
	double dt,          // Time.
	const Cell &sys     // System for which to calculate rate.
        ) const;

    static double CalcRates(
        double t, // Time.
        const Cell &sys, // System for which to calculate rates.
        const Geometry::LocalGeometry1d &local_geom, // Position information
        const PartProcPtrVector &proc, // Vector of processes to calculate.
        fvector &rates, // Return vector for rates.
        unsigned int start = 0 // Start index in output vector.
        );

	// SINGLE PARTICLE RATE CALCULATIONS.

    // Returns the rate of the process for the given particle in
    // the system. Process must be linear in particle number.
    virtual double Rate(
        double t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const = 0;


	// RATE TERM CALCULATIONS.
    //   These routines return the individual rate terms for a 
    //   process, which may have multiple terms (e.g. condensation).

	// PERFORMING THE PROCESS.

    //! Declared in parent class, implemented in dervied classes
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng
        ) const = 0;

    // Performs the process on a given particle in the system.  Particle
    // is given by index.  The process is performed n times.
    virtual int Perform(
        double t,        // Current time (s).
        Cell &sys,     // System to which the particle belongs.
        Particle &sp,  // Particle for which to perform process.
        rng_type &rng,
        unsigned int n // Number of times to perform the process.
        ) const = 0;

	// Surface growth updates for the hybrid particle model (particle-number updates)
	// ==============================================================================
	// Just do gas-phase adjustment for surface growth
	virtual int Perform(
		double t,                        // Current time (s).
		Cell &sys,                       // System to which the particle belongs.
		rng_type &rng,                   // Random number generator
		unsigned int n                   // Number of times to perform the process.)
		) const;

	// Performs the process on a given particle in the system.  Particle
	// is given by index.  The process is performed n times.
	virtual int Perform(
		double t,                        // Current time (s).
		Cell &sys,                       // System to which the particle belongs.
		Particle &sp,                    // Particle for which to perform process.
		rng_type &rng,                   // Random number generator
		unsigned int n,                  // Number of times to perform the process.
		bool isParticleNumberUpdate      // Differentiate this call as a particle-number update
		) const;
	// ==============================================================================

    // READ/WRITE/COPY.
    
    // Creates a copy of the particle process.
    virtual ParticleProcess *const Clone(void) const = 0;

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
}
}

#endif
