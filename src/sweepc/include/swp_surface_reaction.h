/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of a surface reaction process.  Surface reactions are single
    particle processes, hence are deferred by default.

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

#ifndef SWEEP_SURFRXN_H
#define SWEEP_SURFRXN_H

#include "swp_params.h"
#include "swp_particle_process.h"
#include "swp_property_indices.h"

#include "gpc_rate_params.h"

namespace Sweep
{
// Forward declare mechanism class.
class Mechanism;

namespace Processes
{
/*!
 * \brief Arrhenius like reactions on the surface of particles
 */
class SurfaceReaction : public ParticleProcess
{
public:
    // Constructors.
    SurfaceReaction(const Sweep::Mechanism &mech); // Default constructor.
    SurfaceReaction(const SurfaceReaction &copy);  // Copy constructor.
    SurfaceReaction(                 // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism
        );

    // Operators.
    SurfaceReaction &operator=(const SurfaceReaction &rhs);


    // ARRHENIUS COEFFICIENTS.

    // Returns the Arrhenius parameter.
    Sprog::Kinetics::ARRHENIUS &Arrhenius();
    const Sprog::Kinetics::ARRHENIUS &Arrhenius() const;

    // Sets the Arrhenius parameters.
    void SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr);


    // PARTICLE PROPERTY ID.

    // Returns the ID number of the particle property to which
    // the rate of this process is proportional.
    Sweep::PropID PropertyID(void) const;

    //! ID number of the particle property to which the rate of this process is proportional.
    void SetPropertyID(
        Sweep::PropID pid
        );


    // TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

    // Returns rate of the process for the given system.
    virtual double Rate(
        double t,         // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d &local_geom // cell location information
        ) const;

    // Return rate constant and chemistry part for hybrid method
    virtual double Rate(double t, const Cell &sys) const;

	// SINGLE PARTICLE RATE CALCULATIONS.

    // Returns the rate of the process for the given particle in
    // the system. Process must be linear in particle number.
    virtual double Rate(
        double t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const;

    // Returns majorant rate of the process for the given system.
    double MajorantRate(
        double t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const;


	// RATE TERM CALCULATIONS.
    //   These routines return the individual rate terms for a
    //   process, which may have multiple terms (e.g. condensation).

    // Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a double vector. The
    // iterator is advanced to the position after the last term for this
    // process.
    virtual double RateTerms(
        double t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        const Geometry::LocalGeometry1d &local_geom, //Position information
        fvector::iterator &iterm // Iterator to the first term.
        ) const;


    // PERFORMING THE PROCESS.

    //! Perform one surface reaction event
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng
        ) const;

    // Performs the process on a given particle in the system.  Particle
    // is given by index.  The process is performed n times.
    virtual int Perform(
        double t,        // Current time (s).
        Cell &sys,     // System to which the particle belongs.
        Particle &sp,  // Particle for which to perform process.
        rng_type &rng,
        unsigned int n // Number of times to perform the process.
        ) const;

	// Surface growth updates for the hybrid particle model (particle-number updates)
	// ==============================================================================
	// Just do gas-phase adjustment for surface growth
	virtual int Perform(double t,        // Current time (s).
		Cell &sys,                       // System to which the particle belongs.
		rng_type &rng,                   // Random number generator
		unsigned int n) const;           // Number of times to perform the process.)

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
    virtual SurfaceReaction *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    virtual ProcessType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Surface reaction majorant parameter.  The true rate
    // is multiplied by this parameter to get the majorised rate.
    const static double m_majfactor;

    //! Arrhenius rate parameters.
    Sprog::Kinetics::ARRHENIUS m_arr;

    //! Particle property to which the rate of the process is proportional.
    Sweep::PropID m_pid;

    // Default constructor is protected to prevent reactions being
    // defined without knowledge of the parent mechanism.
    SurfaceReaction(void);

    // Adjusts a primary particle according to the rules of the reaction.
    unsigned int adjustPri(
        Sweep::AggModels::Primary &pri, // Primary to adjust.
        rng_type &rng,
        unsigned int n=1     // Number of times to perform adjustment.
        ) const;
};
}
}

#endif
