/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Defines a condensation reaction.  Condensation reactions are treated
    differently from surface reactions as they are modelled as free molecular
    collisions.  The assumptions made in this model are:

    1.  The colliding species is much smaller than the recipient particle.

    Before the condensation process can be used it must be provided the mass and
    diameter of the condensing species in order to calculate the rate terms.  The
    SetCondensingSpecies() function is used for this.

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

#ifndef SWEEP_CONDENSATION_H
#define SWEEP_CONDENSATION_H

#include "swp_particle_process.h"
#include "swp_process_type.h"
#include <iostream>

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Processes
{
class Condensation : public ParticleProcess
{
public: 
    // Constructors.
    Condensation(const Sweep::Mechanism &mech); // Default constructor.
    Condensation(const Condensation &copy);     // Copy constructor.
    Condensation(                    // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    // Destructor.
    ~Condensation(void);

    // Operators.
    Condensation &operator=(const Condensation &rhs);


    // RATE CONSTANT AND PARAMETERS.

    // Sets the coagulation kernel parameters given the mass and
    // collision diameter of the condensing species.
    void SetCondensingSpecies(
        double m, // Mass of the condensing species.
        double d  // Diameter of the condensing species.
        );


    // TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

    // Returns rate of the process for the given system.
    double Rate(
        double t,          // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d &local_geom
        ) const;


    // SINGLE PARTICLE RATE CALCULATIONS.

    // Returns the rate of the process for the given particle in
    // the system. Process must be linear in particle number.
    double Rate(
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
    unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a double vector. The 
    // iterator is advanced to the position after the last term for this
    // process.
    double RateTerms(
        double t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        const Geometry::LocalGeometry1d &local_geom,                  // position information
        fvector::iterator &iterm // Iterator to the first term.
        ) const;


    // PERFORMING THE PROCESS.

    //! Perform one condensation event
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng
        ) const;

    // Performs the process on a given particle in the system.  Particle
    // is given by index.  The process is performed n times.
    int Perform(
        double t,        // Current time (s).
        Cell &sys,     // System to which the particle belongs.
        Particle &sp,  // Particle for which to perform process.
        rng_type &rng, // Random number for leaf node
        unsigned int n // Number of times to perform the process.
        ) const;


    // READ/WRITE/COPY.
    
    // Creates a copy of the particle process.
    virtual Condensation *const Clone(void) const;

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
    // Number of terms in the condensation rate expression.
    static const unsigned int TERM_COUNT;

    // Condensation majorant parameter.  The true rate
    // is multiplied by this parameter to get the majorised rate.
    static const double m_majfactor;

    // Free-molecular enhancement factor.
    const double m_efm;

    double m_kfm1, m_kfm2, m_kfm3; // Free-mol term parameters.

    // Default constructor is protected to prevent condensations being
    // defined without knowledge of the parent mechanism.
    Condensation(void);

    // Adjusts a primary particle according to the rules of the condensation.
    unsigned int adjustPri(
        Sweep::AggModels::Primary &pri, // Primary to adjust.
        rng_type &rng,
        unsigned int n=1     // Number of times to perform adjustment.
        ) const;
};
};
};

#endif
