/*
  Author(s):      Shraddha Shekar (ss663)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2010 Shraddha Shekar.

  File purpose:
    Defines a InterParticle reaction.  InterParticle reactions are treated
    differently from surface reactions as they are modelled as free molecular
    collisions.  The assumptions made in this model are:

    1.  The colliding species is much smaller than the recipient particle.

    Before the InterParticle process can be used it must be provided the mass and
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

#ifndef SWEEP_INTERPART_H
#define SWEEP_INTERPART_H

#include "swp_particle_process.h"
#include "swp_process_type.h"

#include "gpc_rate_params.h"

#include <iostream>


namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Processes
{
class InterParticle : public ParticleProcess
{
public:
    //! Default constructor
    InterParticle(const Sweep::Mechanism &mech,
                  const EnvironmentInterface::SpeciesIndex h4o4si);
    
    //! Copy constructor
    InterParticle(const InterParticle &copy);\
    
    //! Stream-reading constructor
    InterParticle(                    // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    //! Destructor.
    ~InterParticle(void);

    //! Equal operator definition
    InterParticle &operator=(const InterParticle &rhs);



	// PARTICLE PROPERTY ID.

    //! Returns the PropID to which the rate is proportional
    unsigned int PropertyID(void) const;

    //! Sets the PropID to which the rate is proportional
    void SetPropertyID(
        Sweep::PropID pid
        );

    //! Index for this species which appears in various calculations
    void SetH4O4SI_Index(const EnvironmentInterface::SpeciesIndex index) {m_H4O4SI_Index = index;}

	// RATE CONSTANT AND PARAMETERS.

    //! Returns the Arrhenius parameter.
    Sprog::Kinetics::ARRHENIUS &Arrhenius();
    const Sprog::Kinetics::ARRHENIUS &Arrhenius() const;

    //! Sets the Arrhenius parameters.
    void SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr);


    // TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

     //! Returns rate of the process for the given system.
     double Rate(
         double t,          // Time.
         const Cell &sys, // System for which to calculate rate.
         const Geometry::LocalGeometry1d& local_geom // Information regarding surrounding cells and boundaries
         ) const;


     // SINGLE PARTICLE RATE CALCULATIONS.

     //! Returns the rate of the process for the given particle.
     double Rate(
         double t,             // Current time (s).
         const Cell &sys,    // System to which the particle belongs.
         const Particle &sp  // Particle for which to calculate rate.
         ) const;


     //! Returns majorant rate of the process for the given system.
     double MajorantRate(
         double t,             // Current time (s).
         const Cell &sys,    // System to which the particle belongs.
         const Particle &sp  // Particle for which to calculate rate.
         ) const;


 	// RATE TERM CALCULATIONS.
     //   These routines return the individual rate terms for a
     //   process, which may have multiple terms (e.g. InterParticle).

     //! Returns the number of rate terms for this process.
     unsigned int TermCount(void) const;


     //! Calculates the rate terms given an iterator to a double vector.
     double RateTerms(
         double t,                  // Time.
         const Cell &sys,         // System for which to calculate rate terms.
         const Geometry::LocalGeometry1d& local_geom, // Information regarding surrounding cells and boundaries
         fvector::iterator &iterm // Iterator to the first term.
         ) const;


     // PERFORMING THE PROCESS.

     //! Performs one interparticle reaction event
     virtual int Perform(
         double t,
         Cell &sys,
         const Geometry::LocalGeometry1d& local_geom,
         unsigned int iterm,
         rng_type &rng) const;

     //! Performs n interparticle reaction events
     virtual int Perform(
         double t,        // Current time (s).
         Cell &sys,     // System to which the particle belongs.
         Particle &sp,  // Particle for which to perform process.
         rng_type &rng,  // Random generator
         unsigned int n // Number of times to perform the process.
         ) const;


     // READ/WRITE/COPY.

     //! Creates a copy of the particle process.
     virtual InterParticle *const Clone(void) const;

     //! Returns the process type.
     virtual ProcessType ID(void) const;

     //! Writes the object to a binary stream.
     virtual void Serialize(std::ostream &out) const;

     //! Reads the object from a binary stream.
     virtual void Deserialize(
         std::istream &in,            // Input stream.
         const Sweep::Mechanism &mech // Parent mechanism.
         );

protected:
    //! Surface reaction majorant parameter. 
    const static double m_majfactor;

    //! Arrhenius rate parameters.
    Sprog::Kinetics::ARRHENIUS m_arr;

    //! Particle property to which the rate of the process is proportional.
    Sweep::PropID m_pid;

    //! Index for H4O4SI
    EnvironmentInterface::SpeciesIndex m_H4O4SI_Index;

     // Default constructor is protected to prevent InterParticles being
     // defined without knowledge of the parent mechanism.
     InterParticle(void);

}; // InterParticle class
}  // Process namespace
}  // Sweep namespace

#endif

