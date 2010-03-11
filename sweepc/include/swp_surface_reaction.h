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
#include "swp_submodel_type.h"

namespace Sweep
{
// Forward declare mechanism class.
class Mechanism;

namespace Transport
{
    // Forward declaration of unused argument type
    struct TransportOutflow;
}

namespace Processes
{
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

    // Destructor.
    virtual ~SurfaceReaction(void);

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
    unsigned int PropertyID(void) const;

    // Returns the ID number of the particle number for which the
    // PropertyID is valid.  The mechanism should check that this
    // model is enabled.
    SubModels::SubModelType ModelID(void) const;

    // Sets the ID number of the particle property to which
    // the rate of this process is proportional.
    void SetPropertyID(
        unsigned int i,   // ID number of particle property.
        SubModels::SubModelType modelid // The model for which this ID is valid.
          = SubModels::BasicModel_ID //  - Default model is basic particle properties.
        );


    // TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

    // Returns rate of the process for the given system.
    virtual real Rate(
        real t,         // Time.
        const Cell &sys // System for which to calculate rate.
        ) const;


	// SINGLE PARTICLE RATE CALCULATIONS.

    // Returns the rate of the process for the given particle in
    // the system. Process must be linear in particle number.
    virtual real Rate(
        real t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const;

    // Returns majorant rate of the process for the given system.
    real MajorantRate(
        real t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const;


	// RATE TERM CALCULATIONS.
    //   These routines return the individual rate terms for a
    //   process, which may have multiple terms (e.g. condensation).

    // Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a real vector. The
    // iterator is advanced to the position after the last term for this
    // process.
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        fvector::iterator &iterm // Iterator to the first term.
        ) const;


    // PERFORMING THE PROCESS.

    //! Perform one surface reaction event
    virtual int Perform(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        int (*rand_int)(int, int), 
        real(*rand_u01)(),
        Transport::TransportOutflow *out = 0
        ) const;

    // Performs the process on a given particle in the system.  Particle
    // is given by index.  The process is performed n times.
    virtual int Perform(
        real t,        // Current time (s).
        Cell &sys,     // System to which the particle belongs.
        Particle &sp,  // Particle for which to perform process.
        unsigned int n // Number of times to perform the process.
        ) const;


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
    const static real m_majfactor;

    // Arrhenius rate parameters.
    Sprog::Kinetics::ARRHENIUS m_arr;

    // Particle property to which the rate of the process is
    // proportional.
    unsigned int m_pid;

    // Particle model for which the above particle property ID
    // is valid.
    SubModels::SubModelType m_modelid;

    // Default constructor is protected to prevent reactions being
    // defined without knowledge of the parent mechanism.
    SurfaceReaction(void);

    // Adjusts a primary particle according to the rules of the reaction.
    unsigned int adjustPri(
        Sweep::Primary &pri, // Primary to adjust.
        unsigned int n=1     // Number of times to perform adjustment.
        ) const;
};
};
};

#endif
