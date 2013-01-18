/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of a sampled birth process.  Particles are sampled
    from a given distribution and introduced at the defined rate.  If
    no ensemble/particle is defined from which to sample, then an
    empty particle is introduced.

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

#ifndef SWEEP_BIRTH_PROCESS_H
#define SWEEP_BIRTH_PROCESS_H

#include "swp_params.h"
#include "swp_process.h"
#include "swp_process_type.h"
#include "swp_particle.h"

#include <iostream>
#include <vector>

namespace Sweep
{
// Forward declare the Mechanism class.
class Mechanism;
// Forward declare the Cell class.
class Cell;

namespace Processes
{
class BirthProcess : public Process
{
public: 
    // Constructors.
    BirthProcess(const Sweep::Mechanism &mech); // Initialising constructor.
    BirthProcess(const BirthProcess &copy);     // Copy constructor.
    BirthProcess(                               // Stream-reading constructor.
        std::istream &in,                       //  - Input stream.
        const Sweep::Mechanism &mech            //  - Parent mechanism.
        );

    // Destructors.
    ~BirthProcess(void);

    // Operators.
    BirthProcess &operator=(const BirthProcess &rhs);


    // INFORMATION FOR THE SOLVER
    //! Does the Cell inflow have particles present?
    bool HasParticlesInCell() const;

    // SET PROCESS PROPERTIES

    //! Set the Cell from which this process samples.
    void SetCell(Cell* c) {m_cell = c;}

    //! Set the default particle.
    void SetParticle(Particle* p) {m_particle = p;}

    // Returns a pointer to the default particle used for
    // the birth process.

	// TOTAL RATE CALCULATIONS.

    // Returns rate of the process for the given system.
    double Rate(
        double t,          // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d &local_geom
        ) const;

	// RATE TERM CALCULATIONS.

    // Returns the number of rate terms for this process.
    unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a double vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all terms.
    double RateTerms(
        double t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        const Geometry::LocalGeometry1d &local_geom,                  // position information
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

	// PERFORMING THE PROCESS.

    //! Give birth to one particle
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng
        ) const;


    // READ/WRITE/COPY.

    // Creates a copy of the inception.
    BirthProcess *const Clone(void) const;

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

    // Particle sampling.
    Particle *m_particle; // A uniform birth process particle.
    Cell *m_cell; // The cell from which to sample particle.

    // Default constructor is protected to prevent a process being
    // defined without knowledge of the parent mechanism.
    BirthProcess(void);
};
typedef std::vector<BirthProcess*> BirthPtrVector;
};
};

#endif
