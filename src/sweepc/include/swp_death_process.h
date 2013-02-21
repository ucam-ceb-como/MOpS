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

    //! Kill one particle
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng
        ) const;

    //! Perform the process over time dt
    void PerformDT (
            const double t,
            const double dt,
            Sweep::Cell &sys,
            rng_type &rng) const;

    // READ/WRITE/COPY.

    // Creates a copy of the inception.
    DeathProcess *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    ProcessType ID(void) const;

    //! Available death processes
    enum DeathType {
        iDeathDelete,       // Deletes a particle from the ensemble
        iDeathMove,         // Moves a particle downstream
        iDeathRescale       // Rescale the sample volume of the cell
    };

    //! Set the type of process
    void SetDeathType(const DeathType t) {m_dtype = t;}

    //! Set whether adaptive death should be used
    void SetAdaptive(const bool a) {m_adaptive = a;}

    //! Activate changes to cells when the adaptive process is toggled
    void Adapt(Sweep::Cell &sys);

    //! Set the downstream cell
    void SetCell(Cell* c) {m_cell = c;}

protected:

    // Default constructor is protected to prevent a process being
    // defined without knowledge of the parent mechanism.
    DeathProcess(void);

private:
    //! Flag for type of death process
    DeathType m_dtype;

    //! Should the adaptive death process be used?
    bool m_adaptive;

    //! Has the adaptive process been toggled? (True = ON)
    bool m_toggled;

    //! The downstream cell
    Cell *m_cell;

    //! A helper function for doing the process
    void DoParticleDeath(
            const double t,
            const int isp,
            Sweep::Cell &sys,
            rng_type &rng) const;

};
typedef std::vector<DeathProcess*> DeathPtrVector;
};
};

#endif
