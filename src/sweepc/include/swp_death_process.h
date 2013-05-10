 /*!
  * @file   swp_death_process.h
  * @author Matthew Celnik, William Menz
  * @brief  Declaration of a death process
  *
  *   About:
  *      Death processes model the outflow of particles from a reactor.
  *      Several types have been developed in this class. Like birth processes,
  *      they can be stochastic (using the normal jump rate interface) or
  *      continuous (LPDA-like). Each of the process types is explained below.
  *
  *      Delete: just deletes a particle from the ensemble
  *      Move: move a particle downstream (will turn off any inception terms
  *            for the downstream reactor)
  *      Rescale: Model ouflow by rescaling the sample volume
  *      Adaptive: Change the death process depending on the state of the
  *                particle ensemble
  *
  *   Licence:
  *      sweepc is free software; you can redistribute it and/or
  *      modify it under the terms of the GNU Lesser General Public License
  *      as published by the Free Software Foundation; either version 2
  *      of the License, or (at your option) any later version.
  *
  *      This program is distributed in the hope that it will be useful,
  *      but WITHOUT ANY WARRANTY; without even the implied warranty of
  *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *      GNU Lesser General Public License for more details.
  *
  *      You should have received a copy of the GNU Lesser General Public
  *      License along with this program; if not, write to the Free Software
  *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  *      02111-1307, USA.
  *
  *   Contact:
  *      Prof Markus Kraft
  *      Dept of Chemical Engineering
  *      University of Cambridge
  *      New Museums Site
  *      Pembroke Street
  *      Cambridge
  *      CB2 3RA, UK
  *
  *      Email:       mk306@cam.ac.uk
  *      Website:     http://como.cheng.cam.ac.uk
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
            const Geometry::LocalGeometry1d& local_geom,
            rng_type &rng) const;

    // READ/WRITE/COPY.

    // Creates a copy of the inception.
    DeathProcess *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    ProcessType ID(void) const;

    //! Available death processes
    enum DeathType {
        // 'Continuous' processes
        iContDelete,       // Deletes a particle from the ensemble
        iContMove,         // Moves a particle downstream
        iContRescale,      // Rescale the sample volume of the cell
        iContAdaptive,     // Changes the nature of the death process
        // 'Stochastic' processes
        iStochDelete,
        iStochMove
    };

    //! Set the type of process
    void SetDeathType(const DeathType t);

    //! Get the type of the process
    DeathType GetDeathType() const;

    //! Set the downstream cell
    void SetCell(Cell* c);

    //! Activate changes to cells when the adaptive process is toggled
    void Adapt(Sweep::Cell &sys);

protected:

    // Default constructor is protected to prevent a process being
    // defined without knowledge of the parent mechanism.
    DeathProcess(void);

private:
    //! Flag for type of death process
    DeathType m_dtype;

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

    //! Helper function to get the rate.
    double InternalRate(
            double t,
            const Cell &sys,
            const Geometry::LocalGeometry1d &local_geom) const;

    //! Is the process of stochastic type?
    bool IsStochastic() const;

};
typedef std::vector<DeathProcess*> DeathPtrVector;
};
};

#endif
