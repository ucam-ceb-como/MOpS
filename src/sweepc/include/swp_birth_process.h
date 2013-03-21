 /*!
  * @file   swp_birth_process.h
  * @author Matthew Celnik, William Menz
  * @brief  Declaration of a birth process
  *
  *   About:
  *      Birth processes are used to model an inflow of particles from another
  *      Cell. The rate of inflow depends on the number of particles in the
  *      upstream Cell. The process can be chosen to be stochastic (using the
  *      usual jump process infrastructure) or continuous (like LPDA).
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
    //! Initialising constructor
    BirthProcess(const Sweep::Mechanism &mech);

    //! Copy constructor
    BirthProcess(const BirthProcess &copy);

    //! Assignment operator
    BirthProcess &operator=(const BirthProcess &rhs);

    //! Does the Cell inflow have particles present?
    bool HasParticlesInCell() const;

	// TOTAL RATE CALCULATIONS.

    //! Returns rate of the process for the given system.
    double Rate(
        double t,          // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d &local_geom
        ) const;

	// RATE TERM CALCULATIONS.

    //! Returns the number of rate terms for this process.
    unsigned int TermCount(void) const;

    //! Calculates the rate terms given an iterator to a double vector.
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

    //! Perform the process over time dt
    void PerformDT (
            const double t,
            const double dt,
            Sweep::Cell &sys,
            const Geometry::LocalGeometry1d& local_geom,
            rng_type &rng) const;

    // READ/WRITE/COPY.

    //! Creates a copy of the inception.
    BirthProcess *const Clone(void) const;

    //! Returns the process type.
    ProcessType ID(void) const;

    //! Available birth processes
    enum BirthType {
        iStochastic,        // A jump process
        iContinuous         // A continuous process
    };


    // SET PROCESS PROPERTIES

    //! Set the Cell from which this process samples.
    void SetCell(Cell* c);

    //! Sets the birth process type
    void SetBirthType(const BirthType t);

    //! Turn the process on or off
    void SetProcessSwitch(const bool s);

protected:

    //! Particle sampling.
    Cell *m_cell; // The cell from which to sample particle.

    //! Protected default constructor.
    BirthProcess(void);

private:

    //! Returns the scaling factor for moving particles between cells
    double F(const Cell &sys) const;

    //! Do the birth for particle of index isp
    void DoParticleBirth(
            const double t,
            const int isp,
            Sweep::Cell &sys,
            const double wt,
            rng_type &rng) const;

    //! The process birth type
    BirthType m_btype;

    //! Is the process on or off?
    bool m_on;

    //! Helper function to get the rate.
    double InternalRate(
            double t,
            const Cell &sys,
            const Geometry::LocalGeometry1d &local_geom) const;
};
typedef std::vector<BirthProcess*> BirthPtrVector;
};
};

#endif
