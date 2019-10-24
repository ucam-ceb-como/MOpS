/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Solver class defines the stochastic stepping algorithm and
    holds information about the progress of the solution.

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

#ifndef SWEEP_SOLVER_H
#define SWEEP_SOLVER_H

#include "swp_mechanism.h"
#include "swp_cell.h"

#include <vector>
#include <map>

// Forward declaration
namespace Geometry
{
    class LocalGeometry1d;
}

namespace Sweep
{
class Solver
{
public:
    // Constructors.
    Solver(void); // Default constructor.

    //! Copy constructor
    Solver(const Solver &sol);

    // Destructor.
    ~Solver(void);

    // Performs stochastic stepping algorithm up to specified stop time using
    // the given mechanism to define the stochastic processes.  Updates given
    // system accordingly.  On error returns <0, otherwise returns 0.
    int Run(
        double &t,        // Simulation start time.  Will return the stop time.
        double tstop,     // Stop time for simulation.
        Cell &sys,      // System to solve.
        const Mechanism &mech, // Mechanism to use to solve system.
        rng_type &rng
        );

    //! Performs a single stochastic event on the ensemble
    static void timeStep(
        double &t,                // Current solution time.
        double t_stop,            // Steps may not go past this time
        Cell &sys,              // System to update.
        const Geometry::LocalGeometry1d &geom, // Details of cell size
        const Mechanism &mech,  // Mechanism to use.
        const fvector &rates,   // Current process rates as an array.
        double jrate,             // The total jump rate (non-deferred processes).
        rng_type &rng
        );

protected:
    // TIME STEPPING ROUTINES.

    // Calculates the splitting end time after which all particles
    // shallbe updated using LPDA.
    double calcSplitTime(
        double t,         // Current time.
        double tstop,     // Stop time.
        double jrate,     // Sum of all jump process rates.
        unsigned int n  // Number of particles.
        ) const;

    // Selects a process using a DIV algorithm and the process rates
    // as weights.
    static int chooseProcess(const fvector &rates, double (*rand_u01)());

private:
    // Numerical parameters.

    //! Parameter defining number of LPDA updates per particle events.
    double m_splitratio;


};
}

#endif
