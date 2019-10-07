/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This solver only solves the particle population balance by
    using a constant profile gas-phase chemistry profile.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#ifndef SWEEP_FLAME_SOLVER_H
#define SWEEP_FLAME_SOLVER_H

#include "mops_particle_solver.h"
#include "swp_gas_profile.h"
#include "sweep.h"
#include <map>

namespace Sweep
{
class FlameSolver : public Mops::ParticleSolver, public Sweep::Solver
{
public:
    // Constructors.
    FlameSolver(void); // Default constructor.

    //! Copy constructor
    FlameSolver(const FlameSolver &sol);

    //! Clone the object
    FlameSolver *const Clone() const;

    // Destructors.
    virtual ~FlameSolver(void); // Default destructor.

	// Set the end conditions flag
	void SetEndConditions(bool val);

    // PROFILE INPUT.

    // Reads a flame gas-phase profile from a TAB formatted file.
    void LoadGasProfile(
        const std::string &file, // File name.
        Mops::Mechanism &mech    // Mechanism will be update with species.
        );

    // SOLUTION AND POST-PROCESSING.

    // Performs stochastic stepping algorithm up to specified stop time using
    // the given mechanism to define the stochastic processes.  Updates given
    // system accordingly.  On error returns <0, otherwise returns 0.  In this
    // flavour the gas-phase chemistry is interpolated from a vector of
    // IdealGas objects rather than being taken from the given system object.
    // However, the particles in the system object are updated accordingly.
    void Solve(
        Mops::Reactor &r,           // The reactor to solve.
        double tstop,                 // The end time for the step.
        int nsteps,                 // Number of internal steps to take.
        int niter,                  // Number of internal iterations to take.
        rng_type &rng,              // Random number generator
        Mops::Solver::OutFnPtr out, // Output function pointer.
        void *data                  // Custom data object which will be passed as argument to out().
        );

    GasProfile* Gasphase(void);

protected:

    //* The gas-phase chemistry profile.
    GasProfile m_gas_prof;

	//! Stagnation flame correction flag
	bool m_stagnation;

	//! Flag to impose the gas-phase conditions at the end of the time step over the entire time step
	bool m_endconditions;

    // HELPER FUNCTIONS.

    // Uses linear interpolation to return the chemical conditions
    // at a given time using a profile of Idealgas objects.
    double linInterpGas(
        double t,                      // Time.
        Sprog::Thermo::IdealGas &gas // Output gas conditions.
        ) const;
};
}

#endif
