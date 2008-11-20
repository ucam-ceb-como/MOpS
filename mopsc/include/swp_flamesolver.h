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
/*
    // A map of Time/Gas-Phase pairs which describes a gas-phase
    // chemistry profile.
    struct GasPoint {
        real Time;
        Sprog::Thermo::IdealGas Gas;
        GasPoint(const Sprog::SpeciesPtrVector &species);
    };
    typedef std::vector<GasPoint> GasProfile;
    typedef std::map<real, Sprog::Thermo::IdealGas> GasProfile;
    typedef std::pair<real, Sprog::Thermo::IdealGas> GasPoint;
*/

    // Constructors.
    FlameSolver(void); // Default constructor.

    // Destructors.
    virtual ~FlameSolver(void); // Default destructor.


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
        real tstop,                 // The end time for the step.
        int nsteps,                 // Number of internal steps to take.
        int niter,                  // Number of internal iterations to take.
        Mops::Solver::OutFnPtr out, // Output function pointer.
        void *data                  // Custom data object which will be passed as argument to out().
        );

/*
    // Run the solver for the given reactor and the 
    // given time intervals.
    void SolveReactor(
        Mops::Reactor &r,              // Reactor object to solve.
        const Mops::timevector &times, // Vector of time intervals.
        unsigned int nruns = 1         // Number of runs to perform.
        );

    // Post-processes binary output files with the given file name
    // into CSV files.
    void PostProcess(
        const std::string &filename, // Output file name root (no extension).
        unsigned int nruns = 1       // Number of runs to post-process.
        ) const;
*/

private:
    /*
    static const real CONFA;
*/
    // The gas-phase chemistry profile.
    GasProfile m_gasprof;
/*
    // Stats output.
    EnsembleStats *m_stats;

    // CONSOLE OUTPUT.

    // Sets up console output using the given mechanism as a template.
    void setupConsole(const Sweep::Mechanism &mech);

    // Writes current reactor state to the console.
    void consoleOutput(real time, const Sweep::Cell &sys) const;
    // FILE OUTPUT.

    // Sets up the file output by outputting an auxilliary file
    // which stores all the information required to post-process the
    // simulation.
    void beginFileOutput(
        const Mops::Mechanism &mech, // Mechanism which defines the output.
        const Mops::timevector &times // Vector of time intervals.
        );

    // Sets up file output for a new run given the run number.
    void beginRunFileOutput(unsigned int run);

    // Writes the current system state to the output file.
    void fileOutput(real time, const Sweep::Cell &sys);

    // Ends file output by closing all open files.
    void endFileOutput();
*/

    // HELPER FUNCTIONS.

    // Uses linear interpolation to return the chemical conditions
    // at a given time using a profile of Idealgas objects.
    void linInterpGas(
        real t,                      // Time.
        const GasProfile &gasphase,  // Gas-phase profile.
        Sprog::Thermo::IdealGas &gas // Output gas conditions.
        ) const;
};
};

#endif
