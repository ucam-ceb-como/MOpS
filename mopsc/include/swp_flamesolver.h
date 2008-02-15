/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The Solver class holds simulation settings for mops and solves reactors.
    This basic solver only solves gas-phase chemistry equations, there is no
    operator splitting to solve gas-phase chemistry coupled to the particle
    population balance using sweep.
*/

#ifndef SWEEP_FLAME_SOLVER_H
#define SWEEP_FLAME_SOLVER_H

#include "mops_solver.h"
#include "sweep.h"
#include <map>

namespace Sweep
{
class FlameSolver : public Mops::Solver, public Sweep::Solver
{
public:
    // A map of Time/Gas-Phase pairs which describes a gas-phase
    // chemistry profile.
    typedef std::map<real, Sprog::Thermo::IdealGas> GasProfile;
    typedef std::pair<real, Sprog::Thermo::IdealGas> GasPoint;

    // Constructors.
    FlameSolver(void); // Default constructor.

    // Destructors.
    ~FlameSolver(void); // Default destructor.


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
    int Run(
        real &t,        // Simulation start time.  Will return the stop time.
        real tstop,     // Stop time for simulation.
        const GasProfile &gasphase, // Gas-phase profile.
        Cell &sys,      // System to solve.
        const Mechanism &mech // Mechanism to use to solve system.
        );

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

private:
    static const real CONFA;

    // The gas-phase chemistry profile.
    GasProfile m_gasprof;

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
        const Sweep::Mechanism &mech, // Mechanism which defines the output.
        const Mops::timevector &times // Vector of time intervals.
        );

    // Sets up file output for a new run given the run number.
    void beginRunFileOutput(unsigned int run);

    // Writes the current system state to the output file.
    void fileOutput(real time, const Sweep::Cell &sys);

    // Ends file output by closing all open files.
    void endFileOutput();


    // POST-PROCESSING ROUTINES.

    void ppAux(
        const std::string &filename,
        Mechanism &mech, 
        Mops::timevector &times
        ) const;

    void buildOutputVector(
        const Sweep::Cell &sys, // Reactor object holding data to be output.
        fvector &out            // Vector to build.
        ) const;

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
