/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The StrangSolver class holds solves gas-phase chemistry equations coupled
    to the particle population balance (using Sweep) and a Strang operator
    splitting technique (see PP39).
*/

#ifndef MOPS_STRANG_SOLVER_H
#define MOPS_STRANG_SOLVER_H

#include "mops_params.h"
#include "mops_particle_solver.h"
#include "mops_reactor.h"
#include "mops_timeinterval.h"
#include "mops_mechanism.h"
#include "sweep.h"
#include <string>

namespace Mops
{
class StrangSolver : public ParticleSolver, Sweep::Solver
{
public:
    // Constructors.
    StrangSolver(void); // Default constructor.

    // Destructors.
    ~StrangSolver(void); // Default destructor.


    // SOLUTION AND POST-PROCESSING.

    // Run the solver for the given reactor and the 
    // given time intervals.
    void SolveReactor(
        Reactor &r,              // Reactor object to solve.
        const timevector &times, // Vector of time intervals.
        unsigned int nruns = 1   // Number of runs to perform.
        );

    // Post-processes binary output files with the given file name
    // into CSV files.
    void PostProcess(
        const std::string &filename, // Filename to post-process.
        unsigned int nruns = 1       // Number of runs.
        ) const;

private:
    // SIMULATION.

    // Performs n Strang splitting steps.
    void multiStrangStep(
        real dt,        // Splitting step size.
        unsigned int n, // Number of splitting steps.
        Reactor &r      // Reactor to solve.
        );
};
};

#endif
