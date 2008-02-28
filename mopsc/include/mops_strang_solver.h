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
#include "mops_solver.h"
#include "mops_reactor.h"
#include "mops_timeinterval.h"
#include "mops_mechanism.h"
#include "sweep.h"
#include <string>

namespace Mops
{
class StrangSolver : public Solver, Sweep::Solver
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

protected:
    // COMPUTATION TIME.

    double m_swptime;

private:
    // Stats output.
    Sweep::EnsembleStats *m_stats;

    // CONSOLE OUTPUT.

    // Sets up console output using the given mechanism as a template.
    void setupConsole(const Mechanism &mech);

    // Writes current reactor state to the console.
    void consoleOutput(const Reactor &r) const;


    // FILE OUTPUT.

    // Sets up the file output by outputting an auxilliary file
    // which stores all the information required to post-process the
    // simulation and by opening the output file.
    void beginFileOutput(
        const Mechanism &mech,  // Mechanism which defines the output.
        const timevector &times // Vector of time intervals.
        );

    // Sets up file output for a new run given the run number.
    void beginRunFileOutput(unsigned int run);

    // Writes the current reactor state to the output file.
    void fileOutput(const Reactor &r);

    // Ends file output by closing all open files.
    void endFileOutput();


    // SAVE POINTS AND PSL POST-PROCESSING.

    // Creates a simulation save point.  The save points can be
    // used to restart an interrupted simulation, but are primarily
    // used as output points for the particle size distributions.
    void createSavePoint(
        const Reactor &r,  // Reactor to output.
        unsigned int step, // Step number.
        unsigned int run   // Run number.
        ) const;

    // Reads a save point file.
    Reactor *const readSavePoint(
        unsigned int step,    // Step number.
        unsigned int run,     // Run number.
        const Mechanism &mech // Mechanism used to define reactor.
        ) const;

    // Processes the PSLs.
    void postProcessPSLs(
        unsigned int nruns,     // Number of runs to post-process.
        const Mechanism &mech,  // Mechanism use to solve system.
        const timevector &times // Simulation output time intervals.
        ) const;


    // SIMULATION.

    void multiStrangStep(
        real dt,        // Splitting step size.
        unsigned int n, // Number of splitting steps.
        Reactor &r      // Reactor to solve.
        );


    // POST-PROCESSING ROUTINES.

    // Reads a chemistry data point.
    static void readChemDataPoint(
        std::istream &in, // Input stream.
        const Mops::Mechanism &mech, // Chemical mechanism.
        fvector &sum,      // Sums of chemistry data.
        fvector &sumsqr,   // Sums of the squares.
        bool calcsqrs = false // Set =true to also calculate sums of squares.
        );

    // Reads a particle stats data point.
    static void readStatDataPoint(
        std::istream &in, // Input stream.
        const Sweep::Mechanism &mech, // Particle mechanism.
        fvector &sum,      // Sums of stats data.
        fvector &sumsqr,   // Sums of the squares.
        bool calcsqrs = false // Set =true to also calculate sums of squares.
        );

    // Reads a CPU timing data point.
    static void readCPUDataPoint(
        std::istream &in, // Input stream.
        fvector &sum,      // Sums of CPU timing data.
        fvector &sumsqr,   // Sums of the squares.
        bool calcsqrs = false // Set =true to also calculate sums of squares.
        );

    // Builds output vector of particle stats.
    static void buildOutputVector(
        unsigned int step, // Step number.
        real time,         // Step time.
        fvector &avg,      // Averages.
        const fvector &err // Confidence intervals (will be inserted into avg).
        );
};
};

#endif
