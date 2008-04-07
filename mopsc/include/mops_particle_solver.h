/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The ParticleSolver class is an extension of the Solver class which
    also contains routines common to all mops solvers which run Sweep.
*/

#ifndef MOPS_PARTICLE_SOLVER_H
#define MOPS_PARTICLE_SOLVER_H

#include "mops_params.h"
#include "mops_reactor.h"
#include "mops_timeinterval.h"
#include "mops_mechanism.h"
#include "mops_solver.h"
#include "console_io.h"
#include <vector>
#include <string>
#include <fstream>
#include <time.h>

namespace Mops
{
class ParticleSolver : public Solver
{
public:
    // Constructors.
    ParticleSolver(void); // Default constructor.

    // Destructors.
    virtual ~ParticleSolver(void); // Default destructor.

protected:
    // COMPUTATION TIME.

    // Sweep computation time.
    double m_swp_ctime;


    // FILE OUTPUT.

    // Writes the particle stats to the binary output file.
    void outputParticleStats(const Reactor &r) const;

    // Writes the current reactor state to the output file.  This overrides
    // the default routine in order to also output the particle stats.
    void fileOutput(const Reactor &r) const;


    // POST-PROCESSING ROUTINES.

    // Reads a particle stats data point from the binary file.
    // To allow the averages and confidence intervals to be calculated
    // the data point is added to a vector of sums, and the squares are
    // added to the vector sumsqr if necessary.
    static void readParticleDataPoint(
        std::istream &in,             // Input stream.
        const Sweep::Mechanism &mech, // Chemical mechanism.
        fvector &sum,                 // Sums of chemistry data.
        fvector &sumsqr,              // Sums of the squares.
        bool calcsqrs = false         // Set =true to also calculate sums of squares.
        );

    // Writes particle stats profile to a CSV file.
    static void writeParticleStatsCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining particle ensemble.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of gas-phase time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );

    // Writes computation times profile to a CSV file.  This overrides
    // the default function in order to allow output of the Sweep CPU
    // time.
    void writeCT_CSV(
        const std::string &filename,     // Output file name (incl. extension).
        const timevector &times,         // Output time point profile.
        std::vector<fvector> &avg,       // Vector of computation-time time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        ) const;


    // SAVE POINTS AND PSL POST-PROCESSING.

    // Processes the PSLs at each save point into single files.
    void postProcessPSLs(
        unsigned int nruns,     // Number of runs to post-process.
        const Mechanism &mech,  // Mechanism use to solve system.
        const timevector &times // Simulation output time intervals.
        ) const;
};
};

#endif
