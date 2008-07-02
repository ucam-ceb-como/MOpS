/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ParticleSolver class is an extension of the Solver class which
    also contains routines common to all mops solvers which run Sweep.

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

    
    // PARTICLE TRACKING.

    // Number of particles for which to produce tracking output.  Tracked
    // particles will have their PSL properties output at each time-step in
    // a separate CSV file.  If a primary-particle model in implemented, then
    // TEM-style images will also be generated.
    unsigned int m_ptrack_count;


    // FILE OUTPUT.

    // Writes the particle stats to the binary output file.
    void outputParticleStats(const Reactor &r) const;

    // Writes tracked particles to the binary output file.
    void outputPartTrack(const Reactor &r) const;

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

    // Reads the tracked particles from the binary file.  The particles are
    // processed so that only a vector of vectors is returned, which contains
    // the PSL data for each tracked particle at that point.
    void readPartTrackPoint(
        std::istream &in,             // Input stream.
        const Sweep::Mechanism &mech, // Particle mechanism.
        std::vector<fvector> &pdata   // Tracked particle output data for this point.
        ) const;

    // Writes particle stats profile to a CSV file.
    static void writeParticleStatsCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining particle ensemble.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of gas-phase time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );

    // Writes particle tracking for multiple particles to CSV files.
    static void writePartTrackCSV(
        const std::string &filename, // Output file name (excl. extension).
        const Mechanism &mech,       // Mechanism defining particle ensemble.
        const timevector &times,     // Output time profile.
        std::vector<std::vector<fvector> > &track // Vector of tracking data of multiple particles at multiple time points.
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
