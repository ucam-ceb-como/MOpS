/*!
 * \file   simulator.h
 * \author Robert I A Patterson
 *
 * \brief  Top level simulation driver
 *
 *  Copyright (C) 2009 Robert I A Patterson.
 *

 Licence:
    This file is part of "brush".

    brush is free software; you can redistribute it and/or
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
    Prof Markus Kraft
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
#ifndef BRUSH_SIMULATOR_H
#define BRUSH_SIMULATOR_H

#include "brush_params.h"

#include "reactor1d.h"
#include "reset_chemistry.h"
#include "mops_timeinterval.h"

#include <string>

namespace Brush {
    class ResetChemistry;

//! Specify and run a multiple path simulation
class Simulator {
public:
    //! Construct from individual components
    Simulator(const size_t n_paths, const size_t n_corrector_iterations,
              const Mops::timevector &output_times,
              const Reactor1d &initial_reactor,
              const ResetChemistry &reset_chem,
              const std::string output_file)
        : mPaths(n_paths)
        , mCorrectorIterations(n_corrector_iterations)
        , mRtol(0.0)
        , mAtol(0.0)
        , mOutputTimeSteps(output_times)
        , mInitialReactor(initial_reactor)
        , mResetChemistry(reset_chem)
        , mOutputFile(output_file)
    {}

    //! Run the simulation paths and store output
    void runSimulation(const size_t seed_offset);

protected:
    //! Run a single Monte Carlo path
    void runOnePath(const int seed);

    //! Calculate and save statistics in a text file
    void saveParticleStats(const Reactor1d &reac, std::ostream &out);

    //! Save list of particles in a text file
    void saveParticleList(const Reactor1d &reac, std::ostream &out);

    //! Save process rates in a text file
    void saveProcessRates(const Reactor1d &reac, std::ostream &out);

    //! Build the name of the moments file for a path
    std::string buildParticleStatsFileName(const int seed) const;

    //! Build the name of the particle list file for a path
    std::string buildParticleListFileName(const int seed) const;

    //! Build the name of the particle rates log file for a path
    std::string buildParticleRatesFileName(const int seed) const;

private:
    //! Number of paths to run
    size_t mPaths;

    //! Number of corrector iterations to perform per fundamental time step
    size_t mCorrectorIterations;

    //! Relative error tolerance for ODE solver for chemistry
    real mRtol;

    //! Absolute error tolerance for ODE solver for chemistry
    real mAtol;

    //! Use specified output steps
    Mops::timevector mOutputTimeSteps;

    //! Initial condition
    Reactor1d mInitialReactor;

    //! Object for setting reactor chemistry to fixed values
    ResetChemistry mResetChemistry;

    //! Base string for output file names
    std::string mOutputFile;

    //! Default simulator is meaningless
    Simulator();

    //! Starting point for the random number generator seed
    static const size_t sFirstSeed;
};

} //namespace Brush

#endif //BRUSH_SIMULATOR_H
