/*!
 * \file   simulator.cpp
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
#include "simulator.h"

#include "mops_timeinterval.h"
#include "rng.h"
#include "swp_model_stats.h"

#include "reactor1d.h"
#include "pred_corr_solver.h"

#include <sstream>
#include <cassert>
#include <vector>

using namespace Brush;

const size_t Brush::Simulator::sFirstSeed = 123;

/*!
 *@param[in]    n_paths                     Number of independent paths to simulation
 *@param[in]    n_corrector_iterations      Number of corrector iterations per step
 *@param[in]    output_times                Times at which to save output
 *@param[in]    initial_reactor             Initial condition
 *@param[in]    reset_chem                  Object to specify chemical conditions
 *@param[in]    output_file                 Base output file name
 *@param[in]    split_diffusion             Activate split simulation of diffusion
 *@param[in]    split_advection             Activate split simulation of advection
 */
Brush::Simulator::Simulator(const size_t n_paths,
                            const size_t n_corrector_iterations,
                            const Mops::timevector &output_times,
                            const Reactor1d &initial_reactor,
                            const ResetChemistry &reset_chem,
                            const std::string& output_file,
                            const bool split_diffusion,
                            const bool split_advection)
        : mPaths(n_paths)
        , mCorrectorIterations(n_corrector_iterations)
        , mRtol(0.0)
        , mAtol(0.0)
        , mSplitDiffusion(split_diffusion)
        , mSplitAdvection(split_advection)
        , mOutputTimeSteps(output_times)
        , mInitialReactor(initial_reactor)
        , mResetChemistry(reset_chem)
        , mOutputFile(output_file)
{
    for(size_t i = 0; i != initial_reactor.getNumCells(); ++i) {
        assert(initial_reactor.getCell(i).Mixture() != NULL);
        assert(mInitialReactor.getCell(i).Mixture() != NULL);
    }
}


/*!
 *@param[in]    seed_offset     Offset from standard seed for random number generation.
 *
 * Run the simulation paths
 */
void Brush::Simulator::runSimulation(const size_t seed_offset) {
    for(size_t r = 0; r < mPaths; ++r) {
        std::cout << "starting path " << r << '\n';
        // Use sFirstSeed as an offset for the seed to keep it > 0
        runOnePath(sFirstSeed + seed_offset + r);
    }

}

/*!
 * Run one path starting at the specified seed for the random
 * number generator.  This is intended to help distribution
 * of paths across nodes.
 *
 *\param[in]    seed    Seed to set on random number generator
 */
void Brush::Simulator::runOnePath(const int seed) {
    // reset the random number generator
    Sweep::srnd(seed);

    // Create a reactor instance to be modified by the solver
    Reactor1d reac(mInitialReactor);

    // It might be better to change the solver from a class
    // to a namespace
    PredCorrSolver solver(mResetChemistry, mCorrectorIterations, mRtol, mAtol,
                          mSplitDiffusion, mSplitAdvection);

    //==================== File to store the moments for this run
    std::ofstream momentsFile(buildParticleStatsFileName(seed).c_str());

    // Write the column headings to the file, with commas between names, but not after the last name
    momentsFile << "t,x";
    //BOOST_FOREACH(const std::string &colName, Sweep::Stats::EnsembleStats(mInitialReactor.getMechanism().ParticleMech()).Names())
    typedef std::vector<std::string> string_vector;
    {
        const string_vector &colNames = Sweep::Stats::EnsembleStats(mInitialReactor.getMechanism().ParticleMech()).Names();
        const string_vector::const_iterator itEnd = colNames.end();
        for(string_vector::const_iterator it = colNames.begin(); it != itEnd; ++it) {
            momentsFile << ',' << *it;
        }
    }
    momentsFile << '\n';

    //===================== File to store the lists of particles for this run
    std::ofstream particleListFile(buildParticleListFileName(seed).c_str());

    // Output column headings
    // Create a local scope to hide the vector of heading names
    {
        //Build a vector of column headings
        string_vector particleListHeadings;
        particleListHeadings.push_back("t");
        particleListHeadings.push_back("x");

        // Append the remaining (model defined) headings
        Sweep::Stats::EnsembleStats(mInitialReactor.getMechanism().ParticleMech()).PSL_Names(particleListHeadings, 2);

        // Write the column headings to the file, with commas between names, but not after the last name
        string_vector::const_iterator it = particleListHeadings.begin();
        const string_vector::const_iterator itStop = particleListHeadings.end() - 1;
        while(it != itStop) {
            particleListFile << *it++ << ',';
        }

        // Final heading
        particleListFile << particleListHeadings.back() << '\n';
    }

    //==================== File in which to log the process rates for this run
    std::ofstream ratesFile(buildParticleRatesFileName(seed).c_str());

    // Write the column headings to the file, with commas between names, but not after the last name
    ratesFile << "t,x";

    // Create a local scope to hide vector of process names
    {
        string_vector processNames;
        mInitialReactor.getMechanism().ParticleMech().GetProcessNames(processNames);
        string_vector::const_iterator it = processNames.begin();
        const string_vector::const_iterator itEnd = processNames.end();
        while(it != itEnd) {
            ratesFile << ',' << *it++;
        }
    }
    ratesFile << '\n';


    // write initial moments to file
    saveParticleStats(reac, momentsFile);

    // write initial particle list to file
    saveParticleList(reac, particleListFile);

    // log initial process rates
    saveProcessRates(reac, ratesFile);

    // Now loop over the user specified time steps
    for(Mops::timevector::const_iterator it = mOutputTimeSteps.begin(); it != mOutputTimeSteps.end(); ++it) {
        const Mops::TimeInterval &tInt = *it;
        const real dt = tInt.StepSize();

        std::cout << "stepping from " << tInt.StartTime() << ' ' << reac.getTime() << '\n';

        for(unsigned int i = 0; i != tInt.StepCount(); ++i) {
            //std::cout << "solving from " << reac.getTime() << '\n';

            const real dt2 = dt / tInt.SubSplittingStepCount();
            for(unsigned int j = 0; j != tInt.SubSplittingStepCount(); ++j) {
                solver.solve(reac, tInt.StartTime() + i * dt + (j + 1) * dt2,
                             tInt.SplittingStepCount(), mCorrectorIterations);
            }

            //std::cout << "solved up to " << tInt.StartTime() + i * dt  + (j + 1) * dt2 << ' ' << reac.getTime() << '\n';

            // write moments to file
            saveParticleStats(reac, momentsFile);

            // log process rates
            saveProcessRates(reac, ratesFile);
        }

        //std::cout << "stepped up to " << tInt.EndTime() << ' ' << reac.getTime() << '\n';

        // write particle list to file
        saveParticleList(reac, particleListFile);
    }
    std::cout << "finished path\n";
}

/*!
 * Calculate and serialise statistics for each cell in the reactor
 *
 *\param[in]    reac        Reactor for which statistics are to be calculated
 *\param[in]    out         File handle into which to write the moment data
 */
void Brush::Simulator::saveParticleStats(const Reactor1d &reac, std::ostream &out) {
    // Create a stats object to now so there are not a lot of string operations
    // each time one is needed below
    Sweep::Stats::EnsembleStats stats(mInitialReactor.getMechanism().ParticleMech());

    // Reject particles with a collision diameter outside the range (0, 1.0e30)
    Sweep::Stats::IModelStats::StatBound statsBound;
    statsBound.Lower = 0;
    statsBound.Upper = 1.0e30;
    statsBound.PID = Sweep::ParticleCache::iDcol;
    stats.SetStatBoundary(statsBound);

    for(size_t i = 0; i < reac.getNumCells(); ++i) {
        // Collect the statistics
        stats.Calculate(reac.getCell(i).Mixture()->Particles(),
                        1.0 / reac.getCell(i).Mixture()->SampleVolume(),
                        1.0 / reac.getCell(i).Mixture()->SecondarySampleVolume());

        // Output the time and place to which the statistics apply
        out << reac.getTime() << ',' << reac.getCellCentre(i);

        // Put the stats data into the file
        //BOOST_FOREACH(real r, stats.Get())
        {
            const fvector &statsVector = stats.Get();
            const fvector::const_iterator itEnd = statsVector.end();
            for(fvector::const_iterator it = statsVector.begin(); it != itEnd; ++it) {
                out << ',' << *it;
            }
        }
        out << std::endl;
    }

    // Flush the buffer, so that a future crash does not lose data before it actually
    // reaches its destination.
    out.flush();
}

/*!
 * Calculate and serialise a list of particles from the cells in the reactor
 *
 *\param[in]    reac        Reactor for which statistics are to be calculated
 *\param[in]    out         File handle into which to write the moment data
 */
void Brush::Simulator::saveParticleList(const Reactor1d &reac, std::ostream &out) {
    // Create a stats object to now so there are not a lot of string operations
    // each time one is needed below
    Sweep::Stats::EnsembleStats stats(mInitialReactor.getMechanism().ParticleMech());

    // Loop over the cells in the reactor
    for(size_t i = 0; i < reac.getNumCells(); ++i) {
        const Mops::Mixture &mix = *reac.getCell(i).Mixture();

        // Loop over the main particles
        for(unsigned int particleIndex = 0; particleIndex < mix.ParticleCount(); ++particleIndex) {
            // Get a vector containing the particle details
            fvector particleListEntry;
            stats.PSL(*mix.Particles().At(particleIndex), reac.getMechanism().ParticleMech(),
                      reac.getTime(), particleListEntry, 1.0 / mix.SampleVolume());

            // Output the time and place at which the particle is found
            out << reac.getTime() << ',' << mix.Particles().At(particleIndex)->getPosition();

            // Output the particle details
            for(fvector::const_iterator it = particleListEntry.begin();
                it != particleListEntry.end(); ++it) {
                out << ',' << *it;
            }

            // New line ready for next particle
            out << '\n';
        }

        // Loop over the secondary particles
        for(unsigned int particleIndex = 0; particleIndex < mix.Particles().SecondaryCount(); ++particleIndex) {
            // Get a vector containing the particle details
            fvector particleListEntry;
            stats.PSL(*mix.Particles().SecondaryParticleAt(particleIndex), reac.getMechanism().ParticleMech(),
                      reac.getTime(), particleListEntry, 1.0 / mix.SecondarySampleVolume());

            // Output the time and place at which the particle is found
            out << reac.getTime() << ',' << mix.Particles().SecondaryParticleAt(particleIndex)->getPosition();

            // Output the particle details
            for(fvector::const_iterator it = particleListEntry.begin();
                it != particleListEntry.end(); ++it) {
                out << ',' << *it;
            }

            // New line ready for next particle
            out << '\n';
        }
    }

    // Flush the buffer, so that a future crash does not lose data before it actually
    // reaches its destination.
    out.flush();
}

/*!
 * Calculate and log process rates for each cell in the reactor
 *
 *\param[in]    reac        Reactor for which rates are to be calculated
 *\param[in]    out         File handle into which to write the data
 */
void Brush::Simulator::saveProcessRates(const Reactor1d &reac, std::ostream &out) {
    for(size_t i = 0; i < reac.getNumCells(); ++i) {
        // Collect the rates
        fvector rates;
        const Geometry::LocalGeometry1d geom(reac.getGeometry(), i);
        reac.getMechanism().ParticleMech().CalcRateTerms(reac.getTime(),
                                                         *reac.getCell(i).Mixture(),
                                                         geom, rates);

        // Output the time and place to which the rates apply
        out << reac.getTime() << ',' << reac.getCellCentre(i);

        // Put the data into the file
        {
            fvector::const_iterator it = rates.begin();
            const fvector::const_iterator itEnd = rates.end();
            while(it != itEnd) {
                out << ',' << *it++;
            }
        }
        out << '\n';
    }

    // Flush the buffer, so that a future crash does not lose data before it actually
    // reaches its destination.
    out.flush();
}

/*!
 * Construct the filename for the moments file for a run with the
 * specified seed.  This function takes a seed as an argument to match
 * runOnePath().
 *
 *\param[in]    seed        Seed used for the path
 *
 *\return       Name for moments file
 */
std::string Brush::Simulator::buildParticleStatsFileName(const int seed) const {
    std::stringstream name;
    name << mOutputFile << seed << "_partstats.csv";
    return name.str();
}

/*!
 * Construct the filename for the psl file for a run with the
 * specified seed.  This function takes a seed as an argument to match
 * runOnePath().
 *
 *\param[in]    seed        Seed used for the path
 *
 *\return       Name for psl file
 */
std::string Brush::Simulator::buildParticleListFileName(const int seed) const {
    std::stringstream name;
    name << mOutputFile << seed << "_psl.csv";
    return name.str();
}

/*!
 * Construct the filename for particle rates log file for a run with the
 * specified seed.  This function takes a seed as an argument to match
 * runOnePath().
 *
 *\param[in]    seed        Seed used for the path
 *
 *\return       Name for particle rates log file
 */
std::string Brush::Simulator::buildParticleRatesFileName(const int seed) const {
    std::stringstream name;
    name << mOutputFile << seed << "_prates.csv";
    return name.str();
}

