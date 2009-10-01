/*!
 * \file   pred_corr_solver.cpp
 * \author Robert I A Patterson
 *
 * \brief  Routines for solving with predictor corrector coupling
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
#include "pred_corr_solver.h"

#include "reactor1d.h"

#include "rng.h"
#include "swp_transport_outflow.h"

#include "choose_index.hpp"

#include <stdexcept>
#include <limits>
#include <algorithm>

using namespace Brush;

/*!
 * Advance the solution to t_stop by means of n_steps equal length application
 * of the predictor corrector algorithm.
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_stop      Time to which to advance solution
 *\param[in]            n_steps     Number of predictor corector steps
 *\param[in]            n_iter      Number of corrector iterations per step
 */
void Brush::PredCorrSolver::solve(Reactor1d &reac, const real t_stop, const int n_steps, const int n_iter) const {
    const real startTime = reac.getTime();
    const real dt = (t_stop - startTime) / n_steps;

    for(int i = 1; i <= n_steps; ++i) {
        predictorCorrectorStep(reac, startTime + i * dt, n_iter);
    }
}

/*!
 * Advance the solution to t_stop by means of a predictor step followed
 * by zero or more corrector iterations
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_stop      Time to which to advance solution
 *\param[in]            n_iter      Number of corrector iterations per step
 */
void Brush::PredCorrSolver::predictorCorrectorStep(Reactor1d &reac, const real t_stop, const int n_iter) const {
    const real startTime = reac.getTime();

    //=========== Predictor step =====================
    // Advance the chemistry
    solveChemistry(reac, t_stop);

    // reset the time and advance the particles
    reac.setTime(startTime);
    solveParticles(reac, t_stop);

    //=========== Corrector steps ====================
    for(int i = 0; i < n_iter; ++i) {
        ;
    }
}


/*!
 * Advance the chemistry part of the solution without performing any particle
 * events.  In this implementation the stop time is ignored and the chemistry
 * is interpolated from the reset object for the reactor time.
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_stop      Time to which to advance solution
 */
void Brush::PredCorrSolver::solveChemistry(Reactor1d &reac, const real t_stop) const {
    reac.ReplaceChemistry(mResetChemistry, true);
}

/*!
 * Advance the particle part of the solution, which may also affect the chemical
 * species concentrations.
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_stop      Time to which to advance solution
 */
void Brush::PredCorrSolver::solveParticles(Reactor1d &reac, const real t_stop) const {
    

/*    while(reac.getTime() < t_stop) {
            const real dt = particleTimeStep(reac, t_stop);

            if(dt >= 0)
                reac.setTime(reac.getTime() + dt);
            else
                throw std::runtime_error("Failed to take time step");
        }

    // Perform Linear Process Deferment Algorithm to
    // update all deferred processes.
        if (reac.getMechanism().ParticleMech().AnyDeferred()) {
        size_t numCells = reac.getNumCells();
            for(size_t i = 0; i < numCells; ++i) {
            reac.getMechanism().ParticleMech().LPDA(t_stop, *(reac.getCell(i).Mixture()));
            }
        }
*/
    const size_t numCells = reac.getNumCells();
    const Sweep::Mechanism &mech = reac.getMechanism().ParticleMech();

    // Carry out a repeated sequence of jump process simulation (with LPDA updates
    // for particles involved in jumps) followed by LPDA updates for the full
    // population.
    do {
        // Add up the rate of all deferred process in the whole of the reactor,
        // the result is the expected number of simulated events per computational
        // particle per simulated second of physical time.
        //@todo This is a very clumsy average that ignores spatial inhomogeneity
        real deferredRate = 0.0;
        real jumpRate = 0.0;
        for(size_t i = 0; i < numCells; ++i) {
            // Ignore the individual rate term information in the final argument
            fvector dummyVec;
            // Add 1 to the denominator to avoid dividing by 0
            deferredRate += mech.CalcDeferredRateTerms(reac.getTime(), *(reac.getCell(i).Mixture()),
                                                       Geometry::LocalGeometry1d(reac.getGeometry(), i),
                                                       dummyVec)
                            / (reac.getCell(i).Mixture()->ParticleCount() + 1);

            jumpRate     += mech.CalcJumpRateTerms(    reac.getTime(), *(reac.getCell(i).Mixture()),
                                                       Geometry::LocalGeometry1d(reac.getGeometry(), i),
                                                       dummyVec)
                            / (reac.getCell(i).Mixture()->ParticleCount() + 1);
        }

        // Calculate the latest time at which all particles must be updated
        // with deferred events.
        real maxDeferralEnd = t_stop;
        if(deferredRate / jumpRate > mDeferralRatio) {
            // A large number of deferred events are expected between jump events
            // so force updates to take place more frequently by shortening the
            // time until the next update of all particles.
            maxDeferralEnd = std::min(t_stop, mDeferralRatio * deferredRate);
        }

        // Simulate jumps upto maxDeferralEnd
        // Put a very small tolerance on the comparison
        while(reac.getTime() < maxDeferralEnd * (1.0 - std::numeric_limits<real>::epsilon())) {
            // Perform one non-deferred event
            const real dt = particleTimeStep(reac, maxDeferralEnd);

            if(dt >= 0)
                reac.setTime(reac.getTime() + dt);
            else
                throw std::runtime_error("Failed to take time step");
        }

        // Perform all events from deferred processes.
        if (reac.getMechanism().ParticleMech().AnyDeferred()) {
            for(size_t i = 0; i < numCells; ++i) {
                reac.getMechanism().ParticleMech().LPDA(maxDeferralEnd, *(reac.getCell(i).Mixture()));
            }
        }
    } while(reac.getTime() < t_stop * (1.0 - std::numeric_limits<real>::epsilon()));

}

/*!
 * Carry out one jump of the stochastic particle processes, if it occurs before
 * the stop time, otherwise do nothing and just advance time.  Only one stochastic
 * jump event will occur, the active process and active cell are chosen according
 * to the probabilities of the Markov Chain, which are proportional to their
 * rates.
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_stop      Stop time
 *
 *\return       Length of step taken
 */
real Brush::PredCorrSolver::particleTimeStep(Reactor1d &reac, const real t_stop) const {
    //std::cout << "Start of time step " << reac.getTime() << ", max step " << t_stop << '\n';

    const size_t numCells = reac.getNumCells();

    // Each entry in this vector will be a vector of rates, one vector per cell.
    // Within the vector for each cell will be the jump rates, one for each process
    std::vector<fvector> jRates(numCells, fvector(reac.getMechanism().ParticleMech().ProcessCount()));

    // Put total jumps rates for each cell in this array.  This will mean that
    // cellRates[i] == sum(jRates[i]).
    fvector cellRates(numCells);

    // This will be the total rate of all processes in all cells
    real totRate = 0;

    // Get the rate data from the cells
    for(size_t i = 0; i < numCells; ++i) {
        // Get the geometry information for the cell
        Geometry::LocalGeometry1d cellGeom(reac.getGeometry(), i);

        cellRates[i] = reac.getMechanism().ParticleMech().CalcJumpRateTerms(reac.getTime(),
                                                                            *(reac.getCell(i).Mixture()),
                                                                            cellGeom, jRates[i]);
        totRate += cellRates[i];
    }

    // Calculate the random time step length
    real dt;
    if(totRate > 0) {
        // Exponential random variable with mean 1/totRate
        dt = Sweep::rnd();
        dt = -log(dt) / totRate;
    }
    else {
        // Nothing is happening so take a near infinite time step
        dt = 1e30;
    }

    // If the waiting time is too long truncate it, and say that nothing has
    // happened (memoryless property of exponential random variables) by leaving
    // the index of the cell in which to perform an event as -1 (which does not refer
    // to any cell).
    int activeCell = -1;
    if(dt > mMaxDt ) {
        dt = mMaxDt;
        if( reac.getTime() + dt > t_stop) {
            dt = t_stop - reac.getTime();
        }
    }
    else {
        // Check if the waiting time is longer than the time until the end of the
        // step.  In this case step up to the end of the time interval without doing
        // anything by leaving the cell index as -1 (see above)
        if( reac.getTime() + dt > t_stop) {
            dt = t_stop - reac.getTime();
        }
        else {
            // This step is an actual event, choose between the cells with weights
            // proportional to their weights
            activeCell = chooseIndex(cellRates, Sweep::rnd);
            //std::cout << "cell with active process is " << activeCell << '\n';
        }
    }

    int activeProcess = -1;
    if(activeCell >= 0) {
        // Work out which process is responsible for the event in the activeCell
        activeProcess = chooseIndex(jRates[activeCell], Sweep::rnd);
        //std::cout << "active process is " << activeProcess << '\n';

        // Save the initial chemical conditions in sys so that we
        // can restore them at the end of the run.
        //Sprog::Thermo::IdealGas chem = *r.Mixture();

        // Store if chemical conditions are fixed at present, because we
        // shall set them to be fixed during this run, to be restored afterwards.
        //bool fixedchem = r.Mixture()->FixedChem();
        //r.Mixture()->SetFixedChem();

        // Perform the event and find out if a particle was transported out of its cell
        Sweep::Transport::TransportOutflow out;

        // Get the geometry information for the cell
        Geometry::LocalGeometry1d cellGeom(reac.getGeometry(), activeCell);

        reac.getMechanism().ParticleMech().DoProcess(activeProcess, reac.getTime() + dt,
                                                     *(reac.getCell(activeCell).Mixture()), cellGeom, &out);

        if(out.particle) {
            if(out.destination >= 0) {
                // Particle has stayed in the system
                transportIn(reac, out.destination, out);
            }

            // The original copy of the transported particle is no longer needed
            delete out.particle;
            out.particle = NULL;
        }
    }

    //std::cout << ", actual step " << dt << " with process " << activeProcess << " in cell " << activeCell << std::endl;

    return dt;
}

/*!
 * Add a particle that has left one cell into the specified cell.  The particle_details
 * must contain a pointer to a valid particle, ownership of which remains with
 * the caller.
 *
 *\param[in,out]    reac                    1D reactor within which the transport is occuring
 *\param[in]        destination_index       Index of destination cell
 *\param[in]        particle_details        Information about the particle and its weight
 */
void Brush::PredCorrSolver::transportIn(Reactor1d & reac, const size_t destination_index,
                                        const Sweep::Transport::TransportOutflow &particle_details) const {
    real incomingWeight = particle_details.weight;

    // Need to match weights of particle between source and destination
    unsigned int safetyCounter = 0;
    while(true) {
        real destinationWeight = 1.0 / reac.getCell(destination_index).Mixture()->SampleVolume();

        if(incomingWeight >= destinationWeight) {
            // The important thing is that the statistical weight of particles added to the destination cell
            // is equal to the statistical weight of the particle that was removed from the originating cell.
            // Otherwise mass will be lost.

            // Insert one copy of the particle into the destination cell
            reac.getCell(destination_index).Mixture()->Particles().Add(*(new Sweep::Particle(*particle_details.particle)));
            incomingWeight -= destinationWeight;

            // Avoid infinite loops
            if(++safetyCounter > 100000) {
                throw std::runtime_error("Failed to match particle weights in PredCorrSolver::transportIn()");
            }
        }
        else
            break;
    }

    // Unfortunately we cannot quite conserve statistical weight, there will always be a bit
    // left over after the loop above.  This can only be handled in an average sense.
    if(Sweep::rnd() < incomingWeight * reac.getCell(destination_index).Mixture()->SampleVolume()) {
        reac.getCell(destination_index).Mixture()->Particles().Add(*(new Sweep::Particle(*particle_details.particle)));
    }
}

