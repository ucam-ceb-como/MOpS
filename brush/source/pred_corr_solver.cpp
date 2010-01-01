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
#include "swp_cell.h"

#include "choose_index.hpp"

#include <stdexcept>
#include <limits>
#include <algorithm>
#include <vector>
#include <list>
#include <cassert>

using namespace Brush;

/*!
 * This method has a number of anonymous arguments, which are provided as part
 * of the design for a two-way gas-particle coupled solver.  The arguments are
 * anonymous, because their values are not currently used.  Names are given in
 * the header file and should be copied here when they are needed.
 *
 *@param[in]    reset_chem              Chemical species concentrations as function of spatial position
 *@param[in]                            Number of corrector iterations to perform during full coupling
 *@param[in]                            Relative tolerance for ODE solver used for gas phase
 *@param[in]                            Absolute tolerance for ODE solver used for gas phase
 *@param[in]    split_diffusion         True if diffusion is to be simulated via splitting
 *@param[in]    split_advection         True if advection is to be simulated via splitting
 *
 *@todo     Need generalisation to move away from using ResetChemistry
 */
Brush::PredCorrSolver::PredCorrSolver(const ResetChemistry& reset_chem,
                                      const size_t,
                                      const real , const real ,
                                      const bool split_diffusion,
                                      const bool split_advection)
    : mResetChemistry(reset_chem)
    , mMaxDt(3.0e-4)
    , mDeferralRatio(10.0)
    , mSplitDiffusion(split_diffusion)
    , mSplitAdvection(split_advection)
{}

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

    for(size_t i = 0; i < reac.getNumCells(); ++i) {
        assert(reac.getCell(i).Mixture() != NULL);
    }

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

    // Expansion and contraction of the gas phase changes the particle
    // concentration, since the same mass of gas still contains the same
    // number of soot particles.  The change in volume containing a fixed mass
    // of gas is inversely proportional to the change in gas mass density.
    const size_t numCells = reac.getNumCells();
    fvector oldDensity(numCells);
    for(size_t i = 0; i != numCells; ++i) {
        oldDensity[i] = reac.getCell(i).Mixture()->MassDensity();
    }

    // Update the chemistry to the new time
    reac.ReplaceChemistry(mResetChemistry, true);

    // Update the soot particle concentration to account for the expansion /
    // contraction of the gas.
    for(size_t i = 0; i != numCells; ++i) {
        Sweep::Cell * const pMix = reac.getCell(i).Mixture();

        // Collect data for calculating the new value
        const real newDensity = pMix->MassDensity();
        const real oldM0 = pMix->ParticleCount() / pMix->SampleVolume();

        pMix->SetM0(oldM0 * newDensity / oldDensity[i]);
    }
}

/*!
 * Advance the particle part of the solution, which may also affect the chemical
 * species concentrations.
 *
 * Note that the flags activating split simulation of particle transport do not
 * disable any diffusion processes that may have been specified in the particle
 * mechanism.
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_stop      Time to which to advance solution
 */
void Brush::PredCorrSolver::solveParticles(Reactor1d &reac, const real t_stop) const {
    
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

        // Cache to avoid recalculating jump rates for all cells, when only
        // very few are changed in each time step.
        JumpRateCache rateCache(reac);

        // Simulate jumps upto maxDeferralEnd
        // Put a very small tolerance on the comparison
        while(reac.getTime() <= maxDeferralEnd * (1.0 - std::numeric_limits<real>::epsilon())) {
            // Perform one non-deferred event
            const real dt = particleTimeStep(reac, maxDeferralEnd, rateCache);

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
    } while(reac.getTime() <= t_stop * (1.0 - std::numeric_limits<real>::epsilon()));

    // Now do the split particle transport, if there is any
    if(mSplitDiffusion || mSplitAdvection) {
        splitParticleTransport(reac, t_stop);
    }
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
real Brush::PredCorrSolver::particleTimeStep(Reactor1d &reac, const real t_stop,
                                             JumpRateCache &rate_cache) const {
    //std::cout << "Start of time step " << reac.getTime() << ", max step " << t_stop << '\n';

     // This will be the total rate of all processes in all cells
    rate_cache.update();
    real totRate = rate_cache.totalRate();

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
            activeCell = chooseIndex(rate_cache.mCellRates, Sweep::rnd);
            //std::cout << "cell with active process is " << activeCell << '\n';
        }
    }

    int activeProcess = -1;
    if(activeCell >= 0) {
        // Work out which process is responsible for the event in the activeCell
        activeProcess = chooseIndex(rate_cache.mProcessRates[activeCell], Sweep::rnd);
        //std::cout << "active process is " << activeProcess << '\n';

        // Perform the event and find out if a particle was transported out of its cell
        Sweep::Transport::TransportOutflow out;

        // Get the geometry information for the cell
        Geometry::LocalGeometry1d cellGeom(reac.getGeometry(), activeCell);

        reac.getMechanism().ParticleMech().DoProcess(activeProcess, reac.getTime() + dt,
                                                     *(reac.getCell(activeCell).Mixture()), cellGeom, &out);

        // Contents of active cell has changed so mark it for update in the rate cache
        rate_cache.mInvalidCells.push(activeCell);

        if(out.particle) {
            if(out.destination >= 0) {
                // Particle has stayed in the system
                transportIn(reac, out.destination, out);

                // Contents of destination cell has changed so mark it for update in the rate cache
                rate_cache.mInvalidCells.push(out.destination);
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

/*!
 *
 *
 *@param[in]    reac    Reactor in which particles are being transported
 *@param[in]    t_stop  Time upto which transport is to be simulated
 */
void Brush::PredCorrSolver::splitParticleTransport(Reactor1d &reac, const real t_stop) const {
    const size_t numCells = reac.getNumCells();

    // Lists of particles due to be added to cells and their statistical weights
    typedef std::vector<std::list<Sweep::Transport::TransportOutflow> > inflow_lists_vector;
    inflow_lists_vector inflowLists(numCells);

    // Loop over each cell updating particles positions and removing any
    // particles that are moving to new cells.
    for(unsigned int i = 0; i != numCells; ++i) {
        Mops::Mixture &mix = *(reac.getCell(i).Mixture());
        Sweep::Ensemble &particles = mix.Particles();
        
        // Flamelet advection needs some information on adjoining cells to calculate gradients
        const Geometry::LocalGeometry1d geom(reac.getGeometry(), i);
        std::vector<const Sweep::Cell*> neighbouringCells(2, NULL);
        
        // See if there is a left neighbour
        int neighbourIndex = geom.calcDestination(Geometry::left);
        if(neighbourIndex >= 0)
            neighbouringCells[0] = reac.getCell(neighbourIndex).Mixture();

        // See if there is a right neighbour
        neighbourIndex = geom.calcDestination(Geometry::right);
        if(neighbourIndex >= 0)
            neighbouringCells[1] = reac.getCell(neighbourIndex).Mixture();

        // Get rid of any scaling in the binary tree, so that a list of
        // particles can be put into the tree at the end of the transport
        // without needing to update any scaling variables.
        mix.SetM0(particles.Count() / mix.SampleVolume());

        // Get a list of the particles from the ensemble.  New instances of the
        // particles are required as the existing particles will be deleted
        // when a new list is set on the ensemble.
        Sweep::PartPtrList partList;
        Sweep::Ensemble::const_iterator itEns = particles.begin();
        const Sweep::Ensemble::const_iterator itEnsEnd = particles.end();
        while(itEns != itEnsEnd) {
            partList.push_back(new Sweep::Particle(**itEns++));
        }

        // Now go through the particles updating their position
        Sweep::PartPtrList::iterator itPart = partList.begin();
        const Sweep::PartPtrList::iterator itPartEnd = partList.end();
        while(itPart != itPartEnd ) {

            real newPosition = (*itPart)->getPosition();
            //calculate position increment
            {
                // Time since position was last calculated
                const real dt = t_stop - (*itPart)->getPositionTime();

                const real velocity = 
                    reac.getMechanism().ParticleMech().AdvectionVelocity(
                                                        mix,
                                                        **itPart,
                                                        neighbouringCells,
                                                        geom); 
                // Distance covered due to advection
                if(mSplitAdvection)
                    newPosition += dt * velocity;

                const real diffusionCoeff = reac.getMechanism().ParticleMech().DiffusionCoefficient(mix, **itPart);

                if(mSplitDiffusion) {
                    //Sweep::randNorm(0.0, std::sqrt(diffusionCoeff * dt));
                    newPosition += Sweep::randNorm(0.0, std::sqrt(diffusionCoeff * dt));
                }
            }

            // Update the particle position, overwriting its old value
            (*itPart)->setPositionAndTime(newPosition, t_stop);

            // check whether the particle has changed cells
            if(!reac.getGeometry().isInCell(i, newPosition)) {
                //particle is moving between cells
                Sweep::Transport::TransportOutflow out;

                // Dereference iterator to get raw pointer to particle
                out.particle = *itPart;

                // Find the index of the cell into which the particle is moving
                out.destination = reac.getGeometry().containingCell(newPosition);

                if(out.destination >= 0) {
                    // Statistical weight is adjusted by the ratio of physical volumes
                    // of the cells so that the number of physical particles
                    // represented by the computational particle does not change
                    // during transport.  Recall that statistical weight is the
                    // concentration of physical particles represented by a computational
                    // particle and that the number of physical particles in a cell
                    // will be the conentration multiplied by the cell volume.
                    out.weight = reac.getGeometry().cellVolume(i)
                                  / reac.getGeometry().cellVolume(out.destination)
                                  / mix.SampleVolume();

                    // Add the details of the particle to a list ready for inserting
                    // into its destination cell
                    inflowLists[out.destination].push_back(out);

//                    std::cout << "Moved particle from cell " << i << " to "
//                              << out.destination <<'\n';
                }
                else {
                    // Particle has left simulation domain
                    delete out.particle;
                    out.particle = NULL;
                    out.weight = 0.0;
                }

                // Remove the particle from this cell and move the iterator on
                // to the next particle
                itPart = partList.erase(itPart);
            }
            else {
                // particle remains in this cell so move on to next item in list
                ++itPart;
            }
        } // loop over all particles in cell i updating their position

        // Now put the particles that are staying in cell i back into that cell
        particles.SetParticles(partList.begin(), partList.end());
        
    }// loop over all cells and build up lists of particles that are moving cells

    // Add particles that are moving to new cells to their destinations
    for(unsigned int i = 0; i != numCells; ++i) {
        // Process the list waiting to be added to cell i
        for(std::list<Sweep::Transport::TransportOutflow>::iterator itPart = inflowLists[i].begin();
            itPart != inflowLists[i].end();
            ++itPart) {
                // At the moment particles have to be added one by one, a more
                // efficient method could be devised.
                transportIn(reac, i, *itPart);

                // Delete the particles, which have now had clones added to
                // their destination populations
                delete itPart->particle;
                itPart->particle = NULL;
        }

        // Now the population of cell i has been updated doubling can be
        // reactivated.
        reac.getCell(i).Mixture()->Particles().UnfreezeDoubling();
    }
}


/*!
 * Set up the cache and fill it with rates
 * 
 *@param[in]    reac    Reactor with which this cache will be used
 */
Brush::PredCorrSolver::JumpRateCache::JumpRateCache(const Reactor1d& reac)
    : mCellRates(reac.getNumCells())
    , mProcessRates(reac.getNumCells(),
                    fvector(reac.getMechanism().ParticleMech().ProcessCount(), 0.0))
    , mInvalidCells()
    , mReac(reac)
{
    // Update needs to calculate all the cells on construction, so mark them
    // all as invalid
    for(size_t i = 0; i < reac.getNumCells(); ++i) {
        mInvalidCells.push(i);
    }

    update();
}

/*!
 * Return the total rate across all cells assuming all data is correct and does
 * not need updating.
 *
 *@return   Total process rate for entire reactor
 */
Brush::real Brush::PredCorrSolver::JumpRateCache::totalRate() const {
    return std::accumulate(mCellRates.begin(), mCellRates.end(), static_cast<real>(0.0));
}

/*!
 * Work through the list of cells for which recalculations are required
 */
void Brush::PredCorrSolver::JumpRateCache::update() {
    while(!mInvalidCells.empty()) {
        const size_t cellIndex = mInvalidCells.top();

        //std::cout << "JumpRateCache updating cell " << cellIndex << '\n';

         // Get the geometry information for the cell
        Geometry::LocalGeometry1d cellGeom(mReac.getGeometry(), cellIndex);

        // Update the cache data, note that mProcessRates[i] is passed as a
        // non-const reference for overwriting.
        mCellRates[cellIndex] =
            mReac.getMechanism().ParticleMech().CalcJumpRateTerms(
                    mReac.getTime(),
                    *(mReac.getCell(cellIndex).Mixture()),
                    cellGeom, mProcessRates[cellIndex]);
        
        // Remove the index now the update has been carried out
        mInvalidCells.pop();
    }
}
