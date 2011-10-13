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
#include "local_geometry1d.h"

#include "swp_transport_outflow.h"
#include "swp_cell.h"
#include "swp_solver.h"

#include "choose_index.hpp"

#include <stdexcept>
#include <limits>
#include <algorithm>
#include <vector>
#include <list>
#include <cassert>
#include <boost/functional/hash.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>

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
 *@param[in]    weighted_transport      Adjust weights during inter-cell transport to avoid killing/cloning particles
 *
 *@todo     Need generalisation to move away from using ResetChemistry
 */
Brush::PredCorrSolver::PredCorrSolver(const ResetChemistry& reset_chem,
                                      const size_t,
                                      const real , const real ,
                                      const bool split_diffusion,
                                      const bool split_advection,
                                      const bool weighted_transport)
    : mResetChemistry(reset_chem)
    , mDeferralRatio(10.0)
    , mSplitDiffusion(split_diffusion)
    , mSplitAdvection(split_advection)
    , mWeightedTransport(weighted_transport)
{}

/*!
 * Advance the solution to t_stop by means of n_steps equal length application
 * of the predictor corrector algorithm.
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_stop      Time to which to advance solution
 *\param[in]            n_steps     Number of predictor corrector steps (ie number of splits between transport and particle processes)
 *\param[in]            n_iter      Number of corrector iterations per step
 *\param[in]            seed        Value that is unique to this particular time interval and path to use in seeding the RNGs
 */
void Brush::PredCorrSolver::solve(Reactor1d &reac, const real t_stop, const int n_steps,
                                  const int n_iter, size_t seed) const {
    const real startTime = reac.getTime();
    const real dt = (t_stop - startTime) / n_steps;

    // Start building a RNG seed for this part of the calculation.  The
    // result should be different in each call, if this function is called
    // repeatedly in order to simulate forward in time.
    boost::hash_combine(seed, reac.getTime());
    boost::hash_combine(seed, t_stop);

    // Build one RNG for each cell, so that the cells can be simulated
    // independently apart from inter-cell transport.
    std::vector<Sweep::rng_type> cellRNGs(reac.getNumCells());
    for(unsigned i = 0; i < reac.getNumCells(); ++i) {
        // Now make a RNG seed that is unique for each cell
        size_t cellSeed = seed;
        boost::hash_combine(cellSeed, i);
        cellRNGs[i].seed(cellSeed);
    }

    for(int i = 1; i <= n_steps; ++i) {
        predictorCorrectorStep(reac, startTime + i * dt, n_iter, cellRNGs);
    }
}

/*!
 * Advance the solution to t_stop by means of a predictor step followed
 * by zero or more corrector iterations
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_stop      Time to which to advance solution
 *\param[in]            n_iter      Number of corrector iterations per step
 *\param[in]            cell_rngs   Vector of independent RNGs, one for each cell
 *
 *\pre  reac.getNumCells() == cell_rngs.size()
 */
void Brush::PredCorrSolver::predictorCorrectorStep(Reactor1d &reac, const real t_stop,
                                                   const int n_iter,
                                                   std::vector<Sweep::rng_type>& cell_rngs) const {
    const real startTime = reac.getTime();

    for(size_t i = 0; i < reac.getNumCells(); ++i) {
        assert(reac.getCell(i).Mixture() != NULL);
    }

    //=========== Predictor step =====================
    // Advance the chemistry
    solveChemistry(reac, t_stop);

    // reset the time and advance the particles
    reac.setTime(startTime);
    solveParticlesByCell(reac, t_stop, cell_rngs);

//    std::cout << "Particle counts at end of predictorCorrectorStep ";
//    for(size_t i = 0; i < reac.getNumCells(); ++i) {
//        std::cout << reac.getCell(i).Mixture()->ParticleCount() << ' ';
//    }
//    std::cout << std::endl;

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
        oldDensity[i] = reac.getCell(i).Mixture()->GasPhase().MassDensity();
    }

    // Update the chemistry to the new time
    reac.ReplaceChemistry(mResetChemistry, true);

    // Update the soot particle concentration to account for the expansion /
    // contraction of the gas.
    for(size_t i = 0; i != numCells; ++i) {
        Sweep::Cell * const pMix = reac.getCell(i).Mixture();

        // Collect data for calculating the new value
        pMix->AdjustSampleVolume(oldDensity[i] / pMix->GasPhase().MassDensity());
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
 *\param[in]            cell_rngs   Vector of independent RNGs, one for each cell
 */
void Brush::PredCorrSolver::solveParticlesByCell(Reactor1d &reac, const real t_stop,
                                                 std::vector<Sweep::rng_type>& cell_rngs) const {

    const size_t numCells = reac.getNumCells();
    const Sweep::Mechanism &mech = reac.getMechanism().ParticleMech();

    // Need to know the time for the split simulation of the transport processes
    const real t_start = reac.getTime();

    #pragma omp parallel for
    for(size_t i = 0; i < numCells; ++i) {
        // Get details of cell i
        Mops::Reactor& cell = reac.getCell(i);
        Geometry::LocalGeometry1d geom(reac.getGeometry(), i);

        solveParticlesInOneCell(cell, geom, mech, t_stop, cell_rngs[i]);
    }

    // Now do the split particle transport, if there is any
    if(mSplitDiffusion || mSplitAdvection) {
        splitParticleTransport(reac, t_start, t_stop, cell_rngs);
    }
}

/*!
 * Advance the particle part of the solution, which may also affect the chemical
 * species concentrations.
 *
 * @TODO  The choice of deferral length should be integrated with Sweep::Solver
 * with a view to using Sweep::Solver::Run()
 *
 *\param[in,out]        cell        Contents of one grid cell
 *\param[in]            geom        Position information regarding adjoining cells
 *\param[in]            mech        Mechanism defining the particle processes
 *\param[in]            t_stop      Time to which to advance solution
 *\param[in,out]        rng         RNG to use for the simulation
 */
void Brush::PredCorrSolver::solveParticlesInOneCell(Mops::Reactor &cell, const Geometry::LocalGeometry1d &geom,
                                                    const Sweep::Mechanism &mech, const real t_stop,
                                                    Sweep::rng_type &rng) const {
    // Carry out a repeated sequence of jump process simulation (with LPDA updates
    // for particles involved in jumps) followed by LPDA updates for the full
    // population.
    real t = cell.Time();
    do {
        // Ignore the individual rate term information in the final argument
        fvector dummyVec;
        // Add 1 to the denominator to avoid dividing by 0
        const real deferredRate = mech.CalcDeferredRateTerms(t, *(cell.Mixture()) ,geom, dummyVec);

        // Calculate both the rates of the individual jump processes and the total
        fvector jumpRates;
        real jumpRate = mech.CalcJumpRateTerms(t, *(cell.Mixture()), geom, jumpRates);

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
        while(t <= maxDeferralEnd * (1.0 - std::numeric_limits<real>::epsilon())) {
            // Perform one non-deferred event

            jumpRate = mech.CalcJumpRateTerms(t, *(cell.Mixture()), geom, jumpRates);
            Sweep::Solver::timeStep(t, maxDeferralEnd, *(cell.Mixture()), geom, mech,
                                    jumpRates, jumpRate, rng);
        }

        cell.SetTime(t);

        // Perform all events from deferred processes.
        if (mech.AnyDeferred()) {
            mech.LPDA(maxDeferralEnd, *(cell.Mixture()), rng);
        }
    } while(t <= t_stop * (1.0 - std::numeric_limits<real>::epsilon()));

 }

/*!
 * Add a particle that has left one cell into the specified cell.  The particle_details
 * must contain a pointer to a valid particle, ownership of which is relinquished by
 * the caller.
 *
 *\param[in,out]    reac                    1D reactor within which the transport is occurring
 *\param[in]        destination_index       Index of destination cell
 *\param[in]        particle_details        Information about the particle and its weight
 *\param[in,out]    rng                     Random number generator
 */
void Brush::PredCorrSolver::transportIn(Reactor1d & reac, const size_t destination_index,
                                        const Sweep::Transport::TransportOutflow &particle_details,
                                        Sweep::rng_type &rng) const {
    real incomingWeight = particle_details.weight;

    if(mWeightedTransport) {
        // Exploit weights to ensure no particles are created or destroyed

        // This will exactly preserve the physical contribution of the particle under consideration, because
        // statistical weights are always divided by the sample volume.
        particle_details.particle->setStatisticalWeight(particle_details.particle->getStatisticalWeight() *
                                                        incomingWeight * reac.getCell(destination_index).Mixture()->SampleVolume());

        // Ownership of the particle is now taken by the ensemble
        reac.getCell(destination_index).Mixture()->Particles().Add(*particle_details.particle, rng);
    }
    else {
        // Need to match weights of particle between source and destination in a DSA setting
        unsigned int safetyCounter = 0;
        while(true) {
            real destinationWeight = 1.0 / reac.getCell(destination_index).Mixture()->SampleVolume();

            if(incomingWeight >= destinationWeight) {
                // The important thing is that the statistical weight of particles added to the destination cell
                // is equal to the statistical weight of the particle that was removed from the originating cell.
                // Otherwise mass will be lost.

                // Insert one copy of the particle into the destination cell
                reac.getCell(destination_index).Mixture()->Particles().Add(*(new Sweep::Particle(*particle_details.particle)), rng);

                // One unit of destinationWeight has now been added to an ensemble
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
        const real moveParticleProb = incomingWeight * reac.getCell(destination_index).Mixture()->SampleVolume();
        typedef boost::bernoulli_distribution<real> bernoulli_distrib;
        bernoulli_distrib moveDistrib(moveParticleProb);
        boost::variate_generator<Sweep::rng_type&, bernoulli_distrib> moveDecider(rng, moveDistrib);

        if(moveDecider()) {
            reac.getCell(destination_index).Mixture()->Particles().Add(*particle_details.particle, rng);

            // Testing output
            const real extraWeight = (1.0 / reac.getCell(destination_index).Mixture()->SampleVolume()) - incomingWeight;
//            if(std::abs(extraWeight) > 100 * std::numeric_limits<real>::epsilon() / reac.getCell(destination_index).Mixture()->SampleVolume()) {
//                std::cerr << "Transport in added " << extraWeight << " of statistical weight to cell centred at "
//                          << reac.getCellCentre(destination_index) << std::endl;
//            }
        }
        else {
            // Ownership of the particle has not been passed on to an ensemble so the memory must be released
            delete particle_details.particle;

            // Testing output
//            if(std::abs(incomingWeight) > 100 * std::numeric_limits<real>::epsilon() / reac.getCell(destination_index).Mixture()->SampleVolume()) {
//                std::cerr << "Transport in dropped " << incomingWeight
//                          << " of statistical weight from cell centred at "
//                          << reac.getCellCentre(destination_index) << std::endl;
//            }
        }
    }
}

/*!
 *@param[in,out]    reac        Reactor in which particles are being transported
 *@param[in]        t_start     Time at which position was last calculated by splitting
 *@param[in]        t_stop      Time upto which transport is to be simulated
 *@param[in,out]    cell_rngs   Random number generators, one per cell
 */
void Brush::PredCorrSolver::splitParticleTransport(Reactor1d &reac, const real t_start,
                                                   const real t_stop, std::vector<Sweep::rng_type>& cell_rngs) const {
    const size_t numCells = reac.getNumCells();

    // Element i of the outer vector will contain the outflow from cell i of the reactor.
    // Element j of the inner vector will contain particles destined for cell j of the reactor.
    // After this vector of vectors of lists had been populated, it will be necessary to merge
    // the lists of particles being transported into each cell of the reactor.  The origins
    // of the cells will not matter, because statistical weight information is already in
    // the list entries.
    std::vector<inflow_lists_vector> inflowLists(numCells, inflow_lists_vector(numCells));
    std::vector<inflow_lists_vector> secondaryInflowLists(numCells, inflow_lists_vector(numCells));

    // Loop over each cell updating particles positions and removing any
    // particles that are moving to new cells.
#pragma omp parallel for
    for(unsigned int i = 0; i < numCells; ++i) {
        Mops::Mixture &mix = *(reac.getCell(i).Mixture());
        Sweep::Ensemble &particles = mix.Particles();

        // Flamelet advection needs some information on adjoining cells to calculate gradients
        const Geometry::LocalGeometry1d geom(reac.getGeometry(), i);
        //std::vector<const Sweep::Cell*> neighbouringCells(2, NULL);
        std::vector<const Sweep::Cell*> neighbouringCells(2);
        neighbouringCells[0] = neighbouringCells[1] = NULL;

        // See if there is a left neighbour
        int neighbourIndex = geom.calcDestination(Geometry::left);
        if(neighbourIndex >= 0)
            neighbouringCells[0] = reac.getCell(neighbourIndex).Mixture();

        // See if there is a right neighbour
        neighbourIndex = geom.calcDestination(Geometry::right);
        if(neighbourIndex >= 0)
            neighbouringCells[1] = reac.getCell(neighbourIndex).Mixture();

        // The weight will be needed when the particles are put back into the cell
        const real statisticalWeight = 1.0 / mix.SampleVolume();

        // Get a list of the particles from the ensemble.  New instances of the
        // particles are required as the existing particles will be deleted
        // when a new list is set on the ensemble.
        Sweep::PartPtrList partList = particles.TakeParticles();

        // Remove particles that are moving to a new cell from partLists and add them
        // to the inflow list for the appropriate cell
        inflowLists[i] = updateParticleListPositions(t_start, t_stop, mix, i,
                                                     reac.getMechanism().ParticleMech(),
                                                     reac.getGeometry(), neighbouringCells,
                                                     partList, cell_rngs[i]);

        // Now put the particles that are staying in cell i back into that cell
        mix.SetParticles(partList.begin(), partList.end(), statisticalWeight);

    }// loop over all cells and build up lists of particles that are moving cells

    // Join up the lists of particles moving into each cell, by appending them to the
    // first list
    for(size_t i = 0; i != numCells; ++i) {
        // i is the index of the destination cell

        for(size_t j = 1; j != numCells; ++j) {
            // j is the index of the source cell.
            // Data is added to the list coming from the cell j=0
            inflowLists.front()[i].splice(inflowLists.front()[i].end(), inflowLists[j][i],
                                          inflowLists[j][i].begin(), inflowLists[j][i].end());
        }
    }
    // Now put the particles in their destination cells
    moveParticlesToDifferentCells(reac, inflowLists.front(), cell_rngs);
}

/*!
 *@param[in]        t_start             Time at which position was last calculated by splitting
 *@param[in]        t_stop              Time at which new position must be calculated
 *@param[in]        mix                 Mixture in which the particle is moving
 *@param[in]        mech                Mechanism specifying calculation of particle transport properties
 *@param[in]        geom                Information on locations of surrounding cells
 *@param[in]        neighbouringCells   Pointers to the contents of surrounding cells
 *@param[in,out]    sp                  Particle requiring updated position
 *@param[in,out]    rng                 Random number generator
 */
void Brush::PredCorrSolver::updateParticlePosition(const real t_start, const real t_stop, const Mops::Mixture &mix,
                                                   const Sweep::Mechanism &mech,
                                                   const Geometry::LocalGeometry1d & geom,
                                                   const std::vector<const Sweep::Cell*> & neighbouringCells,
                                                   Sweep::Particle& sp, Sweep::rng_type &rng) const
{
    real newPosition = sp.getPosition();
    {
        const real dt = t_stop - t_start; //sp.getPositionTime();
        assert(dt >= 0.0);
        if(mSplitAdvection){
            real velocity = mech.AdvectionVelocity(mix, sp, neighbouringCells, geom);
            newPosition += dt * velocity;
        }
        if(mSplitDiffusion){
            const real diffusionCoeff = mech.DiffusionCoefficient(mix, sp);
            typedef boost::normal_distribution<real> normal_distrib;
            normal_distrib diffusionDistrib(0.0, std::sqrt(diffusionCoeff * dt));
            boost::variate_generator<Sweep::rng_type&, normal_distrib> diffusionGenerator(rng, diffusionDistrib);
            newPosition += diffusionGenerator();
        }
    }
    sp.setPositionAndTime(newPosition, t_stop);
}


/*!
 *@param[in]        t_start             Time at which position was last calculated by splitting
 *@param[in]        t_stop              Time at which new position must be calculated
 *@param[in]        mix                 Mixture in which the particle is moving
 *@param[in]        cell_index          Index of cell containing the particles to be transported
 *@param[in]        mech                Mechanism specifying calculation of particle transport properties
 *@param[in]        geom                Information on locations of surrounding cells
 *@param[in]        neighbouringCells   Pointers to the contents of surrounding cells
 *@param[in,out]    particle_list       List of pointers to the particles that belong in the cell
 *@param[in,out]    rng                 Random number generator
 *
 *@return       Vector of lists of particles to be transported into other cells
 */
Brush::PredCorrSolver::inflow_lists_vector
  Brush::PredCorrSolver::updateParticleListPositions(const real t_start, const real t_stop, const Mops::Mixture &mix,
                                                     const size_t cell_index, const Sweep::Mechanism &mech,
                                                     const Geometry::Geometry1d & geom,
                                                     const std::vector<const Sweep::Cell*> & neighbouringCells,
                                                     Sweep::PartPtrList& particle_list, Sweep::rng_type &rng) const {
    // Build up the return value in this vector of lists
    inflow_lists_vector outflow(geom.numCells());

    // The local geometry will be needed repeatedly for this cell
    Geometry::LocalGeometry1d localGeom(geom, cell_index);

    // Now go through the particles updating their position
    Sweep::PartPtrList::iterator itPart = particle_list.begin();
    const Sweep::PartPtrList::iterator itPartEnd = particle_list.end();
    while(itPart != itPartEnd ) {

        updateParticlePosition(t_start, t_stop, mix, mech, localGeom, neighbouringCells, **itPart, rng);

        if(!geom.isInCell(cell_index, (*itPart)->getPosition())) {
            //particle is moving between cells
            Sweep::Transport::TransportOutflow out;

            // Dereference iterator to get raw pointer to particle
            out.particle = *itPart;

            // Find the index of the cell into which the particle is moving
            out.destination = geom.containingCell((*itPart)->getPosition());

            if(out.destination >= 0) {
                // Statistical weight is adjusted by the ratio of physical volumes
                // of the cells so that the number of physical particles
                // represented by the computational particle does not change
                // during transport.  Recall that statistical weight is the
                // concentration of physical particles represented by a computational
                // particle and that the number of physical particles in a cell
                // will be the concentration multiplied by the cell volume.
                out.weight = geom.cellVolume(cell_index)
                              / geom.cellVolume(out.destination)
                              / mix.SampleVolume();

                // Add the details of the particle to a list ready for inserting
                // into its destination cell
                outflow[out.destination].push_back(out);

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
            itPart = particle_list.erase(itPart);
        }
        else {
            // particle remains in this cell so move on to next item in list
            ++itPart;
        }
    } // loop over all particles in cell i updating their position

    return outflow;
}

/*!
 *@param[in,out]    reac            Reactor in which transport is taking place
 *@param[in]        inflow_lists    Details of particles to be added to different cells
 *@param[in,out]    cell_rngs       Random number generators, one for each cell
 */
void Brush::PredCorrSolver::moveParticlesToDifferentCells(Reactor1d & reac,
                                                          const inflow_lists_vector & inflow_lists,
                                                          std::vector<Sweep::rng_type>& cell_rngs) const {
#pragma omp parallel for
    for(unsigned int i = 0; i < reac.getNumCells(); ++i) {
        // Process the list waiting to be added to cell i
        for(std::list<Sweep::Transport::TransportOutflow>::const_iterator itPart = inflow_lists[i].begin();
            itPart != inflow_lists[i].end();
            ++itPart) {
                // Reset the coagulation count now the particle is moving to a new cell
                itPart->particle->resetCoagCount();
                // At the moment particles have to be added one by one, a more
                // efficient method could be devised.
                transportIn(reac, i, *itPart, cell_rngs[i]);
        }

        // Now the population of cell i has been updated doubling can be
        // reactivated.
        reac.getCell(i).Mixture()->Particles().UnfreezeDoubling();
    }
}

