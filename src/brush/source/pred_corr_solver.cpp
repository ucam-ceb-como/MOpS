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
 *@param[in]    split_diffusion         True if diffusion is to be simulated via splitting
 *@param[in]    drift_adjustment        Only relevant when split_diffusion is true, see \ref mDiffusionDriftAdjustment
 *@param[in]    split_advection         True if advection is to be simulated via splitting
 *@param[in]    strang_advection        Use Strang splitting between transport and particle processes
 *@param[in]    cstr_transport          Transport between cell centres with CSTR random jumps
 *
 *@todo     Need generalisation to move away from using ResetChemistry
 */
/*
 * Following parameters are skipped from the doxygen documentation, becuase they are not
 * names, which is because they are not yet used.
 * param[in]                            Number of corrector iterations to perform during full coupling
 * param[in]                            Relative tolerance for ODE solver used for gas phase
 * param[in]                            Absolute tolerance for ODE solver used for gas phase
 */
Brush::PredCorrSolver::PredCorrSolver(const ResetChemistry& reset_chem,
                                      const size_t,
                                      const real , const real ,
                                      const bool split_diffusion,
                                      const real drift_adjustment,
                                      const bool split_advection,
                                      const bool strang_splitting,
                                      const bool cstr_transport)
    : mResetChemistry(reset_chem)
    , mDeferralRatio(10.0)
    , mSplitDiffusion(split_diffusion)
    , mDiffusionDriftAdjustment(drift_adjustment)
    , mSplitAdvection(split_advection)
    , mStrangTransportSplitting(strang_splitting)
    , mCSTRTransport(cstr_transport)
{}

/*!
 * Advance the solution to t_stop by means of n_steps equal length application
 * of the predictor corrector algorithm.
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_start     Time from which to advance solution
 *\param[in]            t_stop      Time to which to advance solution
 *\param[in]            n_steps     Number of predictor corrector steps (ie number of splits between transport and particle processes)
 *\param[in]            n_iter      Number of corrector iterations per step
 *\param[in]            seed        Value that is unique to this particular time interval and path to use in seeding the RNGs
 */
void Brush::PredCorrSolver::solve(Reactor1d &reac, const real t_start, const real t_stop, const int n_steps,
                                  const int n_iter, size_t seed) const {
    const real dt = (t_stop - t_start) / n_steps;

    // Start building a RNG seed for this part of the calculation.  The
    // result should be different in each call, if this function is called
    // repeatedly in order to simulate forward in time.
    boost::hash_combine(seed, t_start);
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
        predictorCorrectorStep(reac, t_start + (i - 1) * dt, t_start + i * dt, n_iter, cellRNGs);
    }
}

/*!
 * Advance the solution to t_stop by means of a predictor step followed
 * by zero or more corrector iterations
 *
 *\param[in,out]        reac        Reactor describing system state
 *\param[in]            t_start     Time from which to advance solution
 *\param[in]            t_stop      Time to which to advance solution
 *\param[in]            n_iter      Number of corrector iterations per step
 *\param[in]            cell_rngs   Vector of independent RNGs, one for each cell
 *
 *\pre  reac.getNumCells() == cell_rngs.size()
 */
void Brush::PredCorrSolver::predictorCorrectorStep(Reactor1d &reac, const real t_start,
                                                   const real t_stop, const int n_iter,
                                                   std::vector<Sweep::rng_type>& cell_rngs) const {
    //std::cout << "Predictor corrector step from " << t_start << " to " << t_stop << '\n';

    //=========== Predictor step =====================
    // Advance the chemistry
    solveChemistry(reac, t_stop);


    solveParticlesByCell(reac, t_start, t_stop, cell_rngs);

//    std::cout << "Particle counts and sample volume at end of predictorCorrectorStep ";
//    for(size_t i = 0; i < reac.getNumCells(); ++i) {
//        std::cout << reac.getCell(i).ParticleCount() << ' ' << reac.getCell(i).SampleVolume();
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
        oldDensity[i] = reac.getCell(i).GasPhase().MassDensity();
    }

    // Update the chemistry to the new time
    reac.ReplaceChemistry(mResetChemistry, true);

    // Update the soot particle concentration to account for the expansion /
    // contraction of the gas.
    for(size_t i = 0; i != numCells; ++i) {
        Sweep::Cell & mix = reac.getCell(i);

        // Collect data for calculating the new value
        mix.AdjustSampleVolume(oldDensity[i] / mix.GasPhase().MassDensity());
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
 *\param[in]            t_start     Time from which to advance solution
 *\param[in]            t_stop      Time to which to advance solution
 *\param[in]            cell_rngs   Vector of independent RNGs, one for each cell
 */
void Brush::PredCorrSolver::solveParticlesByCell(Reactor1d &reac, const real t_start, const real t_stop,
                                                 std::vector<Sweep::rng_type>& cell_rngs) const {

    const size_t numCells = reac.getNumCells();
    const Sweep::Mechanism &mech = reac.getParticleMechanism();

    // Select time step for Strang or first order splitting
    const real firstStop = mStrangTransportSplitting ? (t_start + t_stop) / 2.0 : t_stop;

#pragma omp parallel for schedule(dynamic) ordered
    for(size_t i = numCells; i > 0; --i) {
        // Get details of cell i
        Sweep::Cell& cell = reac.getCell(i - 1);
        Geometry::LocalGeometry1d geom(reac.getGeometry(), i - 1);

        solveParticlesInOneCell(cell, geom, mech, t_start, firstStop, cell_rngs[i - 1]);
    }

    // Now do the split particle transport, if there is any
    if(mSplitDiffusion || mSplitAdvection) {
        splitParticleTransport(reac, t_start, t_stop, cell_rngs);
    }

    if(mStrangTransportSplitting) {
        // Do the second part of the particle simulation for the Strang splitting
#pragma omp parallel for schedule(dynamic) ordered
        for(size_t i = numCells; i > 0; --i) {
            // Get details of cell i
            Sweep::Cell& cell = reac.getCell(i - 1);
            Geometry::LocalGeometry1d geom(reac.getGeometry(), i - 1);

            solveParticlesInOneCell(cell, geom, mech, firstStop, t_stop, cell_rngs[i - 1]);
        }
    }

}

/*!
 * Advance the particle part of the solution, which may also affect the chemical
 * species concentrations.
 *
 * @todo  The choice of deferral length should be integrated with Sweep::Solver
 * with a view to using Sweep::Solver::Run()
 *
 *\param[in,out]        cell        Contents of one grid cell
 *\param[in]            geom        Position information regarding adjoining cells
 *\param[in]            mech        Mechanism defining the particle processes
 *\param[in]            t           Time from which to advance solution
 *\param[in]            t_stop      Time to which to advance solution
 *\param[in,out]        rng         RNG to use for the simulation
 */
void Brush::PredCorrSolver::solveParticlesInOneCell(Sweep::Cell &cell, const Geometry::LocalGeometry1d &geom,
                                                    const Sweep::Mechanism &mech, real t, real t_stop,
                                                    Sweep::rng_type &rng) const {
    // Carry out a repeated sequence of jump process simulation (with LPDA updates
    // for particles involved in jumps) followed by LPDA updates for the full
    // population.
    do {
        // Ignore the individual rate term information in the final argument
        fvector dummyVec;
        const real deferredRate = mech.CalcDeferredRateTerms(t, cell ,geom, dummyVec);

        // Calculate both the rates of the individual jump processes and the total
        fvector jumpRates;
        real jumpRate = mech.CalcJumpRateTerms(t, cell, geom, jumpRates);

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

            jumpRate = mech.CalcJumpRateTerms(t, cell, geom, jumpRates);
            Sweep::Solver::timeStep(t, maxDeferralEnd, cell, geom, mech,
                                    jumpRates, jumpRate, rng);
        }

        // Perform all events from deferred processes.
        if (mech.AnyDeferred()) {
            mech.LPDA(maxDeferralEnd, cell, rng);
        }
    } while(t <= t_stop * (1.0 - std::numeric_limits<real>::epsilon()));

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
    // the lists of particles being transported into each cell of the reactor.
    std::vector<inflow_lists_vector> inflowLists(numCells, inflow_lists_vector(numCells));
    std::vector<real> statisticalWeights(numCells);
    std::vector<real> cellVolumes(numCells);
    std::vector<real> maxCapacities(numCells);
    std::vector<real> velocities(numCells);

    // Loop over each cell updating particles positions and removing any
    // particles that are moving to new cells.
#pragma omp parallel for schedule(dynamic)
    for(unsigned int i = 0; i < numCells; ++i) {
        Sweep::Cell &mix = reac.getCell(i);

        // Flamelet advection needs some information on adjoining cells to calculate gradients
        const Geometry::LocalGeometry1d geom(reac.getGeometry(), i);
        //std::vector<const Sweep::Cell*> neighbouringCells(2, NULL);
        std::vector<const Sweep::Cell*> neighbouringCells(2);
        neighbouringCells[0] = neighbouringCells[1] = NULL;

        // See if there is a left neighbour
        int neighbourIndex = geom.calcDestination(Geometry::left);
        if(neighbourIndex >= 0)
            neighbouringCells[0] = &(reac.getCell(neighbourIndex));

        // See if there is a right neighbour
        neighbourIndex = geom.calcDestination(Geometry::right);
        if(neighbourIndex >= 0)
            neighbouringCells[1] = &(reac.getCell(neighbourIndex));

        // The weight will be needed when the particles are put back into the cell
        statisticalWeights[i] = 1.0 / mix.SampleVolume();
        cellVolumes[i] = geom.cellVolume();
        maxCapacities[i] = static_cast<real>(mix.Particles().Capacity());
        velocities[i] = mix.GasPhase().Velocity();

        // Remove all particles from the cell and add them
        // to the inflow list for the appropriate cell
        inflowLists[i] = updateParticleListPositions(t_start, t_stop, mix, i,
                                                     reac.getParticleMechanism(),
                                                     reac.getGeometry(), neighbouringCells,
                                                     cell_rngs[i]);
    }// loop over all cells and build up lists of particles that are moving cells

    // There has to be a synchronisation point here, in that parallel replacement of
    // particles in cells cannot begin until the positions and cells of all particles
    // have been calculated.

#pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < numCells; ++i) {
        // i is the index of the destination cell

        // Build up a list of particle pointers to pass into the cell
        // this variable is local to each loop iteration
        Sweep::PartPtrList partList;

        for(size_t j = 0; j != numCells; ++j) {
            // Adjustment factor for particle statistical weight
            const real weightFactor = (maxCapacities[j] * statisticalWeights[j]) /
                                      (maxCapacities[i] * statisticalWeights[i]);

            const real repeatCountConst =  cellVolumes[j] * maxCapacities[i] * velocities[i]
                                         / cellVolumes[i] / maxCapacities[j] / velocities[j];
//            if(i == j + 1) {
//                std::cout << j << ' ' << " max capac " << maxCapacities[j] << " inv sample vol " << statisticalWeights[j] << " cell volume " << cellVolumes[j] <<std::endl;
//                std::cout << i << ' ' << " max capac " << maxCapacities[i] << " inv sample vol " << statisticalWeights[i] << " cell volume " << cellVolumes[i] <<std::endl;
//                std::cout << j << ' ' << i << ' ' << " repeat count " << repeatCountConst << " weight factor " << weightFactor  << std::endl;
//            }

            // Work through the particles coming from cell j
            std::list<Sweep::Particle*>::const_iterator it = inflowLists[j][i].begin();
            const std::list<Sweep::Particle*>::const_iterator itEnd = inflowLists[j][i].end();
            while(it != itEnd) {
                // Repeat count in destination cell (may be less than one), fractional parts are probability
                // of a particle being added to destination.
                real repeatCount = repeatCountConst;

                (*it)->setStatisticalWeight((*it)->getStatisticalWeight() * weightFactor);

                while(repeatCount > 0) {
                    if(repeatCount > 1) {
                        partList.push_back((*it)->Clone());
                    }
                    else {
                        // Final copy of the particle can be the original, but is only added with
                        // probability repeatCount.
                        boost::random::bernoulli_distribution<real> particleDecider(repeatCount);
                        if(particleDecider(cell_rngs[i]))
                            partList.push_back((*it));
                        else {
                            delete (*it);
                        }
                    }
                    repeatCount -= 1.0;
                }
                ++it;
            }
        }
        if(partList.size() > reac.getCell(i).Particles().Capacity())
            std::cout << "Setting " << partList.size() << " particles on cell with center " << reac.getGeometry().cellCentre(i) << std::endl;

        reac.getCell(i).SetParticles(partList.begin(), partList.end(), statisticalWeights[i], cell_rngs[i]);
    }
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
void Brush::PredCorrSolver::updateParticlePosition(const real t_start, const real t_stop, const Sweep::Cell &mix,
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
            const real velocity = mech.AdvectionVelocity(mix, sp, neighbouringCells, geom);
            newPosition += dt * velocity;
        }
        if(mSplitDiffusion){
            const real diffusionCoeff = mech.DiffusionCoefficient(mix, sp);
            typedef boost::normal_distribution<real> normal_distrib;
            normal_distrib diffusionDistrib(0.0, std::sqrt(diffusionCoeff * dt));
            boost::variate_generator<Sweep::rng_type&, normal_distrib> diffusionGenerator(rng, diffusionDistrib);
            newPosition += diffusionGenerator();

            if(mDiffusionDriftAdjustment > 0.0) {
                // Drift correction to allow for different forms of stochastic integral
                const real drift = mDiffusionDriftAdjustment * mech.GradDiffusionCoefficient(mix, sp, neighbouringCells, geom);
                newPosition += dt * drift;
            }
        }
    }
    sp.setPositionAndTime(newPosition, t_stop);
}


/*!
 * CSTR stands for Continuously Stirred Tank Reactor, a chemical engineering term
 * referring to a perfectly mixed reactor.  Because it is perfectly mixed particles
 * cannot be said to flow through the reactor, rather they have a residence time
 * distribution.  In this case the expected residence time is taken as the length
 * of the reactor divided by the velocity and each particle has a probability of
 * leaving the reactor at each time step.  The probability is given by the ratio
 * of the time step length to the expected residence time.
 *
 *@param[in]        t_start             Time at which position was last calculated by splitting
 *@param[in]        t_stop              Time at which new position must be calculated
 *@param[in]        mix                 Mixture in which the particle is moving
 *@param[in]        mech                Mechanism specifying calculation of particle transport properties
 *@param[in]        geom                Information on locations of surrounding cells
 *@param[in]        neighbouringCells   Pointers to the contents of surrounding cells
 *@param[in,out]    sp                  Particle requiring updated position
 *@param[in,out]    rng                 Random number generator
 *
 *@pre      The expected residence time of the particle must be greater than the time step length.
 */
void Brush::PredCorrSolver::updateParticlePositionCSTR(const real t_start, const real t_stop, const Sweep::Cell &mix,
                                                       const Sweep::Mechanism &mech,
                                                       const Geometry::LocalGeometry1d & geom,
                                                       const std::vector<const Sweep::Cell*> & neighbouringCells,
                                                       Sweep::Particle& sp, Sweep::rng_type &rng) const
{
    // Time over which transport is occurring
    const real dt = t_stop - t_start;

    // Probability according to CSTR that particle leaves cell in this time step
    const real p = mech.AdvectionVelocity(mix, sp, neighbouringCells, geom) * dt / geom.cellVolume();
    assert(p <= 1);

    // The particle will either end up at the centre of this cell or the next
    real newPosition = geom.cellCentre();

    // Generate a sample to decide if the particle leaves
    boost::random::bernoulli_distribution<real> leaveCell(p);
    if(leaveCell(rng))
        // Add the distance to the next cell centre
        newPosition += geom.calcSpacing(Geometry::right);

    sp.setPositionAndTime(newPosition, t_stop);
}


/*!
 *@param[in]        t_start             Time at which position was last calculated by splitting
 *@param[in]        t_stop              Time at which new position must be calculated
 *@param[in,out]    mix                 Mixture from which the particles are moving
 *@param[in]        cell_index          Index of cell containing the particles to be transported
 *@param[in]        mech                Mechanism specifying calculation of particle transport properties
 *@param[in]        geom                Information on locations of surrounding cells
 *@param[in]        neighbouringCells   Pointers to the contents of surrounding cells
 *@param[in,out]    rng                 Random number generator
 *
 *@return       Vector of lists of particles to be transported into other cells
 */
Brush::PredCorrSolver::inflow_lists_vector
  Brush::PredCorrSolver::updateParticleListPositions(const real t_start, const real t_stop, Sweep::Cell &mix,
                                                     const size_t cell_index, const Sweep::Mechanism &mech,
                                                     const Geometry::Geometry1d & geom,
                                                     const std::vector<const Sweep::Cell*> & neighbouringCells,
                                                     Sweep::rng_type &rng) const {
    // Build up the return value in this vector of lists
    inflow_lists_vector outflow(geom.numCells());

    // The local geometry will be needed repeatedly for this cell
    Geometry::LocalGeometry1d localGeom(geom, cell_index);

    // Now go through the particles updating their position
    Sweep::Ensemble::iterator itPart = mix.Particles().begin();
    const Sweep::Ensemble::iterator itPartEnd = mix.Particles().end();
    while(itPart != itPartEnd ) {
        // Actually calculate new position of particle
        if(mCSTRTransport)
            updateParticlePositionCSTR(t_start, t_stop, mix, mech, localGeom, neighbouringCells, **itPart, rng);
        else
            updateParticlePosition(t_start, t_stop, mix, mech, localGeom, neighbouringCells, **itPart, rng);

        // Find the index of the cell into which the particle is moving
        const int destination = geom.containingCell((*itPart)->getPosition());

        if(destination >= 0) {
            // On moving to a new cell reset the coagulation count
            if(static_cast<unsigned>(destination) != cell_index)
               (*itPart)->resetCoagCount();

            // Add the details of the particle to a list ready for inserting
            // into its destination cell
            outflow[destination].push_back(*itPart);

        }
        else {
            // Particle has left simulation domain
            delete *itPart;
            *itPart = NULL;
        }

        // Remove the particle from this cell and move the iterator on
        // to the next particle
        ++itPart;
    } // loop over all particles in cell i updating their position

    // Empty the particle ensemble
    mix.Particles().TakeParticles();

    return outflow;
}
