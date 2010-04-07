/*!
 * \file   pred_corr_solver.h
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
#ifndef BRUSH_PRED_CORR_SOLVER_H
#define BRUSH_PRED_CORR_SOLVER_H

#include "brush_params.h"

#include "reset_chemistry.h"

#include "swp_particle.h"

#include <vector>
#include <stack>
#include <list>

// Forward declarations
namespace Brush {
    class Reactor1d;
}
namespace Sprog {
    class Mechanism;
}

namespace Sweep {
    class Cell;
    class Particle;
    class Mechanism;

namespace Transport {
    struct TransportOutflow;
}
}

namespace Mops {
    class Mixture;
}

namespace Geometry {
    class LocalGeometry1d;
    class Geometry1d;
}

namespace Brush {

//! Methods for stepping forward in time.  Consider converting to namespace.
class PredCorrSolver {

public:
    //! Create with an object with details of the chemistry, which is time invariant
    PredCorrSolver(const ResetChemistry& reset_chem,
                   const size_t corrector_iterations,
                   const real rtol, const real atol,
                   const bool split_diffusion, const bool split_advection);

    //! Advance solution to specified time
    void solve(Reactor1d &reac, const real t_stop, const int n_steps, const int n_iter) const;

protected:
        /*!
     *@brief Store jump rates to avoid recalculating rates for all cells
     *
     * Rates for cell i of the reactor used for construction will be stored in
     * position i of the member vectors.
     *
     */
    class JumpRateCache {
        public:
            //! Initialise vectors
            JumpRateCache(const Reactor1d &reac);

            //! Total jump rate by cell
            fvector mCellRates;

            //! Individual process rates by cell
            std::vector<fvector> mProcessRates;

            //! Cell indices for which rates  need recalculating
            std::stack<size_t> mInvalidCells;

            //! Sum up all cell rates
            real totalRate() const;

            //! Recalculate rates for any cells listed in mInvalidCells
            void update();
        private:
            //! Cannot contrstruct without knowing size of reactor
            JumpRateCache();

            //! Reactor with which instance is used
            const Reactor1d& mReac;
    };

    //! Perform one time step using a predictor followed by some corrector iterations
    void predictorCorrectorStep(Reactor1d &reac, const real t_stop, const int n_iter) const;

    //! Advance particle population to specified time
    void solveParticles(Reactor1d &reac, const real t_stop) const;

    //! Perform one stochastic jump
    real particleTimeStep(Reactor1d &reac, const real t_stop, JumpRateCache &rate_cache) const;

    //! Advance chemistry over specified time interval
    void solveChemistry(Reactor1d & reac, const real t_stop) const; //,?some kind of workspace for ODE solver);

    //! Put a particle that has left one cell into its destination
    void transportIn(Reactor1d & reac, const size_t destination_index, const Sweep::Transport::TransportOutflow &particle_details) const;

    //! Carry out split transport on all particles from all cells
    void splitParticleTransport(Reactor1d &reac, const real t_stop) const;

private:
    //! Not possible to have a solver of this type without a reset chemistry object
    PredCorrSolver();

    //! Calculate and set a new position on one particle
    void updateParticlePosition(const real t_stop, const Mops::Mixture &mix,
                                const Sweep::Mechanism &mech,
                                const Geometry::LocalGeometry1d & geom,
                                const std::vector<const Sweep::Cell*> & neighbouring_cells,
                                Sweep::Particle& sp) const;
    
    /*!
     *  \brief Lists of particles due to be added to cells and their statistical weights
     *
     *  Particles in the the lists at index i in the vector are destined to be added to
     *  cell i of a reactor.
     */
    typedef std::vector<std::list<Sweep::Transport::TransportOutflow> > inflow_lists_vector;

    //! Update the positions on a list of particles from one cell
    inflow_lists_vector updateParticleListPositions(const real t_stop, const Mops::Mixture &mix,
                                                    const size_t cell_index, const Sweep::Mechanism &mech,
                                                    const Geometry::Geometry1d &geom,
                                                    const std::vector<const Sweep::Cell*> &neighbouring_cells,
                                                    Sweep::PartPtrList& particle_list) const;

    //! Put particles that are actually moving to new cells into their destination
    void moveParticlesToDifferentCells(Reactor1d & reac, const inflow_lists_vector & inflow_lists) const;

    //! Object for setting reactor chemistry to fixed values (? just use pointer)
    ResetChemistry mResetChemistry;

    //! Length of maximum permitted stochastic jump
    real mMaxDt;

    //! Expected number of events deferred per particle between fixed updates
    real mDeferralRatio;

    //! Indicate if diffusion is to be split from the main particle processes
    bool mSplitDiffusion;

    //! Indicate if advection is to be split from the main particle processes
    bool mSplitAdvection;

}; //class PredCorrSolver
} //namespace Brush

#endif //BRUSH_PRED_CORR_SOLVER_H
