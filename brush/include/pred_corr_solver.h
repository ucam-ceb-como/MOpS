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

#include <vector>

// Forward declarations
namespace Brush {
    class Reactor1d;
}
namespace Sprog {
    class Mechanism;
}

namespace Sweep {
namespace Transport {
    struct TransportOutflow;
}
}

namespace Brush {

//! Methods for stepping forward in time.  Consider converting to namespace.
class PredCorrSolver {

public:
    //! Create with an object to specify the fixed chemistry
    PredCorrSolver(const ResetChemistry& reset_chem)
        : mResetChemistry(reset_chem)
        , mMaxDt(3.0e-4)
        , mDeferralRatio(10.0)
        {}

    //! Advance solution to specified time
    void solve(Reactor1d &reac, const real t_stop, const int n_steps, const int n_iter) const;

protected:
    //! Perform one time step using a predictor followed by some corrector iterations
    void predictorCorrectorStep(Reactor1d &reac, const real t_stop, const int n_iter) const;

    //! Advance particle population to specified time
    void solveParticles(Reactor1d &reac, const real t_stop) const;

    //! Perform one stochastic jump
    real particleTimeStep(Reactor1d &reac, const real t_stop) const;

    //! Advance chemistry over specified time interval
    void solveChemistry(Reactor1d & reac, const real t_stop) const; //,?some kind of workspace for ODE solver);

    //! Put a particle that has left one cell into its destination
    void transportIn(Reactor1d & reac, const size_t destination_index, const Sweep::Transport::TransportOutflow &particle_details) const;

private:
    //! Not possible to have a solver of this type without a reset chemistry object
    PredCorrSolver();

    //! Object for setting reactor chemistry to fixed values (? just use pointer)
    ResetChemistry mResetChemistry;

    //! Length of maximum permitted stochastic jump
    real mMaxDt;

    //! Expected number of events deferred per particle between fixed updates
    real mDeferralRatio;

}; //class PredCorrSolver
} //namespace Brush

#endif //BRUSH_PRED_CORR_SOLVER_H
