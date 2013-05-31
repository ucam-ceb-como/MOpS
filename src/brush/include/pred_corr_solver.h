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

} // namespace Sweep


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
                   const unsigned corrector_iterations,
                   const double rtol, const double atol,
                   const bool split_diffusion, const double drift_correction,
                   const bool split_advection, const bool strang_splitting,
                   const bool cstr_transport);

    //! Advance solution to specified time
    void solve(Reactor1d &reac, const double t_start, const double t_stop, const int n_steps,
               const size_t seed) const;

protected:
    //! Perform one time step using a predictor followed by some corrector iterations
    void predictorCorrectorStep(Reactor1d &reac, const double t_start, const double t_stop,
                                std::vector<Sweep::rng_type>& cell_rngs) const;

    //! Advance particle population to specified time
    void solveParticlesByCell(Reactor1d &reac, const double t_start, const double t_stop,
                              std::vector<Sweep::rng_type>& cell_rngs) const;

    //! Advance particle population to specified time
    void solveParticlesInOneCell(Sweep::Cell &cell, const Geometry::LocalGeometry1d &geom,
                                 const Sweep::Mechanism &mech, double t, double t_stop,
                                 Sweep::rng_type &rng) const;

    //! Advance chemistry over specified time interval
    void solveChemistry(Reactor1d & reac, const double t_start, const double t_stop) const; //,?some kind of workspace for ODE solver);

    //! Carry out split transport on all particles from all cells
    void splitParticleTransport(Reactor1d &reac, const double t_start, const double t_stop,
                                std::vector<Sweep::rng_type>& cell_rngs) const;

private:
    //! Not possible to have a solver of this type without a reset chemistry object
    PredCorrSolver();

    //! Calculate and set a new position on one particle
    void updateParticlePosition(const double t_start, const double t_stop, const Sweep::Cell &mix,
                                const Sweep::Mechanism &mech,
                                const Geometry::LocalGeometry1d & geom,
                                const std::vector<const Sweep::Cell*> & neighbouring_cells,
                                Sweep::Particle& sp, Sweep::rng_type &rng) const;
    
    //! Move particle to centre of this or next cell according to CSTR probabilities
    void updateParticlePositionCSTR(const double t_start, const double t_stop, const Sweep::Cell &mix,
                                    const Sweep::Mechanism &mech,
                                    const Geometry::LocalGeometry1d & geom,
                                    const std::vector<const Sweep::Cell*> & neighbouring_cells,
                                    Sweep::Particle& sp, Sweep::rng_type &rng) const;

    /*!
     *  \brief Lists of particles due to be added to cells
     *
     *  Particles in the the lists at index i in the vector are destined to be added to
     *  cell i of a reactor.
     */
    typedef std::vector<std::list<Sweep::Particle*> > inflow_lists_vector;

    //! Take the particles from a cell and put them in lists according to their new positions
    inflow_lists_vector updateParticleListPositions(const double t_start, const double t_stop, Sweep::Cell &mix,
                                                    const size_t cell_index, const Sweep::Mechanism &mech,
                                                    const Geometry::Geometry1d &geom,
                                                    const std::vector<const Sweep::Cell*> &neighbouring_cells,
                                                    Sweep::rng_type &rng) const;


    //! Object for setting reactor chemistry to fixed values (? just use pointer)
    ResetChemistry mResetChemistry;

    //! Number of corrector iterations (there will always be a predictor)
    unsigned mCorrectorIterations;

    //! Expected number of events deferred per particle between fixed updates
    double mDeferralRatio;

    //! Indicate if diffusion is to be split from the main particle processes
    bool mSplitDiffusion;

    //! This is a stochastic calculus technicality, a value of 1 is probably a safe choice
    /*!
     * The kind of stochastic integration to use when defining the diffusion processes followed by
     * particle is not completely clear.  See for example section II of Schnitzer (93)
     * http://dx.doi.org/10.1103/PhysRevE.48.2553 and also
     * Bach & Duerr (78) http://dx.doi.org/10.1016/0375-9601(78)90001-4, where this parameter
     * is the \f$\lambda\f$ that appears in equation (2).
     *
     * Diffusing particles have trajectories that are modelled by the SDE
     * \f[
     *    dY_t = \mu(Y_t)dt + \sigma(Y_t)dW_t
     *  \f]
     * where \f$W_t\f$ is a standard Brownian Motion.  The question arises as
     * to whether to interpret the \f$\sigma(Y_t)dW_t\f$ in an Ito, Stratonovich
     * or other sense.  This is significant, because the numerical method has to
     * be chosen to match the model.  The different models lead to different
     * equations for the particle concentrations.  The different interpretations
     * of \f$\sigma(Y_t)dW_t\f$ can all be included in an SDE using the Ito
     * interpretation by adjusting the drift \f$\mu\f$ and that is the approach
     * currently taken in the code since Ito diffusion is easiest to simulate.
     * The adjustment is \f$\mu(x) + \lambda \sigma^\prime(x)\sigma(x) \f$.
     *
     * The case \f$\lambda = 0\f$ is the Ito interpretation and gives rise to the
     * following equation for particle concentration (\f$x\f$ is particle type,
     * which can be ignored) during transport:
     * \f[
     *    \partial_t c(t,x,y) = \frac{1}{2} \partial^2_y \left(\sigma(x)^2 c(t,x,y) \right).
     * \f]
     * This is the form that would occur first to a mathematician with an interest
     * in stochastic calculus.
     *
     * The case \f$\lambda = \frac{1}{2}\f$ is the Stratonovich interpretation and gives rise to the
     * following equation for particle concentration (\f$x\f$ is particle type,
     * which can be ignored) during transport:
     * \f[
     *    \partial_t c(t,x,y) = \frac{1}{2} \partial_y \left(\sigma(x) \partial_y\left(\sigma(y) c(t,x,y)\right) \right).
     * \f]
     * The books that I have read that come from the physics direction and do not
     * focus on stochastic calculus as maths tend to claim Stratonovich is the way
     * forward in many physical settings, although one should note that one does
     * not recover Fick's law.
     *
     * The case \f$\lambda = \frac{1}{2}\f$ is the Isothermal interpretation and gives rise to the
     * following equation for particle concentration (\f$x\f$ is particle type,
     * which can be ignored) during transport:
     * \f[
     *    \partial_t c(t,x,y) = \frac{1}{2} \partial_y \left(\sigma(x)^2 \partial_y c(t,x,y) \right).
     * \f]
     * This does recover Fick's law.
     *
     * Note that the diffusion coefficient \f$D(x) = \frac{1}{2}\sigma(x)^2\f$.
     */
    double mDiffusionDriftAdjustment;

    //! Indicate if advection is to be split from the main particle processes or simply ignored
    bool mSplitAdvection;

    //! True if transport  splitting should be Strang (second order), otherwise first order splitting is used
    bool mStrangTransportSplitting;

    //! True if transport should be CSTR like random jumps between cell centres.
    bool mCSTRTransport;

}; //class PredCorrSolver
} //namespace Brush

#endif //BRUSH_PRED_CORR_SOLVER_H
