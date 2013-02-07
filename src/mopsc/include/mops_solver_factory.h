 /*!
  * @file   mops_solver_factory.h
  * @author William Menz
  * @brief  Declaration of the solver factory
  *
  *   About:
  *      A factory class used to create new solver objects, given a supplied
  *      reactor type.
  *
  *   Licence:
  *      mops is free software; you can redistribute it and/or
  *      modify it under the terms of the GNU Lesser General Public License
  *      as published by the Free Software Foundation; either version 2
  *      of the License, or (at your option) any later version.
  *
  *      This program is distributed in the hope that it will be useful,
  *      but WITHOUT ANY WARRANTY; without even the implied warranty of
  *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *      GNU Lesser General Public License for more details.
  *
  *      You should have received a copy of the GNU Lesser General Public
  *      License along with this program; if not, write to the Free Software
  *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  *      02111-1307, USA.
  *
  *   Contact:
  *      Prof Markus Kraft
  *      Dept of Chemical Engineering
  *      University of Cambridge
  *      New Museums Site
  *      Pembroke Street
  *      Cambridge
  *      CB2 3RA, UK
  *
  *      Email:       mk306@cam.ac.uk
  *      Website:     http://como.cheng.cam.ac.uk
  */

#ifndef MOPS_SOLVER_FACTORY_H_
#define MOPS_SOLVER_FACTORY_H_

#include "mops_solver.h"
#include "mops_simplesplit_solver.h"
#include "mops_simulator.h"
#include "mops_strang_solver.h"
#include "mops_predcor_solver.h"

namespace Mops {

enum SolverType {
    GPC,     // Gas-phase chemistry only, default.
    OpSplit, // Use simple operator splitting.
    Strang,  // Strang splitting.
    PredCor, // Split-Predictor---Split-Corrector.
    FlamePP, // Post-process a gas-phase profile (like sweep1).
    MoMIC,   // Method-of-moments for 1D particles.
    Marchenko, // Marchenko's particle transport method - no gas-phase.
};

namespace SolverFactory {
    Mops::Solver* Create(const Mops::SolverType soltype);
};

}

#endif /* MOPS_SOLVER_FACTORY_H_ */
