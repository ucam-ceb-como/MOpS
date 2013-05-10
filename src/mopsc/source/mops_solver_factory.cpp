 /*!
  * @file   mops_solver_factory.cpp
  * @author William Menz
  * @brief  Implementation of solver factory
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

#include "mops_solver_factory.h"

namespace Mops {

Mops::Solver* SolverFactory::Create(const Mops::SolverType soltype) {
    Mops::Solver* solver;

    switch (soltype) {
    case Mops::OpSplit:
        solver = new Mops::SimpleSplitSolver();
        break;
    case Mops::Strang:
        solver = new Mops::StrangSolver();
        break;
    case Mops::PredCor:
        solver = new Mops::PredCorSolver();
        break;
    case Mops::FlamePP:
        solver = new Sweep::FlameSolver();
        break;
    case Mops::GPC:
    default:
        solver = new Solver();
        break;
    }

    return solver;
}

}
