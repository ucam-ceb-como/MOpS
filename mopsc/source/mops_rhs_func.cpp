 /*
  Author(s):      Weerapong Phadungsukana (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    Implementation of Right hand side functions which are used by CVODES.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
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
    Dr Markus Kraft
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

#include "mops_rhs_func.h"
#include "mops_ode_solver.h"

// RHS FUNCTION AND GOVERNING EQUATIONS.

// The right-hand side evaluator.  This function calculates the RHS of
// the reactor differential equations.  CVODE uses a void* pointer to
// allow the calling code to pass whatever information it wants to
// the RHS function.  In this case the void* pointer should be cast
// into a Reactor object.
int rhsFn_CVODE(double t,      // Independent variable.
                            N_Vector y,    // Solution array.
                            N_Vector ydot, // Derivatives of y wrt t.
                            void* solver) // Pointer to ODE solver object.
{
    // Cast the Solver object.
    Mops::ODE_Solver *s = static_cast<Mops::ODE_Solver*>(solver);
    Mops::Reactor *r    = s->GetReactor();

    // Get the RHS from the system model.
    if (r->EnergyEquation() == Mops::Reactor::ConstT) {
        r->RHS_ConstT(t, NV_DATA_S(y), NV_DATA_S(ydot));
    } else {
        r->RHS_Adiabatic(t, NV_DATA_S(y), NV_DATA_S(ydot));
    }

    // Add the source terms, if defined.
    if (s->GetsrcTermsFn() != NULL) 
        s->GetsrcTermsFn()(NV_DATA_S(ydot), s->GetNEquations(), t, *s->GetsrcTerms());

    return 0;
};

// The right-hand side evaluator.  This function calculates the RHS of
// the reactor differential equations. This function work exactly the same
// as rhsFn_CVODE but it allows CVODES to have access to problem parameters.
// This is vital for CVODES to use internal sensitivity Rhs estimator as we
// do not privide CVODES a function to evaluate
int rhsFn_CVODES(double t,      // Current flow time.
                N_Vector y,     // The current solution variables.
                N_Vector ydot,  // Derivatives to return.
                void* solver)   // An ODE_Solver object (to be cast).
{
    // Cast the Solver object.
    Mops::ODE_Solver *s = static_cast<Mops::ODE_Solver*>(solver);
    Mops::Reactor *r    = s->GetReactor();

    s->GetSensitivity().ChangeMechParams();
    int rvalue = rhsFn_CVODE(t, y, ydot, solver);
    //s->GetSensitivity().ResetMechParams();
    return rvalue;
}
