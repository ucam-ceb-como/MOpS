 /*
  Author(s):      Weerapong Phadungsukana (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    Header of Right hand side functions which are used by CVODES.

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

#ifndef MOPS_TEST_H
#define MOPS_TEST_H
// CVODE includes.
#include "nvector/nvector_serial.h"
#include "cvodes/cvodes.h"
#include "cvodes_impl.h"
#include "cvodes/cvodes_dense.h"
#include "cvodes_direct_impl.h"

// The right-hand side evaluator.  This function calculates the RHS of
// the reactor differential equations.  CVODE uses a void* pointer to
// allow the calling code to pass whatever information it wants to
// the RHS function.  In this case the void* pointer should be cast
// into an ODE_Solver object.
int rhsFn_CVODE(
    double t,      // Current flow time.
    N_Vector y,    // The current solution variables.
    N_Vector ydot, // Derivatives to return.
    void* solver   // An ODE_Solver object (to be cast).
    );

// The right-hand side evaluator.  This function calculates the RHS of
// the reactor differential equations. This function work exactly the same
// as rhsFn_CVODE but it allows CVODES to have access to problem parameters.
// This is vital for CVODES to use internal sensitivity Rhs estimator as we
// do not privide CVODES a function to evaluate. Computational time of this
// function will be higher than.
int rhsFn_CVODES(
    double t,      // Current flow time.
    N_Vector y,    // The current solution variables.
    N_Vector ydot, // Derivatives to return.
    void* solver   // An ODE_Solver object (to be cast).
    );

// The Jacobian matrix evaluator.  This function calculates the 
// Jacobian matrix given the current state.  CVODE uses a void* pointer to
// allow the calling code to pass whatever information it wants to
// the function.  In this case the void* pointer should be cast
// into an ODE_Solver object.
int jacFn_CVODE(
    long int N,    // Problem size.
    double t,      // Time.
    N_Vector y,    // Current solution variables.
    N_Vector ydot, // Current value of the vector f(t,y), the RHS.
    DlsMat J,    // Jacobian matrix.
    void* solver,  // An ODE_Solver object (to be cast).
    N_Vector tmp1, // Temporary array available for calculations.
    N_Vector tmp2, // Temporary array available for calculations.
    N_Vector tmp3  // Temporary array available for calculations.
    );

// The Jacobian matrix evaluator.  This function calculates the Jacobian matrix
// given the current state. This function work exactly the same
// as jacFn_CVODE but it allows CVODES to have access to problem parameters.
// This is vital for CVODES to use internal sensitivity Rhs estimator as we
// do not privide CVODES a function to evaluate
int jacFn_CVODES(
    long int N,    // Problem size.
    DlsMat J,    // Jacobian matrix.
    double t,      // Time.
    N_Vector y,    // Current solution variables.
    N_Vector ydot, // Current value of the vector f(t,y), the RHS.
    void* solver,  // An ODE_Solver object (to be cast).
    N_Vector tmp1, // Temporary array available for calculations.
    N_Vector tmp2, // Temporary array available for calculations.
    N_Vector tmp3  // Temporary array available for calculations.
    );
int rhsSensFn_CVODES(
    int Ns, realtype t,
    N_Vector y, N_Vector ydot,
    N_Vector *yS, N_Vector *ySdot,
    void *solver,
    N_Vector tmp1, N_Vector tmp2);

#endif
