/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ODE_Solver class wraps the CVODE ODE solver to solve single
    reactor models.

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

#ifndef MOPS_ODE_SOLVER_H
#define MOPS_ODE_SOLVER_H

#include "mops_params.h"
#include "mops_src_terms.h"
#include "mops_reactor.h"
#include "mops_gpc_sensitivity.h"

// CVODE includes.
#include "nvector/nvector_serial.h"
#include "cvodes_impl.h" // For CVodeMem.
#include "cvodes_dense_impl.h" // For DenseMat.

#include <istream>

#define ZERO  RCONST(0.0)

namespace Mops
{
class ODE_Solver
{
public:
    // Constructors.
    ODE_Solver(); // Default constructor.
    ODE_Solver(const ODE_Solver &copy);// Copy constructor.
    ODE_Solver(          // Stream-reading constructor.
        std::istream &in //   - Input stream.
        );

    // Destructor.
    ~ODE_Solver(void); // Default destructor.

    // Operators.
    ODE_Solver &operator=(const ODE_Solver &rhs);

    // Enumeration of ODE solvers.
    // enum SolverType {CVODE_Solver, RADAU5_Solver};

    
    // SOLVER SETUP.

    // Initialises the solver at the given time.
    void Initialise(const Reactor &reac);

    // Reset the solver.  Need to do this if the the reactor
    // contents has been changed between calls to Solve().
    void ResetSolver(void);

    // Reset the solver.  Need to do this if the the reactor
    // contents has been changed between calls to Solve().
    void ResetSolver(const Reactor &reac);

    // Sets the time in the ODE solver.
//    void SetTime(real time);


    // RUNNING THE SOLVER.

    // Solves the reactor equations up to the given time, assuming
    // that it is in future to the current time.
    void Solve(Reactor &reac, real stop_time);


    // ERROR TOLERANCES.

    // Returns the absolute error tolerance used for ODE
    // calculations.
    real ATOL() const;

    // Sets the absolute error tolerance used for ODE
    // calculations.
    void SetATOL(real atol);

    // Returns the relative error tolerance used for ODE
    // calculations.
    real RTOL() const;

    // Sets the relative error tolerance used for ODE
    // calculations.
    void SetRTOL(real rtol);


    // EXTERNAL SOURCE TERMS.

    // Returns the vector of external source terms.
    const SrcProfile *const ExtSrcTerms(void) const;

    // Sets the external source terms.
    void SetExtSrcTerms(const SrcProfile &src);

    // Returns the external source term function.
    SrcTermFnPtr ExtSrcTermFn(void) const;

    // Sets the source term function pointer.
    void SetExtSrcTermFn(SrcTermFnPtr fn);


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the solver object.
    ODE_Solver* Clone() const;

    // Writes the solver to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the solver data from a binary data stream.
    void Deserialize(std::istream &in);

    // RHS FUNCTION INTERFACE.
    // These function are implemented for being interface with Rhs functions of
    // CVODE and CVODES. Rhs functions are exported to mops_rhs_func.cpp/.h
    // so that the Rhs functions can be made non-static.

    // Return Reactor object.
    Reactor *GetReactor() { return m_reactor;}

    // Return a source term function pointer.
    SrcTermFnPtr GetsrcTermsFn() { return _srcTerms;};
    
    // Return a source term as SrcProfile object.
    const SrcProfile *GetsrcTerms() { return m_srcterms;};

    // Return number of equation.
    unsigned int GetNEquations() { return m_neq; }

    // Set sensitivity object by making a copy of given sensitivity object.
    void SetSensitivity(Mops::SensitivityAnalyzer &sensi) const {m_sensi = sensi;};

    // Get sensitivity object.
    Mops::SensitivityAnalyzer &GetSensitivity() const {return m_sensi;};

protected:
    // ODE solution variables.
    real m_rtol, m_atol;    // Relative and absolute tolerances.
    unsigned int m_neq;     // Number of equations solved.
//    unsigned int m_nsp;     // Number of species in current mechanism.
//    int m_iT;               // Index of temperature in solution vectors.
//    int m_iDens;            // Index of density in solution vectors.

    // Solution variables.
    real m_time;        // Current solution time.
    Reactor *m_reactor; // The reactor being solved.
    real *m_soln;       // Pointer to solution array (comes from Reactor object).
    real *m_deriv;      // Array to hold current solution derivatives.

    // External source terms.
    const SrcProfile *m_srcterms; // Vector of externally defined source terms on  the RHS.
    SrcTermFnPtr _srcTerms; // Source term function pointer.

    // Sensitivity related variables
    mutable Mops::SensitivityAnalyzer m_sensi;
    N_Vector *m_yS;


private:
    // CVODE variables.
    void *m_odewk;     // CVODE workspace.
    N_Vector m_solvec; // Internal solution array for CVODE interface.
    N_Vector m_yvec;   // Internal y work space for CVODE interface.


    // VERY IMPORTANT FOR CVODE INTEGRATION!


    // The Jacobian matrix evaluator.  This function calculates the 
    // Jacobian matrix given the current state.  CVODE uses a void* pointer to
    // allow the calling code to pass whatever information it wants to
    // the function.  In this case the void* pointer should be cast
    // into an ODE_Solver object.
    //static int jacFn_CVODE(
    //    long int N,    // Problem size.
    //    DenseMat J,    // Jacobian matrix.
    //    double t,      // Time.
    //    N_Vector y,    // Current solution variables.
    //    N_Vector ydot, // Current value of the vector f(t,y), the RHS.
    //    void* solver,  // An ODE_Solver object (to be cast).
    //    N_Vector tmp1, // Temporary array available for calculations.
    //    N_Vector tmp2, // Temporary array available for calculations.
    //    N_Vector tmp3  // Temporary array available for calculations.
    //    );


    // INITIALISATION AND DESTRUCTION.
    
    // Initialises the solver to the default state.
    void init(void);

    // Releases all memory used by the solver object.
    void releaseMemory(void);

    // Initialises the CVode ODE solver assuming that the
    // remainder of the the solver has been correctly set up.
    void InitCVode(void);

    // Copies the given CVode workspace into this ODE_Solver object.
    void assignCVMem(const CVodeMemRec &mem); 
};
};

#endif
