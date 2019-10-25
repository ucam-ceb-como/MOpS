/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The PredCorSolver class implements a split-predictor---split-corrector
    algorithm to couple the gas-phase and particle systems.  The principle is
    to include an approximation of the change to the gas-phase due to 
    particle processes in the calculation of the gas-phase ODEs.  In this
    way the ODE solver does not need to be restarted at each splitting step,
    which has a significant advantage in terms of run-time.

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

#ifndef MOPS_PREDCOR_SOLVER_H
#define MOPS_PREDCOR_SOLVER_H

#include "mops_params.h"
#include "swp_flamesolver.h"
#include "swp_gas_profile.h"
#include "mops_reactor.h"
#include "mops_timeinterval.h"
#include "mops_mechanism.h"
#include "mops_src_terms.h"
#include "sweep.h"
#include "sprog.h"
#include <string>

namespace Mops
{
class PredCorSolver : public Sweep::FlameSolver
{
public:
    // Constructors.
    PredCorSolver(void); // Default constructor.

    //! Copy constructor
    PredCorSolver(const PredCorSolver &sol);

    //! Clone the object
    PredCorSolver *const Clone() const;

    // Destructors.
    ~PredCorSolver(void); // Default destructor.


    // SOLVER INITIALISATION.

    // Initialises the solver to solve the given reactor.
    virtual void Initialise(Reactor &r);

    // Resets the solver to solve the given reactor.
    virtual void Reset(Reactor &r);


    // SOLUTION.

    // Solves the coupled reactor using the predictor-corrector splitting
    // algorithm up to the stop time.  Calls the output function after
    // each iteration of the last internal step.
    virtual void Solve(
            Reactor &r,   // The reactor to solve.
            double tstop,   // The end time for the step.
            int nsteps,   // Number of internal steps to take.
            int niter,    // Number of internal iterations to take.
            Sweep::rng_type &rng,  // Random number generator
            OutFnPtr out, // Output function pointer.
            void *data    // Custom data object which will be passed as argument to out().
        );

    /*
    // Run the solver for the given reactor and the 
    // given time intervals.
    void SolveReactor(
        Reactor &r,              // Reactor object to solve.
        const timevector &times, // Vector of time intervals.
        unsigned int nruns = 1   // Number of runs to perform.
        );

    // Post-processes binary output files with the given file name
    // into CSV files.
    void PostProcess(
        const std::string &filename, // Filename to post-process.
        unsigned int nruns = 1       // Number of runs.
        ) const;
*/
/*
    // UNDER-RELAXATION.

    // Returns the under-relaxation coefficient.
    double UnderRelaxCoeff(void) const;

    // Sets the under-relaxation coefficient.
    void SetUnderRelaxCoeff(double relax);
*/

    // SOURCE TERM CALCULATION (REQUIRED FOR REACTOR ODE SOLVER).

    // Adds the source term contribution to the RHS supplied by the
    // reactor class.
    static void AddSourceTerms(
        double *rhs,             // ODE right-hand sides (to be updated)
        unsigned int n,        // Number of values in rhs array.
        double t,                // Current time.
        const SrcProfile &prof // Source term profile.
        );

private:
    // The source terms which describe the effect of the soot processes
    // on the gas-phase conditions.
    SrcProfile m_srcterms, m_srcterms_copy;

    // A copy of the operating reactor.
    Reactor *m_reac_copy;

    // Under-relaxation coefficient.
//    double m_rlx_coeff;

    // A copy of the ODE solver.
    ODE_Solver m_ode_copy;

    // Total number of times Solve() has been called since last
    // reset.
    unsigned int m_ncalls;


    // SIMULATION.

    // Sets up the workspace to begin a new iteration.  This includes
    // generating an initial chemistry profile over the step.  As a first
    // guess the chemistry is assumed to be constant in this interval.  The
    // initial source terms from sweep are then calculated using this chemistry
    // profile.
    void beginIteration(
        Reactor &r,        // The reactor which is being solved.
        unsigned int step, // The step number.
        double dt            // The size of the next time step.
        );

    // Performs a step-wise iteration on the reactor to recalculate
    // the source terms for the gas-phase effect on the particle model.
    void iteration(
        Reactor &r, // Reactor to solve.
        double dt,    // Time step size.
        Sweep::rng_type &rng  // Random number generator
        );

    // Terminates an iteration step.
    void endIteration();

    // Generates a chemistry profile over the required time interval
    // by solving the gas-phase chemistry equations with the current
    // source terms.
    void generateChemProfile(
        Reactor &r, // Reactor to solve.
        double dt     // Time step size over which to calculate profile.
        );

    // Calculates the instantaneous source terms.
    void calcSrcTerms(
        SrcPoint &src,   // Source term point to fill.
        const Reactor &r // Reactor used to get conditions and calculate terms.
        );

    // Calculate the adiabatic temperature change rate due to the
    // species source terms.
    double energySrcTerm(
        const Reactor &r, // Reactor used to calculate terms.
        fvector &src      // Vector of species source terms.
        );

    // Extrapolates the source terms over the given time interval.  Two
    // points are given, from which the gradient is calculated, which 
    // allows the third point to be calculated.
    static void linExSrcTerms(
        SrcPoint &src,          // The point to which to extrapolate.
        const SrcProfile &prof, // The profile of previous points used to find gradient.
        double dt                 // The delta-t from the last profile point to the new point.
        );

    // Applies under-relaxation to the first source point, using the
    // second source point as the initial values.
    static void relaxSrcTerms(
        SrcPoint &src,        // Source terms to relax.
        const SrcPoint &init, // Initial source terms.
        double rcoeff           // Under-relaxation coefficient.
        );

    // Applies under-relaxation to the first source point, using the
    // second source point as the initial values.
    static void copySrcPoint(
        const SrcPoint &from, // The source point from which to copy.
        SrcPoint &to          // The source point to overwrite.
        );

/*
    // POST-PROCESSING.

    // Multiplies all values in a vector by a scaling factor.
    static void multVals(
        fvector &vals, // The values to multiply by the scaling factor.
        double scale     // The scaling factor (numner of runs).
        );
*/
};
};

#endif
