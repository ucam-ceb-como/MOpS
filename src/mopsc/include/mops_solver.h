/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Solver class holds simulation settings for mops and solves reactors.
    This basic solver only solves gas-phase chemistry equations, there is no
    operator splitting to solve gas-phase chemistry coupled to the particle
    population balance using sweep.

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

#ifndef MOPS_SOLVER_H
#define MOPS_SOLVER_H

#include "mops_params.h"
#include "mops_reactor.h"
#include "mops_timeinterval.h"
#include "mops_mechanism.h"
#include "mops_ode_solver.h"
#include "console_io.h"
#include "swp_gas_profile.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>

namespace Mops
{
class Solver
{
public:
    // A function pointer type definition for the solver 
    // output routine.
    typedef void (*OutFnPtr)(
        unsigned int,   // Current internal step number since last Solve() call.
        unsigned int,   // Iteration number for current step.
        const Reactor&, // Reactor being solved.
        const Solver&,  // Reference to current solver.
        void*           // User data object.
        );

    // Constructors.
    Solver(void); // Default constructor.

    //! Copy constructor
    Solver(const Mops::Solver &sol);

    //! Clone the object
    virtual Solver *const Clone() const;

    // Destructors.
    virtual ~Solver(void); // Default destructor.


    // SOLVER INITIALISATION AND RESET.

    // Initialises the solver to solve the given reactor.
    virtual void Initialise(Reactor &r);

    // Resets the solver to solve the given reactor.
    virtual void Reset(Reactor &r);


    // ERROR TOLERANCES FOR ODE SOLVER.

    // Returns the absolute error tolerance used for ODE
    // calculations.
    double ATOL() const;

    // Sets the absolute error tolerance used for ODE
    // calculations.
    void SetATOL(double atol);

    // Returns the relative error tolerance used for ODE
    // calculations.
    double RTOL() const;

    // Sets the relative error tolerance used for ODE
    // calculations.
    void SetRTOL(double rtol);

    // LOI STATUS FOR ODE SOLVER.

    //! Enables LOI status to true.
    void SetLOIStatusTrue();

    //! Sets LOI status to false
    void SetLOIStatusFalse();

    //! Retrieves the LOI status
    bool GetLOIStatus() const;

    //! Sets the cutoff value for LOI comparison
    void SetLOICompValue(double CompValue);

    //! Returns the LOI comparison value
    double ReturnCompValue() const;
    
    // UNDER-RELAXATION.

    // Returns the under-relaxation coefficient.
    double UnderRelaxCoeff(void) const;

    // Sets the under-relaxation coefficient.
    void SetUnderRelaxCoeff(double relax);

    // Calculates and stores various properties used to complete the 
    // energy balance so they can be computed less frequently. 
    void storeTemperatureProperties(
            Reactor &r,           // The reactor to solve.
            Sweep::rng_type &rng  // Random number generator
    );

    // SOLUTION.

    // Runs the solver for the given reactor, advancing it
    // to the given stop time.  The numerical parameters given
    // are the number of internal steps to take, and the number
    // of internal iterations.  Default values of <=0 will use
    // an adaptive method (NOT YET IMPLEMENTED).  Internal solver
    // output is provided after each step/iteration by passing
    // a function pointer.
    virtual void Solve(
            Reactor &r,   // The reactor to solve.
            double tstop,   // The end time for the step.
            int nsteps,   // Number of internal steps to take.
            int niter,    // Number of internal iterations to take.
            Sweep::rng_type &rng,  // Random number generator
            OutFnPtr out, // Output function pointer.
            void *data    // Custom data object which will be passed as argument to out().
        );
    
    //SENSITIVITY.

    //!Returns the double** m_sensitivity matrix with values calculated by CVODES.
    virtual double** GetSensSolution(int n_sensi, int n_species);

    //! Initialises the double** m_sensitivity matrix and sets to identity matrix
    virtual void InitialiseSensMatrix(int n_sensi, int n_species);

    //! Retrieves the number of sensitivity parameters in the problem.
    virtual unsigned int GetNumSens() const;

    //! Destroys the double** m_sensitivity matrix
    virtual void DestroySensMatrix(int n_species);

    //! Retrieves the number of equations solved by the ode solver
    unsigned int GetNEquations() const;

    //LOI KEPT SPECIES

    //! Adds the name of a user-defined kept species to a string vector.
    void AddKeptSpecies(std::string spec_name);

    //!Returns the string vector of user-defined species that need to be kept in the reaction.
    std::vector<std::string> ReturnKeptSpecies();


    // COMPUTATION TIME.
    
    // Returns the number of CT time variables tracked by this
    // solver type.
    virtual unsigned int CT_Count(void) const;

    // Outputs internal computation time data to the given
    // binary stream.
    virtual void OutputCT(std::ostream &out) const;

    // Adds the CT descriptions to a vector of strings.
    virtual void CT_Names(
        std::vector<std::string> &names, // Vector of CT names.
        unsigned int start=0 // Optional start index in vector.
        ) const;

    // Attach sensitivity to ODE_Solver by making copy.
    void AttachSensitivity(SensitivityAnalyzer &sensi) const;

    // Outputs sensitivity results to given file stream.
    void OutputSensitivity(std::fstream &fout, const Mops::Reactor &r, void *sim) const;

protected:
    // ODE SOLVER.
    //! The ODE solver used by the mops solver.
    ODE_Solver m_ode;


    // SOLVER SETTINGS.

    // Default error tolerances for the ODE solver.
    double m_atol, m_rtol;

    // SENSITIVITY SETTINGS

    //! Boolean variable for setting LOI status
    bool m_LOIEnable;

    //! Minimum possible LOI value to keep a species in a mechanism.
    double m_LOIComp;

    // Under-relaxation coefficient.
    double m_rlx_coeff;

    // LOI KEPT SPECIES

    //! String vector containing the names of kept species
    std::vector<std::string> Kept_Spec;

    // COMPUTATION TIME.

    clock_t m_cpu_start, m_cpu_mark;
    double m_tottime, m_chemtime;

    //! Calculates the time duration from a time mark to the current time.
    double calcDeltaCT(double markt) const;

};
};

#endif
