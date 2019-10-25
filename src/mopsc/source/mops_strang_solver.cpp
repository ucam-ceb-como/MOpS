/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the StrangSolver class declared in the
    mops_strang_solver.h header file.

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

#include "mops_strang_solver.h"
#include "mops_reactor_factory.h"
#include "sweep.h"
#include "string_functions.h"
#include "csv_io.h"
#include <stdexcept>

using namespace Mops;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
StrangSolver::StrangSolver(void) {}

// Copy constructor
StrangSolver::StrangSolver(const StrangSolver &sol)
: Mops::ParticleSolver(sol),
  Sweep::Solver(sol) {}

//! Clone the object
StrangSolver *const StrangSolver::Clone() const {
    return new StrangSolver(*this);
}

// Default destructor.
StrangSolver::~StrangSolver(void)
{
}


// SOLVING REACTORS.

// Solves the coupled reactor using a Strang splitting algorithm
// up to the stop time.  calls the output routine once at the
// end of the function.  niter is ignored.
void StrangSolver::Solve(Reactor &r, double tstop, int nsteps, int niter, 
                         Sweep::rng_type &rng,
                         OutFnPtr out, void *data)
{
    // Mark the time at the start of the step, in order to
    // calculate total computation time.
    clock_t totmark = clock();

    // Time counters.
    double t2 = r.Time();
    double dt = (tstop - t2) / (double)nsteps; // Step size.
    double h  = dt * 0.5; // Half step size.

    // Sweep time counters.
    double ts1 = r.Time();
    double ts2 = ts1;

    // Variables required to ensure particle number density is correctly
    // scaled with gas-phase expansion.
    double rho = 0.0;

    // This function stores heat capacity and particle density for the step
    storeTemperatureProperties(r, rng);

    m_cpu_mark = clock();
        // Solve first half-step of gas-phase chemistry.
        rho = r.Mixture()->GasPhase().MassDensity();
        m_ode.Solve(r, t2+=h);
        r.SetTime(t2);
    m_chemtime += calcDeltaCT(m_cpu_mark);

    m_cpu_mark = clock();
	
    // Update heat capacity and particle density for the step
    storeTemperatureProperties(r, rng);

    // Solve one whole step of population balance (Sweep).
    r.Mixture()->AdjustSampleVolume(rho / r.Mixture()->GasPhase().MassDensity());
    Run(ts1, ts2+=dt, *r.Mixture(), r.Mech()->ParticleMech(), rng);

    m_swp_ctime += calcDeltaCT(m_cpu_mark);

    for (int i=1; i<nsteps; ++i) {
        m_cpu_mark = clock();
        // Solve whole step of gas-phase chemistry.
        rho = r.Mixture()->GasPhase().MassDensity();
        m_ode.ResetSolver();
        m_ode.Solve(r, t2+=dt);
        r.SetTime(t2);
        m_chemtime += calcDeltaCT(m_cpu_mark);
        m_cpu_mark = clock();

        // Update heat capacity and particle density for the step
        storeTemperatureProperties(r, rng);

        // Solve whole step of population balance (Sweep).
        r.Mixture()->AdjustSampleVolume(rho / r.Mixture()->GasPhase().MassDensity());
        Run(ts1, ts2+=dt, *r.Mixture(), r.Mech()->ParticleMech(), rng);
        m_swp_ctime += calcDeltaCT(m_cpu_mark);
    }

    m_cpu_mark = clock();
    // Solve last half-step of gas-phase chemistry.  
    rho = r.Mixture()->GasPhase().MassDensity();
    m_ode.ResetSolver();
    m_ode.Solve(r, t2+=h);
    r.Mixture()->AdjustSampleVolume(rho / r.Mixture()->GasPhase().MassDensity());
    r.SetTime(t2);
    m_chemtime += calcDeltaCT(m_cpu_mark);

    // Calculate total computation time.
    m_tottime += calcDeltaCT(totmark);

    // Call the output function.
    if (out) out(nsteps, niter, r, *this, data);
}


// SOLUTION ROUTINES.

void StrangSolver::multiStrangStep(double dt, unsigned int n, Mops::Reactor &r,
                                   Sweep::rng_type &rng)
{
    // Time counters.
    double t2 = r.Time();
    double h  = dt * 0.5; // Half step size.

    // Sweep time counters.
    double ts1 = r.Time();
    double ts2 = ts1;

    // Variables required to ensure particle number density is correctly
    // scaled with gas-phase expansion.
    double rho = 0.0;

    m_cpu_mark = clock();
        // Solve first half-step of gas-phase chemistry.
        rho = r.Mixture()->GasPhase().MassDensity();
        m_ode.Solve(r, t2+=h);
        r.SetTime(t2);
    m_chemtime += calcDeltaCT(m_cpu_mark);

    m_cpu_mark = clock();

    // Solve one whole step of population balance (Sweep).
    r.Mixture()->AdjustSampleVolume(rho / r.Mixture()->GasPhase().MassDensity());
    Run(ts1, ts2+=dt, *r.Mixture(), r.Mech()->ParticleMech(), rng);

        m_swp_ctime += calcDeltaCT(m_cpu_mark);
    
    for (unsigned int i=1; i!=n; ++i) {
        m_cpu_mark = clock();
            // Solve whole step of gas-phase chemistry.
            rho = r.Mixture()->GasPhase().MassDensity();
            m_ode.ResetSolver();
            m_ode.Solve(r, t2+=dt);
            r.SetTime(t2);
        m_chemtime += calcDeltaCT(m_cpu_mark);

        m_cpu_mark = clock();
        // Solve whole step of population balance (Sweep).
        r.Mixture()->AdjustSampleVolume(rho / r.Mixture()->GasPhase().MassDensity());
        Run(ts1, ts2+=dt, *r.Mixture(), r.Mech()->ParticleMech(), rng);
        m_swp_ctime += calcDeltaCT(m_cpu_mark);
        
    }

    m_cpu_mark = clock();
        // Solve last half-step of gas-phase chemistry.    
        m_ode.ResetSolver();
        m_ode.Solve(r, t2+=h);
        r.SetTime(t2);
    m_chemtime += calcDeltaCT(m_cpu_mark);
}
