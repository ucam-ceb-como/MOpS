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



////////////////// aab64 - Debugging purposes only //////////////////
#include <iostream>
#include <fstream>
#include <string>
/////////////////////////////////////////////////////////////////////



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
	// aab64 Make note of energy balance state
	if (r.EnergyEquation() != r.ConstT) {
		r.Mixture()->SetIsAdiabaticFlag(true);
	}
	else {
		r.Mixture()->SetIsAdiabaticFlag(false);
	}

	// Mark the time at the start of the step, in order to
	// calculate total computation time.
	clock_t totmark = clock();

	// Time counters.
	double t2 = r.Time();
	double dt = (tstop - t2) / (double)nsteps; // Step size.
	double h = dt * 0.5; // Half step size.

	// Sweep time counters.
	double ts1 = r.Time();
	double ts2 = ts1;

	// Variables required to ensure particle number density is correctly
	// scaled with gas-phase expansion.
	double rho = 0.0;

	m_cpu_mark = clock();
	// Solve first half-step of gas-phase chemistry.
	rho = r.Mixture()->GasPhase().MassDensity();
	m_ode.Solve(r, t2 += h);
	r.SetTime(t2);
	m_chemtime += calcDeltaCT(m_cpu_mark);

	m_cpu_mark = clock();

	// Solve one whole step of population balance (Sweep).
	if (!r.IsConstV()) r.Mixture()->AdjustSampleVolume(rho / r.Mixture()->GasPhase().MassDensity()); // aab64
	Run(ts1, ts2 += dt, *r.Mixture(), r.Mech()->ParticleMech(), rng);

	m_swp_ctime += calcDeltaCT(m_cpu_mark);

	for (int i = 1; i<nsteps; ++i) {

		m_cpu_mark = clock();
		// Solve whole step of gas-phase chemistry.
		rho = r.Mixture()->GasPhase().MassDensity();
		m_ode.ResetSolver();
		m_ode.Solve(r, t2 += dt);
		r.SetTime(t2);
		m_chemtime += calcDeltaCT(m_cpu_mark);

		m_cpu_mark = clock();
		// Solve whole step of population balance (Sweep).
		if (!r.IsConstV()) r.Mixture()->AdjustSampleVolume(rho / r.Mixture()->GasPhase().MassDensity()); // aab64
		Run(ts1, ts2 += dt, *r.Mixture(), r.Mech()->ParticleMech(), rng);
		m_swp_ctime += calcDeltaCT(m_cpu_mark);
	}

	m_cpu_mark = clock();
	// Solve last half-step of gas-phase chemistry.    
	m_ode.ResetSolver();
	m_ode.Solve(r, t2 += h);
	r.SetTime(t2);
	m_chemtime += calcDeltaCT(m_cpu_mark);

	// Calculate total computation time.
	m_tottime += calcDeltaCT(totmark);

	// Call the output function.
	if (out) out(nsteps, niter, r, *this, data);
}



//////////////////////////////////////////// aab64 ////////////////////////////////////////////
// Should really replace the previous function so that duplicates do not exist
// Writes diagnostics for process events per split

// Solves the coupled reactor using a Strang splitting algorithm
// up to the stop time.  calls the output routine once at the
// end of the function.  niter is ignored.
void StrangSolver::Solve(Reactor &r, double tstop, int nsteps, int niter, 
                         Sweep::rng_type &rng,
                         OutFnPtr out, void *data, bool writediags)
{
	// aab64 Make note of energy balance state
	if (r.EnergyEquation() != r.ConstT) {
		r.Mixture()->SetIsAdiabaticFlag(true);
	}
	else {
		r.Mixture()->SetIsAdiabaticFlag(false);
	}

	//Diagnostics file
    ofstream partProcFile, gasConcFile;
	
    // Diagnostic variables
	double tmpSVin, tmpSVout, tmpWtVarin, tmpWtVarout, tmpWtMassin, tmpWtMassout, tmpIncWeightin, tmpIncWeightout, tmpIncFactorin, tmpIncFactorout;
	double tmpTin, tmpTout;
    unsigned int tmpSPin, tmpSPout, tmpAddin, tmpAddout, tmpInfin, tmpInfout, tmpOutfin, tmpOutfout;
    unsigned int process_iter;
    std::vector<unsigned int> tmpPCin, tmpPCout, tmpFCin, tmpFCout;
    Sprog::fvector tmpGPin, tmpGPout;

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

    m_cpu_mark = clock();
    // Solve first half-step of gas-phase chemistry.
    rho = r.Mixture()->GasPhase().MassDensity();
    m_ode.Solve(r, t2+=h);
    r.SetTime(t2);
    m_chemtime += calcDeltaCT(m_cpu_mark);

    m_cpu_mark = clock();

	if (writediags) {
		// Diagnostics variables at start of split step
		tmpSVin = r.Mixture()->SampleVolume();
		tmpSPin = r.Mixture()->ParticleCount();
		tmpPCin = r.Mech()->ParticleMech().GetProcessUsageCounts();
		tmpFCin = r.Mech()->ParticleMech().GetFictitiousProcessCounts();
		tmpAddin = r.Mech()->ParticleMech().GetDeferredAddCount();
		tmpInfin = r.Mech()->ParticleMech().GetInflowCount();
		tmpOutfin = r.Mech()->ParticleMech().GetOutflowCount();
		tmpWtVarin = r.Mixture()->Particles().GetSum(Sweep::iW);
		tmpWtMassin = r.Mixture()->Particles().GetSum(Sweep::iWM);
		tmpIncFactorin = r.Mixture()->GetInceptionFactor();
		tmpIncWeightin = r.Mixture()->GetInceptingWeight();
		r.Mixture()->GasPhase().GetConcs(tmpGPin);
		tmpTin = r.Mixture()->GasPhase().Temperature();
	}

    // Solve one whole step of population balance (Sweep).
    if (!r.IsConstV()) r.Mixture()->AdjustSampleVolume(rho / r.Mixture()->GasPhase().MassDensity()); //aab64
    Run(ts1, ts2+=dt, *r.Mixture(), r.Mech()->ParticleMech(), rng);

    m_swp_ctime += calcDeltaCT(m_cpu_mark);

	if (writediags) {
		// Diagnostic variables at end of split step
		tmpSVout = r.Mixture()->SampleVolume();
		tmpSPout = r.Mixture()->ParticleCount();
		tmpPCout = r.Mech()->ParticleMech().GetProcessUsageCounts();
		tmpFCout = r.Mech()->ParticleMech().GetFictitiousProcessCounts();
		tmpAddout = r.Mech()->ParticleMech().GetDeferredAddCount();
		tmpInfout = r.Mech()->ParticleMech().GetInflowCount();
		tmpOutfout = r.Mech()->ParticleMech().GetOutflowCount();
		tmpWtVarout = r.Mixture()->Particles().GetSum(Sweep::iW);
		tmpWtMassout = r.Mixture()->Particles().GetSum(Sweep::iWM);
		tmpIncFactorout = r.Mixture()->GetInceptionFactor();
		tmpIncWeightout = r.Mixture()->GetInceptingWeight();
		r.Mixture()->GasPhase().GetConcs(tmpGPout);
		tmpTout = r.Mixture()->GasPhase().Temperature();
		std::string rname(r.GetName());
		std::string partfname, chemfname;
		partfname = "Part-split-diagnostics(" + rname + ").csv";
		chemfname = "Chem-split-diagnostics(" + rname + ").csv";

		// Output particle diagnostics to file
		partProcFile.open(partfname.c_str(), ios::app);
		partProcFile << ts2 << " , " << tstop << " , " << 0 << " , "
			<< tmpSVin << " , " << tmpSVout << " , "
			<< tmpSPin << " , " << tmpSPout << " , "
			<< tmpWtVarin << " , " << tmpWtVarout << " , "
			<< tmpWtMassin << " , " << tmpWtMassout << " , "
			<< tmpIncWeightin << " , " << tmpIncWeightout << " , "
			<< tmpIncFactorin << " , " << tmpIncFactorout << " , ";
		for (process_iter = 0; process_iter < r.Mech()->ParticleMech().Inceptions().size();
			process_iter++) {
			partProcFile << tmpPCout[process_iter] - tmpPCin[process_iter] << " , ";
		}
		if (r.Mech()->ParticleMech().AnyDeferred()) {
			    partProcFile << tmpAddout - tmpAddin << " , ";
		}
		else {
		    partProcFile << tmpPCout[process_iter] - tmpPCin[process_iter] << " , ";
		}
		for (process_iter = r.Mech()->ParticleMech().Inceptions().size() + 1;
			process_iter < tmpPCin.size(); process_iter++) {
			partProcFile << tmpPCout[process_iter] - tmpPCin[process_iter] << " , ";
		}
		for (process_iter = r.Mech()->ParticleMech().Inceptions().size() + 1;
			process_iter < tmpPCin.size(); process_iter++) {
			partProcFile << tmpFCout[process_iter] - tmpFCin[process_iter] << " , ";
		}
		partProcFile << tmpInfout - tmpInfin << " , " << tmpOutfout - tmpOutfin << "\n";
		partProcFile.close();
		
		// Output gasphase diagnostics to file
		gasConcFile.open(chemfname.c_str(), ios::app);
		gasConcFile << ts2 << " , " << tstop << " , " << 0 << " , ";
		for (process_iter = 0; process_iter < tmpGPin.size(); process_iter++) {
			gasConcFile << tmpGPin[process_iter] << " , " << tmpGPout[process_iter] << " , ";
		}
		gasConcFile << tmpTin << " , " << tmpTout << " , ";
		gasConcFile << "\n";
		gasConcFile.close();
	}

	for (int i=1; i<nsteps; ++i) {

	    // Diagnostics variables at start of split step
	    if (writediags) {
		    tmpSVin = r.Mixture()->SampleVolume();
		    tmpSPin = r.Mixture()->ParticleCount();
		    tmpPCin = r.Mech()->ParticleMech().GetProcessUsageCounts();
		    tmpFCin = r.Mech()->ParticleMech().GetFictitiousProcessCounts();
		    tmpAddin = r.Mech()->ParticleMech().GetDeferredAddCount();
		    tmpInfin = r.Mech()->ParticleMech().GetInflowCount();
			tmpOutfin = r.Mech()->ParticleMech().GetOutflowCount();
			tmpWtVarin = r.Mixture()->Particles().GetSum(Sweep::iW);
			tmpWtMassin = r.Mixture()->Particles().GetSum(Sweep::iWM);
			tmpIncFactorin = r.Mixture()->GetInceptionFactor();
			tmpIncWeightin = r.Mixture()->GetInceptingWeight();
			r.Mixture()->GasPhase().GetConcs(tmpGPin);
			tmpTin = r.Mixture()->GasPhase().Temperature();
	    }

	    m_cpu_mark = clock();
        // Solve whole step of gas-phase chemistry.
        rho = r.Mixture()->GasPhase().MassDensity();
        m_ode.ResetSolver();
        m_ode.Solve(r, t2+=dt);
        r.SetTime(t2);
        m_chemtime += calcDeltaCT(m_cpu_mark);

        m_cpu_mark = clock();
        // Solve whole step of population balance (Sweep).
	    if (!r.IsConstV()) r.Mixture()->AdjustSampleVolume(rho / r.Mixture()->GasPhase().MassDensity()); // aab64
	    Run(ts1, ts2+=dt, *r.Mixture(), r.Mech()->ParticleMech(), rng);
        m_swp_ctime += calcDeltaCT(m_cpu_mark);  

		if (writediags) {
			// Diagnostic variables at end of split step
			tmpSVout = r.Mixture()->SampleVolume();
			tmpSPout = r.Mixture()->ParticleCount();
			tmpPCout = r.Mech()->ParticleMech().GetProcessUsageCounts();
			tmpFCout = r.Mech()->ParticleMech().GetFictitiousProcessCounts();
			tmpAddout = r.Mech()->ParticleMech().GetDeferredAddCount();
			tmpInfout = r.Mech()->ParticleMech().GetInflowCount();
			tmpOutfout = r.Mech()->ParticleMech().GetOutflowCount();
			tmpWtVarout = r.Mixture()->Particles().GetSum(Sweep::iW);
			tmpWtMassout = r.Mixture()->Particles().GetSum(Sweep::iWM);
			tmpIncFactorout = r.Mixture()->GetInceptionFactor();
			tmpIncWeightout = r.Mixture()->GetInceptingWeight();
			r.Mixture()->GasPhase().GetConcs(tmpGPout);
			tmpTout = r.Mixture()->GasPhase().Temperature();
			std::string rname(r.GetName());
			std::string partfname, chemfname;
			partfname = "Part-split-diagnostics(" + rname + ").csv";
			chemfname = "Chem-split-diagnostics(" + rname + ").csv";

			// Output particle diagnostics to file
			partProcFile.open(partfname.c_str(), ios::app);
			partProcFile << ts2 << " , " << tstop << " , " << i << " , "
				<< tmpSVin << " , " << tmpSVout << " , "
				<< tmpSPin << " , " << tmpSPout << " , "
			    << tmpWtVarin << " , " << tmpWtVarout << " , "
			    << tmpWtMassin << " , " << tmpWtMassout << " , "
				<< tmpIncWeightin << " , " << tmpIncWeightout << " , "
				<< tmpIncFactorin << " , " << tmpIncFactorout << " , ";
			for (process_iter = 0; process_iter < r.Mech()->ParticleMech().Inceptions().size();
				process_iter++) {
				partProcFile << tmpPCout[process_iter] - tmpPCin[process_iter] << " , ";
			}
			if (r.Mech()->ParticleMech().AnyDeferred()) {
			    partProcFile << tmpAddout - tmpAddin << " , ";
		    }
		    else {
			    partProcFile << tmpPCout[process_iter] - tmpPCin[process_iter] << " , ";
		    }
			for (process_iter = r.Mech()->ParticleMech().Inceptions().size() + 1;
				process_iter < tmpPCin.size(); process_iter++) {
				partProcFile << tmpPCout[process_iter] - tmpPCin[process_iter] << " , ";
			}
			for (process_iter = r.Mech()->ParticleMech().Inceptions().size() + 1;
			    process_iter < tmpPCin.size(); process_iter++) {
			    partProcFile << tmpFCout[process_iter] - tmpFCin[process_iter] << " , ";
			}
			partProcFile << tmpInfout - tmpInfin << " , " << tmpOutfout - tmpOutfin << "\n";
			partProcFile.close();
			
		    // Output gasphase diagnostics to file
		    gasConcFile.open(chemfname.c_str(), ios::app);
		    gasConcFile << ts2 << " , " << tstop << " , " << i << " , ";
		    for (process_iter = 0; process_iter < tmpGPin.size(); process_iter++) {
				gasConcFile << tmpGPin[process_iter] << " , " << tmpGPout[process_iter] << " , ";
			}
			gasConcFile << tmpTin << " , " << tmpTout << " , ";
		    gasConcFile << "\n";
		    gasConcFile.close();
		}
    }

    m_cpu_mark = clock();
    // Solve last half-step of gas-phase chemistry.    
    m_ode.ResetSolver();
    m_ode.Solve(r, t2+=h);
    r.SetTime(t2);
    m_chemtime += calcDeltaCT(m_cpu_mark);

    // Calculate total computation time.
    m_tottime += calcDeltaCT(totmark);

    // Call the output function.
    if (out) out(nsteps, niter, r, *this, data);
}
//////////////////////////////////////////// aab64 ////////////////////////////////////////////



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
