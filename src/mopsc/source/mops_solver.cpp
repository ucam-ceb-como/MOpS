/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Solver class declared in the
    mops_solver.h header file.

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

#include "mops_solver.h"
#include "mops_reactor_factory.h"
#include "csv_io.h"
#include "string_functions.h"
#include <vector>
#include <string>
#include <time.h>
#include <stdexcept>

using namespace Mops;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Solver::Solver(void)
: m_atol(1.0e-3), m_rtol(6.0e-4),
  m_LOIEnable(false), m_rlx_coeff(0.0),
  m_cpu_start((clock_t)0.0), m_cpu_mark((clock_t)0.0),
  m_tottime(0.0), m_chemtime(0.0)
{
}

// Default destructor.
Solver::~Solver(void)
{
}

// Copy constructor
Solver::Solver(const Mops::Solver &sol)
: m_ode(sol.m_ode),
  m_atol(sol.m_atol),
  m_rtol(sol.m_rtol),
  m_LOIEnable(sol.m_LOIEnable),
  m_LOIComp(sol.m_LOIComp),
  m_rlx_coeff(sol.m_rlx_coeff),
  Kept_Spec(sol.Kept_Spec),
  m_cpu_start(sol.m_cpu_start),
  m_cpu_mark(sol.m_cpu_mark),
  m_tottime(sol.m_tottime),
  m_chemtime(sol.m_chemtime)
{}

//! Clone the object
Solver *const Solver::Clone() const {
    return new Solver(*this);
}

// SOLVER INITIALISATION AND RESET.

// Initialises the solver to solve the given reactor.
void Solver::Initialise(Reactor &r)
{
    // Set up the ODE solver.
    m_ode.Initialise(r);
    m_ode.SetATOL(m_atol);
    m_ode.SetRTOL(m_rtol);
}

// Resets the solver to solve the given reactor.
void Solver::Reset(Reactor &r)
{
    // Reset the ODE solver.
    m_ode.ResetSolver(r);
    m_ode.SetATOL(m_atol);
    m_ode.SetRTOL(m_rtol);
}


// ERROR TOLERANCES.

double Solver::ATOL() const
{
    return m_atol;
}

void Solver::SetATOL(double atol)
{
    m_atol = atol;
    m_ode.SetATOL(atol);
}

double Solver::RTOL() const
{
    return m_rtol;
}

void Solver::SetRTOL(double rtol)
{
    m_rtol = rtol;
    m_ode.SetRTOL(rtol);
}

/*!
Sets the solver status to true
*/
void Solver::SetLOIStatusTrue()
{
    m_LOIEnable = true;
}

/*!
Sets the solver status to false
*/
void Solver::SetLOIStatusFalse()
{
    m_LOIEnable = false;
}

/*!
Retrieves the LOI status
@param[out]     LOIEnable    Boolean status of LOI reduction
*/
bool Solver::GetLOIStatus() const
{
    return m_LOIEnable;
}

/*!
Sets the LOI comparison value
@param[in]      CompValue   Minimum value of a species LOI that will be kept in a reduced mechanism
*/
void Solver::SetLOICompValue(double CompValue)
{
    m_LOIComp = CompValue;
}

/*!
Returns the LOI comparison value
return      m_LOIComp       Minimum value of a species LOI that will be kept in a reduced mechanism
*/
double Solver::ReturnCompValue() const
{
    return m_LOIComp;
}

/*!
@param[in]      String name of a species to be kept by the LOI Method
*/
void Solver::AddKeptSpecies(std::string spec_name)
{
    Kept_Spec.push_back(spec_name);
}

/*!
@return     Kept_Spec       The string vector of user-defined species to be kept in the mechanism.
*/
std::vector<std::string> Solver::ReturnKeptSpecies()
{
    return Kept_Spec;
}

// UNDER-RELAXATION.

// Returns the under-relaxation coefficient.
double Solver::UnderRelaxCoeff(void) const
{
    return m_rlx_coeff;
}

// Sets the under-relaxation coefficient.
void Solver::SetUnderRelaxCoeff(double relax) {m_rlx_coeff = relax;}


// SOLVING REACTORS.

// Runs the solver for the given reactor, advancing it
// to the given stop time.  The numerical parameters given
// are the number of internal steps to take, and the number
// of internal iterations.  Default values of <=0 will use
// an adaptive method (NOT YET IMPLEMENTED).  Internal solver
// output is provided after each step/iteration by passing
// a function pointer.
void Solver::Solve(Reactor &r, double tstop, int nsteps, int niter,
                   Sweep::rng_type &rng, OutFnPtr out, void *data)
{
    // Mark the current time.
    m_cpu_mark = clock();
    // Solve reactor.
    m_ode.Solve(r, tstop);
    r.SetTime(tstop);
    // Calculate CPU time.
    double dt = calcDeltaCT(m_cpu_mark);
    m_tottime += dt;
    m_chemtime += dt;
    // Perform output.
    if (out) out(nsteps, niter, r, *this, data);
}

// Calculates and stores various properties used to complete the 
// energy balance with particles so they can be computed less frequently. 
void Solver::storeTemperatureProperties(Reactor &r, Sweep::rng_type &rng)
{
	// Check if particle terms are to be included in the energy balance
	if (r.IncludeParticles())
	{
		Sprog::Thermo::IdealGas *tmpGasPhase = (&r.Mixture()->GasPhase());
		double mw = r.Mixture()->ParticleModel()->Components()[0]->MolWt();
		fvector Hs, Cs;
		double bulkCg;
		double rhop = (r.Mixture()->Particles().GetSum(Sweep::iM) +
			r.Mixture()->Particles().GetTotalMass())
			/ (mw * r.Mixture()->SampleVolume());
		if (r.IsConstV()) {
			// Constant volume reactor: Use Cv, Us.
			tmpGasPhase->Us(Hs);
			tmpGasPhase->Cvs(Cs);
			bulkCg = tmpGasPhase->BulkCv();
		}
		else {
			// Constant pressure reactor: Use Cp, Hs.
			tmpGasPhase->Hs(Hs);
			tmpGasPhase->Cps(Cs);
			bulkCg = tmpGasPhase->BulkCp();
		}
		int pindex = r.Mech()->ParticleMech().GetParticleSpeciesIndex();
		r.Mixture()->setGasPhaseProperties(bulkCg, Cs[pindex], rhop * Cs[pindex], Hs);
	}
}

/*!
Retrieves the solution vector if the basic gpc ode solver is called
@param[in] n_sensi         Number of sensitivities computed
@param[in] n_species       Number of species in solution
@see       GetSensSolution Retrieves the pointer to sensitivity solution vector
@return    double**       Sensitivity solution array
*/
double** Solver::GetSensSolution(int n_sensi, int n_species)
{
    return m_ode.GetSensSolution(n_sensi, n_species);

}

/*!
Retrieves the solution vector if the basic gpc ode solver is called
@param[in]      n_sensi     Number of sensitivites computed
@param[in]      n_species   Number of species in solution
*/
void Solver::InitialiseSensMatrix(int n_sensi, int n_species)
{
    m_ode.InitialiseSensArray(n_sensi, n_species);
}

/*!
Retrieves the solution vector if the basic gpc ode solver is called
@param[in]      n_species   Number of species in solution
*/
void Solver::DestroySensMatrix(int n_species)
{
    m_ode.DestroySensArray(n_species);
}

/*!
Retrieves the solution vector if the basic gpc ode solver is called
@return     m_NS       Number of sensitivities initialised for computation
*/
unsigned int Solver::GetNumSens() const
{
    return m_ode.GetNSensitivities();
}

/*!
Retrieves the solution vector if the basic gpc ode solver is called
@return     m_eq       Number of equations computed by ODE solver
*/
unsigned int Solver::GetNEquations() const
{
    return m_ode.GetNEquations();
}

// COMPUTATION TIME.

// Returns the number of CT time variables tracked by this
// solver type.
unsigned int Mops::Solver::CT_Count(void) const {return 2;}

// Outputs internal computation time data to the given
// binary stream.
void Mops::Solver::OutputCT(std::ostream &out) const
{
    out.write((char*)&m_tottime, sizeof(m_tottime));
    out.write((char*)&m_chemtime, sizeof(m_chemtime));
}

// Attach sensitivity to ODE_Solver by making copy.
void Solver::AttachSensitivity(SensitivityAnalyzer &sensi) const
{
    m_ode.SetSensitivity(sensi);
}

// Outputs sensitivity results to given file stream.
void Solver::OutputSensitivity(std::fstream &fout, const Mops::Reactor &r, void *sim) const
{
    m_ode.GetSensitivity().OutputSens(fout, r, sim);
}

// Adds the CT descriptions to a vector of strings.
void Solver::CT_Names(vector<string> &names, unsigned int start) const
{
    // Resize output vector to hold names, and get iterator
    // to first insertion point.
    if (start+CT_Count() > names.size()) names.resize(start+CT_Count(), "");
    std::vector<std::string>::iterator i = names.begin()+start;

    // Add names to output array.
    *(i++) = "Total Comput. Time (s)";
    *(i++) = "ODE (gas-phase) Comput. Time (s)";
}


// COMPUTATION TIME CALCULATION.

// Calculates the time duration from a time mark to the
// current time.
double Solver::calcDeltaCT(double markt) const
{
    return (double)(clock() - markt) / (double)CLOCKS_PER_SEC;
}
