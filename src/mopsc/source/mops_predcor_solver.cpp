/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the PredCorSolver class declared in the
    mops_predcor_solver.h header file.

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

#include "mops_predcor_solver.h"
#include "mops_reactor_factory.h"
#include "mops_psr.h"

#include "sweep.h"
#include "string_functions.h"
#include "csv_io.h"
#include "local_geometry1d.h"

#include <stdexcept>

//#define CHECK_PTR

using namespace Mops;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
PredCorSolver::PredCorSolver(void)
: Sweep::FlameSolver(), m_reac_copy(NULL), m_ncalls(0)
{
}

//! Copy constructor
PredCorSolver::PredCorSolver(const PredCorSolver &sol)
: FlameSolver(sol),
  m_srcterms(sol.m_srcterms),
  m_srcterms_copy(sol.m_srcterms_copy),
  m_reac_copy(sol.m_reac_copy),
  m_ode_copy(sol.m_ode_copy),
  m_ncalls(sol.m_ncalls) {}

//! Clone the object
PredCorSolver *const PredCorSolver::Clone() const {
    return new PredCorSolver(*this);
}


// Default destructor.
PredCorSolver::~PredCorSolver(void)
{
    delete m_reac_copy;
}


// SOLVER INITIALISATION AND RESET.

// Initialises the solver to solve the given reactor.
void PredCorSolver::Initialise(Reactor &r)
{
    // Ensure 'fixed-chemistry' is set so no stochastic adjustments are made
    r.Mixture()->SetFixedChem(true);
    if (r.IsConstV())
    {
        r.Mixture()->setConstV(true);
    }

    // Set up ODE solver.
    FlameSolver::Initialise(r);
    m_ode.Initialise(r);
    m_ode.SetExtSrcTermFn(AddSourceTerms);
    m_ode_copy.Initialise(r);
    m_ode_copy.SetExtSrcTermFn(AddSourceTerms);

    // Reset number of Solve() calls.
    m_ncalls = 0;

    // Set up internal solver settings.
    m_gas_prof.resize(2, Sweep::GasPoint(r.Mech()->GasMech().Species()));

    // Set up source terms.
    m_srcterms.resize(2, SrcPoint(r.Mech()->GasMech().SpeciesCount()+2));
    m_srcterms_copy.resize(2, SrcPoint(r.Mech()->GasMech().SpeciesCount()+2));
    m_ode.SetExtSrcTerms(m_srcterms);
    m_ode_copy.SetExtSrcTerms(m_srcterms);

    // Clone the reactor.
    delete m_reac_copy;
    m_reac_copy = r.Clone();
}

// Resets the solver to solve the given reactor.
void PredCorSolver::Reset(Reactor &r)
{
    // Set up ODE solver.
    FlameSolver::Reset(r);
    m_ode.SetExtSrcTermFn(AddSourceTerms);
    m_ode_copy.ResetSolver(r);
    m_ode_copy.SetExtSrcTermFn(AddSourceTerms);

    // Reset number of Solve() calls.
    m_ncalls = 0;

    // Set up internal solver settings.
    m_gas_prof.resize(2, Sweep::GasPoint(r.Mech()->GasMech().Species()));

    // Set up source terms.
    m_srcterms.resize(2, SrcPoint(r.Mech()->GasMech().SpeciesCount()+2));
    m_srcterms_copy.resize(2, SrcPoint(r.Mech()->GasMech().SpeciesCount()+2));
    m_ode.SetExtSrcTerms(m_srcterms);
    m_ode_copy.SetExtSrcTerms(m_srcterms);

    // Ensure 'fixed-chemistry' is set so no stochastic adjustments are made
    r.Mixture()->SetFixedChem(true);
    if (r.IsConstV())
    {
        r.Mixture()->setConstV(true);
    }

    // Clone the reactor.
    delete m_reac_copy;
    m_reac_copy = r.Clone();
}


// SOLVING REACTORS.


/* aab64 - warning on using predcor with mops network simulator:
 The predcor solver appears to be incompatible with the network simulator 
 This is due to the way that the reactor class and PSR derived class are 
 set up. Essentially, because the PSR class has inflow and outflow pointers
 but the setup in the solver is for a generic reactor object, the information
 about the mixtures pointed to is lost when the reactor is copied for the 
 predcor iterations e.g. in the assignment r=*m_reac_copy. This has not currently
 been fixed. One option would be to restructure things so that the reactor, r, is 
 passed by reference to the iteration() function etc. The clone function r.Clone() 
 works correctly and could be used for the assignment. 
 To investigate enable the ifdef CHECK_PTR checks */


// Solves the coupled reactor using the predictor-corrector splitting
// algorithm up to the stop time.  Calls the output function after
// each iteration of the last internal step.
void PredCorSolver::Solve(Reactor &r, double tstop, int nsteps, int niter, 
                          Sweep::rng_type &rng, OutFnPtr out, void *data)
{
    int step=0, iter=0;

    // Calculate step size.
    double dt = (tstop - r.Time()) / (double)nsteps;

    // Internal splits without file output.
    for (step=0; step<nsteps-1; ++step) {
        // Start the iteration procedure.
        beginIteration(r, step, dt);

        // Iterate this step for the required number of runs.
        for (iter=0; iter<niter; ++iter) {
            if (m_ncalls==0) m_ode.ResetSolver(*m_reac_copy);
            iteration(r, dt, rng);
        }

        // Wind up the iteration algorithm.
        endIteration();
    }

    // Internal split with file output.
    beginIteration(r, step, dt);
    for (iter=0; iter!=niter; ++iter) {
        if (m_ncalls==0) m_ode.ResetSolver(*m_reac_copy);
        iteration(r, dt, rng);
        out(step+1, iter+1, r, *this, data);
    }
    endIteration();

    // Increment the total call count (number of times
    // the solver has been called).
    ++m_ncalls;
}

/*
// Solves the given reactor for the given time intervals.
void PredCorSolver::SolveReactor(Mops::Reactor &r, 
                                 const timevector &times, 
                                 unsigned int nruns)
{
    unsigned int icon;
    double t1;     // Current time.
    double dt, t2; // Stop time for each step.

    // Start the CPU timing clock.
    m_cpu_start = clock();
    m_chemtime  = 0.0;
    m_swp_ctime = 0.0;

    // Store the initial conditions.
    Mixture initmix(r.Mech()->ParticleMech());
    initmix.SetFracs(r.Mixture()->GasPhase().MoleFractions());
    initmix.SetTemperature(r.Mixture()->GasPhase().Temperature());
    initmix.SetDensity(r.Mixture()->GasPhase().Density());

    // Initialise the reactor with the start time.
    t1 = times[0].StartTime();
    t2 = t1;
    r.Mixture()->GasPhase().Particles().Initialise(m_pcount, r.Mech()->ParticleMech());
    r.Mixture()->GasPhase().SetMaxM0(m_maxm0);

    // Initialise the ODE solver.
    m_ode.Initialise(r);
    m_ode.SetExtSrcTermFn(AddSourceTerms);
    m_ode_copy.Initialise(r);
    m_ode_copy.SetExtSrcTermFn(AddSourceTerms);

    // Set up file output.
    writeAux(m_output_filename, *r.Mech(), times);
    openOutputFile(0); // Only one output file for this algorithm.
    
    // Write the initial conditions to file.
    fileOutput(r);

    // Set up the console output.
    icon = m_console_interval;
    setupConsole(*r.Mech());
    consoleOutput(r);

    // Set up internal solver settings.
    m_gas_prof.resize(2, Sweep::GasPoint(r.Mech()->GasMech().Species()));
    m_srcterms.resize(2, SrcPoint(r.Mech()->GasMech().SpeciesCount()+2));
    m_srcterms_copy.resize(2, SrcPoint(r.Mech()->GasMech().SpeciesCount()+2));
    m_ode.SetExtSrcTerms(m_srcterms);
    m_ode_copy.SetExtSrcTerms(m_srcterms);
    m_reac_copy = r.Clone();

    // Loop over all time intervals.
    unsigned int global_step = 0;
    timevector::const_iterator iint;
    for (iint=times.begin(); iint!=times.end(); ++iint) {
        // Print interval info to console.
        printf("Time interval from %f s to %f s (%d steps).\n", 
               iint->StartTime(), iint->EndTime(), iint->StepCount());

        // Get the step size for this interval.
        dt = (*iint).StepSize();

        // Loop over the steps in this interval.  Note that there are
        // a set number of internal split between file outputs, which
        // ensures that unnecessarily large files are not generated when
        // using short splitting steps.  Therefore we need to perform N-1
        // splits before performing the last split with file output.  This
        // is a bit messed up here because of the iterative nature of the
        // algorithm compared to Strang splitting.
        for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++global_step) {
            // Internal splits without file output.
            for (unsigned int isplit=1; isplit!=iint->SplittingStepCount(); ++isplit) {
                // Start the iteration procedure.
                beginIteration(r, global_step, iint->SplitStepSize());

                // Iterate this step for the required number of runs.
                for (unsigned int irun=0; irun!=nruns; ++irun) {
                    if (global_step==0) m_ode.ResetSolver(*m_reac_copy);
                    iteration(r, iint->SplitStepSize());
                }

                // Wind up the iteration algorithm.
                endIteration();
            }

            // Internal split with file output (and save points if 
            // required).
            beginIteration(r, global_step, iint->SplitStepSize());
            for (unsigned int irun=0; irun!=nruns; ++irun) {
                if (global_step==0) m_ode.ResetSolver(*m_reac_copy);
                iteration(r, iint->SplitStepSize());
                fileOutput(r);
                if (istep == iint->StepCount()-1) {
                    createSavePoint(r, global_step+1, irun);
                }
            }
            endIteration();

            // Generate console output (last iteration only).
            if (--icon == 0) {
                consoleOutput(r);
                icon = m_console_interval;
            }
        } // End of output step loop.
    } // End of time interval loop.

    // Clear up internal solver settings.
    m_gas_prof.clear();
    m_srcterms.clear();
    m_srcterms_copy.clear();
    delete m_reac_copy;

    // Close the output files.
    closeOutputFile();
}
*/

// Sets up the workspace to begin a new iteration.  This includes
// generating an initial chemistry profile over the step.  As a first
// guess the chemistry is assumed to be constant in this interval.  The
// initial source terms from sweep are then calculated using this chemistry
// profile.
void PredCorSolver::beginIteration(Reactor &r, unsigned int step, double dt)
{
    // Calculate gas profile time interval based on number of profile
    // points.
    double h = dt / (m_gas_prof.size()-1);

    // Generate the initial chemistry profile over this step.  The
    // first point should have already been set before this routine
    // was called, either by the previous step or by the SolveReactor()
    // function, if this is the first step.  The chemistry profile is
    // initially assumed to be constant, so we just copy the first point
    // to the other points.
    for (unsigned int i=1; i!=m_gas_prof.size(); ++i) {
        m_gas_prof[i].Time = m_gas_prof[i-1].Time + h;
        m_gas_prof[i].Gas = m_gas_prof[0].Gas;
    }

    // The source terms at the beginning of this step are those from the
    // end of the previous step.  It is assumed that for the first step
    // the source terms in m_srcterms_copy are appropriately assigned
    // before this routine is called.
    m_srcterms[0] = m_srcterms_copy[1];

    // Now extrapolate the source terms for the end of the step.
    if (m_ncalls == 0) {
        // This is the first call.  Just use a constant source-term
        // extrapolation.
        m_srcterms[1] = m_srcterms[0];
        m_srcterms[1].Time = m_srcterms[0].Time + dt;
        //for (unsigned int i=0; i!=r.Mech()->GasMech().SpeciesCount()+2; ++i) {
        //    m_srcterms[1].Terms[i] = m_srcterms[0].Terms[i];
        //}
    } else {
        // Use linear extrapolation from previous step to get source
        // terms for this step.
        linExSrcTerms(m_srcterms[1], m_srcterms_copy, dt);
    }

    // Now make the copy of the source terms.
    m_srcterms_copy[0] = m_srcterms[0];
    m_srcterms_copy[1] = m_srcterms[1];

    // Copy the reactor object in order to reset the solver, to iterate
    // over the step.
    *m_reac_copy = r;

    // Make a copy of the ODE workspace.
    m_ode_copy = m_ode;
}

// Performs a step-wise iteration on the reactor to recalculate
// the source terms for the gas-phase effect on the particle model.
void PredCorSolver::iteration(Reactor &r, double dt, Sweep::rng_type &rng)
{
	
    // aab64: check what happens to flow pointers during assignment
#ifdef CHECK_PTR
    cout << "in iteration check 1\n";
    // Cast the reactor to a PSR reactor and check outflow pointers
    Mops::PSR *psr = dynamic_cast<Mops::PSR*>(&r);
    if (psr != NULL)
    {
	Mops::FlowPtrVector::const_iterator itbeg = psr->Mops::PSR::Outflows().begin();
	Mops::FlowPtrVector::const_iterator itend = psr->Mops::PSR::Outflows().end();
	while (itbeg != itend) {
	    assert(&(*itbeg)->Mixture()->GasPhase());
	    ++itbeg;
	}
    }
    
    cout << "in iteration check 2\n";
    // Cast the reactor copy to a PSR reactor and check outflow pointers
    Mops::PSR *psr_copy = dynamic_cast<Mops::PSR*>(&(*m_reac_copy));
    if (psr_copy != NULL)
    {
	Mops::FlowPtrVector::const_iterator itbeg = psr_copy->Mops::PSR::Outflows().begin();
	Mops::FlowPtrVector::const_iterator itend = psr_copy->Mops::PSR::Outflows().end();
	while (itbeg != itend) {
	    assert(&(*itbeg)->Mixture()->GasPhase());
	    ++itbeg;
	}
    }
#endif

    // Reset reactor and solver for another iteration.
    r = *m_reac_copy;
	
    // aab64: check what happens to flow pointers during assignment
#ifdef CHECK_PTR
    cout << "in iteration check 3, post assignment\n";
    // Cast the reactor to a PSR reactor and check outflow pointers
    //Mops::PSR *psr = dynamic_cast<Mops::PSR*>(&r);
    if (psr != NULL)
    {
	Mops::FlowPtrVector::const_iterator itbeg = psr->Mops::PSR::Outflows().begin();
	Mops::FlowPtrVector::const_iterator itend = psr->Mops::PSR::Outflows().end();
	while (itbeg != itend) {
	    assert(&(*itbeg)->Mixture()->GasPhase());
            ++itbeg;
	}
    }
    cout << "in iteration checks complete\n";
#endif

    // Note the start time.
    double ts1=r.Time();

    // Refresh from the copy of the ODE solver.
    m_ode = m_ode_copy;
    //m_ode.ResetSolver();

    // Set source terms in ODE solver object.
    //m_srcterms[0] = m_srcterms_copy[0];
    m_srcterms[1] = m_srcterms_copy[1];
    m_ode.SetExtSrcTerms(m_srcterms);

    // Generate chemistry profile over the step using the current
    // source terms.
    m_cpu_mark = clock();
        generateChemProfile(r, dt);
    m_chemtime += calcDeltaCT(m_cpu_mark);

    // Solve step using Sweep in order to update the
    // source terms.
    m_cpu_mark = clock();
        // Set the stop time.    
        double ts2 = ts1+dt;

        // Scale M0 according to gas-phase expansion.
       r.Mixture()->AdjustSampleVolume(m_reac_copy->Mixture()->GasPhase().MassDensity()		
                / r.Mixture()->GasPhase().MassDensity());									

        // Run Sweep for this time step.
        Run(ts1, ts2, *r.Mixture(), r.Mech()->ParticleMech(), rng);
    m_swp_ctime += calcDeltaCT(m_cpu_mark);

    // Now update the source terms at the end of the step.
    calcSrcTerms(m_srcterms_copy[1], r);

    // Apply under-relaxation to source terms.
    relaxSrcTerms(m_srcterms_copy[1], m_srcterms[1], m_rlx_coeff);
}

// Terminates an iteration step.
void PredCorSolver::endIteration()
{
    // Copy last profile point to the first point.
    m_gas_prof[0].Time = (m_gas_prof.end()-1)->Time;
    m_gas_prof[0].Gas  = (m_gas_prof.end()-1)->Gas;
}

// Generates a chemistry profile over the required time interval
// by solving the gas-phase chemistry equations with the current
// source terms.
void PredCorSolver::generateChemProfile(Reactor &r, double dt)
{
    // Calculate time interval between gas profile points.
    double h = dt / (double)(m_gas_prof.size()-1);
    double t1 = r.Time();

    // Save initial profile point.
    m_gas_prof[0].Time = t1;
    m_gas_prof[0].Gas = r.Mixture()->GasPhase();

    // Solve the reactor for the remaining points.
    for (unsigned int i=1; i!=m_gas_prof.size(); ++i) {
        m_ode.Solve(r, t1+=h);
        r.SetTime(t1);
        m_gas_prof[i].Time = t1;
        m_gas_prof[i].Gas = r.Mixture()->GasPhase();
    }
}

// Calculates the instantaneous source terms.
void PredCorSolver::calcSrcTerms(SrcPoint &src, const Reactor &r)
{
    // Get rates-of-change using Sweep mechanism.
    src.Time = r.Time();
    r.Mech()->ParticleMech().CalcGasChangeRates(r.Time(), *r.Mixture(), 
                                                Geometry::LocalGeometry1d(), src.Terms);

    // Convert to concentrations to compute change in enthalpy
    Mops::fvector csrc;
    csrc.resize(r.Mech()->GasMech().SpeciesCount() + 2, 0.0);
    fvector::iterator rhodot = src.Terms.begin() + r.Mech()->GasMech().SpeciesCount() + 1;
    for (unsigned int k=0; k!=r.Mech()->GasMech().SpeciesCount(); ++k) {
        csrc[k] = (r.Mixture()->GasPhase().Density() * src.Terms[k]) + 
            (r.Mixture()->GasPhase().MoleFraction(k) * (*rhodot));
    }

    // Calculate the enthalpy-based temperature change 
    // rate using the species change rates and an adiabatic
    // assumption.
    double dT_dt = 0.0;
    if (r.EnergyEquation() == Reactor::ConstT) {
        src.Terms[r.Mech()->GasMech().SpeciesCount()] = dT_dt;
    } else {
        // Note this does not account for particle formation so it should only be
        // used for systems where enthalpy of particle formation is insignificant
        dT_dt = energySrcTerm(r, csrc); 
        src.Terms[r.Mech()->GasMech().SpeciesCount()] += dT_dt;
    }

    // Calculate density change based on whether reactor is constant
    // volume or constant pressure.
    if (r.IsConstP()) {
        // Constant pressure: density change due to change in temperature.
        double rho_T = r.Mixture()->GasPhase().Density() / r.Mixture()->GasPhase().Temperature();
        src.Terms[r.Mech()->GasMech().SpeciesCount() + 1] = dT_dt * (-1.0 * rho_T);
    } else {
        // Constant volume: add gdot to wdot for density derivative.
        // This was already done in calcGasChangeRates - no action required. 
        src.Terms[r.Mech()->GasMech().SpeciesCount()+1] += 0.0;
    }
}

// Calculate the adiabatic temperature change rate due to the
// species source terms.
double PredCorSolver::energySrcTerm(const Reactor &r, fvector &src)
{
    if (r.EnergyEquation() == Reactor::Adiabatic) {
        // Adiabatic temperature model.
        fvector Hs;
        double C;

        // Calculate heat capacity and enthalpy/internal energy.
        if (r.IsConstV()) {
            // Constant volume reactor: Use Cv, Us.
            C = r.Mixture()->GasPhase().BulkCv();
            r.Mixture()->GasPhase().Us(Hs);
        } else {
            // Constant pressure reactor: Use Cp, Hs.
            C = r.Mixture()->GasPhase().BulkCp();
            r.Mixture()->GasPhase().Hs(Hs);
        }

        // Calculate and return temperature source term.
        double Tdot = 0.0;
        for (unsigned int i=0; i!=r.Mech()->GasMech().SpeciesCount(); ++i) {
            Tdot -= Hs[i] * src[i];
        }
        return Tdot / (C * r.Mixture()->GasPhase().Density());

    } else if (r.EnergyEquation() == Reactor::ConstT) {
        // Constant temperature model.
        return 0.0;
    }

    // Assume constant temperature by default.
    return 0.0;
}

// Extrapolates the source terms over the given time interval.  Two
// points are given, from which the gradient is calculated, which 
// allows the third point to be calculated.
void PredCorSolver::linExSrcTerms(SrcPoint &src, 
                                  const SrcProfile &prof, 
                                  double dt)
{
    if (prof.size() > 1) {
        // Set the time.
        src.Time = prof[prof.size()-1].Time + dt;

        // Pre-calculate gradient factors (use only last 2 profile points).
        SrcProfile::const_iterator j1 = prof.begin()+(prof.size()-2);
        SrcProfile::const_iterator j2 = prof.begin()+(prof.size()-1);
        double t_ratio = dt / (j2->Time - j1->Time);

        // Find new source terms with linear extrapolation.
        for (unsigned int i=0; i!=j2->Terms.size(); ++i) {
            src.Terms[i] = (j2->Terms[i] - j1->Terms[i]) * t_ratio;
        }
    } else if (prof.size() > 0) {
        // Set the time.
        src.Time = prof[0].Time + dt;

        // Assume the source terms are constant.
        for (unsigned int i=0; i!=prof[0].Terms.size(); ++i) {
            src.Terms[i] = prof[0].Terms[i];
        }
    } else {
        // The profile is empty, this is an error.  
        // Set all source terms to zero.
        src.Time = dt;
        for (unsigned int i=0; i!=src.Terms.size(); ++i) {
            src.Terms[i] = 0.0;
        }
    }
}

// Applies under-relaxation to the first source point, using the
// second source point as the initial values.
void PredCorSolver::relaxSrcTerms(SrcPoint &src, const SrcPoint &init, 
                                  double rcoeff)
{
    for (unsigned int i=0; i!=src.Terms.size(); ++i) {
        src.Terms[i] = ((1.0-rcoeff)*src.Terms[i]) + (rcoeff*init.Terms[i]);
//        src.Terms[i] += rcoeff * (init.Terms[i] - src.Terms[i]);
    }
}

// Applies under-relaxation to the first source point, using the
// second source point as the initial values.
void PredCorSolver::copySrcPoint(const SrcPoint &from, SrcPoint &to)
{
    to.Time = from.Time;
    to.Terms.resize(from.Terms.size());
    for (unsigned int i=0; i!=from.Terms.size(); ++i) {
        to.Terms[i] = from.Terms[i];
    }
}


// POST-PROCESSING.
/*
void PredCorSolver::PostProcess(const std::string &filename, 
                               unsigned int nruns) const
{
    // READ AUXILLIARY INFORMATION.

    // Read auxilliary information about the simulation (mechanism
    // and time intervals).
    Mops::Mechanism mech;
    Mops::timevector times;
    readAux(filename, mech, times);

    // SETUP OUTPUT DATA STRUCTURES.

    // Get reference to particle mechanism.
    Sweep::Mechanism &pmech = mech.ParticleMech();

    // Calculate number of output points.
    unsigned int npoints = 1; // 1 for initial conditions.
    for(Mops::timevector::const_iterator i=times.begin(); i!=times.end(); ++i) {
        npoints += i->StepCount();
    }

    // Declare chemistry outputs (averages and errors).
    vector<fvector> achem(npoints), echem(npoints);

    // Declare particle stats outputs (averages and errors).
    vector<fvector> astat(npoints), estat(npoints);
    Sweep::Stats::EnsembleStats stats(pmech);

    // Declare CPU time outputs (averages and errors).
    vector<vector<double> > acpu(npoints), ecpu(npoints);


    // OPEN SIMULATION OUTPUT FILE.

    // Build the simulation input file name.
    string fname = filename + "(0).dat";

    // Open the simulation output file.
    fstream fin(fname.c_str(), ios_base::in | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!fin.good()) {
        throw runtime_error("Failed to open file for post-processing "
                            "(Mops, PredCorSolver::PostProcess).");
    }

    // READ INITIAL CONDITIONS.

    // The initial conditions were only written to the file
    // once, as they are the same for all runs.
    readGasPhaseDataPoint(fin, mech, achem[0], echem[0], nruns>1);
    readParticleDataPoint(fin, pmech, astat[0], estat[0], nruns>1);
    readCTDataPoint(fin, 3, acpu[0], ecpu[0], nruns>1);
    
    // The initial conditions need to be multiplied by the number
    // of runs in order for the calcAvgConf to give the correct
    // result.
    multVals(achem[0], nruns);
    multVals(astat[0], nruns);
    multVals(acpu[0], nruns);

    // READ ALL TIME POINTS FOR ALL RUNS.

    // Loop over all time intervals.
    unsigned int step = 1;
    for (Mops::timevector::const_iterator iint=times.begin(); 
         iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
            // Loop over all runs.
            for(unsigned int irun=0; irun!=nruns; ++irun) {
                readGasPhaseDataPoint(fin, mech, achem[step], echem[step], nruns>1);
                readParticleDataPoint(fin, pmech, astat[step], estat[step], nruns>1);
                readCTDataPoint(fin, 3, acpu[step], ecpu[step], nruns>1);
            }
        }

        // Close the input file.
        fin.close();
    }

    // CALCULATE AVERAGES AND CONFIDENCE INTERVALS.
    
    calcAvgConf(achem, echem, nruns);
    calcAvgConf(astat, estat, nruns);
    calcAvgConf(acpu, ecpu, nruns);

    // CREATE CSV FILES.

    writeGasPhaseCSV(filename+"-chem.csv", mech, times, achem, echem);
    writeParticleStatsCSV(filename+"-part.csv", mech, times, astat, estat);
    writeCT_CSV(filename+"-cpu.csv", times, acpu, ecpu);

    // POST-PROCESS PSLs.

    // Now post-process the PSLs.
    postProcessPSLs(nruns, mech, times);
}

// Multiplies all values in a vector by a scaling factor.
void PredCorSolver::multVals(fvector &vals, double scale)
{
    for (fvector::iterator i=vals.begin(); i!=vals.end(); ++i) {
        (*i) *= scale;
    }
}
*/
/*
// UNDER-RELAXATION.

// Returns the under-relaxation coefficient.
double PredCorSolver::UnderRelaxCoeff(void) const
{
    return m_rlx_coeff;
}

// Sets the under-relaxation coefficient.
void PredCorSolver::SetUnderRelaxCoeff(double relax) {m_rlx_coeff = relax;}
*/

// SOURCE TERM CALCULATION (REQUIRED FOR REACTOR ODE SOLVER).

// Adds the source term contribution to the RHS supplied by the
// reactor class.
void PredCorSolver::AddSourceTerms(double *rhs, unsigned int n, 
                                   double t, const SrcProfile &src)
{
    // Locate the first src point after the given time t.
    SrcProfile::const_iterator j = LocateSrcPoint(src, t);

    if (j == src.begin()) {
        // This time is before the beginning of the profile, just
        // use the first point and assume constant.
        for (unsigned int k=0; k!=min(n,(unsigned int)j->Terms.size()); ++k) {
            rhs[k] += j->Terms[k];
        }
    } else if (j == src.end()) {
        // This time is after the end of the profile.  Just assume
        // constant in this range.
        --j;
        for (unsigned int k=0; k!=min(n,(unsigned int)j->Terms.size()); ++k) {
            rhs[k] += j->Terms[k];
        }
    } else {
        // Use linear interpolation (trapezium rule) to calculate the 
        // source terms.
        SrcProfile::const_iterator i = j; --i; // Point before j;
        double tterm = (t - i->Time) / (j->Time - i->Time);
        for (unsigned int k=0; k!=min(n,(unsigned int)min(i->Terms.size(),j->Terms.size())); ++k) {
            rhs[k] += i->Terms[k] + (tterm * (j->Terms[k] - i->Terms[k]));
        }
    }
}
