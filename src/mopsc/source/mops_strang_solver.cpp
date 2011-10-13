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
StrangSolver::StrangSolver(void)
{
}

// Default destructor.
StrangSolver::~StrangSolver(void)
{
}


// SOLVING REACTORS.

// Solves the coupled reactor using a Strang splitting algorithm
// up to the stop time.  calls the output routine once at the
// end of the function.  niter is ignored.
void StrangSolver::Solve(Reactor &r, real tstop, int nsteps, int niter, 
                         Sweep::rng_type &rng,
                         OutFnPtr out, void *data)
{
    // Mark the time at the start of the step, in order to
    // calculate total computation time.
    clock_t totmark = clock();

    // Time counters.
    real t2 = r.Time();
    real dt = (tstop - t2) / (real)nsteps; // Step size.
    real h  = dt * 0.5; // Half step size.

    // Sweep time counters.
    real ts1 = r.Time();
    real ts2 = ts1;

    // Variables required to ensure particle number density is correctly
    // scaled with gas-phase expansion.
    real rho = 0.0;

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
    
    for (int i=1; i<nsteps; ++i) {
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

    // Calculate total computation time.
    m_tottime += calcDeltaCT(totmark);

    // Call the output function.
    if (out) out(nsteps, niter, r, *this, data);
}

/*
// Solves the given reactor for the given time intervals.
void StrangSolver::SolveReactor(Mops::Reactor &r, 
                                const timevector &times, 
                                unsigned int nruns)
{
    unsigned int icon;
    real t1;     // Current time.
    real dt, t2; // Stop time for each step.

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

    // Set up file output.
    writeAux(m_output_filename, *r.Mech(), times);

    // Set up the console output.
    icon = m_console_interval;
    setupConsole(*r.Mech());

    for (unsigned int irun=0; irun!=nruns; ++irun) {
        // Start the CPU timing clock.
        m_cpu_start = clock();
        m_chemtime  = 0.0;
        m_swp_ctime = 0.0;

        // Initialise the reactor with the start time.
        t2 = times[0].StartTime();
        r.Mixture()->GasPhase().SetFracs(initmix.MoleFractions());
        r.Mixture()->GasPhase().SetTemperature(initmix.Temperature());
        r.Mixture()->GasPhase().SetDensity(initmix.Density());
        r.Mixture()->GasPhase().Reset(m_maxm0);
        r.SetTime(t2);

        // Set up the ODE solver for this run.
        m_ode.ResetSolver(r);
        m_ode.SetATOL(m_atol);
        m_ode.SetRTOL(m_rtol);

        // Begin file output for this run.
        openOutputFile(irun);
        printf("Run Number %d of %d.\n", irun+1, nruns);
        m_console.PrintDivider();

        // Output initial conditions.
        fileOutput(r);
        consoleOutput(r);

        // Loop over the time intervals.
        unsigned int global_step = 0;
        timevector::const_iterator iint;
        for (iint=times.begin(); iint!=times.end(); ++iint) {
            // Get the step size for this interval.
            dt = (*iint).StepSize();

            // Loop over the steps in this interval.
            unsigned int istep;
            for (istep=0; istep!=(*iint).StepCount(); ++istep, ++global_step) {
                // Run the reactor solver for this step (timed).
                multiStrangStep(iint->SplitStepSize(), iint->SplittingStepCount(), r);

                // Generate file output.
                fileOutput(r);

                // Generate console output.
                if (--icon == 0) {
                    consoleOutput(r);
                    icon = m_console_interval;
                }
            }

            createSavePoint(r, global_step, irun);
        }

        // Close the output files.
        closeOutputFile();
    }
}
*/

// SOLUTION ROUTINES.

void StrangSolver::multiStrangStep(Mops::real dt, unsigned int n, Mops::Reactor &r,
                                   Sweep::rng_type &rng)
{
    // Time counters.
    real t2 = r.Time();
    real h  = dt * 0.5; // Half step size.

    // Sweep time counters.
    real ts1 = r.Time();
    real ts2 = ts1;

    // Variables required to ensure particle number density is correctly
    // scaled with gas-phase expansion.
    real rho = 0.0;

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


// POST-PROCESSING.
/*
void StrangSolver::PostProcess(const std::string &filename, 
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

    // Declare data structure for particle tracking data.
    // Runs -> time steps -> particles -> coordinate.
    vector<vector<vector<fvector> > > ptrack(nruns);

    // READ ALL RUNS.

    // Now we must read the reactor conditions and all time points
    // and all runs.
    for(unsigned int irun=0; irun!=nruns; ++irun) {
        // Build the simulation input file name.
        string fname = filename + "(" + cstr(irun) + ").dat";

        // Open the simulation output file.
        fstream fin(fname.c_str(), ios_base::in | ios_base::binary);

        // Throw error if the output file failed to open.
        if (!fin.good()) {
            throw runtime_error("Failed to open file for post-processing "
                                "(Mops, StrangSolver::PostProcess).");
        }

        // Read initial stats and computation times.
        readGasPhaseDataPoint(fin, mech, achem[0], echem[0], nruns>1);
        readParticleDataPoint(fin, pmech, astat[0], estat[0], nruns>1);
        readCTDataPoint(fin, 3, acpu[0], ecpu[0], nruns>1);
        ptrack[irun].resize(npoints); // Resize particle tracking vector.
        readPartTrackPoint(fin, pmech, ptrack[irun][0]);

        // Loop over all time intervals.
        unsigned int step = 1;
        for (Mops::timevector::const_iterator iint=times.begin(); 
             iint!=times.end(); ++iint) {
            // Loop over all time steps in this interval.
            for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
                readGasPhaseDataPoint(fin, mech, achem[step], echem[step], nruns>1);
                readParticleDataPoint(fin, pmech, astat[step], estat[step], nruns>1);
                readCTDataPoint(fin, 3, acpu[step], ecpu[step], nruns>1);
                readPartTrackPoint(fin, pmech, ptrack[irun][step]);
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
    for(unsigned int irun=0; irun!=nruns; ++irun) {
        writePartTrackCSV(filename+"("+cstr(irun)+")-track", mech, times, ptrack[irun]);
    }

    // POST-PROCESS PSLs.

    // Now post-process the PSLs.
    postProcessPSLs(nruns, mech, times);
}
*/
