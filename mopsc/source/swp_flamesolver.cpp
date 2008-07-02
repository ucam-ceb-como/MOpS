/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the FlameSolver class declared in the
    swp_flamesolver.h header file.

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

#include "swp_flamesolver.h"
#include "swp_gas_profile.h"
#include "mops_timeinterval.h"
#include "sweep.h"
#include "string_functions.h"
#include "csv_io.h"
#include <fstream>
#include <stdexcept>
#include <string>

using namespace Sweep;
using namespace Sweep::ActSites;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
FlameSolver::FlameSolver()
{
}

// Default destructor.
FlameSolver::~FlameSolver()
{
}

// PROFILE INPUT.

// Reads a flame gas-phase profile from a TAB formatted file.
void FlameSolver::LoadGasProfile(const std::string &file, Mops::Mechanism &mech)
{
    map<real,real> alpha_prof;

    // Clear the current gas-phase profile.
    m_gasprof.clear();

    // Open the file to read.
    ifstream fin(file.c_str(), ios::in);
    if (!fin.good()) {
        throw runtime_error("Unable to open gas profile input "
                            "file (Mops, Sweep::FlameSolver::LoadGasProfile).");
    }

    // Variables read from file.
    vector<string> subs;
    string delim = ",\t \r"; // Possible file delimiters (comma, tab and space).
    string line;

    // Get the first line (should be the header defining the columns).
    if (!getline(fin, line).eof()) {

        // Split the line to get the column headings.
        split(line, subs, delim);

        // Get important column indices (time, temperature and pressure).
        int tcol=-1, Tcol=-1, Pcol=-1, Acol = -1;
        tcol = findinlist(string("Time"), subs);
        Tcol = findinlist(string("T"), subs);
        Pcol = findinlist(string("P"), subs);
        Acol = findinlist(string("Alpha"), subs);

        // Check that the file contains required columns.
        if (tcol < 0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no Time "
                                "column (Mops, Sweep::FlameSolver::LoadGasProfile).");
        }
        if (Tcol < 0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no temperature "
                                "column (Mops, Sweep::FlameSolver::LoadGasProfile).");
        }
        if (Pcol < 0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no pressure "
                                "column (Mops, Sweep::FlameSolver::LoadGasProfile).");
        }


        // All other columns are chemical species.  Add them, and note
        // their columns.
        map<unsigned int,int> spcols;
        for (int i=0; (unsigned)i!=subs.size(); ++i) {
            if ((i!=tcol) && (i!=Tcol) && (i!=Pcol) && (i!=Acol)) {
                Sprog::Species *sp = mech.AddSpecies();
                sp->SetName(subs[i]);
                spcols[i] = mech.FindSpecies(subs[i]);
            }
        }

        // Now we can read the profile.
        while(!getline(fin, line).eof()) {
            // Set t=0 and create a new IdealGas object.
            real t = 0.0;
            real T = 0.0;
            real P = 0.0;
            real alpha = 0.0;
            GasPoint gpoint(mech.Species());
//            Sprog::Thermo::IdealGas gas(mech.Species());

            // Split the line by columns.
            split(line, subs, delim);

            // Loop over all the elements in the line and save them
            // to the correct gas-phase variable.
            for (int i=0; (unsigned)i!=subs.size(); ++i) {
                if (i==tcol) {
                    // This is the Time column.
                    t = cdble(subs[i]);
                    gpoint.Time = t;
                } else if (i==Tcol) {
                    // This is the temperature column.
                    T = cdble(subs[i]);
                } else if (i==Pcol) {
                    // This is the pressure column.
                    P = cdble(subs[i]);
                } else if (i==Acol) {
                    alpha = cdble(subs[i]);
                } else {
                    // This is a gas-phase species column.
                    map<unsigned int,int>::iterator isp = spcols.find(i);
                    if (isp != spcols.end()) {
                        gpoint.Gas.RawData()[isp->second] = cdble(subs[i]);
                    }
                }
            }

            // Set up the gas-phase by setting temperature, pressure and
            // normalising the mixture fractions.
            // TODO:  This will give the wrong component densities
            //        unless all species are specified!
            gpoint.Gas.SetTemperature(T);
            gpoint.Gas.SetPressure(P*1.0e5);
            gpoint.Gas.Normalise();

            // Add the profile point.
            alpha_prof[t] = alpha;
            m_gasprof.push_back(gpoint);
        }

        // Set up ABF model to use alpha profile.
        ABFModel::Instance().Initialise(mech.ParticleMech());
        ABFModel::Instance().SetAlphaProfile(alpha_prof);

        // Close the input file.
        fin.close();

        // Sort the profile by time.
        SortGasProfile(m_gasprof);
    } else {
        // There was no data in the file.
        fin.close();
        throw runtime_error("Input file contains no data "
                            "(Mops, Sweep::FlameSolver::LoadGasProfile).");
    }
}


// SOLUTION AND POST-PROCESSING.

// Performs stochastic stepping algorithm up to specified stop time using
// the given mechanism to define the stochastic processes.  Updates given
// system accordingly.  On error returns <0, otherwise returns 0.  In this
// flavour the gas-phase chemistry is interpolated from a vector of
// IdealGas objects rather than being taken from the given system object.
// However, the particles in the system object are updated accordingly.
int FlameSolver::Run(real &t, real tstop, const GasProfile &gasphase, 
                     Cell &sys, const Mechanism &mech)
{
    int err = 0;
    real tsplit, dtg, dt, jrate;
    fvector rates(mech.TermCount(), 0.0);
    
    // Save the initial chemical conditions in sys so that we
    // can restore them at the end of the run.
    Sprog::Thermo::IdealGas chem = sys;

    // Store if chemical conditions are fixed at present, because we
    // shall set them to be fixed during this run, to be restored afterwards.
    bool fixedchem = sys.FixedChem();
    sys.SetFixedChem();

    // Global maximum time step.
    dtg     = tstop - t;
    m_maxdt = dtg / 3.0;
    m_tstop = tstop;

    // Loop over time until we reach the stop time.
    while (t < tstop)
    {
        // Calculate LPDA splitting time step.
        if (mech.AnyDeferred() && (sys.ParticleCount() > 0)) {
            // Calculate the chemical conditions.
            linInterpGas(t, gasphase, sys);

            // Get the process jump rates (and the total rate).
            jrate = mech.CalcJumpRates(t, sys, sys, rates);

            // Calculate the splitting end time.
            tsplit = calcSplitTime(t, tstop, jrate, sys.ParticleCount(), dtg);
        } else {
            // There are no deferred processes, therefore there
            // is no need to perform LPDA splitting steps.
            tsplit = tstop;
        }

        // Perform stochastic jump processes.
        while (t < tsplit) {
            // Calculate the chemical conditions.
            linInterpGas(t, gasphase, sys);

            // Calculate jump rates.
            jrate = mech.CalcJumpRates(t, sys, sys, rates);

            // Perform time step.
            dt = timeStep(t, sys, mech, rates, jrate);
            if (dt >= 0.0) {
                t += dt; t = min(t, tstop);
            } else {
                return -1;
            }
        }

        // Perform Linear Process Deferment Algorithm to
        // update all deferred processes.
        if (mech.AnyDeferred()) {
            linInterpGas(t, gasphase, sys);
            mech.LPDA(t, sys);
        }
    }

    // Restore initial chemical conditions to sys.
    sys = chem;
    sys.SetFixedChem(fixedchem);

    return err;
}

// Run the solver for the given reactor and the 
// given time intervals.
void FlameSolver::SolveReactor(Mops::Reactor &r, 
                               const Mops::timevector &times,
                               unsigned int nruns)
{
    unsigned int icon;
    real t1;     // Current time.
    real dt, t2; // Stop time for each step.

    // Initialise the reactor with the start time.
    t1 = times[0].StartTime();
    t2 = t1;
    dynamic_cast<Sweep::Cell&>(*r.Mixture()) = m_gasprof[0].Gas;
    r.Mixture()->Particles().Initialise(m_pcount, r.Mech()->ParticleMech());
    r.Mixture()->SetMaxM0(m_maxm0);

    // Initialise ODE solver.
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

        // Re-initialise the solution.
        r.Mixture()->Reset(m_maxm0);
        t1 = times[0].StartTime();
        t2 = t1;
        r.SetTime(t1);

        // Begin file output for this run.
        openOutputFile(irun);
        printf("Run Number %d of %d.\n", irun+1, nruns);
        m_console.PrintDivider();

        // Output initial conditions.
        fileOutput(r);
        consoleOutput(r);

        // Loop over the time intervals.
        unsigned int global_step = 0;
        Mops::timevector::const_iterator iint;
        for (iint=times.begin(); iint!=times.end(); ++iint) {
            // Get the step size for this interval.
            dt = (*iint).StepSize();

            // Loop over the steps in this interval.
            for (unsigned int istep=0; istep<(*iint).StepCount(); 
                 ++istep, ++global_step) {
                // Run the reactor solver for this step (timed).
                m_cpu_mark = clock();
                    Run(t1, t2+=dt, m_gasprof, *r.Mixture(), 
                        r.Mech()->ParticleMech());
                    r.SetTime(t1);
                m_swp_ctime += calcDeltaCT(m_cpu_mark);

                // Generate file output.
                r.SetTime(t1);
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

// Post-processes binary output files with the given file name
// into CSV files.
void FlameSolver::PostProcess(const std::string &filename, unsigned int nruns) const
{
    // Read auxilliary information about the simulation (mechanism
    // and time intervals).
    Mops::Mechanism mech;
    Mops::timevector times;
    readAux(filename, mech, times);
    Sweep::Mechanism &pmech = mech.ParticleMech();

    // Calculate number of output points.
    unsigned int npoints = 1; // 1 for initial conditions.
    for(Mops::timevector::const_iterator i=times.begin(); i!=times.end(); ++i) {
        npoints += i->StepCount();
    }
    // Declare gas-phase output variables (temporary, not processed).
    fvector achem, echem;
    // Declare stats outputs (averages and errors).
    vector<fvector> astat(npoints), estat(npoints);
    Stats::EnsembleStats stats(pmech);

    // Declare CPU time outputs (averages and errors).
    vector<fvector> acpu(npoints), ecpu(npoints);


    // Now we must read the mixture conditions at all time points for
    // all runs:

    for(unsigned int irun=0; irun!=nruns; ++irun) {
        // Build the simulation input file name.
        string fname = filename + "(" + cstr(irun) + ").dat";

        // Open the simulation output file.
        fstream fin(fname.c_str(), ios_base::in | ios_base::binary);

        // Throw error if the output file failed to open.
        if (!fin.good()) {
            throw runtime_error("Failed to open file for post-processing "
                                "(Mops, Sweep::FlameSolver::PostProcess).");
        }

        // Read initial stats and computation times.
        {
            // Read the gas-phase.
            readGasPhaseDataPoint(fin, mech, achem, echem, nruns>1);
            // Read the stats.
            readParticleDataPoint(fin, pmech, astat[0], estat[0], nruns>1);
            // Read the computation time.
            readCTDataPoint(fin, 3, acpu[0], ecpu[0], nruns>1);
        }

        // Loop over all time intervals.
        unsigned int step = 1;
        for (Mops::timevector::const_iterator iint=times.begin(); 
             iint!=times.end(); ++iint) {
            // Loop over all time steps in this interval.
            for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
                readGasPhaseDataPoint(fin, mech, achem, echem, nruns>1);
                readParticleDataPoint(fin, pmech, astat[step], estat[step], nruns>1);
                readCTDataPoint(fin, 3, acpu[step], ecpu[step], nruns>1);
            }
        }

        // Close the input file.
        fin.close();
    }
    
    // CALCULATE AVERAGES AND CONFIDENCE INTERVALS.

    calcAvgConf(astat, estat, nruns);
    calcAvgConf(acpu, ecpu, nruns);

    // OUTPUT TO CSV FILE.
    
    writeParticleStatsCSV(filename+"-part.csv", mech, times, astat, estat);
    writeCT_CSV(filename+"-cpu.csv", times, acpu, ecpu);

    // POST-PROCESS PSLs.

    // Now post-process the PSLs.
    postProcessPSLs(nruns, mech, times);
}


// HELPER ROUTINES.

void FlameSolver::linInterpGas(Sweep::real t, 
                               const GasProfile &gasphase, 
                               Sprog::Thermo::IdealGas &gas) const
{
    // Get the time point after the required time.
    GasProfile::const_iterator j = LocateGasPoint(gasphase, t); //gasphase.upper_bound(t);
    
    if (j == gasphase.begin()) {
        // This time is before the beginning of the profile.  Return
        // the first time point.
        gas = j->Gas;
    } else if (j == gasphase.end()) {
        // This time is after the profile.  Return the last time
        // point
        --j;
        gas = j->Gas;
    } else {
        // Get the time point before the required time.
        GasProfile::const_iterator i = j; --i;

        // Assign the conditions to this point.
        gas = i->Gas;
        
        // Calculate time interval between points i and j.
        real dt_pro = j->Time - i->Time;

        // Calculate time interval between point i and current time.
        real dt = t - i->Time;

        // Calculate the intermediate gas-phase mole fractions by linear
        // interpolation of the molar concentrations.
        real dens = 0.0;
        for (unsigned int k=0; k<gas.Species()->size(); ++k) {
            real dc = (j->Gas.MolarConc(k) - i->Gas.MolarConc(k)) * dt / dt_pro;
            gas.RawData()[k] = gas.MolarConc(k) + dc;
            dens += gas.RawData()[k];
        }
        gas.Normalise();

        // Now use linear interpolation to calculate the temperature.
        real dT = (j->Gas.Temperature() - i->Gas.Temperature()) * dt / dt_pro;
        gas.SetTemperature(gas.Temperature()+dT);

        // Now set the gas density, calculated using the values above.
        gas.SetDensity(dens);
    }
}
