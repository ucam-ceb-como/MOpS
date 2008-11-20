/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Simulator class declared in the
    mops_simulator.h header file.

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

#include "mops_simulator.h"
#include "mops_ode_solver.h"
#include "mops_flux_postprocessor.h"
#include "mops_reactor_factory.h"
#include "string_functions.h"
#include "csv_io.h"
#include <stdexcept>

using namespace Mops;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Simulator::Simulator(void)
: m_nruns(1), m_niter(1), m_pcount(0), m_maxm0(0.0), 
  m_cpu_start((clock_t)0.0), m_cpu_mark((clock_t)0.0), m_runtime(0.0), 
  m_console_interval(1), m_console_msgs(true), 
  m_output_filename("mops-out"), m_output_every_iter(false),
  m_output_step(0), m_output_iter(0), m_ptrack_count(0)
{
}

// Default destructor.
Simulator::~Simulator(void)
{
}


// SIMULATION SETTINGS.

// Returns the number of runs to perform
unsigned int Simulator::RunCount(void) const {return m_nruns;}

// Sets the number of runs to peform.
void Simulator::SetRunCount(unsigned int n) {m_nruns = n;}

// Returns the number of iteration to perform per step.
unsigned int Simulator::IterCount(void) const {return m_niter;}

// Sets the number of iterations to perform per step.
void Simulator::SetIterCount(unsigned int n) {m_niter = n;}

    
// These particle count related settings should be part of the
// Reactor class!

// Returns the max. stochastic particle count.
unsigned int Simulator::MaxPartCount(void) const {return m_pcount;}

// Sets the max. stochastic particle count.
void Simulator::SetMaxPartCount(unsigned int n) {m_pcount = n;}

// Returns the max. M0, for initial ensemble scaling.
real Simulator::MaxM0(void) const {return m_maxm0;}

// Sets max. M0.
void Simulator::SetMaxM0(real m0) {m_maxm0 = m0;}


// CONSOLE INTERVAL.

unsigned int Simulator::ConsoleInterval() const
{
    return m_console_interval;
}

void Simulator::SetConsoleInterval(unsigned int cint)
{
    m_console_interval = cint;
}


// CONSOLE VARIABLE NAMES.

const std::vector<std::string> &Simulator::ConsoleVariables() const
{
    return m_console_vars;
}

const std::string Simulator::ConsoleVariable(unsigned int i) const
{
    if (i < m_console_vars.size()) {
        return m_console_vars[i];
    } else {
        // Returns the first variable name if index is invalid.     
        return "";
    }
}

void Simulator::AddConsoleVariable(const std::string &var)
{
    m_console_vars.push_back(var);
}

void Simulator::RemoveConsoleVariable(const std::string &var)
{
    vector<string>::iterator i;
    for (i=m_console_vars.begin(); i!=m_console_vars.end(); ++i) {
        if ((*i).compare(var) == 0) {
            m_console_vars.erase(i);
            return;
        }
    }
}

void Simulator::RemoveConsoleVariable(unsigned int i)
{
    if (i < m_console_vars.size()) {
        m_console_vars.erase(m_console_vars.begin()+i);
    }
}


// CONSOLE MESSAGES.

bool Simulator::UseConsoleMsgs() const
{
    return m_console_msgs;
}

void Simulator::SetUseConsoleMsgs(bool msgs)
{
    m_console_msgs = msgs;
}


// OUTPUT FILE NAME.

const std::string &Simulator::OutputFile() const
{
    return m_output_filename;
}

void Simulator::SetOutputFile(const std::string &name)
{
    m_output_filename = name;
}

// Set simulator to output every iteration.
void Simulator::SetOutputEveryIter(bool fout) {m_output_every_iter=fout;}

// STATISTICAL BOUNDS OUTPUT

void Simulator::SetOutputStatBoundary(Sweep::ParticleCache::PropID pid, double lower, double upper) 
{
    m_statbound.Lower = lower;
    m_statbound.Upper = upper;
    m_statbound.PID   = pid;
}

// POVRAY OUTPUT.

void Simulator::SetParticleTrackCount(unsigned int ptcount) {
    m_ptrack_count = ptcount;
}

// SOLVING REACTORS.

// Solves the given reactor for the given time intervals.
void Simulator::RunSimulation(Mops::Reactor &r, 
                              const timevector &times, 
                              Solver &s)
{
    unsigned int icon;
    real dt, t2; // Stop time for each step.

    // Store the initial conditions.
    Mixture initmix(r.Mech()->ParticleMech());
    initmix.SetFracs(r.Mixture()->MoleFractions());
    initmix.SetTemperature(r.Mixture()->Temperature());
    initmix.SetDensity(r.Mixture()->Density());

    // Initialise the reactor with the start time.
    t2 = times[0].StartTime();
    r.SetTime(t2);
    r.Mixture()->Particles().Initialise(m_pcount, r.Mech()->ParticleMech());
    r.Mixture()->SetMaxM0(m_maxm0);

    // Set up the solver.
    s.Initialise(r);

    // Set up file output.
    writeAux(*r.Mech(), times, s);
    openOutputFile();

    // Output initial conditions.
    fileOutput(m_output_step, m_output_iter, r, s, this);

    // Set up the console output.
    icon = m_console_interval;
    setupConsole(*r.Mech());

    // Loop over runs.
    for (unsigned int irun=0; irun!=m_nruns; ++irun) {
        // Start the CPU timing clock.
        m_cpu_start = clock();
        m_runtime  = 0.0;

        // Initialise the reactor with the start time.
        t2 = times[0].StartTime();
        r.Mixture()->SetFracs(initmix.MoleFractions());
        r.Mixture()->SetTemperature(initmix.Temperature());
        r.Mixture()->SetDensity(initmix.Density());
        r.Mixture()->Reset(m_maxm0);
        r.SetTime(t2);

        // Set up the ODE solver for this run.
        s.Reset(r);

        // Print initial conditions to the console.
        printf("mops: Run number %d of %d.\n", irun+1, m_nruns);
        m_console.PrintDivider();
        consoleOutput(r);

        // Loop over the time intervals.
        unsigned int global_step = 0;
        timevector::const_iterator iint;
        for (iint=times.begin(); iint!=times.end(); ++iint) {
            // Get the step size for this interval.
            dt = (*iint).StepSize();

            // Set output parameters for this time interval.
            m_output_step = max((int)iint->SplittingStepCount(), 0);
            m_output_iter = max((int)m_niter, 0);

            // Loop over the steps in this interval.
            unsigned int istep;
            for (istep=0; istep<iint->StepCount(); ++istep, ++global_step) {
                // Run the solver for this step (timed).
                m_cpu_mark = clock();
                    s.Solve(r, t2+=dt, iint->SplittingStepCount(), 
                            m_niter, &fileOutput, (void*)this);
                m_runtime += calcDeltaCT(m_cpu_mark);

                // Generate console output.
                if (--icon == 0) {
                    consoleOutput(r);
                    icon = m_console_interval;
                }
            }

            // Create a save point at the end of this time
            // interval.
            createSavePoint(r, global_step, irun);
        }

        // Print run time to the console.
        printf("mops: Run number %d completed in %.1f s.\n", irun+1, m_runtime);
    }    

    // Close the output files.
    closeOutputFile();
}


// POST-PROCESSING.

void Simulator::PostProcess()
{
    // READ AUXILLIARY INFORMATION.

    // Read auxilliary information about the simulation (mechanism
    // and time intervals).
    Mops::Mechanism mech;
    Mops::timevector times;
    unsigned int ncput = 0;
    vector<string> cput_head; // CPU time headings.
    readAux(mech, times, ncput, cput_head);

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

    // Declare gas-phase rate outputs (averages and errors).
    vector<fvector> agprates(npoints), egprates(npoints);
    vector<fvector> agpfwdrates(npoints), egpfwdrates(npoints);
    vector<fvector> agprevrates(npoints), egprevrates(npoints);
    vector<fvector> agpwdot(npoints), egpwdot(npoints);

    // Declare particle-phase rate outputs (averages and errors).
    vector<fvector> apprates(npoints), epprates(npoints);
    vector<fvector> appwdot(npoints), eppwdot(npoints);

    // Declare data structure for particle tracking data.
    // Runs -> time steps -> particles -> coordinate.
    vector<vector<vector<fvector> > > ptrack(m_nruns);

    // OPEN SIMULATION OUTPUT FILE.

    // Build the simulation input file name.
    string fname = m_output_filename + ".sim";

    // Open the simulation output file.
    fstream fin(fname.c_str(), ios_base::in | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!fin.good()) {
        throw std::runtime_error("Failed to open file for post-processing "
                                 "(Mops, Simulator::PostProcess).");
    }

    // READ INITIAL CONDITIONS.

    // The initial conditions were only written to the file
    // once, as they are the same for all runs.
    readGasPhaseDataPoint(fin, mech, achem[0], echem[0], true);
    readParticleDataPoint(fin, pmech, astat[0], estat[0], true);
    readGasRxnDataPoint(fin, mech,
                        agprates[0], egprates[0],
                        agpfwdrates[0], egpfwdrates[0],
                        agprevrates[0], egprevrates[0],
                        agpwdot[0], egpwdot[0],
                        true);
    readPartRxnDataPoint(fin, mech.ParticleMech(),
                         apprates[0], epprates[0],
                         appwdot[0], eppwdot[0],
                         true);
    readCTDataPoint(fin, ncput, acpu[0], ecpu[0], true);
    for(unsigned int irun=0; irun!=m_nruns; ++irun) {
        ptrack[irun].resize(npoints); // Resize particle tracking vector.
    }
    readPartTrackPoint(fin, pmech, ptrack[0][0]);
    for(unsigned int irun=1; irun!=m_nruns; ++irun) {
        ptrack[irun][0].assign(ptrack[0][0].begin(), ptrack[0][0].end());
    }

    // The initial conditions need to be multiplied by the number
    // of runs and iterations in order for the calcAvgConf to give 
    // the correct result.
    if (m_output_every_iter) {
        multVals(achem[0], m_nruns*m_niter);
        multVals(astat[0], m_nruns*m_niter);
        multVals(agprates[0], m_nruns*m_niter);
        multVals(agpfwdrates[0], m_nruns*m_niter);
        multVals(agprevrates[0], m_nruns*m_niter);
        multVals(agpwdot[0], m_nruns*m_niter);
        multVals(apprates[0], m_nruns*m_niter);
        multVals(appwdot[0], m_nruns*m_niter);
        multVals(acpu[0], m_nruns*m_niter);
        multVals(echem[0], m_nruns*m_niter);
        multVals(estat[0], m_nruns*m_niter);
        multVals(egprates[0], m_nruns*m_niter);
        multVals(egpfwdrates[0], m_nruns*m_niter);
        multVals(egprevrates[0], m_nruns*m_niter);
        multVals(egpwdot[0], m_nruns*m_niter);
        multVals(epprates[0], m_nruns*m_niter);
        multVals(eppwdot[0], m_nruns*m_niter);
        multVals(ecpu[0], m_nruns*m_niter);
    } else {
        multVals(achem[0], m_nruns);
        multVals(astat[0], m_nruns);
        multVals(agprates[0], m_nruns);
        multVals(agpfwdrates[0], m_nruns);
        multVals(agprevrates[0], m_nruns);
        multVals(agpwdot[0], m_nruns);
        multVals(apprates[0], m_nruns);
        multVals(appwdot[0], m_nruns);
        multVals(acpu[0], m_nruns);
        multVals(echem[0], m_nruns);
        multVals(estat[0], m_nruns);
        multVals(egprates[0], m_nruns);
        multVals(egpfwdrates[0], m_nruns);
        multVals(egprevrates[0], m_nruns);
        multVals(egpwdot[0], m_nruns);
        multVals(epprates[0], m_nruns);
        multVals(eppwdot[0], m_nruns);
        multVals(ecpu[0], m_nruns);
    }

    // READ ALL OUTPUT POINTS.

    // Now we must read the reactor conditions at all time points
    // and all runs.
    for(unsigned int irun=0; irun!=m_nruns; ++irun) {
        // Loop over all time intervals.
        unsigned int step = 1;
        for (Mops::timevector::const_iterator iint=times.begin(); 
             iint!=times.end(); ++iint) {
            // Loop over all time steps in this interval.
            for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
                if (m_output_every_iter) {
                    // Loop over all iterations.
                    for(unsigned int iiter=0; iiter!=m_niter; ++iiter) {
                        // Read output point for all iterations for step.
                        readGasPhaseDataPoint(fin, mech, achem[step], echem[step], true);
                        readParticleDataPoint(fin, pmech, astat[step], estat[step], true);
                        readGasRxnDataPoint(fin, mech,
                                            agprates[step], egprates[step],
                                            agpfwdrates[step], egpfwdrates[step],
                                            agprevrates[step], egprevrates[step],
                                            agpwdot[step], egpwdot[step],
                                            true);
                        readPartRxnDataPoint(fin, mech.ParticleMech(), apprates[step], epprates[step], appwdot[step], eppwdot[step], true);
                        readCTDataPoint(fin, ncput, acpu[step], ecpu[step], true);
                        readPartTrackPoint(fin, pmech, ptrack[irun][step]);
                    }
                } else {
                    // Read single output point for step.
                    readGasPhaseDataPoint(fin, mech, achem[step], echem[step], true);
                    readParticleDataPoint(fin, pmech, astat[step], estat[step], true);
                    readGasRxnDataPoint(fin, mech,
                                        agprates[step], egprates[step],
                                        agpfwdrates[step], egpfwdrates[step],
                                        agprevrates[step], egprevrates[step],
                                        agpwdot[step], egpwdot[step],
                                        true);
                    readPartRxnDataPoint(fin, mech.ParticleMech(), apprates[step], epprates[step], appwdot[step], eppwdot[step], true);
                    readCTDataPoint(fin, ncput, acpu[step], ecpu[step], true);
                    readPartTrackPoint(fin, pmech, ptrack[irun][step]);
                }
            }
        }

    }

    // Close the simulation output file.
    fin.close();

    // CALCULATE AVERAGES AND CONFIDENCE INTERVALS.
    if (m_output_every_iter) {
        calcAvgConf(achem, echem, m_nruns*m_niter);
        calcAvgConf(astat, estat, m_nruns*m_niter);
        calcAvgConf(agprates, egprates, m_nruns*m_niter);
        calcAvgConf(agpfwdrates, egpfwdrates, m_nruns*m_niter);
        calcAvgConf(agprevrates, egprevrates, m_nruns*m_niter);
        calcAvgConf(agpwdot, egpwdot, m_nruns*m_niter);
        calcAvgConf(apprates, epprates, m_nruns*m_niter);
        calcAvgConf(appwdot, eppwdot, m_nruns*m_niter);
        calcAvgConf(acpu, ecpu, m_nruns*m_niter);
    } else {
        calcAvgConf(achem, echem, m_nruns);
        calcAvgConf(astat, estat, m_nruns);
        calcAvgConf(agprates, egprates, m_nruns);
        calcAvgConf(agpfwdrates, egpfwdrates, m_nruns);
        calcAvgConf(agprevrates, egprevrates, m_nruns);
        calcAvgConf(agpwdot, egpwdot, m_nruns);
        calcAvgConf(apprates, epprates, m_nruns);
        calcAvgConf(appwdot, eppwdot, m_nruns);
        calcAvgConf(acpu, ecpu, m_nruns);
    }

    // POST-PROCESS ELEMENT FLUX
    // Element flux must be output before CSV files since writeXXXCSV will change the contents of what it output afterwards
    writeElementFluxOutput(m_output_filename, mech, times, agpfwdrates, agprevrates, achem);

    // OUTPUT TO CSV FILES.

    writeGasPhaseCSV(m_output_filename+"-chem.csv", mech, times, achem, echem);
    writeParticleStatsCSV(m_output_filename+"-part.csv", mech, times, astat, estat);
    writeGasRxnCSV(m_output_filename+"-gp-rates.csv", mech, times, agprates, egprates);
    //writeGasRxnCSV(m_output_filename+"-gp-fwd-rates.csv", mech, times, agpfwdrates, egpfwdrates);
    //writeGasRxnCSV(m_output_filename+"-gp-rev-rates.csv", mech, times, agprevrates, egprevrates);
    writeProdRatesCSV(m_output_filename+"-gp-wdot.csv", mech, times, agpwdot, egpwdot);
    writePartProcCSV(m_output_filename+"-part-rates.csv", mech.ParticleMech(), times, apprates, epprates);
    writeProdRatesCSV(m_output_filename+"-part-wdot.csv", mech, times, appwdot, eppwdot);
    writeCT_CSV(m_output_filename+"-cput.csv", times, acpu, ecpu, cput_head);
    for(unsigned int irun=0; irun!=m_nruns; ++irun) {
        writePartTrackCSV(m_output_filename+"("+cstr(irun)+")-track", mech, 
                          times, ptrack[irun]);
    }

    // POST-PROCESS PSLs.

    // Now post-process the PSLs.
    postProcessPSLs(mech, times);
}


// FILE OUTPUT (protected).

// Opens the simulation output file.
void Simulator::openOutputFile() const
{
    // SIMULATOR OUTPUT FILE
    // Build the simulation output file name.
    string fname = m_output_filename + ".sim";

    // Open the simulation output file.
    m_file.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!m_file.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "output (Mops, Simulator::openOutputFile).");
    }

    // SIMULATOR SENSITIVITY OUTPUT FILE
    // Build the simulation output file name.
    string fsenname = m_output_filename + ".sen";

    // Open the sensitivity simulation output file.
    m_senfile.open(fsenname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!m_senfile.good()) {
        throw runtime_error("Failed to open file for sensitivity simulation "
                            "output (Mops, Simulator::openOutputFile).");
    }
}

// Closes the output file.
void Simulator::closeOutputFile() const
{
    // Close the simulation output file.
    m_file.close();
    // Close the sensitivity output file.
    m_senfile.close();
}

// Writes the gas-phase conditions of the given reactor to
// the binary output file.
void Simulator::outputGasPhase(const Reactor &r) const
{
    // Write gas-phase conditions to file.
    m_file.write((char*)&r.Mixture()->RawData()[0], 
                 sizeof(r.Mixture()->RawData()[0]) *
                 r.Mech()->SpeciesCount());
    real T = r.Mixture()->Temperature();
    m_file.write((char*)&T, sizeof(T));
    real D = r.Mixture()->Density();
    m_file.write((char*)&D, sizeof(D));
    real P = r.Mixture()->Pressure();
    m_file.write((char*)&P, sizeof(P));
}

// Writes the particle stats to the binary output file.
void Simulator::outputParticleStats(const Reactor &r) const
{
    // Write particle stats to file.
    static Sweep::Stats::EnsembleStats stats(r.Mech()->ParticleMech());
    stats.SetStatBoundary(m_statbound);
    r.Mixture()->GetVitalStats(stats);
    stats.Serialize(m_file);
}

// Writes tracked particles to the binary output file.
void Simulator::outputPartTrack(const Reactor &r) const
{
    // Write the number of tracked particles.
    unsigned int n = min(r.Mixture()->ParticleCount(), m_ptrack_count);
    m_file.write((char*)&n, sizeof(n));
    
    // Output the current time.
    double t = (double)r.Time();
    m_file.write((char*)&t, sizeof(t));
    
    if (n > 0) {
        // Serialize the particles.
        for (unsigned int i=0; i!=n; ++i) {
            r.Mixture()->Particles().At(i)->Serialize(m_file);
        }
    }
}
// Write sensitivity output to the binary file.
void Simulator::outputSensitivity(const Reactor &r) const
{
}

// Writes the gas-phase reaction rates-of-progress and the
// species molar production rates to the binary output file.
void Simulator::outputGasRxnRates(const Reactor &r) const
{
    // Calculate the rates-of-progress.
    static fvector rop, rfwd, rrev;
    r.Mech()->Reactions().GetRatesOfProgress(*r.Mixture(), rop, rfwd, rrev);

    // Calculate the molar production rates.
    static fvector wdot;
    r.Mech()->Reactions().GetMolarProdRates(rop, wdot);

    // Write rates to the file.
    m_file.write((char*)&rop[0], sizeof(rop[0]) * r.Mech()->ReactionCount());
    m_file.write((char*)&rfwd[0], sizeof(rfwd[0]) * r.Mech()->ReactionCount());
    m_file.write((char*)&rrev[0], sizeof(rrev[0]) * r.Mech()->ReactionCount());
    m_file.write((char*)&wdot[0], sizeof(wdot[0]) * r.Mech()->SpeciesCount());
}

// Writes the particle process rates and the
// species molar production rates to the binary output file.
void Simulator::outputPartRxnRates(const Reactor &r) const
{
    if (r.Mech()->ParticleMech().ProcessCount() != 0) {
        // Calculate the process rates.
        static fvector rates;
        r.Mech()->ParticleMech().CalcRates(r.Time(), *r.Mixture(), rates);

        // Calculate the molar production rates (mol/mol).
        static fvector wdot;
        r.Mech()->ParticleMech().CalcGasChangeRates(r.Time(), *r.Mixture(), wdot);

        // Now convert from mol/mol to mol/m3.
        fvector::iterator rhodot = wdot.begin()+r.Mech()->SpeciesCount()+1;
        for (unsigned int k=0; k!=r.Mech()->SpeciesCount(); ++k) {
            wdot[k] = (r.Mixture()->Density() * wdot[k]) + 
                      (r.Mixture()->MoleFraction(k) * (*rhodot));
        }

        // Write rates to the file.
        m_file.write((char*)&rates[0], sizeof(rates[0]) * r.Mech()->ParticleMech().ProcessCount());
        m_file.write((char*)&wdot[0], sizeof(wdot[0]) * r.Mech()->SpeciesCount());
    }
}


// FILE OUTPUT.

// Writes a reactor and solver variables to the output file.
void Simulator::fileOutput(unsigned int step, unsigned int iter, 
                           const Mops::Reactor &r, const Solver &s,
                           void *sim)
{
    // Cast the void pointer to a Simulator object.
    Simulator *me = static_cast<Simulator*>(sim);

    if (step == me->m_output_step) {
        if (me->m_output_every_iter || (iter == me->m_output_iter)) {
            // Write the gas-phase conditions to the output file.
            me->outputGasPhase(r);

            // Write particle stats to file.
            me->outputParticleStats(r);

            // Write gas-phase and particle reaction rates
            // to the file.
            me->outputGasRxnRates(r);
            me->outputPartRxnRates(r);

            // Write CPU times to file.
            s.OutputCT(me->m_file);

            // Do particle tracking output.
            me->outputPartTrack(r);

            // Write sensitivityto file.
            //s.GetODE_Solver().GetSensitivity().outputTo();
        }
    }
}


// CONSOLE OUTPUT (protected).

// Sets up console output.
void Simulator::setupConsole(const Mops::Mechanism &mech)
{
    vector<string> header;
    m_console_mask.clear();

    // Loop over the column variable names.
    vector<string>::const_iterator i;
    for (i=m_console_vars.begin(); i!=m_console_vars.end(); ++i) {
        if ((*i).compare("RHO")==0 || (*i).compare("rho")==0 || 
            (*i).compare("Rho")==0) {
            // This variable is the mixture density.
            header.push_back("Density");
            m_console_mask.push_back(mech.Species().size()+1);
        } else if ((*i).compare("T")==0) {
            // This variable is the mixture temperature.
            header.push_back("T (K)");
            m_console_mask.push_back(mech.Species().size());
        } else if ((*i).compare("TIME")==0 || (*i).compare("time")==0 || 
                   (*i).compare("Time")==0) {
            // This variable is the flow time.
            header.push_back("Time (s)");
            m_console_mask.push_back(mech.Species().size()+2);
        } else if ((*i).compare("#SP")==0 || (*i).compare("#sp")==0) {
            // Number of stochastic particles.
            header.push_back("#SP");
            m_console_mask.push_back(mech.Species().size()+3);
        } else if ((*i).compare("M0")==0 || (*i).compare("m0")==0) {
            // Particle number density.
            header.push_back("M0 (cm-3)");
            m_console_mask.push_back(mech.Species().size()+4);
        } else if ((*i).compare("FV")==0 || (*i).compare("fv")==0 || 
                   (*i).compare("Fv")==0) {
            // Particle volume fraction.
            header.push_back("Fv");
            m_console_mask.push_back(mech.Species().size()+5);
        } else if ((*i).compare("CT")==0 || (*i).compare("ct")==0) {
            // Computation time.
            header.push_back("CPU (s)");
            m_console_mask.push_back(mech.Species().size()+6);
        } else {
            // Check for a species name.
            int isp = mech.FindSpecies((*i));
            if (isp >= 0) {
                // This is a valid species name.
                header.push_back((*i) + " (mol/m3)");
                m_console_mask.push_back(isp);
            } else {
                // This is an invalid variable.  Just print the first species.
                header.push_back(mech.Species()[0]->Name());
                m_console_mask.push_back(0);
            }
        }
    }

    // Print the first header.
    m_console.PrintDivider();
    m_console.PrintRow(header);
    m_console.PrintDivider();

    // Set up the console output class.
    m_console.SetAutoHeader(header);
    m_console.EnableAutoDividers();
}

// Writes output to the console.
void Simulator::consoleOutput(const Mops::Reactor &r) const
{
    // Get output data from gas-phase.
    static vector<real> out;
    r.Mixture()->GetConcs(out);
    out.push_back(r.Mixture()->Temperature());
    out.push_back(r.Mixture()->Density());
    out.push_back(r.Time());

    // Get output data from particles.
    Sweep::Stats::EnsembleStats stats(r.Mech()->ParticleMech());
    stats.SetStatBoundary(m_statbound);
    r.Mixture()->GetVitalStats(stats);
    out.push_back(stats.BasicStats().PCount());
    out.push_back(stats.BasicStats().M0());
    out.push_back(stats.BasicStats().Fv());

    // Get output CPU time.
    double cputime = calcDeltaCT(m_cpu_start);
    out.push_back(cputime);

    // Print data to console.
    m_console.PrintRow(out, m_console_mask);
}


// POST-PROCESSING ROUTINES (protected).

// Writes the auxilliary post-processing information using
// the given file name.  This information is the chemical mechanism
// and the output time intervals.
void Simulator::writeAux(const Mops::Mechanism &mech, 
                         const Mops::timevector &times,
                         const Mops::Solver &solv) const
{
    // Build output file name.
    string fname = m_output_filename + ".aux";

    // Output a file for post-processing information (mechanism and
    // time intervals).
    fstream fout;
    fout.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!fout.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "output (Mops, Simulator::writeAux).");
    }

    // Write the mechanism to the output file.
    mech.Serialize(fout);

    // Write the time intervals to the output file.
    unsigned int n = times.size();
    fout.write((char*)&n, sizeof(n));
    for (timevector::const_iterator i=times.begin(); i!=times.end(); i++) {
        (*i).Serialize(fout);
    }

    // Write simulator settings.
    Serialize(fout);

    // Write number of CPU times generated by solver.
    unsigned int ncput = solv.CT_Count();
    fout.write((char*)&ncput, sizeof(ncput));

    // Write CPU time headings.
    vector<string> head;
    solv.CT_Names(head);
    for (unsigned int i=0; i!=ncput; ++i) {
        // Write length of heading string.
        unsigned int n = head[i].length();
        fout.write((char*)&n, sizeof(n));
        // write heading string.
        fout.write(head[i].c_str(), n);
    }

    // Close the post-processing info file.
    fout.close();
}
    
// Reads auxilliary post-processing information using the
// given file name.  This information is the chemical mechanism
// and the output time intervals.
void Simulator::readAux(Mops::Mechanism &mech, 
                        Mops::timevector &times,
                        unsigned int &ncput,
                        std::vector<std::string> &cput_head)
{
    // Build input file name.
    string fname = m_output_filename + ".aux";

    // Output the file which contains the simulation mechanism and
    // time intervals.
    fstream fin;
    fin.open(fname.c_str(), ios_base::in | ios_base::binary);

    // Throw error if the file failed to open.
    if (!fin.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "post-processing (Mops, Simulator::readAux).");
    }

    // Read the mechanism from the file.
    mech.Deserialize(fin);

    // Read the time intervals from the file.
    times.clear();
    unsigned int ntimes;
    fin.read(reinterpret_cast<char*>(&ntimes), sizeof(ntimes));
    for (unsigned int i=0; i<ntimes; i++) {
        times.push_back(Mops::TimeInterval(fin));
    }

    // Read simulator settings.
    Deserialize(fin);

    // Read number of CPU times generated by solver.
    fin.read(reinterpret_cast<char*>(&ncput), sizeof(ncput));

    // Read CPU time headings.
    cput_head.resize(ncput, "");
    char *str = NULL;
    unsigned int n = 0;
    for (unsigned int i=0; i!=ncput; ++i) {
        // Read length of heading string.
        fin.read(reinterpret_cast<char*>(&n), sizeof(n));
        // Resize temporary character array.
        str = new char[n];
        // Read heading string.
        fin.read(str, n);
        cput_head[i] = string(str, n);
        // Delete temporary character array.
        delete [] str;
    }

    // Close the simulation settings file.
    fin.close();
}

// Reads a gas-phase chemistry data point from the binary file.
// To allow the averages and confidence intervals to be calculated
// the data point is added to a vector of sums, and the squares are
// added to the vector sumsqr if necessary.
void Simulator::readGasPhaseDataPoint(std::istream &in, const Mops::Mechanism &mech,
                                      fvector &sum, fvector &sumsqr, bool calcsqrs)
{
    // Check for valid stream.
    if (in.good()) {
        unsigned int N = mech.SpeciesCount();

        // Read the gas-phase conditions.
        fvector y(N, 0.0);
        real T=0.0, D=0.0, P=0.0;

        in.read(reinterpret_cast<char*>(&y[0]), sizeof(y[0])*N);
        in.read(reinterpret_cast<char*>(&T), sizeof(T));
        in.read(reinterpret_cast<char*>(&D), sizeof(D));
        D *= 1.0e-6; // Convert density from m^3 to cm^3.
        in.read(reinterpret_cast<char*>(&P), sizeof(P));

        // Resize vectors.
        sum.resize(N+3, 0.0);
        sumsqr.resize(N+3, 0.0);

        // Calculate sums and sums of square (for average and
        // error calculation).
        for (unsigned int i=0; i!=N; ++i) {
            // Note conversion to cm-3 from m-3.
            real conc = D * y[i];
            sum[i] += conc;
            if (calcsqrs) sumsqr[i] += conc * conc;
        }

        // Calculates sums and sums of squares of temperature, density
        // and pressure.
        sum[N]   += T;
        sum[N+1] += D; // Note conversion to cm-3 from m-3.
        sum[N+2] += P;
        if (calcsqrs) {
            sumsqr[N]   += (T*T);
            sumsqr[N+1] += (D*D);
            sumsqr[N+2] += (P*P);
        }
    }
}

// Reads a CPU timing data from the binary file.
// To allow the averages and confidence intervals to be calculated
// the data point is added to a vector of sums, and the squares are
// added to the vector sumsqr if necessary.
void Simulator::readCTDataPoint(std::istream &in, unsigned int N,
                                fvector &sum, fvector &sumsqr, 
                                bool calcsqrs)
{
    // Check for valid stream.
    if (in.good()) {
        // Read the computation times.
        real *cpu = new real[N];
        in.read(reinterpret_cast<char*>(&cpu[0]), sizeof(cpu[0])*N);

        // Resize output vectors.
        sum.resize(N, 0.0);
        sumsqr.resize(N, 0.0);

        // Calculate sums and sums of square (for average and
        // error calculation).
        for (unsigned int i=0; i!=N; ++i) {
            sum[i] += cpu[i];
            if (calcsqrs) sumsqr[i] += (cpu[i] * cpu[i]);
        }

        delete [] cpu;
    }
}

// Reads a particle stats data point from the binary file.
// To allow the averages and confidence intervals to be calculated
// the data point is added to a vector of sums, and the squares are
// added to the vector sumsqr if necessary.
void Simulator::readParticleDataPoint(std::istream &in, 
                                      const Sweep::Mechanism &mech, 
                                      fvector &sum, fvector &sumsqr, 
                                      bool calcsqrs)
{
    // Check for valid stream.
    if (in.good()) {
        // Read the stats.
        Sweep::Stats::EnsembleStats stats(mech);
        stats.Deserialize(in, mech);

        // Get the stats vector.
        fvector s;
        stats.Get(s);

        // Resize vectors.
        sum.resize(s.size(), 0.0);
        sumsqr.resize(s.size(), 0.0);

        // Calculate sums and sums of squares (for average and
        // error calculation).
        for (unsigned int i=0; i!=s.size(); ++i) {
            sum[i] += s[i];
            if (calcsqrs) sumsqr[i] += (s[i] * s[i]);
        }
    }
}

// Reads a gasp-phase reaction rates stats data point from the binary file.
// To allow the averages and confidence intervals to be calculated
// the data point is added to a vector of sums, and the squares are
// added to the vector sumsqr if necessary.
void Simulator::readGasRxnDataPoint(std::istream &in, const Mops::Mechanism &mech,
                                    fvector &rates_sum, fvector &rates_sumsqr,
                                    fvector &fwd_rates_sum, fvector &fwd_rates_sumsqr,
                                    fvector &rev_rates_sum, fvector &rev_rates_sumsqr,
                                    fvector &wdot_sum, fvector &wdot_sumsqr,
                                    bool calcsqrs)
{
    // Check for valid stream.
    if (in.good()) {
        // Get the reaction rates-of-progress vector.
        fvector rop(mech.ReactionCount());
        in.read(reinterpret_cast<char*>(&rop[0]), sizeof(rop[0])*mech.ReactionCount());

        // Get the forward reaction rates vector.
        fvector rfwd(mech.ReactionCount());
        in.read(reinterpret_cast<char*>(&rfwd[0]), sizeof(rop[0])*mech.ReactionCount());

        // Get the reverse reaction rates vector.
        fvector rrev(mech.ReactionCount());
        in.read(reinterpret_cast<char*>(&rrev[0]), sizeof(rop[0])*mech.ReactionCount());

        // Get the species molar production rates.
        fvector wdot(mech.SpeciesCount());
        in.read(reinterpret_cast<char*>(&wdot[0]), sizeof(wdot[0])*mech.SpeciesCount());

        // Resize vectors.
        rates_sum.resize(rop.size(), 0.0);
        rates_sumsqr.resize(rop.size(), 0.0);
        fwd_rates_sum.resize(rfwd.size(), 0.0);
        fwd_rates_sumsqr.resize(rfwd.size(), 0.0);
        rev_rates_sum.resize(rrev.size(), 0.0);
        rev_rates_sumsqr.resize(rrev.size(), 0.0);
        wdot_sum.resize(wdot.size(), 0.0);
        wdot_sumsqr.resize(wdot.size(), 0.0);

        // Calculate sums and sums of squares (for average and
        // error calculation).
        for (unsigned int i=0; i!=rop.size(); ++i) {
            rates_sum[i] += rop[i];
            if (calcsqrs) rates_sumsqr[i] += (rop[i] * rop[i]);
        }
        for (unsigned int i=0; i!=rfwd.size(); ++i) {
            fwd_rates_sum[i] += rfwd[i];
            if (calcsqrs) fwd_rates_sumsqr[i] += (rfwd[i] * rfwd[i]);
        }
        for (unsigned int i=0; i!=rrev.size(); ++i) {
            rev_rates_sum[i] += rrev[i];
            if (calcsqrs) rev_rates_sumsqr[i] += (rrev[i] * rrev[i]);
        }
        for (unsigned int i=0; i!=wdot.size(); ++i) {
            wdot_sum[i] += wdot[i];
            if (calcsqrs) wdot_sumsqr[i] += (wdot[i] * wdot[i]);
        }
    }
}

// Reads a particle process rates stats data point from the binary file.
// To allow the averages and confidence intervals to be calculated
// the data point is added to a vector of sums, and the squares are
// added to the vector sumsqr if necessary.
void Simulator::readPartRxnDataPoint(std::istream &in, const Sweep::Mechanism &mech,
                                     fvector &rates_sum, fvector &rates_sumsqr,
                                     fvector &wdot_sum, fvector &wdot_sumsqr,
                                     bool calcsqrs)
{
    // Check for valid stream.
    if (in.good() && (mech.ProcessCount() > 0)) {
        // Get the process rates vector.
        fvector rates(mech.ProcessCount());
        in.read(reinterpret_cast<char*>(&rates[0]),
                sizeof(rates[0])*mech.ProcessCount());

        // Get the species molar production rates.
        fvector wdot(mech.Species()->size());
        in.read(reinterpret_cast<char*>(&wdot[0]),
                sizeof(wdot[0])*mech.Species()->size());

        // Resize vectors.
        rates_sum.resize(rates.size(), 0.0);
        rates_sumsqr.resize(rates.size(), 0.0);
        wdot_sum.resize(wdot.size(), 0.0);
        wdot_sumsqr.resize(wdot.size(), 0.0);

        // Calculate sums and sums of squares (for average and
        // error calculation).
        for (unsigned int i=0; i!=rates.size(); ++i) {
            rates_sum[i] += rates[i];
            if (calcsqrs) rates_sumsqr[i] += (rates[i] * rates[i]);
        }
        for (unsigned int i=0; i!=wdot.size(); ++i) {
            wdot_sum[i] += wdot[i];
            if (calcsqrs) wdot_sumsqr[i] += (wdot[i] * wdot[i]);
        }
    }
}

// Reads the tracked particles from the binary file.  The particles are
// processed so that only a vector of vectors is returned, which contains
// the PSL data for each tracked particle at that point.
void Simulator::readPartTrackPoint(std::istream &in, 
                                   const Sweep::Mechanism &mech,
                                   std::vector<fvector> &pdata) const
{
    // Check for valid stream.
    if (in.good()) {
        vector<fvector>::iterator i;

        // Read the number of tracked particles.
        unsigned int n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(n));

        // Read the output time.
        double t = 0.0;
        in.read(reinterpret_cast<char*>(&t), sizeof(t));

        // Create a particle stats object.
        Sweep::Stats::EnsembleStats stats(mech);

        // Resize output vectors to hold particle data.  Also reset
        // all entries to 0.0.
        if (pdata.size() < m_ptrack_count) pdata.resize(m_ptrack_count);
        for (i=pdata.begin(); i!=pdata.end(); ++i) {
            fill(i->begin(), i->end(), 0.0);
            i->resize(stats.PSL_Count(), 0.0);
        }
        
        // Read the tracked particles and retrieve the PSL stats
        // using the EnsembleStats class.
        i = pdata.begin();
        for (unsigned int j=0; j!=n; ++j) {
            Sweep::Particle sp(in, mech);
            stats.PSL(sp, (real)t, *(i++));
        }
        Sweep::Particle empty(t, mech);
        for (unsigned int j=n; j!=m_ptrack_count; ++j)
        {
            stats.PSL(empty, (real)t, *(i++));
        }
    }
}

// Multiplies all values in a vector by a scaling factor.
void Simulator::multVals(fvector &vals, real scale)
{
    for (fvector::iterator i=vals.begin(); i!=vals.end(); ++i) {
        (*i) *= scale;
    }
}

// Takes vectors of vectors of variable sums and sums of squares, which
// are converted into the average values and the confidence intervals.
void Simulator::calcAvgConf(std::vector<fvector> &avg, 
                            std::vector<fvector> &err, 
                            unsigned int nruns)
{
    const real CONFA = 3.29; // for 99.9% confidence interval.

    // Pre-calc some useful values.
    real invruns = 1.0 / (real)nruns;
//    real invruns_1 = 1.0 / (real)(nruns-1);
    unsigned int npoints = avg.size();

    // Loop over all steps and all variables.
    for (unsigned int step=0; step!=npoints; ++step) {
        for (unsigned int i=0; i!=avg[step].size(); ++i) {
            // Calculate average over all runs.
            avg[step][i] *= invruns;

            if (nruns > 1) {
                // Calculate standard deviation over all runs.
                (err[step][i] *= invruns) -= (avg[step][i] * avg[step][i]);
                // Calculate confidence interval.
                err[step][i] = CONFA * sqrt(abs(err[step][i] * invruns));
            } else {
                err[step][i] = 0.0;
            }
        }
    }
}

// Takes a vector of average values and a vector with the confidence
// bound of each variable and combines them into a single vector:
// BEFORE:
// avg = (a1, a2, a3, ..., aN)
// err = (e1, e2, e3, ..., eN)
// AFTER:
// avg = (a1, e1, a2, e2, a3, e3, ..., aN, eN)
// The step number and time are insert at the beginning of the avg
// vector.
void Simulator::buildOutputVector(unsigned int step, real time, 
                                  fvector &avg, const fvector &err)
{
    for (unsigned int i=avg.size(); i!=0; --i) {
        if (i < err.size()) {
            avg.insert(avg.begin()+i, *(err.begin()+i-1));
        } else {
            avg.insert(avg.begin()+i, 0.0);
        }
    }
    avg.insert(avg.begin(), time);
    avg.insert(avg.begin(), step);
}

// Writes gas-phase conditions profile to a CSV file.
void Simulator::writeGasPhaseCSV(const std::string &filename, 
                                 const Mechanism &mech, 
                                 const timevector &times, 
                                 std::vector<fvector> &avg, 
                                 const std::vector<fvector> &err)
{
    // Open file for the CSV results.
    CSV_IO csv(filename, true);

    // Write the header row to the gas-phase chemistry CSV file.
    vector<string> head;
    head.push_back("Step");
    head.push_back("Time (s)");
    for (unsigned int isp=0; isp<mech.SpeciesCount(); ++isp) {
        head.push_back(mech.Species(isp)->Name() + " (mol/cm3)");
    }
    head.push_back("T (K)");
    head.push_back("Density (mol/cm3)");
    head.push_back("Pressure (Pa)");
    for (unsigned int i=head.size(); i!=2; --i) {
        head.insert(head.begin()+i, "Err");
    }
    csv.Write(head);

    // Output initial conditions.
    buildOutputVector(0, times[0].StartTime(), avg[0], err[0]);
    csv.Write(avg[0]);

    // Loop over all points, performing output.
    unsigned int step = 1;
    for (timevector::const_iterator iint=times.begin(); iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        real t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            buildOutputVector(step, t, avg[step], err[step]);
            csv.Write(avg[step]);
        }
    }

    // Close the CSV files.
    csv.Close();
}

// Writes computation times profile to a CSV file.
void Simulator::writeCT_CSV(const std::string &filename, 
                            const timevector &times, 
                            std::vector<fvector> &avg,
                            const std::vector<fvector> &err,
                            const std::vector<std::string> &head) const
{
    // Open file for the CSV results.
    CSV_IO csv(filename, true);

    // Write the header row to the CPU time CSV file.
    vector<string> hdr;
    hdr.assign(head.begin(), head.end());
    hdr.insert(hdr.begin(), "Time (s)");
    hdr.insert(hdr.begin(), "Step");
    for (unsigned int i=hdr.size(); i!=2; --i) {
        hdr.insert(hdr.begin()+i, "Err");
    }
    csv.Write(hdr);

    // Output initial conditions.
    buildOutputVector(0, times[0].StartTime(), avg[0], err[0]);
    csv.Write(avg[0]);

    // Loop over all points, performing output.
    unsigned int step = 1;
    for (timevector::const_iterator iint=times.begin(); iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        real t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            buildOutputVector(step, t, avg[step], err[step]);
            csv.Write(avg[step]);
        }
    }

    // Close the CSV files.
    csv.Close();
}

// Writes particle stats profile to a CSV file.
void Simulator::writeParticleStatsCSV(const std::string &filename, 
                                      const Mechanism &mech, 
                                      const timevector &times, 
                                      std::vector<fvector> &avg, 
                                      const std::vector<fvector> &err)
{
    // Open file for the CSV results.
    CSV_IO csv(filename, true);

    Sweep::Stats::EnsembleStats stats(mech.ParticleMech());

    // Write the header row to the particle stats CSV file.
    vector<string> head;
    head.push_back("Step");
    head.push_back("Time (s)");
    stats.Names(head, 2);
    for (unsigned int i=head.size(); i!=2; --i) {
        head.insert(head.begin()+i, "Err");
    }
    csv.Write(head);

    // Output initial conditions.
    buildOutputVector(0, times[0].StartTime(), avg[0], err[0]);
    csv.Write(avg[0]);

    // Loop over all points, performing output.
    unsigned int step = 1;
    for (timevector::const_iterator iint=times.begin(); iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        real t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            buildOutputVector(step, t, avg[step], err[step]);
            csv.Write(avg[step]);
        }
    }

    // Close the CSV files.
    csv.Close();
}

// Writes gas-phase reaction rates profile to a CSV file.
void Simulator::writeGasRxnCSV(const std::string &filename, 
                               const Mechanism &mech, 
                               const timevector &times, 
                               std::vector<fvector> &avg, 
                               const std::vector<fvector> &err)
{
    // Open file for the CSV results.
    CSV_IO csv(filename, true);

    // Write the header row to the gas-phase chemistry CSV file.
    vector<string> head;
    head.push_back("Step");
    head.push_back("Time (s)");
    for (unsigned int i=0; i<mech.ReactionCount(); ++i) {
        head.push_back("Rxn " + cstr(i) + " (mol/cm3s)");
    }
    for (unsigned int i=head.size(); i!=2; --i) {
        head.insert(head.begin()+i, "Err");
    }
    csv.Write(head);

    // Output initial conditions.
    buildOutputVector(0, times[0].StartTime(), avg[0], err[0]);
    csv.Write(avg[0]);

    // Loop over all points, performing output.
    unsigned int step = 1;
    for (timevector::const_iterator iint=times.begin(); iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        real t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            buildOutputVector(step, t, avg[step], err[step]);
            csv.Write(avg[step]);
        }
    }

    // Close the CSV files.
    csv.Close();
}

// Writes gas-phase reaction rates profile to a CSV file.
void Simulator::writeProdRatesCSV(const std::string &filename, 
                                  const Mechanism &mech, 
                                  const timevector &times, 
                                  std::vector<fvector> &avg, 
                                  const std::vector<fvector> &err)
{
    // Open file for the CSV results.
    CSV_IO csv(filename, true);

    // Write the header row to the gas-phase chemistry CSV file.
    vector<string> head;
    head.push_back("Step");
    head.push_back("Time (s)");
    for (unsigned int isp=0; isp<mech.SpeciesCount(); ++isp) {
        head.push_back(mech.Species(isp)->Name() + " (mol/cm3s)");
    }
    for (unsigned int i=head.size(); i!=2; --i) {
        head.insert(head.begin()+i, "Err");
    }
    csv.Write(head);

    // Output initial conditions.
    buildOutputVector(0, times[0].StartTime(), avg[0], err[0]);
    csv.Write(avg[0]);

    // Loop over all points, performing output.
    unsigned int step = 1;
    for (timevector::const_iterator iint=times.begin(); iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        real t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            buildOutputVector(step, t, avg[step], err[step]);
            csv.Write(avg[step]);
        }
    }

    // Close the CSV files.
    csv.Close();
}

// Writes particle process rates profile to a CSV file.
void Simulator::writePartProcCSV(const std::string &filename, 
                                 const Sweep::Mechanism &mech, 
                                 const timevector &times, 
                                 std::vector<fvector> &avg, 
                                 const std::vector<fvector> &err)
{
    // Open file for the CSV results.
    CSV_IO csv(filename, true);

    // Write the header row to the gas-phase chemistry CSV file.
    vector<string> head;
    head.push_back("Step");
    head.push_back("Time (s)");
    mech.GetProcessNames(head, 2);
    // Add units.
    for (unsigned int i=2; i!=head.size(); ++i) {
        head[i] = head[i] + " (1/cm3s)";
    }
    // Add error columns.
    for (unsigned int i=head.size(); i!=2; --i) {
        head.insert(head.begin()+i, "Err");
    }
    csv.Write(head);

    // Output initial conditions.
    buildOutputVector(0, times[0].StartTime(), avg[0], err[0]);
    csv.Write(avg[0]);

    // Loop over all points, performing output.
    unsigned int step = 1;
    for (timevector::const_iterator iint=times.begin(); iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        real t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            buildOutputVector(step, t, avg[step], err[step]);
            csv.Write(avg[step]);
        }
    }

    // Close the CSV files.
    csv.Close();
}

// Writes particle tracking for multiple particles to CSV files.
void Simulator::writePartTrackCSV(const std::string &filename,
                                  const Mechanism &mech,
                                  const timevector &times,
                                  std::vector<std::vector<fvector> > &track)
{
    // The track vector<vector<fvector>> should be arranged thus:
    // time-steps / particles / PSL variables.

    // Create an EnsembleStats object to get the PSL variable
    // names vector.
    Sweep::Stats::EnsembleStats stats(mech.ParticleMech());

    // Build the header row vector
    vector<string> head;
    head.push_back("Step");
    head.push_back("Time (s)");
    stats.PSL_Names(head, 2);

    // Determine max. number of particles tracked.
    unsigned int np = 0;
    for (unsigned int i=0; i!=track.size(); ++i) {
        np = max(np, (unsigned int)track[i].size());
    }

    // Open sufficient CSV files for all tracked particles.  Also
    // write the header row and the initial conditions.
    vector<CSV_IO*> csv(track[0].size());
    for (unsigned int i=0; i!=np; ++i) {
        // Open the CSV file for particle i.
        csv[i] = new CSV_IO(filename+"-p"+cstr(i)+".csv", true);
        // Write header row.
        csv[i]->Write(head);
        // Output particle initial conditions.        
        track[0][i].insert(track[0][i].begin(), times[0].StartTime());
        track[0][i].insert(track[0][i].begin(), 0.0);
        csv[i]->Write(track[0][i]);
    }

    // Loop over all points, performing output.
    unsigned int step = 1;
    for (timevector::const_iterator iint=times.begin(); iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        real t = iint->StartTime();
        for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            for (unsigned int i=0; i!=track[step].size(); ++i) {
                track[step][i].insert(track[step][i].begin(), t);
                track[step][i].insert(track[step][i].begin(), (real)step);
                csv[i]->Write(track[step][i]);
            }
        }
    }

    // Close the CSV files.
    for (unsigned int i=0; i!=track[0].size(); ++i) {
        csv[i]->Close();
        delete csv[i];
    }
}


// SAVE POINTS.

// Creates a simulation save point.  The save points can be
// used to restart an interrupted simulation, but are primarily
// used as output points for the particle size distributions.
void Simulator::createSavePoint(const Reactor &r, unsigned int step, 
                                unsigned int run) const
{
    // Build the save point file name.
    string fname = m_output_filename + "(" + cstr(run) + ")-SP(" + 
                   cstr(step) + ").sav";

    // Open the save point file.
    ofstream fout;
    fout.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    if (fout.good()) {
        fout.write((char*)&step, sizeof(step));
        fout.write((char*)&run, sizeof(run));
        ReactorFactory::Write(r, fout);
        fout.close();
    } else {
        // Throw error if the output file failed to open.
        throw runtime_error("Failed to open file for save point "
                            "output (Mops, Simulator::createSavePoint).");
    }
}

// Reads a save point file.
Reactor *const Simulator::readSavePoint(unsigned int step, 
                                        unsigned int run, 
                                        const Mechanism &mech) const
{
    Reactor *r = NULL;

    // Build the save posint file name.
    string fname = m_output_filename + "(" + cstr(run) + ")-SP(" + 
                   cstr(step) + ").sav";

    // Open the save point file.
    ifstream fin;
    fin.open(fname.c_str(), ios_base::in | ios_base::binary);

    if (m_file.good()) {
        // Read the step and run number (for file validation).
        unsigned int fstep=0, frun=0;
        fin.read(reinterpret_cast<char*>(&fstep), sizeof(fstep));
        fin.read(reinterpret_cast<char*>(&frun), sizeof(frun));

        // Check that the input file is valid.
        if (fstep != step) {
            // The file step number does not match!
            throw runtime_error("File step number does not match file name "
                                "(Mops, Simulator::readSavePoint).");
        }
        if (frun != run) {
            // The file run number does not match!
            throw runtime_error("File run number does not match file name "
                                "(Mops, Simulator::readSavePoint).");
        }

        // Deserialize the reactor from the file.
        r = ReactorFactory::Read(fin, mech);

        // Close the input file.
        fin.close();

        return r;
    } else {
        // Throw error if the output file failed to open.
        throw runtime_error("Failed to open file for save point "
                            "input (Mops, Simulator::readSavePoint).");
    }
    return NULL;
}

// Processes the PSLs at each save point into single files.
void Simulator::postProcessPSLs(const Mechanism &mech, 
                                const timevector &times) const
{
    Reactor *r = NULL;
    unsigned int step = 0;
    fvector psl;
    vector<fvector> ppsl;

    // Get reference to the particle mechanism.
    const Sweep::Mechanism &pmech = mech.ParticleMech();

    // Create an ensemble stats object.
    Sweep::Stats::EnsembleStats stats(pmech);

    // Build header row for PSL CSV output files.
    vector<string> header;
    stats.PSL_Names(header);

    // Build header row for primary-PSL CSV output files.
    vector<string> priheader;
    if (stats.GeneratesPPSL()) {
        stats.PPSL_Names(priheader);
    }

    // Open output files for all PSL save points.  Remember to
    // write the header row as well.
    vector<CSV_IO*> out(times.size(), NULL);
    for (unsigned int i=0; i!=times.size(); ++i) {
        real t = times[i].EndTime();
        out[i] = new CSV_IO();
        out[i]->Open(m_output_filename + "-psl(" +
                    cstr(t) + "s).csv", true);
        out[i]->Write(header);
    }

    // Open output files for all primary-PSL save points.  Remember to
    // write the header row as well.
    vector<CSV_IO*> priout(times.size(), NULL);
    if (stats.GeneratesPPSL()) {
        for (unsigned int i=0; i!=times.size(); ++i) {
            real t = times[i].EndTime();
            priout[i] = new CSV_IO();
            priout[i]->Open(m_output_filename + "-pri-psl(" +
                            cstr(t) + "s).csv", true);
            priout[i]->Write(priheader);
        }
    }

    // Loop over all time intervals.
    for (unsigned int i=0; i!=times.size(); ++i) {
        // Calculate the total step count after this interval.
        step += times[i].StepCount();

        // Loop over all runs.
        for (unsigned int irun=0; irun!=m_nruns; ++irun) {
            // Read the save point for this step and run.
            r = readSavePoint(step, irun, mech);

            if (r != NULL) {
                real scale = (real)m_nruns;
                if (m_output_every_iter) scale *= (real)m_niter;

                // Get PSL for all particles.
                for (unsigned int j=0; j!=r->Mixture()->ParticleCount(); ++j) {
                    // Get PSL.
                    stats.PSL(r->Mixture()->Particles(), j, times[i].EndTime(), 
                              psl, 1.0/(r->Mixture()->SampleVolume()*scale));
                    // Output particle PSL to CSV file.
                    out[i]->Write(psl);

                    // Get primary-PSL.
                    if (stats.GeneratesPPSL()) {
                        stats.PPSL(r->Mixture()->Particles(), j, times[i].EndTime(), 
                                   ppsl, 1.0/(r->Mixture()->SampleVolume()*scale));
                        // Output primary-PSL to CSV file.
                        for (unsigned int k=0; k!=ppsl.size(); ++k) {
                            priout[i]->Write(ppsl[k]);
                        }
                    }
                }

                // Draw particle images for tracked particles.
                unsigned int n = min(m_ptrack_count,r->Mixture()->ParticleCount());
                for (unsigned int j=0; j!=n; ++j) {
                    Sweep::Particle *sp = r->Mixture()->Particles().At(j);
                    if (sp != NULL) {
                        real t = times[i].EndTime();
                        string fname = m_output_filename + "-tem(" + cstr(t) + 
                                       "s, " + cstr(j) + ").pov";
                        Sweep::Imaging::ParticleImage img;
                        img.Construct(*sp);
//                        img.ConstructRandom(1.0, 5.0, 10001);
                        ofstream file; file.open(fname.c_str());
                        img.WritePOVRAY(file);
                        file.close();
                    }
                }

                delete r;
            } else {
                // Throw error if the reactor was not read.
                throw runtime_error("Unable to read reactor from save point "
                                    "(Mops, ParticleSolver::postProcessPSLs).");
            }
        }
    }

    // Close output CSV files.
    for (unsigned int i=0; i!=times.size(); ++i) {
        out[i]->Close();
        delete out[i];
    }
}

// Writes element fluxes to FluxViewer format.
void Simulator::writeElementFluxOutput(const std::string &filename, 
                               const Mechanism &mech, 
                               const timevector &times, 
                               const std::vector<fvector> &agpfwdrates,
                               const std::vector<fvector> &agprevrates,
                               const std::vector<fvector> &achem)
{
    Mops::fvector atemperatures;
    for (unsigned int i = 0; i < achem.size(); i++) {
        atemperatures.push_back(achem.at(i).at(achem.at(i).size() - 3));
    }
    FluxAnalyser fa(mech, times, agpfwdrates, agprevrates, atemperatures);
    for (unsigned int i = 0; i < mech.ElementCount(); i++) {
        for (unsigned int j = 0; j < m_flux_elements.size(); j++) {
            if (mech.Elements(i)->Name().compare(Strings::convertToCaps(m_flux_elements.at(j))) == 0) {
                fa.addElement(*mech.Elements(i));
            }
        }
    }
    fa.writeFluxes(filename, true);
}

// Add element for flux analysis postprocessor
void Simulator::AddFluxElement(const std::string &elem_str) {
    m_flux_elements.push_back(elem_str);
}

// Clear element list for flux analysis postprocessor
void Simulator::ClearFluxElements() {
    m_flux_elements.clear();
}

// COMPUTATION TIME CALCULATION.

// Calculates the time duration from a time mark to the
// current time.
double Simulator::calcDeltaCT(double markt) const
{
    return (double)(clock() - markt) / (double)CLOCKS_PER_SEC;
}


// READ/WRITE/COPY FUNCTIONS.

// Writes the simulator to a binary data stream.
void Simulator::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;
    
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write number of runs.
        unsigned int n = (unsigned int)m_nruns;
        out.write((char*)&n, sizeof(n));

        // Number of internal solver iterations to perform.
        n = (unsigned int)m_niter;
        out.write((char*)&n, sizeof(n));

        // Max. number of stochastic particles in sweep.
        n = (unsigned int)m_pcount;
        out.write((char*)&n, sizeof(n));

        // Max. M0 value, for initial scaling of ensemble.
        double val = (double)m_maxm0;
        out.write((char*)&val, sizeof(val));

        // Computation time.
        int i = (int)m_cpu_start;
        out.write((char*)&i, sizeof(i));
        i = (int)m_cpu_mark;
        out.write((char*)&i, sizeof(i));
        val = (double)m_runtime;
        out.write((char*)&val, sizeof(val));

        // Interval of console output data (in terms of time steps).
        n = (unsigned int)m_console_interval;
        out.write((char*)&n, sizeof(n));
        
        // Console column variable names.
        n = (unsigned int)m_console_vars.size();
        out.write((char*)&n, sizeof(n));
        for (unsigned int i=0; i!=n; ++i) {
            // Write the length of the name to the stream.
            unsigned int m = (unsigned int)m_console_vars[i].length();
            out.write((char*)&m, sizeof(m));
            // Write the name to the stream.
            out.write(m_console_vars[i].c_str(), m);
        }

        // Console variable mask.
        n = (unsigned int)m_console_mask.size();
        out.write((char*)&n, sizeof(n));
        for (unsigned int i=0; i!=n; ++i) {
            // Write the mask value to the stream.
            unsigned int m = (unsigned int)m_console_mask[i];
            out.write((char*)&m, sizeof(m));
        }

        // Console messages.
        if (m_console_msgs) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Name of output file.
        n = (unsigned int)m_output_filename.length();
        out.write((char*)&n, sizeof(n));
        out.write(m_output_filename.c_str(), n);


        // Flag controlling iteration output.
        if (m_output_every_iter) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Sub-step number at which output will occur.
        n = (unsigned int)m_output_step;
        out.write((char*)&n, sizeof(n));

        // Sub-step iteration number at which output will occur.
        n = (unsigned int)m_output_iter;
        out.write((char*)&n, sizeof(n));

        // Number of particles for which to produce tracking output.
        n = (unsigned int)m_ptrack_count;
        out.write((char*)&n, sizeof(n));
    } else {
        throw invalid_argument("Output stream not ready (Mops, Simulator::Serialize).");
    }
}

// Reads the Simulator data from a binary data stream.
void Simulator::Deserialize(std::istream &in)
{
    const unsigned int trueval  = 1;

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;
        unsigned int n = 0, m=0;
        int i = 0;
        char *name = NULL;

        switch (version) {
            case 0:
                // Read number of runs.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_nruns = n;

                // Number of internal solver iterations to perform.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_niter = n;

                // Max. number of stochastic particles in sweep.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_pcount = n;

                // Max. M0 value, for initial scaling of ensemble.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_maxm0 = (real)val;

                // Computation time.
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_cpu_start = (clock_t)i;
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_cpu_mark = (clock_t)i;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_runtime = (real)val;

                // Interval of console output data (in terms of time steps).
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_console_interval = n;
                
                // Console column variable names.
                m_console_vars.clear();
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                for (unsigned int k=0; k!=n; ++k) {
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));
                    name = new char[m];
                    in.read(name, m);
                    m_console_vars.push_back(string(name, m));
                    delete [] name;
                }

                // Console variable mask.
                m_console_mask.clear();
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                for (unsigned int k=0; k!=n; ++k) {
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));
                    m_console_mask.push_back(m);
                }

                // Console messages.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_console_msgs = (n==trueval);

                // Name of output file.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                name = new char[n];
                in.read(name, n);
                m_output_filename = string(name, n);
                delete [] name;

                // Flag controlling iteration output.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_output_every_iter = (n==trueval);

                // Sub-step number at which output will occur.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_output_step = n;

                // Sub-step iteration number at which output will occur.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_output_iter = n;

                // Number of particles for which to produce tracking output.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_ptrack_count = n;

                break;
            default:
                throw runtime_error("Simulator serialized version number "
                                    "is invalid (Mops, Simulator::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Mops, Simulator::Deserialize).");
    }
}
