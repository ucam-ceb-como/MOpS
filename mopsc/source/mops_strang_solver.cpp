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

// Solves the given reactor for the given time intervals.
void StrangSolver::SolveReactor(Mops::Reactor &r, 
                                const timevector &times, 
                                unsigned int nruns)
{
    unsigned int icon;
    real t1;     // Current time.
    real dt, t2; // Stop time for each step.

    // Store the initial conditions.
    Mixture initmix(r.Mech()->Species());
    initmix.SetFracs(r.Mixture()->MoleFractions());
    initmix.SetTemperature(r.Mixture()->Temperature());
    initmix.SetDensity(r.Mixture()->Density());

    // Initialise the reactor with the start time.
    t1 = times[0].StartTime();
    t2 = t1;
    r.Initialise(t1);
    r.Mixture()->Particles().Initialise(m_pcount);
    r.Mixture()->SetMaxM0(m_maxm0);

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

        // Set the error tolerances in the reactor.
        r.SetATOL(m_atol);
        r.SetRTOL(m_rtol);

        // Initialise the reactor with the start time.
        t2 = times[0].StartTime();
        r.Mixture()->SetFracs(initmix.MoleFractions());
        r.Mixture()->SetTemperature(initmix.Temperature());
        r.Mixture()->SetDensity(initmix.Density());
        r.Mixture()->Reset(m_maxm0);
        r.SetTime(t2);
        r.ResetSolver();

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
                m_cpu_mark = clock();
                    multiStrangStep(iint->SplitStepSize(), iint->SplittingStepCount(), r);
                m_chemtime += (double)(clock() - m_cpu_mark) / CLOCKS_PER_SEC;

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


// SOLUTION ROUTINES.

void StrangSolver::multiStrangStep(Mops::real dt, unsigned int n, Mops::Reactor &r)
{
    // Time counters.
    real t2 = r.Time();
    real h  = dt * 0.5; // Half step size.

    // Sweep time counters.
    real ts1 = r.Time();
    real ts2 = ts1;

    // Variables required to ensure particle number density is correctly
    // scaled with gas-phase expansion.
    real rho = 0.0, m0 = 0.0;

    m_cpu_mark = clock();
        // Solve first half-step of gas-phase chemistry.
        rho = r.Mixture()->Density();
        r.Solve(t2+=h);
    m_chemtime += calcDeltaCT(m_cpu_mark);

    m_cpu_mark = clock();
        // Solve one whole step of population balance (Sweep).
        m0 = r.Mixture()->ParticleCount()/r.Mixture()->SampleVolume();
        r.Mixture()->SetM0(r.Mixture()->Density() * m0 / rho);
        Run(ts1, ts2+=dt, *r.Mixture(), r.Mech()->ParticleMech());
    m_swp_ctime += calcDeltaCT(m_cpu_mark);
    
    for (unsigned int i=1; i!=n; ++i) {
        m_cpu_mark = clock();
            // Solve whole step of gas-phase chemistry.
            rho = r.Mixture()->Density();
            r.ResetSolver();
            r.Solve(t2+=dt);
        m_chemtime += calcDeltaCT(m_cpu_mark);

        m_cpu_mark = clock();
            // Solve whole step of population balance (Sweep).
            m0 = r.Mixture()->ParticleCount()/r.Mixture()->SampleVolume();
            r.Mixture()->SetM0(r.Mixture()->Density() * m0 / rho);
            Run(ts1, ts2+=dt, *r.Mixture(), r.Mech()->ParticleMech());
        m_swp_ctime += calcDeltaCT(m_cpu_mark);
        
    }

    m_cpu_mark = clock();
        // Solve last half-step of gas-phase chemistry.    
        r.ResetSolver();
        r.Solve(t2+=h);
    m_chemtime += calcDeltaCT(m_cpu_mark);
}


// POST-PROCESSING.

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
    Sweep::EnsembleStats stats(pmech);

    // Declare CPU time outputs (averages and errors).
    vector<vector<double> > acpu(npoints), ecpu(npoints);

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

        // Loop over all time intervals.
        unsigned int step = 1;
        for (Mops::timevector::const_iterator iint=times.begin(); 
             iint!=times.end(); ++iint) {
            // Loop over all time steps in this interval.
            for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
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
