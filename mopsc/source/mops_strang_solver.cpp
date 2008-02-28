#include "mops_strang_solver.h"
#include "mops_reactor_factory.h"
#include "sweep.h"
#include "string_functions.h"
#include "csv_io.h"

using namespace Mops;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
StrangSolver::StrangSolver(void)
: m_stats(NULL), m_swptime(0.0)
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
    Mixture initmix(r.Mechanism()->Species());
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
    beginFileOutput(*r.Mechanism(), times);

    // Set up the console output.
    icon = m_console_interval;
    setupConsole(*r.Mechanism());

    for (unsigned int irun=0; irun!=nruns; ++irun) {
        // Start the CPU timing clock.
        m_cpu_start = clock();
        m_chemtime  = 0.0;
        m_swptime   = 0.0;

        // Initialise the reactor with the start time.
        t2 = times[0].StartTime();
        r.Mixture()->SetFracs(initmix.MoleFractions());
        r.Mixture()->SetTemperature(initmix.Temperature());
        r.Mixture()->SetDensity(initmix.Density());
        r.Mixture()->Reset(m_maxm0);
        r.SetTime(t2);
        r.ResetSolver();

        // Set the error tolerances in the reactor.
        r.SetATOL(m_atol);
        r.SetRTOL(m_rtol);

        // Begin file output for this run.
        beginRunFileOutput(irun);
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
        endFileOutput();
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
    m_chemtime += (double)(clock() - m_cpu_mark) / CLOCKS_PER_SEC;

    m_cpu_mark = clock();
        // Solve one whole step of population balance (Sweep).
        m0 = r.Mixture()->ParticleCount()/r.Mixture()->SampleVolume();
        r.Mixture()->SetM0(r.Mixture()->Density() * m0 / rho);
        Run(ts1, ts2+=dt, *r.Mixture(), r.Mechanism()->ParticleMech());
    m_swptime += (double)(clock() - m_cpu_mark) / CLOCKS_PER_SEC;
    
    for (unsigned int i=1; i!=n; ++i) {
        m_cpu_mark = clock();
            // Solve whole step of gas-phase chemistry.
            rho = r.Mixture()->Density();
            r.ResetSolver();
            r.Solve(t2+=dt);
        m_chemtime += (double)(clock() - m_cpu_mark) / CLOCKS_PER_SEC;

        m_cpu_mark = clock();
            // Solve whole step of population balance (Sweep).
            m0 = r.Mixture()->ParticleCount()/r.Mixture()->SampleVolume();
            r.Mixture()->SetM0(r.Mixture()->Density() * m0 / rho);
            Run(ts1, ts2+=dt, *r.Mixture(), r.Mechanism()->ParticleMech());
        m_swptime += (double)(clock() - m_cpu_mark) / CLOCKS_PER_SEC;
        
    }

    m_cpu_mark = clock();
        // Solve last half-step of gas-phase chemistry.    
        r.ResetSolver();
        r.Solve(t2+=h);
    m_chemtime += (double)(clock() - m_cpu_mark) / CLOCKS_PER_SEC;
}


// POST-PROCESSING.

void StrangSolver::PostProcess(const std::string &filename, unsigned int nruns) const
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
    vector<vector<double>> acpu(npoints), ecpu(npoints);

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
        readChemDataPoint(fin, mech, achem[0], echem[0], nruns>1);
        readStatDataPoint(fin, pmech, astat[0], estat[0], nruns>1);
        readCPUDataPoint(fin, acpu[0], ecpu[0], nruns>1);

        // Loop over all time intervals.
        unsigned int step = 1;
        for (Mops::timevector::const_iterator iint=times.begin(); 
             iint!=times.end(); ++iint) {
            // Loop over all time steps in this interval.
            for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
                readChemDataPoint(fin, mech, achem[step], echem[step], nruns>1);
                readStatDataPoint(fin, pmech, astat[step], estat[step], nruns>1);
                readCPUDataPoint(fin, acpu[step], ecpu[step], nruns>1);
            }
        }

        // Close the input file.
        fin.close();
    }

    // CALCULATE AVERAGES AND CONFIDENCE INTERVALS.
    
    const real CONFA = 3.29;
    real invruns = 1.0 / (real)nruns;

    for (unsigned int step=0; step!=npoints; ++step) {
        // Gas-phase conditions.
        for (unsigned int i=0; i!=achem[0].size(); ++i) {
            // Calculate average over all runs.
            achem[step][i] *= invruns;

            if (nruns > 1) {
                // Calculate standard deviation over all runs.
                (echem[step][i] *= invruns) -= (achem[step][i] * achem[step][i]);
                // Calculate confidence interval.
                echem[step][i] = CONFA * sqrt(abs(echem[step][i] * invruns));
            }
        }

        // Particle statistics.
        for (unsigned int i=0; i!=astat[0].size(); ++i) {
            // Calculate average over all runs.
            astat[step][i] *= invruns;

            if (nruns > 1) {
                // Calculate standard deviation over all runs.
                (estat[step][i] *= invruns) -= (astat[step][i] * astat[step][i]);
                // Calculate confidence interval.
                estat[step][i] = CONFA * sqrt(abs(estat[step][i] * invruns));
            }
        }

        // CPU time.
        for (unsigned int i=0; i!=acpu[0].size(); ++i) {
            // Calculate average over all runs.
            acpu[step][i] *= invruns;

            if (nruns > 1) {
                // Calculate standard deviation over all runs.
                (ecpu[step][i] *= invruns) -= (acpu[step][i] * acpu[step][i]);
                // Calculate confidence interval.
                ecpu[step][i] = CONFA * sqrt(abs(ecpu[step][i] * invruns));
            }
        }
    }

    // CREATE CSV FILES.

    // Now open files for the CSV results.
    CSV_IO csvc(string(filename).append("-chem.csv"), true);
    CSV_IO csvs(string(filename).append("-part.csv"), true);
    CSV_IO csvt(string(filename).append("-cput.csv"), true);

    // Write the header row to the gas-phase chemistry CSV file.
    vector<string> headc;
    headc.push_back("Step");
    headc.push_back("Time (s)");
    for (unsigned int isp=0; isp<mech.SpeciesCount(); isp++) {
        headc.push_back(string(mech.Species(isp)->Name()).append(" (mol/cm3)"));
    }
    headc.push_back("T (K)");
    headc.push_back("Density (mol/cm3)");
    headc.push_back("Pressure (Pa)");
    for (unsigned int i=headc.size(); i!=2; --i) {
        headc.insert(headc.begin()+i, "Err");
    }
    csvc.Write(headc);

    // Write header for particle CSV file.
    vector<string> heads;
    heads.push_back("Step");
    heads.push_back("Time (s)");
    stats.Names(heads, 2);
    for (unsigned int i=heads.size(); i!=2; --i) {
        heads.insert(heads.begin()+i, "Err");
    }
    csvs.Write(heads);

    // Write header for CPU timing CSV file.
    vector<string> headt;
    headt.push_back("Step");
    headt.push_back("Time (s)");
    headt.push_back("CPU Time (s)");
    headt.push_back("Err");
    headt.push_back("Chem CPU Time (s)");
    headt.push_back("Err");
    headt.push_back("Sweep CPU Time (s)");
    headt.push_back("Err");
    csvt.Write(headt);

    // OUTPUT TO CSV FILES.

    // Output initial conditions.
    {
        // Gas-phase conditions.
        buildOutputVector(0, times[0].StartTime(), achem[0], echem[0]);
        csvc.Write(achem[0]);

        // Stats.
        buildOutputVector(0, times[0].StartTime(), astat[0], estat[0]);
        csvs.Write(astat[0]);

        // CPU times.
        buildOutputVector(0, times[0].StartTime(), acpu[0], ecpu[0]);
        csvt.Write(acpu[0]);
    }

    // Loop over all points, performing output.
    unsigned int step = 1;
    for (timevector::const_iterator iint=times.begin(); iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        real t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();

            // Gas-phase conditions.
            buildOutputVector(step, t, achem[step], echem[step]);
            csvc.Write(achem[step]);

            // Stats.
            buildOutputVector(step, t, astat[step], estat[step]);
            csvs.Write(astat[step]);

            // CPU times.
            buildOutputVector(step, t, acpu[step], ecpu[step]);
            csvt.Write(acpu[step]);
        }
    }

    // Close the CSV files.
    csvc.Close();
    csvs.Close();
    csvt.Close();

    // Now post-process the PSLs.
    postProcessPSLs(nruns, mech, times);
}


// POST-PROCESSING ROUTINES.

// Reads a chemistry data point.
void StrangSolver::readChemDataPoint(std::istream &in, const Mops::Mechanism &mech,
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
        in.read(reinterpret_cast<char*>(&P), sizeof(P));

        // Resize vectors.
        sum.resize(N+3, 0.0);
        if (calcsqrs) sumsqr.resize(N+3, 0.0);

        // Calculate sums and sums of square (for average and
        // error calculation).
        for (unsigned int i=0; i!=N; ++i) {
            sum[i] += (D * y[i] * 1.0e-6);
            if (calcsqrs) sumsqr[i] += (D * D * y[i] * y[i] * 1.0e-12);
        }

        // Calculates sums and sums of squares of temperature, density
        // and pressure.
        sum[N]   += T;
        sum[N+1] += (D * 1.0e-6);
        sum[N+2] += P;
        if (calcsqrs) {
            sumsqr[N]   += (T*T);
            sumsqr[N+1] += (D*D*1.0e-12);
            sumsqr[N+2] += (P*P);
        }
    }
}

// Reads a particle stats data point.
void StrangSolver::readStatDataPoint(std::istream &in, const Sweep::Mechanism &mech, 
                                     fvector &sum, fvector &sumsqr, bool calcsqrs)
{
    // Check for valid stream.
    if (in.good()) {
        // Read the stats.
        Sweep::EnsembleStats stats(mech);
        stats.Deserialize(in);

        // Get the stats vector.
        fvector s;
        stats.Get(s);

        // Resize vectors.
        sum.resize(s.size(), 0.0);
        if (calcsqrs) sumsqr.resize(s.size(), 0.0);

        // Calculate sums and sums of square (for average and
        // error calculation).
        for (unsigned int i=0; i!=s.size(); ++i) {
            sum[i] += s[i];
            if (calcsqrs) sumsqr[i] += (s[i] * s[i]);
        }
    }
}

// Reads a CPU timing data point.
void StrangSolver::readCPUDataPoint(std::istream &in, fvector &sum, 
                                    fvector &sumsqr, bool calcsqrs)
{
    const unsigned int N = 3;

    // Check for valid stream.
    if (in.good()) {

        // Read the computation times.
        real cpu[N];
        in.read(reinterpret_cast<char*>(&cpu[0]), sizeof(cpu[0])*N);

        // Resize output vectors.
        sum.resize(N, 0.0);
        if (calcsqrs) sumsqr.resize(N, 0.0);

        // Calculate sums and sums of square (for average and
        // error calculation).
        for (unsigned int i=0; i!=N; ++i) {
            sum[i] += cpu[i];
            if (calcsqrs) sumsqr[i] += (cpu[i] * cpu[i]);
        }
    }
}

// Builds output vector of particle stats.
void StrangSolver::buildOutputVector(unsigned int step, real time, 
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


// CONSOLE OUTPUT.

// Sets up console output using the given mechanism as a template.
void StrangSolver::setupConsole(const Mechanism &mech)
{
    vector<string> header;
    m_console_mask.clear();

    // Loop over the column variable names.
    vector<string>::const_iterator i;
    for (i=m_console_vars.begin(); i!=m_console_vars.end(); i++) {
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

// Writes current reactor state to the console.
void StrangSolver::consoleOutput(const Reactor &r) const
{
    // Get output data from gas-phase.
    vector<real> out;
    r.Mixture()->GetConcs(out);
    out.push_back(r.Mixture()->Temperature());
    out.push_back(r.Mixture()->Density());
    out.push_back(r.Time());

    // Get output data from particles.
    Sweep::EnsembleStats stats(r.Mechanism()->ParticleMech());
    r.Mixture()->GetVitalStats(stats);
    out.push_back(stats.BasicStats().PCount());
    out.push_back(stats.BasicStats().M0());
    out.push_back(stats.BasicStats().Fv());

    // Get output CPU time.
    double cputime = (double)(clock() - m_cpu_start) / (double)CLOCKS_PER_SEC;
    out.push_back(cputime);

    // Print data to console.
    m_console.PrintRow(out, m_console_mask);
}

// FILE OUTPUT.

// Sets up the file output by outputting an auxilliary file
// which stores all the information required to post-process the
// simulation and by opening the output file.
void StrangSolver::beginFileOutput(const Mops::Mechanism &mech, 
                                   const Mops::timevector &times)
{
    // Build output file name.
    string fname(m_output_filename); fname.append(".dat");

    // Output a file for post-processing information (mechanism and
    // time intervals).
    fstream fout;
    fout.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!fout.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "output (Mops, StrangSolver::beginFileOutput).");
    }

    // Write the mechanism to the output file.
    mech.Serialize(fout);

    // Write the time intervals to the output file.
    unsigned int n = times.size();
    fout.write((char*)&n, sizeof(n));
    for (Mops::timevector::const_iterator i=times.begin(); i!=times.end(); i++) {
        (*i).Serialize(fout);
    }

    // Close the post-processing info file.
    fout.close();

    // Set up stats output.
    delete m_stats;
    m_stats = new Sweep::EnsembleStats(mech.ParticleMech());
}

// Sets up file output for a new run given the run number.
void StrangSolver::beginRunFileOutput(unsigned int run)
{
    // Build the simulation output file name.
    string fname = m_output_filename + "(" + cstr(run) + ").dat";

    // Open the simulation output file.
    m_file.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!m_file.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "output (Mops, StrangSolver::beginRunFileOutput).");
    }
}

// Writes the current system state to the output file.
void StrangSolver::fileOutput(const Reactor &r)
{
    // Write gas-phase conditions to file.
    m_file.write((char*)&r.Mixture()->RawData()[0], 
                 sizeof(r.Mixture()->RawData()[0])*r.Mechanism()->SpeciesCount());
    real T = r.Mixture()->Temperature();
    m_file.write((char*)&T, sizeof(T));
    real D = r.Mixture()->Density();
    m_file.write((char*)&D, sizeof(D));
    real P = r.Mixture()->Pressure();
    m_file.write((char*)&P, sizeof(P));

    // Write particle stats to file.
    r.Mixture()->GetVitalStats(*m_stats);
    m_stats->Serialize(m_file);

    // Write CPU times to file.
    double cputime = (double)(clock() - m_cpu_start) / (double)CLOCKS_PER_SEC;
    m_file.write((char*)&cputime, sizeof(cputime));
    m_file.write((char*)&m_chemtime, sizeof(m_chemtime));
    m_file.write((char*)&m_swptime, sizeof(m_swptime));
}

// Ends file output by closing all open files.
void StrangSolver::endFileOutput()
{
    // Close the output file.
    m_file.close();
}


// SAVE POINTS AND PSL POST-PROCESSING.

// Creates a simulation save point.  The save points can be
// used to restart an interrupted simulation, but are primarily
// used as output points for the particle size distributions.
void StrangSolver::createSavePoint(const Reactor &r, unsigned int step, 
                                   unsigned int run) const
{
    // Build the save point file name.
    string fname = m_output_filename + "(" + cstr(run) + ")-SP(" + 
                   cstr(step) + ").mops";

    // Open the save point file.
    ofstream fout;
    fout.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    if (m_file.good()) {
        fout.write((char*)&step, sizeof(step));
        fout.write((char*)&run, sizeof(run));
        ReactorFactory::Write(r, fout);
        fout.close();
    } else {
        // Throw error if the output file failed to open.
        throw runtime_error("Failed to open file for save point "
                            "output (Mops, StrangSolver::createSavePoint).");
    }
}

// Reads a save point file.
Reactor *const StrangSolver::readSavePoint(unsigned int step, 
                                           unsigned int run, 
                                           const Mechanism &mech) const
{
    Reactor *r = NULL;

    // Build the save posint file name.
    string fname = m_output_filename + "(" + cstr(run) + ")-SP(" + 
                   cstr(step) + ").mops";

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
                                "(Mops, StrangSolver::readSavePoint).");
        }
        if (frun != run) {
            // The file run number does not match!
            throw runtime_error("File run number does not match file name "
                                "(Mops, StrangSolver::readSavePoint).");
        }

        // Deserialize the reactor from the file.
        r = ReactorFactory::Read(fin, mech);

        // Close the input file.
        fin.close();

        return r;
    } else {
        // Throw error if the output file failed to open.
        throw runtime_error("Failed to open file for save point "
                            "input (Mops, StrangSolver::readSavePoint).");
    }
}

// Processes the PSLs.
void StrangSolver::postProcessPSLs(unsigned int nruns, const Mechanism &mech, 
                                   const timevector &times) const
{
    Reactor *r = NULL;
    unsigned int step = 0;
    Sweep::EnsembleStats stats(mech.ParticleMech());
    fvector psl;

    // Build header row for CSV output files.
    vector<string> header;
    stats.PSL_Names(header);

    // Open output files for all PSL save points.  Remember to
    // write the header row as well.
    vector<CSV_IO*> out(times.size());
    for (unsigned int i=0; i!=times.size(); ++i) {
        real t = times[i].EndTime();
        out[i] = new CSV_IO();
        out[i]->Open(m_output_filename + "-psl(" +
                    cstr(t) + "s).csv", true);
        out[i]->Write(header);
    }

    // Loop over all time intervals.
    for (unsigned int i=0; i!=times.size(); ++i) {
        // Calculate the total step count after this interval.
        step += times[i].StepCount();

        // Loop over all runs.
        for (unsigned int irun=0; irun!=nruns; ++irun) {
            // Read the save point for this step and run.
            r = readSavePoint(step, irun, mech);

            if (r != NULL) {
                // Get PSL for all particles.
                for (unsigned int j=0; j!=r->Mixture()->ParticleCount(); ++j) {
                    // Get PSL.
                    stats.PSL(r->Mixture()->Particles(), j, times[i].EndTime(), 
                              psl, 1.0/(r->Mixture()->SampleVolume()*nruns));
                    // Output particle PSL to CSV file.
                    out[i]->Write(psl);
                }

                delete r;
            } else {
                // Throw error if the reactor was not read.
                throw runtime_error("Unable to read reactor from save point "
                                    "(Mops, StrangSolver::postProcessPSLs).");
            }
        }
    }

    // Close output CSV files.
    for (unsigned int i=0; i!=times.size(); ++i) {
        out[i]->Close();
        delete out[i];
    }
}
