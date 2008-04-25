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
: m_atol(1.0e-3), m_rtol(6.0e-4), m_nruns(1), m_pcount(0), m_maxm0(0.0), 
  m_cpu_start((clock_t)0.0), m_cpu_mark((clock_t)0.0), m_chemtime(0.0), 
  m_console_interval(1),
  m_console_msgs(true), m_output_filename("mops-out")
{
}

// Default destructor.
Solver::~Solver(void)
{
}


// ERROR TOLERANCES.

real Solver::ATOL() const
{
    return m_atol;
}

void Solver::SetATOL(real atol)
{
    m_atol = atol;
    m_ode.SetATOL(atol);
}

real Solver::RTOL() const
{
    return m_rtol;
}

void Solver::SetRTOL(real rtol)
{
    m_rtol = rtol;
    m_ode.SetRTOL(rtol);
}


// SWEEP SETTINGS.

// Returns the number of runs to perform
unsigned int Solver::RunCount(void) const {return m_nruns;}

// Sets the number of runs to peform.
void Solver::SetRunCount(unsigned int n) {m_nruns = n;}

// Returns the max. stochastic particle count.
unsigned int Solver::MaxPartCount(void) const {return m_pcount;}

// Sets the max. stochastic particle count.
void Solver::SetMaxPartCount(unsigned int n) {m_pcount = n;}

// Returns the max. M0, for initial ensemble scaling.
real Solver::MaxM0(void) const {return m_maxm0;}

// Sets max. M0.
void Solver::SetMaxM0(real m0) {m_maxm0 = m0;}


// CONSOLE INTERVAL.

unsigned int Solver::ConsoleInterval() const
{
    return m_console_interval;
}

void Solver::SetConsoleInterval(unsigned int cint)
{
    m_console_interval = cint;
}


// CONSOLE VARIABLE NAMES.

const std::vector<std::string> &Solver::ConsoleVariables() const
{
    return m_console_vars;
}

const std::string Solver::ConsoleVariable(unsigned int i) const
{
    if (i < m_console_vars.size()) {
        return m_console_vars[i];
    } else {
        // Returns the first variable name if index is invalid.     
        return "";
    }
}

void Solver::AddConsoleVariable(const std::string &var)
{
    m_console_vars.push_back(var);
}

void Solver::RemoveConsoleVariable(const std::string &var)
{
    vector<string>::iterator i;
    for (i=m_console_vars.begin(); i!=m_console_vars.end(); ++i) {
        if ((*i).compare(var) == 0) {
            m_console_vars.erase(i);
            return;
        }
    }
}

void Solver::RemoveConsoleVariable(unsigned int i)
{
    if (i < m_console_vars.size()) {
        m_console_vars.erase(m_console_vars.begin()+i);
    }
}


// CONSOLE MESSAGES.

bool Solver::UseConsoleMsgs() const
{
    return m_console_msgs;
}

void Solver::SetUseConsoleMsgs(bool msgs)
{
    m_console_msgs = msgs;
}


// OUTPUT FILE NAME.

const std::string &Solver::OutputFile() const
{
    return m_output_filename;
}

void Solver::SetOutputFile(const std::string &name)
{
    m_output_filename = name;
}


// SOLVING REACTORS.

// Solves the given reactor for the given time intervals.
void Solver::SolveReactor(Mops::Reactor &r, 
                          const timevector &times, unsigned int nruns)
{
    // Start the CPU timing clock.
    m_cpu_start = clock();
    m_chemtime  = 0.0;

    unsigned int icon;
    real dt, t2; // Stop time for each step.

    // Initialise the reactor with the start time.
    t2 = times[0].StartTime();
    r.SetTime(t2);

    // Set up the ODE solver.
    m_ode.Initialise(r);
    m_ode.SetATOL(m_atol);
    m_ode.SetRTOL(m_rtol);

    // Set up file output.
    writeAux(m_output_filename, *r.Mech(), times);
    openOutputFile(0);

    // Set up the console output.
    icon = m_console_interval;
    setupConsole(*r.Mech());

    // Output initial conditions.
    fileOutput(r);
    consoleOutput(r);

    // Loop over the time intervals.
    timevector::const_iterator iint;
    for (iint=times.begin(); iint!=times.end(); ++iint) {
        // Get the step size for this interval.
        dt = (*iint).StepSize();

        // Loop over the steps in this interval.
        unsigned int istep;
        for (istep=0; istep<(*iint).StepCount(); ++istep) {
            // Run the reactor solver for this step (timed).
            m_cpu_mark = clock();
                m_ode.Solve(r, (t2+=dt));
                r.SetTime(t2);
            m_chemtime += calcDeltaCT(m_cpu_mark);

            // Generate file output.
            fileOutput(r);

            // Generate console output.
            if (--icon == 0) {
                consoleOutput(r);
                icon = m_console_interval;
            }
        }
    }

    // Close the output files.
    closeOutputFile();
}


// POST-PROCESSING.

void Solver::PostProcess(const std::string &filename, 
                         unsigned int nruns) const
{
    // READ AUXILLIARY INFORMATION.

    // Read auxilliary information about the simulation (mechanism
    // and time intervals).
    Mops::Mechanism mech;
    Mops::timevector times;
    readAux(filename, mech, times);

    // SETUP OUTPUT DATA STRUCTURES.

    // Calculate number of output points.
    unsigned int npoints = 1; // 1 for initial conditions.
    for(Mops::timevector::const_iterator i=times.begin(); i!=times.end(); ++i) {
        npoints += i->StepCount();
    }

    // Declare chemistry outputs (averages and errors).
    vector<fvector> achem(npoints), echem(npoints);

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
                                "(Mops, Solver::PostProcess).");
        }

        // Read initial conditions and computation times.
        readGasPhaseDataPoint(fin, mech, achem[0], echem[0], nruns>1);
        readCTDataPoint(fin, 2, acpu[0], ecpu[0], nruns>1);

        // Loop over all time intervals.
        unsigned int step = 1;
        for (Mops::timevector::const_iterator iint=times.begin(); 
             iint!=times.end(); ++iint) {
            // Loop over all time steps in this interval.
            for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
                readGasPhaseDataPoint(fin, mech, achem[step], echem[step], nruns>1);
                readCTDataPoint(fin, 2, acpu[step], ecpu[step], nruns>1);
            }
        }

        // Close the input file.
        fin.close();
    }

    // CALCULATE AVERAGES AND CONFIDENCE INTERVALS.
    
    calcAvgConf(achem, echem, nruns);
    calcAvgConf(acpu, ecpu, nruns);

    // OUTPUT TO CSV FILES.

    writeGasPhaseCSV(filename+"-chem.csv", mech, times, achem, echem);
    writeCT_CSV(filename+"-cpu.csv", times, acpu, ecpu);
}


// FILE OUTPUT (protected).

// Opens an output file for the given run number.
void Solver::openOutputFile(unsigned int run) const
{
    // Build the simulation output file name.
    string fname = m_output_filename + "(" + cstr(run) + ").dat";

    // Open the simulation output file.
    m_file.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!m_file.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "output (Mops, Solver::openOutputFile).");
    }
}

// Closes the output file.
void Solver::closeOutputFile() const
{
    // Close the output file.
    m_file.close();
}

// Writes the gas-phase conditions of the given reactor to
// the binary output file.
void Solver::outputGasPhase(const Reactor &r) const
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


// FILE OUTPUT (private).

void Solver::fileOutput(const Mops::Reactor &r) const
{
    // Write the gas-phase conditions to the output file.
    outputGasPhase(r);

    // Write CPU times to file.
    double cputime = calcDeltaCT(m_cpu_start);
    m_file.write((char*)&cputime, sizeof(cputime));
    m_file.write((char*)&m_chemtime, sizeof(m_chemtime));
}


// CONSOLE OUTPUT (protected).

// Sets up console output.
void Solver::setupConsole(const Mops::Mechanism &mech)
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

// Writes output to the console.
void Solver::consoleOutput(const Mops::Reactor &r) const
{
    // Get output data from gas-phase.
    static vector<real> out;
    r.Mixture()->GetConcs(out);
    out.push_back(r.Mixture()->Temperature());
    out.push_back(r.Mixture()->Density());
    out.push_back(r.Time());

    // Get output data from particles.
    Sweep::EnsembleStats stats(r.Mech()->ParticleMech());
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
void Solver::writeAux(const std::string &filename, 
                      const Mops::Mechanism &mech, 
                      const Mops::timevector &times)
{
    // Build output file name.
    string fname = filename + ".dat";

    // Output a file for post-processing information (mechanism and
    // time intervals).
    fstream fout;
    fout.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!fout.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "output (Mops, Solver::writeAux).");
    }

    // Write the mechanism to the output file.
    mech.Serialize(fout);

    // Write the time intervals to the output file.
    unsigned int n = times.size();
    fout.write((char*)&n, sizeof(n));
    for (timevector::const_iterator i=times.begin(); i!=times.end(); i++) {
        (*i).Serialize(fout);
    }

    // Close the post-processing info file.
    fout.close();
}
    
// Reads auxilliary post-processing information using the
// given file name.  This information is the chemical mechanism
// and the output time intervals.
void Solver::readAux(const std::string &filename, 
                     Mops::Mechanism &mech, 
                     Mops::timevector &times)
{
    // Build input file name.
    string fname = filename + ".dat";

    // Output the file which contains the simulation mechanism and
    // time intervals.
    fstream fin;
    fin.open(fname.c_str(), ios_base::in | ios_base::binary);

    // Throw error if the file failed to open.
    if (!fin.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "post-processing (Mops, Solver::readAux).");
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

    // Close the simulation settings file.
    fin.close();
}

// Reads a gas-phase chemistry data point from the binary file.
// To allow the averages and confidence intervals to be calculated
// the data point is added to a vector of sums, and the squares are
// added to the vector sumsqr if necessary.
void Solver::readGasPhaseDataPoint(std::istream &in, const Mops::Mechanism &mech,
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
        sumsqr.resize(N+3, 0.0);

        // Calculate sums and sums of square (for average and
        // error calculation).
        for (unsigned int i=0; i!=N; ++i) {
            // Note conversion to cm-3 from m-3.
            sum[i] += (D * y[i] * 1.0e-6);
            if (calcsqrs) sumsqr[i] += (D * D * y[i] * y[i] * 1.0e-12);
        }

        // Calculates sums and sums of squares of temperature, density
        // and pressure.
        sum[N]   += T;
        sum[N+1] += (D * 1.0e-6); // Note conversion to cm-3 from m-3.
        sum[N+2] += P;
        if (calcsqrs) {
            sumsqr[N]   += (T*T);
            sumsqr[N+1] += (D*D*1.0e-12);
            sumsqr[N+2] += (P*P);
        }
    }
}

// Reads a CPU timing data from the binary file.
// To allow the averages and confidence intervals to be calculated
// the data point is added to a vector of sums, and the squares are
// added to the vector sumsqr if necessary.
void Solver::readCTDataPoint(std::istream &in, unsigned int N,
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

// Takes vectors of vectors of variable sums and sums of squares, which
// are converted into the average values and the confidence intervals.
void Solver::calcAvgConf(std::vector<fvector> &avg, 
                         std::vector<fvector> &err, 
                         unsigned int nruns)
{
    const real CONFA = 3.29; // for 99.9% confidence interval.

    // Pre-calc some useful values.
    real invruns = 1.0 / (real)nruns;
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
void Solver::buildOutputVector(unsigned int step, real time, 
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
void Solver::writeGasPhaseCSV(const std::string &filename, 
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
void Solver::writeCT_CSV(const std::string &filename, 
                         const timevector &times, 
                         std::vector<fvector> &avg,
                         const std::vector<fvector> &err) const
{
    // Open file for the CSV results.
    CSV_IO csv(filename, true);

    // Write the header row to the CPU time CSV file.
    vector<string> head;
    head.push_back("Step");
    head.push_back("Time (s)");
    head.push_back("CPU Time (s)");
    head.push_back("Err");
    head.push_back("Chem CPU Time (s)");
    head.push_back("Err");
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


// SAVE POINTS.

// Creates a simulation save point.  The save points can be
// used to restart an interrupted simulation, but are primarily
// used as output points for the particle size distributions.
void Solver::createSavePoint(const Reactor &r, unsigned int step, 
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
                            "output (Mops, Solver::createSavePoint).");
    }
}

// Reads a save point file.
Reactor *const Solver::readSavePoint(unsigned int step, 
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
                                "(Mops, Solver::readSavePoint).");
        }
        if (frun != run) {
            // The file run number does not match!
            throw runtime_error("File run number does not match file name "
                                "(Mops, Solver::readSavePoint).");
        }

        // Deserialize the reactor from the file.
        r = ReactorFactory::Read(fin, mech);

        // Close the input file.
        fin.close();

        return r;
    } else {
        // Throw error if the output file failed to open.
        throw runtime_error("Failed to open file for save point "
                            "input (Mops, Solver::readSavePoint).");
    }
    return NULL;
}


// COMPUTATION TIME CALCULATION.

// Calculates the time duration from a time mark to the
// current time.
double Solver::calcDeltaCT(double markt) const
{
    return (double)(clock() - markt) / (double)CLOCKS_PER_SEC;
}
