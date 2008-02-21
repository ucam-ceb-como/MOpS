#include "mops_solver.h"
#include "mops_reactor_factory.h"
#include "csv_io.h"
#include <vector>
#include <string>
#include <time.h>

using namespace Mops;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Solver::Solver(void)
{
    m_atol = 1.0e-3;
    m_rtol = 6.0e-4;
    m_console_interval = 1;
    m_console_vars.clear();
    m_console_msgs = true;
    m_output_filename = "mops-out";
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
}

real Solver::RTOL() const
{
    return m_rtol;
}

void Solver::SetRTOL(real rtol)
{
    m_rtol = rtol;
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
    vector<string>::const_iterator i;
    for (i=m_console_vars.begin(); i!=m_console_vars.end(); i++) {
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
    r.Initialise(t2);

    // Set the error tolerances in the reactor.
    r.SetATOL(m_atol);
    r.SetRTOL(m_rtol);

    // Set up file output.
    beginFileOutput(*r.Mechanism(), times);

    // Set up the console output.
    icon = m_console_interval;
    setupConsole(*r.Mechanism());

    // Output initial conditions.
    fileOutput(r);
    consoleOutput(r);

    // Loop over the time intervals.
    timevector::const_iterator iint;
    for (iint=times.begin(); iint!=times.end(); iint++) {
        // Get the step size for this interval.
        dt = (*iint).StepSize();

        // Loop over the steps in this interval.
        unsigned int istep;
        for (istep=0; istep<(*iint).StepCount(); istep++) {
            // Run the reactor solver for this step (timed).
            m_cpu_mark = clock();
                r.Solve((t2+=dt));
            m_chemtime += (double)(clock() - m_cpu_mark) / CLOCKS_PER_SEC;

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
    endFileOutput();
}


// POST-PROCESSING.

void Solver::PostProcess(const std::string &filename, unsigned int nruns) const
{
    Mechanism mech;
    timevector times;
    Reactor *r;
    fvector out;
    vector<double> cpu(4);

    // Build input file name.
    string fname(filename); fname.append(".dat");

    // Output the file which contains the simulation mechanism and
    // time intervals.
    fstream fin;
    fin.open(fname.c_str(), ios_base::in | ios_base::binary);

    // Throw error if the file failed to open.
    if (!fin.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "post-processing (Mops, Solver::PostProcess).");
    }

    // Read the mechanism from the file.
    mech.Deserialize(fin);

    // Read the time intervals from the file.
    unsigned int ntimes;
    fin.read(reinterpret_cast<char*>(&ntimes), sizeof(ntimes));
    for (unsigned int i=0; i<ntimes; i++) {
        times.push_back(TimeInterval(fin));
    }

    // Close the simulation settings file.
    fin.close();

    // Now we must read the reactor conditions and all time points:

    // Build the simulation input file name.
    fname = string(filename).append("(1).dat");

    // Open the simulation output file.
    fin.open(fname.c_str(), ios_base::in | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!fin.good()) {
        throw runtime_error("Failed to open file for post-processing "
                            "(Mops, Solver::PostProcess).");
    }

    // Now open a file for the CSV results.
    CSV_IO csv(string(filename).append("-chem.csv"), true);

    // Write the header row to the CSV file.
    vector<string> header;
    header.push_back("Step");
    header.push_back("Time (s)");
    for (unsigned int isp=0; isp<mech.SpeciesCount(); isp++) {
        header.push_back(string(mech.Species(isp)->Name()).append(" (mol/cm3)"));
    }
    header.push_back("T (K)");
    header.push_back("Density (mol/cm3)");
    header.push_back("Pressure (Pa)");
    csv.Write(header);

    // Open a CSV file for the computation time results.
    CSV_IO csvcpu(string(filename).append("-ct.csv"), true);
    vector<string> cpuheader;
    cpuheader.push_back("Step");
    cpuheader.push_back("Time (s)");
    cpuheader.push_back("CPU Time (s)");
    cpuheader.push_back("Chem CPU Time (s)");
    csvcpu.Write(cpuheader);

    // Read initial reactor conditions.
    {
        // REACTOR.

        // Read the reactor object (temporary).
        r = ReactorFactory::Read(fin, mech);

        // Get the data for output.
        buildOutputVector(*r, out);

        // Write output data to the CSV file.
        csv.Write(out);

        // Delete temporary reactor object.
        delete r;

        // COMPUTATION TIMES.

        // Read the computation times.
        cpu[0] = out[0]; cpu[1] = out[1];
        fin.read(reinterpret_cast<char*>(&cpu[2]), sizeof(cpu[2]));
        fin.read(reinterpret_cast<char*>(&cpu[3]), sizeof(cpu[3]));

        // Write CPU times to file.
        csvcpu.Write(cpu);
    }

    // Loop over all time intervals.
    for (timevector::const_iterator iint=times.begin(); iint!=times.end(); iint++) {
        // Loop over all time steps in this interval.
        unsigned int istep;
        for (istep=0; istep<(*iint).StepCount(); istep++) {
            // REACTOR.

            // Read the reactor object (temporary).
            r = ReactorFactory::Read(fin, mech);

            // Get the data for output.
            buildOutputVector(*r, out);
            out[0] = (real)(istep + 1); // Step number.

            // Write output data to the CSV file.
            csv.Write(out);

            // Delete temporary reactor object.
            delete r;

            // COMPUTATION TIMES.

            // Read the computation times.
            cpu[0] = out[0]; cpu[1] = out[1];
            fin.read(reinterpret_cast<char*>(&cpu[2]), sizeof(cpu[2]));
            fin.read(reinterpret_cast<char*>(&cpu[3]), sizeof(cpu[3]));
            
            // Write CPU times to file.
            csvcpu.Write(cpu);
        }
    }

    // Close the CSV file and the input file.
    fin.close();
    csv.Close();
}


// FILE OUTPUT (protected).

// Sets up the file output by outputting an auxilliary file
// which stores all the information required to post-process the
// simulation and by opening the output file.
void Solver::beginFileOutput(const Mops::Mechanism &mech, 
                             const timevector &times)
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
                            "output (Mops, Solver::beginFileOutput).");
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

    // Build the simulation output file name.
    fname = string(m_output_filename).append("(1).dat");

    // Open the simulation output file.
    m_file.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!m_file.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "output (Mops, Solver::beginFileOutput).");
    }
}

void Solver::fileOutput(const Mops::Reactor &r)
{
    // Write the reactor to the output file.
    ReactorFactory::Write(r, m_file);

    // Write CPU times to file.
    double cputime = (double)(clock() - m_cpu_start) / (double)CLOCKS_PER_SEC;
    m_file.write((char*)&cputime, sizeof(cputime));
    m_file.write((char*)&m_chemtime, sizeof(m_chemtime));
}

void Solver::endFileOutput()
{
    // Close the output file.
    m_file.close();
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
            // This variable is the mixture temperature.
            header.push_back("Time (s)");
            m_console_mask.push_back(mech.Species().size()+2);
        } else {
            // Check for a species name.
            int isp = mech.FindSpecies((*i));
            if (isp >= 0) {
                // This is a valid species name.
                header.push_back((*i));
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
    // Get output data.
    vector<real> out;
    r.Mixture()->GetConcs(out);
    out.push_back(r.Mixture()->Temperature());
    out.push_back(r.Mixture()->Density());
    out.push_back(r.Time());

    // Print data to console.
    m_console.PrintRow(out, m_console_mask);
}


// POST-PROCESSING ROUTINES.

void Solver::buildOutputVector(const Mops::Reactor &r, fvector &out) const
{
    fvector concs;
    unsigned int n = r.Mechanism()->SpeciesCount();

    // Get the data for output.
    r.Mixture()->GetConcs(concs);
    out.resize(n+5);
    out[0] = 0;
    out[1] = r.Time();
    for (int isp=0; isp<n; isp++) {
        out[isp+2] = concs[isp] * 1.0e-6; // Convert to mol/cm3.
    }
    out[n+2] = r.Mixture()->Temperature();
    out[n+3] = r.Mixture()->Density() * 1.0e-6; // Convert to mol/cm3.
    out[n+4] = r.Mixture()->Pressure();
}


void Solver::readAux(const std::string &filename, Mops::Mechanism &mech, 
                     Mops::timevector &times) const
{
    // Build input file name.
    string fname(filename); fname.append(".dat");

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
