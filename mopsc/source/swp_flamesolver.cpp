#include "swp_flamesolver.h"
#include "mops_timeinterval.h"
#include "sweep.h"
#include "string_functions.h"
#include "csv_io.h"
#include <fstream>
#include <stdexcept>
#include <string>

using namespace Sweep;
using namespace std;
using namespace Strings;

const real FlameSolver::CONFA = 3.29;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
FlameSolver::FlameSolver()
: m_stats(NULL)
{
}

// Default destructor.
FlameSolver::~FlameSolver()
{
    delete m_stats;
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
    string delim = ",\t "; // Possible file delimiters (comma, tab and space).
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
        for (unsigned int i=0; i!=subs.size(); ++i) {
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
            Sprog::Thermo::IdealGas gas(mech.Species());

            // Split the line by columns.
            split(line, subs, delim);

            // Loop over all the elements in the line and save them
            // to the correct gas-phase variable.
            for (unsigned int i=0; i!=subs.size(); ++i) {
                if (i==tcol) {
                    // This is the Time column.
                    t = cdble(subs[i]);                    
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
                        gas.RawData()[isp->second] = cdble(subs[i]);
                    }
                }
            }

            // Set up the gas-phase by setting temperature, pressure and
            // normalising the mixture fractions.
            // TODO:  This will give the wrong component densities
            //        unless all species are specified!
            gas.SetTemperature(T);
            gas.SetPressure(P*1.0e5);
            gas.Normalise();

            // Add the profile point.
            alpha_prof[t] = alpha;
            m_gasprof.insert(GasPoint(t, gas));
        }

        // Set up ABF model to use alpha profile.
        Sweep::ABFModel::Instance().Initialise(mech.ParticleMech());
        Sweep::ABFModel::Instance().SetAlphaProfile(alpha_prof);

        // Close the input file.
        fin.close();
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
//    Sprog::Thermo::IdealGas gas(*mech.Species());

    // Ensure the process counters contain sufficient 
    // entries to track all rate terms.
    m_processcounter.resize(mech.TermCount(), 0);
    m_ficticiouscounter.resize(mech.TermCount(), 0);

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
                t += dt;
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
    r.Initialise(t1);
    r.Mixture()->Particles().Initialise(m_pcount);
    r.Mixture()->SetMaxM0(m_maxm0);

    // Set up file output.
    beginFileOutput(*r.Mechanism(), times);

    // Set up the console output.
    icon = m_console_interval;
    setupConsole(r.Mechanism()->ParticleMech());

    for (unsigned int irun=0; irun!=nruns; ++irun) {
        // Start the CPU timing clock.
        m_cpu_start = clock();
        m_chemtime  = 0.0;

        // Re-initialise the solution.
        r.Mixture()->Reset(m_maxm0);
        t1 = times[0].StartTime();
        t2 = t1;

        // Begin file output for this run.
        beginRunFileOutput(irun);
        printf("Run Number %d of %d.\n", irun+1, nruns);

        // Output initial conditions.
        fileOutput(t1, *r.Mixture());
        consoleOutput(t1, *r.Mixture());

        // Loop over the time intervals.
        Mops::timevector::const_iterator iint;
        for (iint=times.begin(); iint!=times.end(); ++iint) {
            // Get the step size for this interval.
            dt = (*iint).StepSize();

            // Loop over the steps in this interval.
            for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep) {
                // Run the reactor solver for this step (timed).
                m_cpu_mark = clock();
                    Run(t1, t2+=dt, m_gasprof, *r.Mixture(), 
                        r.Mechanism()->ParticleMech());
                    r.SetTime(t1);
                m_chemtime += (double)(clock() - m_cpu_mark) / 
                              (double)CLOCKS_PER_SEC;

                // Generate file output.
                fileOutput(t1, *r.Mixture());

                // Generate console output.
                if (--icon == 0) {
                    consoleOutput(t1, *r.Mixture());
                    icon = m_console_interval;
                }
            }
        }

        // Close the output files.
        endFileOutput();
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

    // Declare stats outputs (averages and errors).
    fvector s;
    vector<fvector> as(npoints), es(npoints);
    EnsembleStats stats(pmech);

    // Declare CPU time outputs (averages and errors).
    double cpu; 
    vector<double> acpu(npoints), ecpu(npoints);


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
            // PARTICLE STATS.

            // Read the stats.
            stats.Deserialize(fin);
            stats.Get(s);

            // Calculate sums and sums of square (for average and
            // error calculation).
            as[0].resize(s.size(), 0.0);
            es[0].resize(s.size(), 0.0);
            for (unsigned int i=0; i!=s.size(); ++i) {
                as[0][i] += s[i];
                if (nruns>1) es[0][i] += s[i] * s[i];
            }

            // COMPUTATION TIMES.

            // Read the computation time.
            fin.read(reinterpret_cast<char*>(&cpu), sizeof(cpu));

            // Calculate sums and sums of square (for average and
            // error calculation).
            acpu[0] += cpu;
            if (nruns>1) ecpu[0] += cpu * cpu;
        }

        // Loop over all time intervals.
        unsigned int step = 1;
        for (Mops::timevector::const_iterator iint=times.begin(); 
             iint!=times.end(); ++iint) {
            // Loop over all time steps in this interval.
            for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
                // PARTICLE STATS.

                // Read the stats.
                stats.Deserialize(fin);
                stats.Get(s);

                // Calculate sums and sums of square (for average and
                // error calculation).
                as[step].resize(s.size(), 0.0);
                es[step].resize(s.size(), 0.0);
                for (unsigned int i=0; i!=s.size(); ++i) {
                    as[step][i] += s[i];
                    if (nruns>1) es[step][i] += s[i] * s[i];
                }

                // COMPUTATION TIMES.

                // Read the computation time.
                fin.read(reinterpret_cast<char*>(&cpu), sizeof(cpu));

                // Calculate sums and sums of square (for average and
                // error calculation).
                acpu[step] += cpu;
                if (nruns>1) ecpu[step] += cpu * cpu;
            }
        }

        // Close the input file.
        fin.close();
    }
    

    // CALCULATE AVERAGES AND CONFIDENCE INTERVALS.

    real invruns = 1.0 / (real)nruns;
    for (unsigned int step=0; step!=npoints; ++step) {
        // Particle statistics.
        for (unsigned int i=0; i!=s.size(); ++i) {
            // Calculate average over all runs.
            as[step][i] *= invruns;

            if (nruns > 1) {
                // Calculate standard deviation over all runs.
                (es[step][i] *= invruns) -= (as[step][i] * as[step][i]);
                // Calculate confidence interval.
                es[step][i] = CONFA * sqrt(abs(es[step][i] * invruns));
            }
        }

        // CPU time.
        acpu[step] *= invruns;
        if (nruns > 1) {
            (ecpu[step] *= invruns) -= (acpu[step] * acpu[step]);
            ecpu[step] = CONFA * sqrt(abs(ecpu[step] * invruns));
        }
    }

    // OUTPUT TO CSV FILE.
    
    // Now open a file for the CSV results and write header.
    CSV_IO csv(string(filename).append(".csv"), true);
    vector<string> header;
    header.push_back("Step");
    header.push_back("Time (s)");
    stats.Names(header, 2);
    for (unsigned int i=header.size(); i!=2; --i) {
        header.insert(header.begin()+i, "Err");
    }
    csv.Write(header);

    // Open a CSV file for the computation time results.
    CSV_IO csvcpu(string(filename).append("-ct.csv"), true);
    vector<string> cpuheader;
    cpuheader.push_back("Step");
    cpuheader.push_back("Time (s)");
    cpuheader.push_back("CPU Time (s)");
    cpuheader.push_back("Err");
    csvcpu.Write(cpuheader);

    // Output initial conditions.
    {
        // Stats.
        for (unsigned int i=as[0].size(); i!=0; --i) {
            as[0].insert(as[0].begin()+i, *(es[0].begin()+i-1));
        }
        as[0].insert(as[0].begin(), 0.0);
        as[0].insert(as[0].begin(), times[0].StartTime());
        csv.Write(as[0]);

        // CPU times.
        vector<double> cpuout(4);
        cpuout[0] = 0.0;
        cpuout[1] = times[0].StartTime();
        cpuout[2] = acpu[0];
        cpuout[3] = ecpu[0];
        csvcpu.Write(cpuout);
    }

    // Loop over all points, performing output.
    unsigned int step = 1;
    for (Mops::timevector::const_iterator iint=times.begin(); iint!=times.end(); ++iint) {
        // Loop over all time steps in this interval.
        real t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();

            // Stats.
            for (unsigned int i=as[step].size(); i!=0; --i) {
                as[step].insert(as[step].begin()+i, *(es[step].begin()+i-1));
            }
            as[step].insert(as[step].begin(), t);
            as[step].insert(as[step].begin(), step);
            csv.Write(as[step]);

            // CPU times.
            vector<double> cpuout(4);
            cpuout[0] = 0.0;
            cpuout[1] = t;
            cpuout[2] = acpu[step];
            cpuout[3] = ecpu[step];
            csvcpu.Write(cpuout);
        }
    }

    // Close the CSV files.
    csv.Close();
    csvcpu.Close();
}


// CONSOLE OUTPUT.

// Sets up console output using the given mechanism as a template.
void FlameSolver::setupConsole(const Sweep::Mechanism &mech)
{
    // Build header. 
    vector<string> header;
    header.push_back("Time (s)");
    header.push_back("#SP");
    header.push_back("M0");
    header.push_back("Fv");
    header.push_back("dcol (nm)");
    header.push_back("S (cm2/cm3)");

    // Print the first header.
    m_console.PrintDivider();
    m_console.PrintRow(header);
    m_console.PrintDivider();

    // Set up the console output class.
    m_console.SetAutoHeader(header);
    m_console.EnableAutoDividers();
}

// Writes current reactor state to the console.
void FlameSolver::consoleOutput(real time, const Sweep::Cell &sys) const
{
    // Get output data.
    fvector out;
    out.push_back(time);
    out.push_back(m_stats->BasicStats().PCount());
    out.push_back(m_stats->BasicStats().M0());
    out.push_back(m_stats->BasicStats().Fv());
    out.push_back(m_stats->BasicStats().AvgCollDiam());
    out.push_back(m_stats->BasicStats().SurfaceArea());

    // Print data to console.
    m_console.PrintRow(out);
}

// FILE OUTPUT.

// Sets up the file output by outputting an auxilliary file
// which stores all the information required to post-process the
// simulation and by opening the output file.
void FlameSolver::beginFileOutput(const Mops::Mechanism &mech, 
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
                            "output (Mops, Solver::beginFileOutput).");
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
    m_stats = new EnsembleStats(mech.ParticleMech());
}

// Sets up file output for a new run given the run number.
void FlameSolver::beginRunFileOutput(unsigned int run)
{
    // Build the simulation output file name.
    string fname = m_output_filename + "(" + cstr(run) + ").dat";

    // Open the simulation output file.
    m_file.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    // Throw error if the output file failed to open.
    if (!m_file.good()) {
        throw runtime_error("Failed to open file for simulation "
                            "output (Mops, Sweep::FlameSolver::beginRunFileOutput).");
    }
}

// Writes the current system state to the output file.
void FlameSolver::fileOutput(real time, const Sweep::Cell &sys)
{
    // Write particle stats to file.
    sys.GetVitalStats(*m_stats);
    m_stats->Serialize(m_file);

    // Write CPU times to file.
    double cputime = (double)(clock() - m_cpu_start) / (double)CLOCKS_PER_SEC;
    m_file.write((char*)&cputime, sizeof(cputime));
}

// Ends file output by closing all open files.
void FlameSolver::endFileOutput()
{
    // Close the output file.
    m_file.close();
}


// HELPER ROUTINES.

void FlameSolver::linInterpGas(Sweep::real t, 
                               const GasProfile &gasphase, 
                               Sprog::Thermo::IdealGas &gas) const
{
    // Get the time point after the required time.
    GasProfile::const_iterator j = gasphase.upper_bound(t);
    
    if (j == gasphase.begin()) {
        // This time is before the beginning of the profile.  Return
        // the first time point.
        gas = j->second;
    } else {       
        // Get the time point before the required time.
        GasProfile::const_iterator i = j; --i;

        // Assign the conditions to this point.
        gas = i->second;
        
        // Calculate time interval between points i and j.
        real dt_pro = j->first - i->first;

        // Calculate time interval between point i and current time.
        real dt = t - i->first;

        // Calculate the intermediate gas-phase mole fractions by linear
        // interpolation of the molar concentrations.
        real dens = 0.0;
        for (unsigned int k=0; k<gas.Species()->size(); ++k) {
            real dc = (j->second.MolarConc(k) - i->second.MolarConc(k)) * dt / dt_pro;
            gas.RawData()[k] = gas.MolarConc(k) + dc;
            dens += gas.RawData()[k];
        }
        gas.Normalise();

        // Now use linear interpolation to calculate the temperature.
        real dT = (j->second.Temperature() - i->second.Temperature()) * dt / dt_pro;
        gas.SetTemperature(gas.Temperature()+dT);

        // Now set the gas density, calculated using the values above.
        gas.SetDensity(dens);
    }
}
