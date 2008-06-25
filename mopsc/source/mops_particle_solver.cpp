#include "mops_particle_solver.h"
#include "mops_reactor_factory.h"
#include "csv_io.h"
#include "string_functions.h"
#include "sweep.h"
#include <vector>
#include <string>
#include <time.h>
#include <stdexcept>

using namespace Mops;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ParticleSolver::ParticleSolver(void)
: m_swp_ctime(0.0), m_ptrack_count(10)
{
}

// Default destructor.
ParticleSolver::~ParticleSolver(void)
{
}


// FILE OUTPUT.

// Writes the particle stats to the binary output file.
void ParticleSolver::outputParticleStats(const Reactor &r) const
{
    // Write particle stats to file.
    static Sweep::EnsembleStats stats(r.Mech()->ParticleMech());
    r.Mixture()->GetVitalStats(stats);
    stats.Serialize(m_file);
}

// Writes tracked particles to the binary output file.
void ParticleSolver::outputPartTrack(const Reactor &r) const
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

// Writes the current reactor state to the output file.  This overrides
// the default routine in order to also output the particle stats.
void ParticleSolver::fileOutput(const Reactor &r) const
{
    // Write gas-phase conditions to file.
    outputGasPhase(r);

    // Write particle stats to file.
    outputParticleStats(r);

    // Write CPU times to file.
    double cputime = calcDeltaCT(m_cpu_start);
    m_file.write((char*)&cputime, sizeof(cputime));
    m_file.write((char*)&m_chemtime, sizeof(m_chemtime));
    m_file.write((char*)&m_swp_ctime, sizeof(m_swp_ctime));

    // Do particle tracking output.
    outputPartTrack(r);
}


// POST-PROCESSING ROUTINES.

// Reads a particle stats data point from the binary file.
// To allow the averages and confidence intervals to be calculated
// the data point is added to a vector of sums, and the squares are
// added to the vector sumsqr if necessary.
void ParticleSolver::readParticleDataPoint(std::istream &in, 
                                           const Sweep::Mechanism &mech, 
                                           fvector &sum, fvector &sumsqr, 
                                           bool calcsqrs)
{
    // Check for valid stream.
    if (in.good()) {
        // Read the stats.
        Sweep::EnsembleStats stats(mech);
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

// Reads the tracked particles from the binary file.  The particles are
// processed so that only a vector of vectors is returned, which contains
// the PSL data for each tracked particle at that point.
void ParticleSolver::readPartTrackPoint(std::istream &in, 
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
        Sweep::EnsembleStats stats(mech);

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

// Writes particle stats profile to a CSV file.
void ParticleSolver::writeParticleStatsCSV(const std::string &filename, 
                                           const Mechanism &mech, 
                                           const timevector &times, 
                                           std::vector<fvector> &avg, 
                                           const std::vector<fvector> &err)
{
    // Open file for the CSV results.
    CSV_IO csv(filename, true);

    Sweep::EnsembleStats stats(mech.ParticleMech());

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

// Writes particle tracking for multiple particles to CSV files.
void ParticleSolver::writePartTrackCSV(const std::string &filename,
                                       const Mechanism &mech,
                                       const timevector &times,
                                       std::vector<std::vector<fvector> > &track)
{
    // The track vector<vector<fvector>> should be arranged thus:
    // time-steps / particles / PSL variables.

    // Create an EnsembleStats object to get the PSL variable
    // names vector.
    Sweep::EnsembleStats stats(mech.ParticleMech());

    // Build the header row vector
    vector<string> head;
    head.push_back("Step");
    head.push_back("Time (s)");
    stats.PSL_Names(head, 2);

    // Determine max. number of particles tracked.
    unsigned int np = 0;
    for (unsigned int i=0; i!=track.size(); ++i) {
        np = max(np, track[i].size());
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

// Writes computation times profile to a CSV file.  This overrides
// the default function in order to allow output of the Sweep CPU
// time.
void ParticleSolver::writeCT_CSV(const std::string &filename, 
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
    head.push_back("Sweep CPU Time (s)");
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


// SAVE POINTS AND PSL POST-PROCESSING.

// Processes the PSLs at each save point into single files.
void ParticleSolver::postProcessPSLs(unsigned int nruns, const Mechanism &mech, 
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
