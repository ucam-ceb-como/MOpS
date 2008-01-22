#include "swpoutput.h"
#include "swpensemble.h"
#include <fstream>
#include "morestrings.h"

using namespace Sweep;
using namespace std;

const real SweepOutput::CONFA = 3.29; // Gives 99% confidence interval.

SweepOutput::SweepOutput(void)
{
    m_filename = "";
}

SweepOutput::~SweepOutput(void)
{
    // Close the file if it is still open.
    if (m_file.is_open()) m_file.close();
}

void SweepOutput::CompileStatistics(const System &sys, const Mechanism &mech, vector<real> &stats)
{
    // Reset all stats to 0.
    stats.clear();
    stats.resize(NSTATS+mech.ComponentCount()+mech.ValueCount(), 0.0);

    // Get particle count.
    stats[0] = sys.ParticleCount();

    if (stats[0] > 0) {
        // stats[1] will be the number density once volume scaling has occured.
        stats[1] = stats[0];

        // Loop over all particles and compile stats.
        Ensemble::const_iterator i;
        vector<real>::const_iterator j;
        unsigned int k;
        for (i=sys.ConstEnsemble().begin(); i!=sys.ConstEnsemble().end(); i++) {
            // Stats.
            stats[2] += (*i)->Mass();
            stats[3] += (*i)->Volume();
            stats[4] += (*i)->SurfaceArea();
            stats[8] += (*i)->CollisionDiameter() * 1.0e7;

            // Component values.
            for (j=(*i)->Composition().begin(),k=NSTATS; j!=(*i)->Composition().end(); j++,k++) {
                stats[k] += *j;
            }

            // Value values.
            for (j=(*i)->Values().begin(); j!=(*i)->Values().end(); j++,k++) {
                stats[k] += *j;
            }
        }

        // Calculate stats which are particle averages.
        stats[5]  = stats[2] / stats[0];
        stats[6]  = stats[3] / stats[0];
        stats[7]  = stats[4] / stats[0];
        stats[8] /= stats[0];

        // Calculate stats which are bulk properties (per unit volume).
        real vol = sys.SampleVolume();
        stats[1] /= vol;
        stats[2] /= vol;
        stats[3] /= vol;
        stats[4] /= vol;
        vector<real>::iterator istat;
        for (istat=stats.begin()+NSTATS; istat!=stats.end(); istat++) {
            *istat /= vol;
        }

        // The sample volume (scaling).
        stats[9] = vol;
    }
}

void SweepOutput::CompilePSL(const Sweep::System &sys, std::vector<vector<real> > &psl)
{
    // Resize PSL to hold all particles.
    psl.clear();
    psl.resize(sys.ParticleCount());

    // Calculate the weight of each particle, scaled to give the true number
    // density on summation.
    real wt = 1.0 / sys.SampleVolume();

    // Loop over all particles, adding an entry to the psl for each.
    Ensemble::const_iterator isp;
    vector<vector<real>>::iterator ipsl;
    vector<real>::const_iterator ic;
    for (isp=sys.ConstEnsemble().begin(),ipsl=psl.begin(); isp!=sys.ConstEnsemble().end(); isp++,ipsl++) {
        ipsl->push_back(wt);
        // Calculated properties.
        ipsl->push_back((*isp)->Volume());
        ipsl->push_back((*isp)->SurfaceArea());
        ipsl->push_back((*isp)->CollisionDiameter() * 1.0e7);
        // Particle composition.
        for (ic=(*isp)->Composition().begin(); ic!=(*isp)->Composition().end(); ic++) {
            ipsl->push_back(*ic);
        }
        // particle values.
        for (ic=(*isp)->Values().begin(); ic!=(*isp)->Values().end(); ic++) {
            ipsl->push_back(*ic);
        }
    }
}

int SweepOutput::Open(const std::string &filename, const unsigned int n)
{
    // Make sure file is closed.
    if (m_file.is_open()) m_file.close();

    // Save file name with appended number.
    m_filename = BuildFilename(filename, n);

    // Open the file.
    m_file.open(m_filename.c_str(), ios_base::binary);
    if (m_file.good()) {
        return 0;
    } else {
        return -1;
    }
}

void SweepOutput::Close()
{
    // Close the file and wipe the file name.
    if (m_file.is_open()) m_file.close();
    m_filename = "";
}

int SweepOutput::Write(const Sweep::real t, const Sweep::System &sys, const Sweep::Mechanism &mech)
{
    if (m_file.is_open()) {
        // Gather output information.
        vector<real> stats; CompileStatistics(sys, mech, stats);
        vector<real>::iterator iter;

        // Write output to file.
        m_file.write(reinterpret_cast<char const*>(&t), sizeof t);
        for (iter=stats.begin(); iter!=stats.end(); iter++) {
            m_file.write(reinterpret_cast<char const*>(&(*iter)), sizeof *iter);
        }
        return 0;
    } else {
        // Can't write because file stream is closed.
        return -1;
    }
}

int SweepOutput::PostProcess(const std::string &filename, const unsigned int i1, 
                             const unsigned int i2, const unsigned int npoints,
                             const Sweep::Mechanism &mech) const
{
    // This function takes the binary file outputs from a sweep run and post-processes
    // them to give a CSV file with averages and confidence intervals.
    string fname;
    ifstream fin;
    unsigned int i, nstat=NSTATS+mech.ComponentCount()+mech.ValueCount();
    real val;
    vector<real> tavg(npoints,0.0); vector<vector<real>> savg(npoints);
    vector<real> terr(npoints,0.0); vector<vector<real>> serr(npoints);
    vector<real>::iterator iavg, ierr, itavg, iterr;
    vector<vector<real>>::iterator isavg, iserr;

    // Size data structures appropriately.
    for (isavg=savg.begin(),iserr=serr.begin(); 
         (isavg!=savg.end()) && (iserr!=serr.end()); 
         isavg++, iserr++)
    {
        (*isavg).resize(nstat, 0.0);
        (*iserr).resize(nstat, 0.0);
    }

    // Loop through all output files.
    for (i=i1; i<=i2; i++) {
        // Open file i for reading.
        fname = BuildFilename(filename, i);
        fin.open(fname.c_str(), ios_base::binary);
        if (!fin.good()) return -1;

        // Read rows from file.
        for (isavg=savg.begin(),iserr=serr.begin(),itavg=tavg.begin(),iterr=terr.begin(); 
             (isavg!=savg.end()) && (iserr!=serr.end()); 
             isavg++,iserr++,itavg++,iterr++)
        {
            // Read time.
            fin.read((char*)&val, sizeof val);
            *itavg += val; *iterr += (val*val); // Sums and sums of squares.

            // Read stats.
            for (iavg=isavg->begin(),ierr=iserr->begin(); 
                 (iavg!=isavg->end()) && (ierr!=iserr->end()); 
                 iavg++,ierr++)
            {
                fin.read((char*)&val, sizeof val);
                *iavg += val; *ierr += (val*val); // Sums and sums of squares.
            }
        }

        // Close file.
        fin.close();
    }

    // Open output CSV file.
    fname = filename.substr(0, filename.find_first_of("."));
    ofstream fout((fname+=".csv").c_str());

    // Write header to output file.
    vector<string> header; GetColumnNames(mech, header);
    vector<string>::iterator ihead;
    fout << "Time (s),error";
    for (ihead=header.begin(); ihead!=header.end(); ihead++) {
        fout << ',' << *ihead << ",error";
    }
    fout << '\n';

    // Now calculate averages and statistical errors, then write output file.
    real nruns = (real)(i2-i1+1);
    for (isavg=savg.begin(),iserr=serr.begin(),itavg=tavg.begin(),iterr=terr.begin(); 
         (isavg!=savg.end()) && (iserr!=serr.end()); 
         isavg++,iserr++,itavg++,iterr++)
    {
        // Time.
        *itavg /= nruns;
        *iterr = (*iterr / nruns) - (*itavg * *itavg);
        *iterr = CONFA * sqrt(*iterr / nruns);
        fout << *itavg << ',' << *iterr;

        // Stats.
        for (iavg=isavg->begin(),ierr=iserr->begin(); 
             (iavg!=isavg->end()) && (ierr!=iserr->end()); 
             iavg++,ierr++)
        {
            *iavg /= nruns;
            *ierr = (*ierr / nruns) - (*iavg * *iavg);
            *ierr = CONFA * sqrt(*ierr / nruns);
            fout << ',' << *iavg << ',' << *ierr;
        }
        // End-of-line.
        fout << '\n';
    }

    // Close the output, and we're done!
    fout.close();
    return 0;
}

int SweepOutput::WritePSL(const Sweep::real t, const Sweep::System &sys, const unsigned int i)
{
    if (m_filename != "") {
        return WritePSL(m_filename, t, sys, i);
    } else {
        // Unknown file name.
        return -1;
    }
}

int SweepOutput::WritePSL(const string &filename, const Sweep::real t, const Sweep::System &sys, 
                          const unsigned int i)
{
    // Build psl file name.
    string name = BuildPSLFilename(filename, i);

    // Open PSl output file.
    ofstream fpsl(name.c_str(), ios::binary);

    if (fpsl.good()) {
        // Compile psl.
        unsigned int n = sys.ParticleCount();
        vector<vector<real>> psl; CompilePSL(sys, psl);

        // Write PSL to file.
        fpsl.write(reinterpret_cast<char const*>(&t), sizeof t);
        fpsl.write(reinterpret_cast<char*>(&n), sizeof n);
        vector<vector<real>>::iterator ipsl;
        vector<real>::iterator irow;
        n=0; int nn=0;
        for (ipsl=psl.begin(); ipsl!=psl.end(); ipsl++) {
            n++;
            for (irow=ipsl->begin(); irow!=ipsl->end(); irow++) {
                fpsl.write(reinterpret_cast<char*>(&(*irow)), sizeof *irow);   
                nn++;
            }
        }
;
        // Close file and return.
        fpsl.close();
        return 0;
    } else {
        // Could not open file.
        return -1;
    }
}

int SweepOutput::PostProcessPSL(const std::string &filename, const unsigned int r1, 
                                const unsigned int r2, const std::vector<unsigned int> psln, 
                                const Sweep::Mechanism &mech)
{
    // Open output file.
    string name = filename.substr(0,filename.find_last_of(".",filename.length()))+="-psl.csv";
    ofstream fout(name.c_str());

    if (fout.good()) {
        // Write header to output file.
        vector<string> header; GetPSLColumnNames(mech, header);
        vector<string>::iterator ihead;
        fout << header[0];
        for (ihead=header.begin()+1; ihead!=header.end(); ihead++)
            fout << ',' << *ihead;
        fout << '\n';

        // Variables read from file.
        string pslfilename;
        ifstream fin;
        real t;
        unsigned int n;

        real nruns = (real)(r2-r1+1);

        // Counters and iterators.
        unsigned int r, i, ipsl;
        vector<unsigned int>::const_iterator ipsln;

        // Loop over all PSLs.
        for (ipsln=psln.begin(); ipsln!=psln.end(); ipsln++) {
            // Loop over all runs.
            for (r=r1; r<=r2; r++) {
                // Open PSL file for reading.
                pslfilename = BuildPSLFilename(BuildFilename(filename, r), *ipsln);
                fin.open(pslfilename.c_str(), ios::binary);

                // Read time and number of particles.
                fin.read((char*)&t, sizeof t);
                fin.read((char*)&n, sizeof n);

                // Loop over all particles.
                for (ipsl=0; ipsl<n; ipsl++) {
                    // Read weight and rescale it by number of runs.
                    fin.read((char*)&t, sizeof t);  fout << (t/nruns);

                    // Read particle properties, composition and values.
                    for (i=0; i<(NPSL+mech.ComponentCount()+mech.ValueCount()-1); i++) {
                        fin.read((char*)&t, sizeof t); fout << ',' << t;
                    }

                    // End-of-line.
                    fout << '\n';
                }

                // Close input file.
                fin.close();
            }
        }

        // Close output file and return.
        fout.close();
        return 0;
    } else {
        // Could not open output file.
        return -1;
    }
}

int SweepOutput::PostProcessPSL(const std::string &filename, const unsigned int r1, 
                                const unsigned int r2, const unsigned int psln, 
                                const Sweep::Mechanism &mech)
{
    vector<unsigned int> ns; ns.push_back(psln);
    return PostProcessPSL(filename, r1, r2, ns, mech);
}

void SweepOutput::PrintToConsole(const Sweep::real t, const Sweep::System &sys, const Sweep::Mechanism &mech)
{
    vector<real> stats; CompileStatistics(sys, mech, stats);
    vector<real> rates; mech.GetRates(rates, t, sys);
    //real rcoag = const_cast<Mechanism&>(mech).GetProcess(0).Rate(t, sys);
    printf("%f\t%d\t%e\t%e\t%e\n", t , (int)stats[0], stats[1], stats[NSTATS], sys.GetSpeciesConc("A4", t));
}