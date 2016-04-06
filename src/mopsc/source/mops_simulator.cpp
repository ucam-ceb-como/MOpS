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
//#define USE_MPI

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "mops_simulator.h"
#include "mops_ode_solver.h"
#include "mops_flux_postprocessor.h"
#include "mops_reactor_factory.h"
#include "string_functions.h"
#include "csv_io.h"
#include "geometry1d.h"
#include <stdexcept>
#include <memory>
#include "gpc_reaction_set.h"
#include "loi_reduction.h"
#include "mops_gpc_sensitivity.h"
#include "gpc_species.h"
#include "swp_particle_image.h"
#include <algorithm>
#include <time.h>

#include <boost/functional/hash.hpp>
#include <boost/random/mersenne_twister.hpp>

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
  m_output_step(0), m_output_iter(0), m_write_jumps(false),
//////////////////////////////////////////// aab64 ////////////////////////////////////////////
  m_write_diags(false), 
//////////////////////////////////////////// aab64 ////////////////////////////////////////////
  m_write_ensemble_file(false),
  m_write_PAH(false), m_mass_spectra(true), m_mass_spectra_ensemble(true), 
  m_mass_spectra_xmer(1), m_mass_spectra_frag(false), 
  m_ptrack_count(0)
{
}

// Default destructor.
Simulator::~Simulator(void)
{
}

// Copy constructor.
Simulator::Simulator(const Mops::Simulator &copy)
{
    *this = copy;
}

// Assignment operator
Simulator &Simulator::operator=(const Mops::Simulator &rhs) {
    if (this != &rhs) {
        m_nruns = rhs.m_nruns;
        m_niter = rhs.m_niter;
        m_pcount = rhs.m_pcount;
        m_maxm0 = rhs.m_maxm0;
        m_cpu_start = rhs.m_cpu_start;
        m_cpu_mark = rhs.m_cpu_mark;
        m_runtime = rhs.m_runtime;
        m_console_interval = rhs.m_console_interval;
        m_console_msgs = rhs.m_console_msgs;
        m_output_filename = rhs.m_output_filename;
        m_output_every_iter = rhs.m_output_every_iter;
        m_output_step = rhs.m_output_step;
        m_output_iter = rhs.m_output_iter;
        m_write_jumps = rhs.m_write_jumps;
//////////////////////////////////////////// aab64 ////////////////////////////////////////////
		m_write_diags = rhs.m_write_diags;
//////////////////////////////////////////// aab64 ////////////////////////////////////////////
        m_write_ensemble_file = rhs.m_write_ensemble_file;
        m_write_PAH = rhs.m_write_PAH;
        m_mass_spectra = rhs.m_mass_spectra;
        m_mass_spectra_ensemble = rhs.m_mass_spectra_ensemble;
        m_mass_spectra_xmer = rhs.m_mass_spectra_xmer;
        m_mass_spectra_frag = rhs.m_mass_spectra_frag;
        m_ptrack_count = rhs.m_ptrack_count;
    }
    return *this;
}


// SIMULATION SETTINGS.

// Returns the number of runs to perform
unsigned int Simulator::RunCount(void) const {return m_nruns;}

// Sets the number of runs to peform.
void Simulator::SetRunCount(unsigned int n) {m_nruns = n;}

// Returns the number of iteration to perform per step.
unsigned int Simulator::IterCount(void) const {return m_niter;}

// Return time vector.
const timevector &Simulator::TimeVector() const { return m_times; }

// Return the number of time steps.
unsigned int Simulator::TimeStepCount() const
{
    unsigned int n_time_steps = 0;
    for (unsigned int i =0; i < m_times.size(); ++i) {
        n_time_steps += m_times.at(i).StepCount();
    }
    return n_time_steps;
}

// Sets the number of iterations to perform per step.
void Simulator::SetIterCount(unsigned int n) {m_niter = n;}

// Sets the time vector.
void Simulator::SetTimeVector(const timevector &times) {m_times = times;}


// These particle count related settings should be part of the
// Reactor class!

// Returns the max. stochastic particle count.
unsigned int Simulator::MaxPartCount(void) const {return m_pcount;}

// Sets the max. stochastic particle count.
void Simulator::SetMaxPartCount(unsigned int n) {m_pcount = n;}

/*!
 *@return       Value which particle number density is not expected to exceed \f$\mathrm{m}^{-3}\f$
 */
double Simulator::MaxM0(void) const {return m_maxm0;}

/*!
 *@param[in]    m0      Value which particle number density is not expected to exceed \f$\mathrm{m}^{-3}\f$
 */
void Simulator::SetMaxM0(double m0) {m_maxm0 = m0;}

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

//! Set simulator to write the jumps CSV file.
void Simulator::SetWriteJumpFile(bool writejumps) {m_write_jumps=writejumps;}



//////////////////////////////////////////// aab64 ////////////////////////////////////////////
//! Set simulator to write the diagnosticss CSV file.
void Simulator::SetWriteDiagsFile(bool writediags) {m_write_diags = writediags;}
//////////////////////////////////////////// aab64 ////////////////////////////////////////////



//! Set simulator to write the jumps CSV file.
void Simulator::SetWriteEnsembleFile(bool writeensemble) {m_write_ensemble_file=writeensemble;}

//! Set simulator to write the detailed PAH info in the psl file.
void Simulator::SetWritePAH(bool postpocessPAH) {m_write_PAH=postpocessPAH;}

// STATISTICAL BOUNDS OUTPUT

void Simulator::SetOutputStatBoundary(Sweep::PropID pid, double lower, double upper)
{
    m_statbound.Lower = lower;
    m_statbound.Upper = upper;
    m_statbound.PID   = pid;
}

// POVRAY OUTPUT.

void Simulator::SetParticleTrackCount(unsigned int ptcount) {
    m_ptrack_count = ptcount;
}
// options for Postprocess (only for PAH-PP model)
const bool Simulator::MassSpectra() const
{
    return m_mass_spectra;
}

void Simulator::SetMassSpectra(const bool val)
{
    m_mass_spectra = val;
}

const bool Simulator::MassSpectraEnsemble() const
{
    return m_mass_spectra_ensemble;
}

void Simulator::SetMassSpectraEnsemble(const bool val)
{
    m_mass_spectra_ensemble = val;
}

const int Simulator::MassSpectraXmer() const
{
    return m_mass_spectra_xmer;
}

void Simulator::SetMassSpectraXmer(const int val)
{
    m_mass_spectra_xmer = val;
}

const bool Simulator::MassSpectraFrag() const
{
    return m_mass_spectra_frag;
}

void Simulator::SetMassSpectraFrag(const bool val)
{
    m_mass_spectra_frag = val;
}

// SOLVING REACTORS.

// Solves the given reactor for the given time intervals.
void Simulator::RunSimulation(Mops::Reactor &r,
                              Solver &s, size_t seed)
{
    unsigned int icon;
    double dt, t2; // Stop time for each step.

    // Make a copy of the initial mixture and store in an auto pointer
    // so that it will be deleted when we leave this scope.
    std::auto_ptr<Mixture> initmix(r.Mixture()->Clone());

    // Initialise the reactor with the start time.
    t2 = m_times[0].StartTime();
    r.SetTime(t2);
    //r.Mixture()->GasPhase().SetMaxM0(m_maxm0);

    // Set up the solver.
    s.Initialise(r);

    if (r.Mech()->GasMech().ReactionCount() == 0)
    {
        s.SetLOIStatusFalse();
    }

    // Set up file output.
	#ifdef USE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{                    //ms785 only write aux file if this is the master nodes
	#endif

    writeAux(*r.Mech(), m_times, s);
    openOutputFile(); //.sim

    // Output initial conditions.
    // - Note: added sensitivity output will write system initial condition, not initial values.
	#ifdef USE_MPI
	}
	if (rank==0)						//ms785
	#endif

    fileOutput(m_output_step, m_output_iter, r, s, this);

	#ifdef USE_MPI
	closeOutputFile();						//ms785
	#endif

    // Set up the console output.
    icon = m_console_interval;
    setupConsole(*r.Mech());
	string m_output_filename_base=m_output_filename;		//ms785

	// Loop over runs.
	#ifdef USE_MPI
	#else
    for (unsigned int irun=0; irun!=m_nruns; ++irun) {
	#endif

		#ifdef USE_MPI

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		int irun=rank;
		m_output_filename=m_output_filename_base+cstr(rank);				//ms785
		openOutputFile();

		#endif

		size_t runSeed = seed;
		boost::hash_combine(runSeed, irun);
		boost::mt19937 rng(runSeed);

        // Start the CPU timing clock.
        m_cpu_start = clock();
        m_runtime  = 0.0;

        // Initialise the reactor with the start time.
        t2 = m_times[0].StartTime();
        r.SetTime(t2);
        // also reset the contents of the reactor
        r.Fill(*(initmix->Clone()), true);

        // Set up the ODE solver for this run.
        s.Reset(r);

        // Print initial conditions to the console.
        printf("mops: Run number %d of %d.\n", irun+1, m_nruns);
        m_console.PrintDivider();
        consoleOutput(r);

        unsigned int istep;

        // Initialise some LOI stuff
        if (s.GetLOIStatus() == true) setupLOI(r, s);
		
		

//////////////////////////////////////////// aab64 ////////////////////////////////////////////
		if (m_write_diags) {
			/* Create partProc diagnostics csv file with pre/post split SV, #SPs, #events in split
			step including additions in LPDA, create gasConcFile with pre/post split concs */
			ofstream partProcFile, gasConcFile;
			int process_iter;
			std::vector<std::string> tmpPNames;

			r.Mech()->ParticleMech().GetProcessNames(tmpPNames, 0);

			// Add headers to partProc diagnostics file
			partProcFile.open("Part-split-diagnostics.csv");
			partProcFile << "Time (s)" << " , " << "Time out (s)" << " , " << "Step number (-)" << " , "
				<< "SV in (-)" << " , " << "SV out (-)" << " , "
				<< "SP in (-)" << " , " << "SP out (-)" << " , ";
			for (process_iter = 0; process_iter < tmpPNames.size() - 1; process_iter++) {
				partProcFile << tmpPNames[process_iter] << " , ";
			}
			partProcFile << "TransitionRegimeCoagulationTerms (kernel specific)";
			for (process_iter = tmpPNames.size()+1; process_iter < r.Mech()->ParticleMech().GetTermCount(); 
				process_iter++) {
				partProcFile << " , ";
			}
			partProcFile << "FictitiousCoagulationTerms (kernel specific)" << "\n";
			partProcFile.close();

			// Add headers to gasConc diagnostics file
			gasConcFile.open("Chem-split-diagnostics.csv");
			gasConcFile << "Time (s)" << " , " << "Time out (s)" << " , " << "Step number (-)" << " , ";
			for (process_iter = 0; process_iter < r.Mech()->GasMech().Species().size() - 1; process_iter++) {
				gasConcFile << r.Mech()->GasMech().Species(process_iter)->Name() << " pre-split (mol/m3)" << " , "
					<< r.Mech()->GasMech().Species(process_iter)->Name() << " post-split (mol/m3)" << " , ";
			}
			gasConcFile << "TiO2 pre-split (mol/m3)" << " , " << "TiO2 post-split (mol/m3)" << "\n";
			gasConcFile.close();
		}
//////////////////////////////////////////// aab64 ////////////////////////////////////////////



        // Loop over the time intervals.
        unsigned int global_step = 0;
        timevector::const_iterator iint;
        for (iint=m_times.begin(); iint!=m_times.end(); ++iint) {
            // Get the step size for this interval.
            dt = (*iint).StepSize();

            // Set output parameters for this time interval.
            m_output_step = max((int)iint->SplittingStepCount(), 0);
            m_output_iter = max((int)m_niter, 0);

            // Loop over the steps in this interval.
            for (istep=0; istep<iint->StepCount(); ++istep, ++global_step) {
                // Run the solver for this step (timed).
                m_cpu_mark = clock();
                


//////////////////////////////////////////// aab64 ////////////////////////////////////////////
/* If the solve function with diagnostics capacity replaces the original solve function, this 
   if statement is no longer necessary */
				if (m_write_diags) {
					s.Solve(r, t2 += dt, iint->SplittingStepCount(), m_niter,
						rng, &fileOutput, (void*)this, m_write_diags);
				} else {
					s.Solve(r, t2 += dt, iint->SplittingStepCount(), m_niter,
						rng, &fileOutput, (void*)this);
				}
//////////////////////////////////////////// aab64 ////////////////////////////////////////////

				

                //Set up and solve Jacobian here
                if (s.GetLOIStatus() == true)
                    solveLOIJacobian(r, s, istep, t2);

                m_runtime += calcDeltaCT(m_cpu_mark);
                // Generate console output.
				#ifdef USE_MPI					//ms785
				if (rank==0)
				#endif
                if (--icon == 0) {
                    consoleOutput(r);
                    icon = m_console_interval;
                }
            } // number of steps

            // Create a save point at the end of this time
            // interval.

			m_output_filename=m_output_filename_base;

            //@todo Reinstate fractal dimension calculations for
            // the PAH-PP model

            createSavePoint(r, global_step, irun);
            if (s.GetLOIStatus() == true){
                r.DestroyJac(m_loi_J, r.Mech()->GasMech().SpeciesCount());
            }

            // Write the ensemble or gas-phase files
            if (m_write_ensemble_file) createEnsembleFile(r, global_step, irun);

        } // number of time intervals
        if (s.GetLOIStatus() == true) {
            std::vector<std::string> rejects;
            LOIReduction::RejectSpecies(m_loi_data, s.ReturnCompValue(), r.Mech(), rejects, s.ReturnKeptSpecies());
            r.Mech()->GasMech().WriteReducedMech(OutputFile() + std::string("-kept.inp"), rejects);
        }

        // Print run time to the console.
        printf("mops: Run number %d completed in %.1f s.\n", irun+1, m_runtime);

        // Reset the process jump count
        r.Mech()->ParticleMech().ResetJumpCount();

		// currently this function is limited to PAH-PP model
		// Produce a file named "primary" which stores information of target (criteria are hard-coded) primary particle
		//r.Mech()->ParticleMech().Mass_pah(r.Mixture()->Particles());

	#ifdef USE_MPI

		 closeOutputFile();			//ms785
	#else
	} // (number of runs)
	#endif

	#ifdef USE_MPI
	#else
    // Close the output files.
    closeOutputFile();
	#endif

    // If we have a PSR, clear any stream memory.
    if (r.SerialType() == Mops::Serial_PSR) {
        Mops::PSR* psr = dynamic_cast<Mops::PSR *>(&r);
        psr->ClearStreamMemory();
    }
}

/*!
 * Set up the LOI details before a run
 *
 * @param r         Reactor object
 * @param s         Solver object
 */
void Simulator::setupLOI(const Mops::Reactor &r, Mops::Solver &s) {

    m_loi_file.open(LOIReduction::buildLOIFileName(OutputFile()).c_str());
    s.InitialiseSensMatrix(s.GetNumSens(),r.Mech()->GasMech().SpeciesCount());
    LOIReduction::CreateLOIFile(m_loi_file, r.Mech());

    m_loi_data.clear();
    fvector z = Mops::fvector(r.Mech()->GasMech().SpeciesCount());
    m_loi_data.resize(s.GetNumSens(), z);

}


/*!
 * Solve the Jacobian matrix for LOI calculations
 *
 * @param r         Reactor object
 * @param s         Solver object
 * @param istep     Step number
 * @param t2        Time
 */
void Simulator::solveLOIJacobian(
        Mops::Reactor &r,
        Mops::Solver &s,
        const unsigned int istep,
        const double t2) {
    if (istep == 0) {
        m_loi_J = r.CreateJac(r.Mech()->GasMech().SpeciesCount());
    }
    // wjm34: 1.0e-7 is hard-coded, not sure why?
    r.RateJacobian(t2, r.Mixture()->GasPhase().RawData(), m_loi_J, 1.0e-7);
    m_loi_data = LOIReduction::CalcLOI(
            m_loi_J,
            s.GetSensSolution(
                    s.GetNumSens(),
                    r.Mech()->GasMech().SpeciesCount()),
            m_loi_data,
            r.Mech()->GasMech().SpeciesCount(),
            s.GetNumSens()
            );
    LOIReduction::SaveLOI(m_loi_data, t2, m_loi_file, r.Mech());
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
	
    // Calculate number of output points.
    unsigned int npoints = 1; // 1 for initial conditions.
    for(Mops::timevector::const_iterator i=times.begin(); i!=times.end(); ++i) {
        npoints += i->StepCount();
    }

    Simulator::postProcessSimulation(mech, times, npoints, ncput, cput_head);

}

/*!
 * Post process the results from a binary file.
 *
 * @param mech      Gas/particle mechanism
 * @param times     Time vector describing simulation output points
 * @param npoints   Number of points
 * @param ncput     Number of CPU time poinst
 * @param cput_head CPU time headers
 */
void Simulator::postProcessSimulation(
    Mops::Mechanism &mech,
    Mops::timevector &times,
    unsigned int npoints,
    unsigned int ncput,
    std::vector<std::string> cput_head
    ) const {

    // Get reference to particle mechanism.
    Sweep::Mechanism &pmech = mech.ParticleMech();

    // Declare chemistry outputs (averages and errors).
    vector<fvector> achem(npoints), echem(npoints);

    // Declare particle stats outputs (averages and errors).
    vector<fvector> astat(npoints), estat(npoints);
    Sweep::Stats::EnsembleStats stats(pmech);

    // Declare particle-phase number of jumps output (averages and errors).
    vector<fvector> appjumps(npoints), eppjumps(npoints);

    // Declare CPU time outputs (averages and errors).
    vector<vector<double> > acpu(npoints), ecpu(npoints);

    // Declare gas-phase rate outputs (averages and errors).
    vector<fvector> agprates(npoints), egprates(npoints);
    vector<fvector> agpfwdrates(npoints), egpfwdrates(npoints);
    vector<fvector> agprevrates(npoints), egprevrates(npoints);
    vector<fvector> agpwdot(npoints), egpwdot(npoints);
    vector<fvector> agpsdot(npoints), egpsdot(npoints); // Added by mm864 

    // Declare particle-phase rate outputs (averages and errors). (At the moment only gas species goes into particle mech)
    vector<fvector> apprates(npoints), epprates(npoints);
    vector<fvector> appwdot(npoints), eppwdot(npoints);

    // Declare data structure for particle tracking data.
    // Runs -> time steps -> particles -> coordinate.
    vector<vector<vector<fvector> > > ptrack(m_nruns);

    // OPEN SIMULATION OUTPUT FILES.

    // Build the simulation input file name.
    string fname = m_output_filename + ".sim";
    std::ostringstream ranstream;
//  ranstream << getpid();
//  string fname = "/scratch/ms785/"+m_output_filename+ranstream.str()+".sim";

    // Open the simulation input file.
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
                        agpwdot[0], egpwdot[0], agpsdot[0], egpsdot[0], // added by mm864
                        true);
    readPartRxnDataPoint(fin, mech.ParticleMech(),
                         apprates[0], epprates[0],
                         appwdot[0], eppwdot[0],
                         appjumps[0], eppjumps[0],
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
    multVals(agpsdot[0], m_nruns*m_niter); // added by mm864
        multVals(apprates[0], m_nruns*m_niter);
        multVals(appjumps[0], m_nruns*m_niter);
        multVals(appwdot[0], m_nruns*m_niter);
        multVals(acpu[0], m_nruns*m_niter);
        multVals(echem[0], m_nruns*m_niter);
        multVals(estat[0], m_nruns*m_niter);
        multVals(egprates[0], m_nruns*m_niter);
        multVals(egpfwdrates[0], m_nruns*m_niter);
        multVals(egprevrates[0], m_nruns*m_niter);
        multVals(egpwdot[0], m_nruns*m_niter);
    multVals(egpsdot[0], m_nruns*m_niter); // added by mm864
        multVals(epprates[0], m_nruns*m_niter);
        multVals(eppjumps[0], m_nruns*m_niter);
        multVals(eppwdot[0], m_nruns*m_niter);
        multVals(ecpu[0], m_nruns*m_niter);
    } else {
        multVals(achem[0], m_nruns);
        multVals(astat[0], m_nruns);
        multVals(agprates[0], m_nruns);
        multVals(agpfwdrates[0], m_nruns);
        multVals(agprevrates[0], m_nruns);
        multVals(agpwdot[0], m_nruns);
    multVals(agpsdot[0], m_nruns); // added by mm864
        multVals(apprates[0], m_nruns);
        multVals(appjumps[0], m_nruns);
        multVals(appwdot[0], m_nruns);
        multVals(acpu[0], m_nruns);
        multVals(echem[0], m_nruns); 
        multVals(estat[0], m_nruns);
        multVals(egprates[0], m_nruns);
        multVals(egpfwdrates[0], m_nruns);
        multVals(egprevrates[0], m_nruns);
        multVals(egpwdot[0], m_nruns);
    multVals(egpsdot[0], m_nruns); // added by mm864
        multVals(epprates[0], m_nruns);
        multVals(eppjumps[0], m_nruns);
        multVals(eppwdot[0], m_nruns);
        multVals(ecpu[0], m_nruns);
    }

    // READ ALL OUTPUT POINTS.

    // Now we must read the reactor conditions at all time points
    // and all runs.
    cout << "mops: postprocessing "<<m_nruns<<" runs"<<endl;
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
                                            agpwdot[step], egpwdot[step], agpsdot[step], egpsdot[step], // modified by mm864
                                            true);

                        readPartRxnDataPoint(fin, mech.ParticleMech(), apprates[step], epprates[step], appwdot[step], eppwdot[step], appjumps[step], eppjumps[step], true);
                        readCTDataPoint(fin, ncput, acpu[step], ecpu[step], true);
                        readPartTrackPoint(fin, pmech, ptrack[irun][step]);
                    }
                }
                else  {
                    // Read single output point for step.
                    readGasPhaseDataPoint(fin, mech, achem[step], echem[step], true);
                    readParticleDataPoint(fin, pmech, astat[step], estat[step], true);

                    readGasRxnDataPoint(fin, mech,
                                        agprates[step], egprates[step],
                                        agpfwdrates[step], egpfwdrates[step],
                                        agprevrates[step], egprevrates[step],
                                        agpwdot[step], egpwdot[step], agpsdot[step], egpsdot[step], // modified by mm864
                                        true);

                    readPartRxnDataPoint(fin, mech.ParticleMech(), apprates[step], epprates[step], appwdot[step], eppwdot[step], appjumps[step], eppjumps[step], true);
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
    calcAvgConf(agpsdot, egpsdot, m_nruns*m_niter); // added by mm864
        calcAvgConf(apprates, epprates, m_nruns*m_niter);
        calcAvgConf(appjumps, eppjumps, m_nruns*m_niter);
        calcAvgConf(appwdot, eppwdot, m_nruns*m_niter);
        calcAvgConf(acpu, ecpu, m_nruns*m_niter);
    } else {
        calcAvgConf(achem, echem, m_nruns);
        calcAvgConf(astat, estat, m_nruns);
        calcAvgConf(agprates, egprates, m_nruns);
        calcAvgConf(agpfwdrates, egpfwdrates, m_nruns);
        calcAvgConf(agprevrates, egprevrates, m_nruns);
        calcAvgConf(agpwdot, egpwdot, m_nruns);
    calcAvgConf(agpsdot, egpsdot, m_nruns); // added by mm864
        calcAvgConf(apprates, epprates, m_nruns);
        calcAvgConf(appjumps, eppjumps, m_nruns);
        calcAvgConf(appwdot, eppwdot, m_nruns);
        calcAvgConf(acpu, ecpu, m_nruns);
    }

    // POST-PROCESS ELEMENT FLUX
    // Element flux must be output before CSV files since writeXXXCSV will change the contents of what it output afterwards
    writeElementFluxOutput(m_output_filename, mech, times, agpfwdrates, agprevrates, achem);

    // Sensitivity output.
    Mops::SensitivityAnalyzer::PostProcess(m_output_filename);

    // OUTPUT TO CSV FILES.

    writeGasPhaseCSV(m_output_filename+"-chem.csv", mech, times, achem, echem);
    writeParticleStatsCSV(m_output_filename+"-part.csv", mech, times, astat, estat);
    writeGasRxnCSV(m_output_filename+"-gp-rates.csv", mech, times, agprates, egprates); // this is for q_i (rate of progress) 
    //writeGasRxnCSV(m_output_filename+"-gp-fwd-rates.csv", mech, times, agpfwdrates, egpfwdrates);
    //writeGasRxnCSV(m_output_filename+"-gp-rev-rates.csv", mech, times, agprevrates, egprevrates);
    writeProdRatesCSV(m_output_filename+"-gp-wdot.csv", mech, times, agpwdot, egpwdot);
    writeSurfaceProdRatesCSV(m_output_filename+"-gp-sdot.csv", mech, times, agpsdot, egpsdot); // added by mm864
    writePartProcCSV(m_output_filename+"-part-rates.csv", mech.ParticleMech(), times, apprates, epprates);
    writeProdRatesCSV(m_output_filename+"-part-wdot.csv", mech, times, appwdot, eppwdot);
    writeCT_CSV(m_output_filename+"-cput.csv", times, acpu, ecpu, cput_head);
    for(unsigned int irun=0; irun!=m_nruns; ++irun) {
        writePartTrackCSV(m_output_filename+"("+cstr(irun)+")-track", mech,
                          times, ptrack[irun]);
    }

    // Only write jump file if the flag has been set.
    if (m_write_jumps) {
        writePartJumpCSV(m_output_filename+"-part-jumps.csv", mech.ParticleMech(), times, appjumps, eppjumps);
    }

    // POST-PROCESS PSLs.

    // Now post-process the ensemble to find interested information, in this case, mass of Xmer
    if (MassSpectra())
        postProcessXmer(mech, times);

    // Now post-process the PSLs.
    postProcessPSLs(mech, times);

    // Now post-process the ensemble to find interested information, in this case, PAH mass distribution of the largest soot aggregate in the ensemnble
    if (m_write_PAH && pmech.WriteBinaryTrees()) {
        // more potential functionality can be added in the below function
        postProcessPAHinfo(mech, times);
        // then output details about the whole particle ensemble
        postProcessPAHPSLs(mech, times);
    }
}


// FILE OUTPUT (protected).

// Opens the simulation output file.
void Simulator::openOutputFile() const
{
    // SIMULATOR OUTPUT FILE
    // Build the simulation output file name.
     string fname = m_output_filename + ".sim";
//	std::ostringstream ranstream;
//	ranstream <<getpid();
	//string fname = "/scratch/ms785/"+m_output_filename+ranstream.str()+".sim";
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
// the binary output file. (No need modification, modification when read) 
void Simulator::outputGasPhase(const Reactor &r) const
{
    // Write gas-phase conditions to file.
    m_file.write(reinterpret_cast<const char*>(&r.Mixture()->GasPhase().RawData()[0]),
                 sizeof(r.Mixture()->GasPhase().RawData()[0]) *
                 r.Mech()->GasMech().SpeciesCount());
    double T = r.Mixture()->GasPhase().Temperature();
    m_file.write((char*)&T, sizeof(T));
    double D = r.Mixture()->GasPhase().Density();
    m_file.write((char*)&D, sizeof(D));
    double P = r.Mixture()->GasPhase().Pressure();
    m_file.write((char*)&P, sizeof(P));
}

// Writes the particle stats to the binary output file.
void Simulator::outputParticleStats(const Reactor &r) const
{
    // Write particle stats to file.
    Sweep::Stats::EnsembleStats stats(r.Mech()->ParticleMech());
    stats.SetStatBoundary(m_statbound);
    r.Mixture()->GetVitalStats(stats);
    stats.Serialize(m_file);
}

/*
 * @brief Writes tracked particles to the binary output file.
 *
 * @param[in]    r    Reactor to output
 */
void Simulator::outputPartTrack(const Reactor &r) const
{
    // Write the number of tracked particles.
    unsigned int n = min(r.Mixture()->ParticleCount(), m_ptrack_count);
    m_file.write((char*)&n, sizeof(n));

    // Output the current time.
    double t = (double)r.Time();
    m_file.write((char*)&t, sizeof(t));

    if (n > 0) {
        // Serialize the particles for tracking
		// Provide a way to detect multiple instances of PAHs
        std::map<void*, boost::shared_ptr<Sweep::AggModels::PAHPrimary> > duplicates;
        for (unsigned int i=0; i!=n; ++i) {
            r.Mixture()->Particles().At(i)->Serialize(m_file, &duplicates);
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
    if(r.Mech()->GasMech().ReactionCount() > 0) {
        // Calculate the rates-of-progress.
        static fvector rop, rfwd, rrev;
        r.Mech()->GasMech().Reactions().GetRatesOfProgress(r.Mixture()->GasPhase(), rop, rfwd, rrev); // GetRatesOfProgress 6

        // Calculate the molar production rates.
        static fvector wdot, sdot;
        r.Mech()->GasMech().Reactions().GetMolarProdRates(rop, wdot);
	r.Mech()->GasMech().Reactions().GetSurfaceMolarProdRates(rop, sdot); // added by mm864
        // Write rates to the file.
        m_file.write((char*)&rop[0], sizeof(rop[0]) * r.Mech()->GasMech().ReactionCount());
        m_file.write((char*)&rfwd[0], sizeof(rfwd[0]) * r.Mech()->GasMech().ReactionCount());
        m_file.write((char*)&rrev[0], sizeof(rrev[0]) * r.Mech()->GasMech().ReactionCount());
        m_file.write((char*)&wdot[0], sizeof(wdot[0]) * r.Mech()->GasMech().SpeciesCount());
	m_file.write((char*)&sdot[0], sizeof(sdot[0]) * r.Mech()->GasMech().SpeciesCount());
    }
}

// Writes the particle process rates and the
// species molar production rates to the binary output file.
void Simulator::outputPartRxnRates(const Reactor &r) const
{
    if (r.Mech()->ParticleMech().ProcessCount() != 0) {
        // Calculate the process rates.
        static fvector rates;
        r.Mech()->ParticleMech().CalcRates(r.Time(), *r.Mixture(), Geometry::LocalGeometry1d(), rates);

        // Calculate the molar production rates (mol/mol).
        static fvector wdot;
        r.Mech()->ParticleMech().CalcGasChangeRates(r.Time(), *r.Mixture(), Geometry::LocalGeometry1d(), wdot);

        // Calculate the number of jumps (-).
        static fvector jumps;
        r.Mech()->ParticleMech().CalcJumps(r.Time(), *r.Mixture(), Geometry::LocalGeometry1d(), jumps);

        // Now convert from mol/mol to mol/m3.
        fvector::iterator rhodot = wdot.begin()+r.Mech()->GasMech().SpeciesCount()+1;
        for (unsigned int k=0; k!=r.Mech()->GasMech().SpeciesCount(); ++k) {
            wdot[k] = (r.Mixture()->GasPhase().Density() * wdot[k]) +
                      (r.Mixture()->GasPhase().MoleFraction(k) * (*rhodot));
        }

        // Write rates to the file.
        m_file.write((char*)&rates[0], sizeof(rates[0]) * r.Mech()->ParticleMech().ProcessCount());
        m_file.write((char*)&wdot[0], sizeof(wdot[0]) * r.Mech()->GasMech().SpeciesCount());
        m_file.write((char*)&jumps[0], sizeof(jumps[0]) * r.Mech()->ParticleMech().ProcessCount());
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

            // Write gas-phase reaction rates
            me->outputGasRxnRates(r);

            me->outputPartRxnRates(r);

            // Write CPU times to file.
            s.OutputCT(me->m_file);

            // Do particle tracking output.
            me->outputPartTrack(r);

            // Write sensitivityto file.
            s.OutputSensitivity(me->m_senfile, r, me);
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
            m_console_mask.push_back(mech.GasMech().Species().size()+1);
        } else if ((*i).compare("T")==0) {
            // This variable is the mixture temperature.
            header.push_back("T (K)");
            m_console_mask.push_back(mech.GasMech().Species().size());
        } else if ((*i).compare("TIME")==0 || (*i).compare("time")==0 ||
                   (*i).compare("Time")==0) {
            // This variable is the flow time.
            header.push_back("Time (s)");
            m_console_mask.push_back(mech.GasMech().Species().size()+2);
        } else if ((*i).compare("#SP")==0 || (*i).compare("#sp")==0) {
            // Number of stochastic particles.
            header.push_back("#SP");
            m_console_mask.push_back(mech.GasMech().Species().size()+3);
        } else if ((*i).compare("M0")==0 || (*i).compare("m0")==0) {
            // Particle number density.
            header.push_back("M0 (m-3)");
            m_console_mask.push_back(mech.GasMech().Species().size()+4);
        } else if ((*i).compare("FV")==0 || (*i).compare("fv")==0 ||
                   (*i).compare("Fv")==0) {
            // Particle volume fraction.
            header.push_back("Fv");
            m_console_mask.push_back(mech.GasMech().Species().size()+5);
        } else if ((*i).compare("CT")==0 || (*i).compare("ct")==0) {
            // Computation time.
            header.push_back("CPU (s)");
            m_console_mask.push_back(mech.GasMech().Species().size()+6);
        } else {
            // Check for a species name.
            int isp = mech.GasMech().FindSpecies((*i));
            if (isp >= 0) {
                // This is a valid species name.
                header.push_back((*i) + " (mol/m3)");
                m_console_mask.push_back(isp);
            } else {
                // This is an invalid variable.  Just print the first species.
                header.push_back(mech.GasMech().Species()[0]->Name());
                m_console_mask.push_back(0);
            }
        }
    }
	#ifdef USE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank==0)
	{
	#endif
    // Print the first header.
    m_console.PrintDivider();
    m_console.PrintRow(header);
    m_console.PrintDivider();

    // Set up the console output class.
    m_console.SetAutoHeader(header);
    m_console.EnableAutoDividers();
	#ifdef USE_MPI
	}
	#endif
}

// Writes output to the console.
void Simulator::consoleOutput(const Mops::Reactor &r) const
{
    // Get output data from gas-phase.
    static vector<double> out;
    r.Mixture()->GasPhase().GetConcs(out);
    out.push_back(r.Mixture()->GasPhase().Temperature());
    out.push_back(r.Mixture()->GasPhase().Density());
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

// Reads a gas-phase and surface chemistry data point from the binary file.
// To allow the averages and confidence intervals to be calculated
// the data point is added to a vector of sums, and the squares are
// added to the vector sumsqr if necessary.
void Simulator::readGasPhaseDataPoint(std::istream &in, const Mops::Mechanism &mech,
                                      fvector &sum, fvector &sumsqr, bool calcsqrs)
{
    // Check for valid stream.
    if (in.good()) {
        unsigned int N = mech.GasMech().SpeciesCount();
	unsigned int N_gas = mech.GasMech().GasSpeciesCount();
        // Read the gas-phase conditions.
        fvector y(N, 0.0);
        double T=0.0, D=0.0, P=0.0;

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

	// Separated by mm864 

	// Gas species 
        for (unsigned int i=0; i!=N_gas; ++i) {
            // Note conversion to cm-3 from m-3.
            double conc = D * y[i];
            sum[i] += conc;
            if (calcsqrs) sumsqr[i] += conc * conc;
        }

	// Surface species 
	for (unsigned int i=N_gas; i!=N; ++i) {

	string spName = mech.GasMech().GetSpecies(i)->Name(); 
	string phName = mech.GasMech().GetSpecies(spName)->PhaseName();
	double site_d =  mech.GasMech().FindSiteDensity(phName);
	int sp_occ = mech.GasMech().FindSiteOccup(spName);
	// Note conversion to cm-2 from m-2.
            double conc = y[i] * site_d/sp_occ * 1.0e-4;	
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
        double *cpu = new double[N];
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
                                    fvector &wdot_sum, fvector &wdot_sumsqr, fvector &sdot_sum, fvector &sdot_sumsqr, 
                                    bool calcsqrs)
{
    // Check for valid stream.
    if (in.good()) {
        // Create vectors into which to read the file data
        const unsigned int rxnCount = mech.GasMech().ReactionCount();
        fvector rop(rxnCount);
        fvector rfwd(rxnCount);
        fvector rrev(rxnCount);
        fvector wdot(mech.GasMech().SpeciesCount());
	fvector sdot(mech.GasMech().SpeciesCount());	

        // Resize vectors passed in as reference arguments.
        rates_sum.resize(rxnCount, 0.0);
        rates_sumsqr.resize(rxnCount, 0.0);
        fwd_rates_sum.resize(rxnCount, 0.0);
        fwd_rates_sumsqr.resize(rxnCount, 0.0);
        rev_rates_sum.resize(rxnCount, 0.0);
        rev_rates_sumsqr.resize(rxnCount, 0.0);
        wdot_sum.resize(mech.GasMech().SpeciesCount(), 0.0);
        wdot_sumsqr.resize(mech.GasMech().SpeciesCount(), 0.0);
	sdot_sum.resize(mech.GasMech().SpeciesCount(), 0.0);
        sdot_sumsqr.resize(mech.GasMech().SpeciesCount(), 0.0);

        // Nothing to read or calculate if there are no reactions
        if(rxnCount > 0) {
            // Get the reaction rates-of-progress vector.
            in.read(reinterpret_cast<char*>(&rop[0]), sizeof(rop[0])*mech.GasMech().ReactionCount());

            // Get the forward reaction rates vector.
            in.read(reinterpret_cast<char*>(&rfwd[0]), sizeof(rop[0])*mech.GasMech().ReactionCount());

            // Get the reverse reaction rates vector.
            in.read(reinterpret_cast<char*>(&rrev[0]), sizeof(rop[0])*mech.GasMech().ReactionCount());

            // Get the species molar production rates.
            in.read(reinterpret_cast<char*>(&wdot[0]), sizeof(wdot[0])*mech.GasMech().SpeciesCount());

	     // Get the species surface molar production rates.
            in.read(reinterpret_cast<char*>(&sdot[0]), sizeof(sdot[0])*mech.GasMech().SpeciesCount());

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
            for (unsigned int i=0; i!=sdot.size(); ++i) {
                sdot_sum[i] += sdot[i];
                if (calcsqrs) sdot_sumsqr[i] += (sdot[i] * sdot[i]);
            }
		
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
                                     fvector &jumps_sum, fvector &jumps_sumsqr,
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

        // Get the number of jumps vector
        fvector jumps(mech.ProcessCount());
        in.read(reinterpret_cast<char*>(&jumps[0]),
                sizeof(jumps[0])*mech.ProcessCount());


        // Resize vectors.
        rates_sum.resize(rates.size(), 0.0);
        rates_sumsqr.resize(rates.size(), 0.0);
        wdot_sum.resize(wdot.size(), 0.0);
        wdot_sumsqr.resize(wdot.size(), 0.0);
        jumps_sum.resize(jumps.size(), 0.0);
        jumps_sumsqr.resize(jumps.size(), 0.0);

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
        for (unsigned int i=0; i!=jumps.size(); ++i) {
            jumps_sum[i] += jumps[i];
            if (calcsqrs) jumps_sumsqr[i] += (jumps[i] * jumps[i]);
        }
    }
}

/*
 * @brief Reads the tracked particles from the binary file
 *
 * @param[in]        in       Input binary stream from which data is read
 * @param[in]        mech     Mechanism defining the particle processes
 * @param[in,out]    pdata    A vector (number of tracked particles) of vectors (PSL data) at that time point
 */
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
		// Provide a way to detect multiple instances of PAHs
        std::map<void*, boost::shared_ptr<Sweep::AggModels::PAHPrimary> > duplicates;
        i = pdata.begin();
        for (unsigned int j=0; j!=n; ++j) {
            Sweep::Particle sp(in, mech, &duplicates);
            stats.PSL(sp, mech, (double)t, *(i++), 1.0);
        }
        Sweep::Particle empty(t, mech);
        for (unsigned int j=n; j!=m_ptrack_count; ++j)
        {
            stats.PSL(empty, mech, (double)t, *(i++), 1.0);
        }
    }
}

// Multiplies all values in a vector by a scaling factor.
void Simulator::multVals(fvector &vals, double scale)
{
    for (fvector::iterator i=vals.begin(); i!=vals.end(); ++i) {
        (*i) *= scale;
    }
}

// Takes vectors of vectors of variable sums and sums of squares, which
// are converted into the average values and the CLT confidence intervals
// for these averages.
void Simulator::calcAvgConf(std::vector<fvector> &avg,
                            std::vector<fvector> &err,
                            unsigned int nruns)
{
    const double CONFA = 3.29; // for 99.9% confidence interval.

    // Pre-calc some useful values.
    double invruns = 1.0 / (double)nruns;
//    double invruns_1 = 1.0 / (double)(nruns-1);
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
void Simulator::buildOutputVector(unsigned int step, double time,
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
    for (unsigned int isp=0; isp<mech.GasMech().GasSpeciesCount(); ++isp) {
        head.push_back(mech.GasMech().Species(isp)->Name() + " (mol/cm3)");
    }
     for (unsigned int isp=mech.GasMech().GasSpeciesCount(); isp<mech.GasMech().SpeciesCount(); ++isp) {
        head.push_back(mech.GasMech().Species(isp)->Name() + " (mol/cm2)"); // added by mm864
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
        double t = iint->StartTime();
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
        double t = iint->StartTime();
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
        double t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            buildOutputVector(step, t, avg[step], err[step]);
            csv.Write(avg[step]);
        }
    }

    // Close the CSV files.
    csv.Close();
}

// Writes number of jump events to a CSV file.
void Simulator::writePartJumpCSV(const std::string &filename,
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
        head[i] = head[i] + " (-)";
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
        double t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            buildOutputVector(step, t, avg[step], err[step]);
            csv.Write(avg[step]);
        }
    }

    // Close the CSV files.
    csv.Close();
}

// Writes gas-phase and surface phase reaction rates profile to a CSV file. (C)
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

	// separated by mm864

	unsigned int gas_react = 0;

   for (unsigned int i =0; i < mech.GasMech().ReactionCount(); ++i){

	if (mech.GasMech().Reactions(i)->IsSURF() == false) // if it is not a surface 
        { gas_react = gas_react + 1; 
        }
    }

    for (unsigned int i=0; i<gas_react; ++i) {
        head.push_back("Rxn " + cstr(i) + " (mol/cm3s)"); // q will be in cm since readData point converts everything to cgs
    }

   for (unsigned int i=gas_react; i<mech.GasMech().ReactionCount(); ++i) {
        head.push_back("Rxn " + cstr(i) + " (mol/cm2s)");
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
        double t = iint->StartTime();
        for (unsigned int istep=0; istep<(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            buildOutputVector(step, t, avg[step], err[step]);
            csv.Write(avg[step]);
        }
    }

    // Close the CSV files.
    csv.Close();
}

// Writes gas-phase reaction rates profile to a CSV file. // modified by mm864
void Simulator::writeSurfaceProdRatesCSV(const std::string &filename,
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
	/*
	* Separated by mm864
	*/

    for (unsigned int isp=0; isp<mech.GasMech().GasSpeciesCount(); ++isp) {
        head.push_back(mech.GasMech().Species(isp)->Name() + " (mol/cm3s)");
    }

    for (unsigned int isp=mech.GasMech().GasSpeciesCount(); isp<mech.GasMech().SpeciesCount(); ++isp) {
        head.push_back(mech.GasMech().Species(isp)->Name() + " (mol/cm2s)"); 
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
        double t = iint->StartTime();
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
	/*
	* Separated by mm864
	*/

    for (unsigned int isp=0; isp<mech.GasMech().SpeciesCount(); ++isp) {
        head.push_back(mech.GasMech().Species(isp)->Name() + " (mol/cm3s)");
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
        double t = iint->StartTime();
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
        head[i] = head[i] + " (1/m3s)";
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
        double t = iint->StartTime();
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
        double t = iint->StartTime();
        for (unsigned int istep=0; istep!=(*iint).StepCount(); ++istep, ++step) {
            t += iint->StepSize();
            for (unsigned int i=0; i!=track[step].size(); ++i) {
                track[step][i].insert(track[step][i].begin(), t);
                track[step][i].insert(track[step][i].begin(), (double)step);
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

/*!
 * @brief           Writes a binary file containing the ensemble
 *
 * This function is activated with the -ensemble argument when running MOPS.
 * It is to be used to write a binary file containing the ensemble and other
 * useful metadata; which may be re-used in other simulations as an input.
 *
 * The particle model MUST be written to ensure that particles are regenerated
 * with the appropriate deserialisation function. The coagulation kernel is
 * also written, as particles with non-unity weight should not be used in a
 * simulation with DSA coagulation.
 *
 * @param r         Reactor pointer
 * @param step      Step number
 * @param run       Run number
 */
void Simulator::createEnsembleFile(const Reactor &r, unsigned int step,
                                unsigned int run) const
{
    // Build the save point file name.
    string fname = m_output_filename + "(" + cstr(run) + ")-SP(" +
                   cstr(step) + ").ens";

    // Open the save point file.
    ofstream fout;
    fout.open(fname.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

    if (fout.good()) {

        // First write the particle model
        r.Mixture()->ParticleModel()->Serialize(fout);

        // Now write the coagulation process, assuming only one process
        Sweep::Processes::CoagPtrVector coaglist;
        coaglist = r.Mech()->ParticleMech().Coagulations();
        int id(0);
        id = coaglist[0]->ID();     // Get the ID of the coag process
        fout.write((char*)&id, sizeof(id));

        // Now write the ensemble
        r.Mixture()->Particles().Serialize(fout);

        fout.close();
    } else {
        // Throw error if the output file failed to open.
        throw runtime_error("Failed to open file for ensemble "
                            "output (Mops, Simulator::createEnsembleFile).");
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

    // Open output files for all PSL save points.  Remember to
    // write the header row as well.
    vector<CSV_IO*> out(times.size(), NULL);
    for (unsigned int i=0; i!=times.size(); ++i) {
        double t = times[i].EndTime();
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
        for (unsigned int irun=0; irun!=m_nruns; ++irun) {
            // Read the save point for this step and run.
            r = readSavePoint(step, irun, mech);

            if (r != NULL) {
                double scale = (double)m_nruns;
                if (m_output_every_iter) scale *= (double)m_niter;

                // Get PSL for all particles.
                for (unsigned int j=0; j!=r->Mixture()->ParticleCount(); ++j) {
                    // Get PSL.
                    stats.PSL(*(r->Mixture()->Particles().At(j)), mech.ParticleMech(),
                              times[i].EndTime(), psl,
                              1.0/(r->Mixture()->SampleVolume()*scale));
                    // Output particle PSL to CSV file.
                    out[i]->Write(psl);
                }

                // Draw particle images for tracked particles.
                unsigned int n = min(m_ptrack_count,r->Mixture()->ParticleCount());
                for (unsigned int j=0; j!=n; ++j) {
                    double t = times[i].EndTime();
                    string fname = m_output_filename + "-tem(" + cstr(t) +
                                   "s, " + cstr(j) + ").pov";
                    std::ofstream file;
                    file.open(fname.c_str());

                    r->Mixture()->Particles().At(j)->writeParticlePOVRAY(file);

                    file.close();
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

// Processes the PAH-PSLs at each save point into single files.
void Simulator::postProcessPAHPSLs(const Mechanism &mech,
                                const timevector &times) const
{
    Reactor *r = NULL;
    unsigned int step = 0;
    std::vector<std::vector<double> > temp_PAH;

    // Get reference to the particle mechanism.
    const Sweep::Mechanism &pmech = mech.ParticleMech();
    double den = pmech.Components(0)->Density();

        // postProcessXmer is only designed for PAH-PP model
    if (pmech.AggModel() == Sweep::AggModels::PAH_KMC_ID){
        // Build header row for CSV output files.
        vector<string> header;
        vector<string> separator;

        // add the name for the columns
        header.push_back("Index");
        header.push_back("#C");   // #C represent the num of Carbon
        header.push_back("#H");
        header.push_back("#Rings6");
        header.push_back("#Rings5");
        header.push_back("#EdgeC");
        header.push_back("Mass(u)");
        header.push_back("Mass(kg)");
        header.push_back("PAHCollDiameter (m)");
        header.push_back("PAH denbsity (kg/m3)");
        header.push_back("PAH volume (m3)");
        header.push_back("diameter (m)");
        header.push_back("collision diameter (m)");
        header.push_back("time created (s)");
        header.push_back("PAH_ID");


        // Open output files for all PSL save points.  Remember to
        // write the header row as well.
        vector<CSV_IO*> out(times.size(), NULL);
        for (unsigned int i=0; i!=times.size(); ++i) {
            double t = times[i].EndTime();
            out[i] = new CSV_IO();
            out[i]->Open(m_output_filename + "-postprocess-PAH(" +
                        cstr(t) + "s).csv", true);
            out[i]->Write(header);
        }

        // Loop over all time intervals.
        for (unsigned int i=0; i!=times.size(); ++i) {
            // Calculate the total step count after this interval.
            step += times[i].StepCount();

            // Loop over all runs.
            for (unsigned int irun=0; irun!=m_nruns; ++irun) {
                // Read the save point for this step and run.
                r = readSavePoint(step, irun, mech);
                //create separator to distinguish the results of independant runs
                separator.push_back(cstr(irun+1)+"runs");
                out[i]->Write(separator);
                separator.clear();

                if (r != NULL) {
                    double scale = (double)m_nruns;
                    if (m_output_every_iter) scale *= (double)m_niter;

                    // Get PSL for all particles.
                    for (unsigned int j=0; j!=r->Mixture()->ParticleCount(); ++j) {

                        Sweep::Particle* sp=r->Mixture()->Particles().At(j);
                        Sweep::AggModels::PAHPrimary *pah = dynamic_cast<Sweep::AggModels::PAHPrimary*>(sp->Primary());
                        // store the mass of individual PAH within the selected particle in the vector temp_PAH
                        pah->OutputPAHPSL(temp_PAH, j, den);
                        // check whether it is a gasphase PAH, yes set index to -1
                        if (pah->Numprimary()==1 && pah->NumPAH()==1)
                            temp_PAH[0][0]=-1;

                        // Output particle PSL to CSV file.
                        for (size_t ii = 0; ii!=temp_PAH.size(); ++ii) 
                            out[i]->Write(temp_PAH[ii]);
                        // temp_PAH must be cleared before next output
                        temp_PAH.clear();
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
}

template<template <typename> class P = std::greater >
struct compare_pair_second {
    template<class T1, class T2> bool operator()(const std::pair<T1, T2>& left, const std::pair<T1, T2>& right) {
        return P<T2>()(left.second, right.second);
    }
};

void Simulator::postProcessPAHinfo(const Mechanism &mech,
                                const timevector &times) const
{
    Reactor *r = NULL;
    fvector temp;//use to hold all the psl information
    vector<unsigned int> temp_max; 
    fvector temp_PAH_mass;  // store number density of each Xmer
    pair<unsigned int, double> m_index_mass;
    vector<pair<unsigned int, double> > index_mass;
    // Get reference to the particle mechanism.
    const Sweep::Mechanism &pmech = mech.ParticleMech();

    // postProcessXmer is only designed for PAH-PP model
    if (pmech.AggModel() == Sweep::AggModels::PAH_KMC_ID){

        // Create an ensemble stats object.
        Sweep::Stats::EnsembleStats stats(pmech);

        // Open output files for all PSL save points.  Remember to
        // write the header row as well.
        vector<CSV_IO*> out(times.size(), NULL);


        unsigned int step = 0;
        for (unsigned int i=0; i!=times.size(); ++i) {
            double t = times[i].EndTime();
            out[i] = new CSV_IO();
            out[i]->Open(m_output_filename + "-PAH-mass-within-largest-agg(" +
                        cstr(t) + "s).csv", true);
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
                    double scale = (double)m_nruns;
                    if (m_output_every_iter) scale *= (double)m_niter;
                    // Get PSL for all particles.
                    for (unsigned int j = 0; j != r->Mixture()->ParticleCount(); ++j) 
                    {
                        // Get PSL.
                        stats.PSL(*(r->Mixture()->Particles().At(j)), mech.ParticleMech(),
                                    times[i].EndTime(), temp,
                                    1.0/(r->Mixture()->SampleVolume()*scale));
                        //temp[11]=>num of PAH
                        //temp[13]=>num of C, temp[14]=>num of H
                        double mass = 12 * temp[13] + temp[14];
                        // store the index and the mass in a vector of pair
                        m_index_mass = make_pair(j, mass);
                        index_mass.push_back(m_index_mass);

                    }
                    // find the largest soot aggregate by sorting the vecotr
                    std::sort(index_mass.begin(), index_mass.end(), compare_pair_second<std::greater>());
                    // the last element in the temp_max shoule be the index of the largest particle in the ensemble
                    // before finding the largest particle, we have to make sure there are particles in the ensemble.
                    if (index_mass.size()!=0)
                    {
                        // count represents that the number of structures of particles (individual PAHs) will be serialized. 
                        // By default, the 10 largest particle is selected. 
                        // But if the num of particle in the ensemble is less than 10. count is modifed accordingly.
                        unsigned int count(0);
                        if (index_mass.size() > 10)
                            count = 10;
                        else count = index_mass.size();

                        for (size_t pp = 0;pp != count; ++pp) {
                            Sweep::Particle* sp=r->Mixture()->Particles().At(index_mass[pp].first);
                        Sweep::AggModels::PAHPrimary *pah = dynamic_cast<Sweep::AggModels::PAHPrimary*>(sp->Primary());
                        // store the mass of individual PAH within the selected particle in the vector temp_PAH_mass
                        pah->mass_PAH(temp_PAH_mass);
                        // Output particle info to CSV file.
                        out[i]->Write(temp_PAH_mass);
                        temp_PAH_mass.clear();
                        }
                    }
                    // clear the temp varialbe
                    index_mass.clear();
                    //max_mass=0;
                    temp_PAH_mass.clear();
                    temp_max.clear();
                    delete r;
                } else {
                    // Throw error if the reactor was not read.
                    throw runtime_error("Unable to read reactor from save point "
                                        "(Mops, ParticleSolver::postProcessXmer).");
                }
            }
        }
        // Close output CSV files.
        for (unsigned int i=0; i!=times.size(); ++i) {
            out[i]->Close();
            delete out[i];
        }
    }
}

//void Simulator::postProcessXmer(const Mechanism &mech,
//                                const timevector &times) const
//{
//    Reactor *r = NULL;
//    fvector temp;//use to hold all the psl information
//    fvector psl_xmer; // store information of Xmer in ensemble and then dump them to a file
//    fvector m0_xmer;  // store number density of each Xmer
//    fvector m_mass;
//    fvector m_m0;
//    std::vector<std::vector<double> > m_allmass;
//    std::vector<std::vector<double> > m_allm0;
//    // Get reference to the particle mechanism.
//    const Sweep::Mechanism &pmech = mech.ParticleMech();
//
//    // postProcessXmer is only designed for PAH-PP model
//    if (pmech.AggModel() == Sweep::AggModels::PAH_KMC_ID){
//
//        // Create an ensemble stats object.
//        Sweep::Stats::EnsembleStats stats(pmech);
//
//        // Open output files for all PSL save points.  Remember to
//        // write the header row as well.
//        vector<CSV_IO*> out(times.size(), NULL);
//
//        // start with monomer and end with trimer currently
//        for (int k = 0; k!=3 ;++k) {
//            unsigned int step = 0;
//            for (unsigned int i=0; i!=times.size(); ++i) {
//                double t = times[i].EndTime();
//                out[i] = new CSV_IO();
//                out[i]->Open(m_output_filename + "-" + cstr(k+1) +"mer-(" +
//                            cstr(t) + "s).csv", true);
//            }
//
//            // Loop over all time intervals.
//            for (unsigned int i=0; i!=times.size(); ++i) {
//                // Calculate the total step count after this interval.
//                step += times[i].StepCount();
//
//                // Loop over all runs.
//                for (unsigned int irun=0; irun!=m_nruns; ++irun) {
//                    // Read the save point for this step and run.
//                    r = readSavePoint(step, irun, mech);
//
//                    //particle count and number density of ensemble
//                    double Pcount = r->Mixture()->ParticleCount();
//                    double PM0 = r->Mixture()->ParticleCount() / r->Mixture()->SampleVolume();
//
//                    if (r != NULL) {
//                        double scale = (double)m_nruns;
//                        if (m_output_every_iter) scale *= (double)m_niter;
//                        // Get PSL for all particles.
//                        for (unsigned int j = 0; j != Pcount; ++j) {
//                            // Get PSL.
//                            stats.PSL(*(r->Mixture()->Particles().At(j)), mech.ParticleMech(),
//                                        times[i].EndTime(), temp,
//                                        1.0/(r->Mixture()->SampleVolume()*scale));
//                            if ( (k+1) == temp[11]){//temp[11]=>num of PAH
//                                //temp[13]=>num of C, temp[14]=>num of H
//                                double mass = 12 * temp[13] + temp[14];
//                                // currently only mass <= 2000 are interested in mass spectra.
//                                if (mass<=2000) psl_xmer.push_back(mass);
//                            }
//                        }
//                        delete r;
//                        // calculate m0 for each xmer, used for single run
//                        calculateM0(psl_xmer,m0_xmer,Pcount,PM0);
//                    } else {
//                        // Throw error if the reactor was not read.
//                        throw runtime_error("Unable to read reactor from save point "
//                                            "(Mops, ParticleSolver::postProcessXmer).");
//                    }
//                    if (psl_xmer.size()  == m0_xmer.size()) {
//                        m_mass.insert(m_mass.end(),psl_xmer.begin(), psl_xmer.end());
//                        m_allmass.push_back(psl_xmer);
//                        m_allm0.push_back(m0_xmer);
//                        psl_xmer.clear();
//                        m0_xmer.clear();
//                    } else {
//                        // Throw error if the size of vector differ.
//                        throw runtime_error("size of psl_xmer is different form that of m0_xmer, the analysis of xmer information fails"
//                                            "(Mops, ParticleSolver::postProcessXmer).");
//                    }
//
//                }
//                // calculate m0 for xmer, used for multiple runs
//                calculateM0(m_mass, m_m0, m_allmass, m_allm0);
//                // Output particle PSL to CSV file.
//                if (m_mass.size() != m_m0.size()) 
//                    throw runtime_error("size of m_mass is different form that of m_m0, the analysis of xmer information fails" 
//                                        "when trying to calculate mo for multiple runs (Mops, ParticleSolver::postProcessXmer).");
//                out[i]->Write(m_mass);
//                out[i]->Write(m_m0);
//                m_allm0.clear();
//                m_allmass.clear();
//                m_mass.clear();
//                m_m0.clear();
//            }
//
//            // Close output CSV files.
//            for (unsigned int i=0; i!=times.size(); ++i) {
//                out[i]->Close();
//                delete out[i];
//            }
//        }
//    }
//}

void Simulator::postProcessXmer(const Mechanism &mech,
                                const timevector &times) const
{
    Reactor *r = NULL;
    fvector temp;//use to hold all the psl information
    fvector psl_xmer; // store information of Xmer in ensemble and then dump them to a file
    fvector m0_xmer;  // store number density of each Xmer
    fvector m_mass;
    fvector m_m0;
    std::vector<std::vector<double> > m_allmass;
    std::vector<std::vector<double> > m_allm0;
    // Get reference to the particle mechanism.
    const Sweep::Mechanism &pmech = mech.ParticleMech();

    // postProcessXmer is only designed for PAH-PP model
    if (pmech.AggModel() == Sweep::AggModels::PAH_KMC_ID){

        // Create an ensemble stats object.
        Sweep::Stats::EnsembleStats stats(pmech);

        // Open output files for all PSL save points.  Remember to
        // write the header row as well.
        vector<CSV_IO*> out(times.size(), NULL);

        // start with monomer and end with trimer currently
        for (int k = 0; k!=MassSpectraXmer() ;++k) {
            unsigned int step = 0;
            if (MassSpectraEnsemble()){
                for (unsigned int i=0; i!=times.size(); ++i) {
                    double t = times[i].EndTime();
                    out[i] = new CSV_IO();
                    out[i]->Open(m_output_filename + "-mass-spectra-(" +
                                cstr(t) + "s).csv", true);
                }
            }
            else {
                for (unsigned int i=0; i!=times.size(); ++i) {
                    double t = times[i].EndTime();
                    out[i] = new CSV_IO();
                    out[i]->Open(m_output_filename + "-" + cstr(k+1) +"mer-(" +
                                cstr(t) + "s).csv", true);
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

                    //particle count and number density of ensemble
                    double Pcount = r->Mixture()->ParticleCount();
                    double PM0 = r->Mixture()->ParticleCount() / r->Mixture()->SampleVolume();

                    if (r != NULL) {
                        double scale = (double)m_nruns;
                        if (m_output_every_iter) scale *= (double)m_niter;
                        // Get PSL for all particles.
                        for (unsigned int j = 0; j != Pcount; ++j) {
                            // Get PSL.
                            stats.PSL(*(r->Mixture()->Particles().At(j)), mech.ParticleMech(),
                                        times[i].EndTime(), temp,
                                        1.0/(r->Mixture()->SampleVolume()*scale));
                            //temp[11]=>num of PAH, temp[15]=>num of primary particles
                            //temp[13]=>num of C, temp[14]=>num of H
                            double mass = 12 * temp[13] + temp[14];
                            // generate Mass spectra for the whole ensemble
                            if (MassSpectraEnsemble())    psl_xmer.push_back(mass);
                            else {
                                if ((k+1) == temp[11] && 1 == temp[15])    psl_xmer.push_back(mass);
                                else {
                                    Sweep::Particle* sp=r->Mixture()->Particles().At(j);
                                    Sweep::AggModels::PAHPrimary *pah = dynamic_cast<Sweep::AggModels::PAHPrimary*>(sp->Primary());
                                    // if including the xmers in the soot aggregate
                                    //if (MassSpectraFrag())
                                    //pah->FindXmer(psl_xmer,k+1);
                                    if (MassSpectraFrag()&& k<=1 && temp[11]>3 && 1 == temp[15])
                                        pah->Fragtest(psl_xmer, k, pmech.Mode(),pmech.Threshold());
                                }
                            }
                        }
                        delete r;
                        // calculate m0 for each xmer, used for single run
                        calculateM0(psl_xmer,m0_xmer,Pcount,PM0);
                    } else {
                        // Throw error if the reactor was not read.
                        throw runtime_error("Unable to read reactor from save point "
                                            "(Mops, ParticleSolver::postProcessXmer).");
                    }
                    if (psl_xmer.size()  == m0_xmer.size()) {
                        m_mass.insert(m_mass.end(),psl_xmer.begin(), psl_xmer.end());
                        m_allmass.push_back(psl_xmer);
                        m_allm0.push_back(m0_xmer);
                        psl_xmer.clear();
                        m0_xmer.clear();
                    } else {
                        // Throw error if the size of vector differ.
                        throw runtime_error("size of psl_xmer is different form that of m0_xmer, the analysis of xmer information fails"
                                            "(Mops, ParticleSolver::postProcessXmer).");
                    }

                }
                // calculate m0 for xmer, used for multiple runs
                calculateM0(m_mass, m_m0, m_allmass, m_allm0);
                // Output particle PSL to CSV file.
                if (m_mass.size() != m_m0.size()) 
                    throw runtime_error("size of m_mass is different form that of m_m0, the analysis of xmer information fails" 
                                        "when trying to calculate mo for multiple runs (Mops, ParticleSolver::postProcessXmer).");
                out[i]->Write(m_mass);
                out[i]->Write(m_m0);
                m_allm0.clear();
                m_allmass.clear();
                m_mass.clear();
                m_m0.clear();
            }

            // Close output CSV files.
            for (unsigned int i=0; i!=times.size(); ++i) {
                out[i]->Close();
                delete out[i];
            }
        }
    }
}

void Mops::calculateM0(fvector &m_xmer, fvector &m_M0, double Pcount, double PM0) 
{
    fvector temp(m_xmer.begin(), m_xmer.end());

    sort(m_xmer.begin(), m_xmer.end());
    fvector::iterator it;
    it = unique (m_xmer.begin(), m_xmer.end());
    m_xmer.resize(it-m_xmer.begin());
    
    for (unsigned i = 0 ; i != m_xmer.size() ; ++i)
        m_M0.push_back(count(temp.begin(),temp.end(),m_xmer[i])* PM0 / Pcount);
}

void Mops::calculateM0(fvector &m_mass, fvector &m_m0, std::vector<std::vector<double> > &m_allmass, std::vector<std::vector<double> > &m_allm0)
{
    // remove redundant value, only unique ones are left
    sort(m_mass.begin(), m_mass.end());
    fvector::iterator it;
    it = unique (m_mass.begin(), m_mass.end());
    m_mass.resize(it-m_mass.begin());

    for (size_t i = 0; i != m_mass.size(); ++i){
        double sum = 0.0;
        for (size_t j = 0; j !=m_allmass.size(); ++j){
            it = find(m_allmass[j].begin(),m_allmass[j].end(),m_mass[i]);
            size_t m_index = distance(m_allmass[j].begin(),it);
            if (it != m_allmass[j].end()) sum+=m_allm0[j][m_index];
        }
        m_m0.push_back(sum/m_allmass.size());
    }
}

// Writes element fluxes to FluxViewer format.
void Simulator::writeElementFluxOutput(const std::string &filename,
                               const Mechanism &mech,
                               const timevector &times,
                               const std::vector<fvector> &agpfwdrates,
                               const std::vector<fvector> &agprevrates,
                               const std::vector<fvector> &achem) const
{
    if(mech.GasMech().ReactionCount() > 0) {
        Mops::fvector atemperatures;
        for (unsigned int i = 0; i < achem.size(); i++) {
            atemperatures.push_back(achem.at(i).at(achem.at(i).size() - 3));
        }
        FluxAnalyser fa(mech, times, agpfwdrates, agprevrates, atemperatures);
        for (unsigned int i = 0; i < mech.GasMech().ElementCount(); i++) {
            for (unsigned int j = 0; j < m_flux_elements.size(); j++) {
                if (mech.GasMech().Elements(i)->Name().compare(Strings::convertToCaps(m_flux_elements.at(j))) == 0) {
                    fa.addElement(*mech.GasMech().Elements(i));
                }
            }
        }
        fa.writeFluxes(filename, true);
    }
}

// Add element for flux analysis postprocessor
void Simulator::AddFluxElement(const std::string &elem_str)
{
    m_flux_elements.push_back(elem_str);
}

// Clear element list for flux analysis postprocessor
void Simulator::ClearFluxElements()
{
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

		#ifdef USE_MPI							//ms785
		int numprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		n = numprocs;
		#endif

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

        // Flag controlling jump file output.
        if (m_write_jumps) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

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
                m_maxm0 = (double)val;

                // Computation time.
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_cpu_start = (clock_t)i;
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_cpu_mark = (clock_t)i;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_runtime = (double)val;

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

                // Flag controlling jump file output.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_write_jumps = (n==trueval);

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
