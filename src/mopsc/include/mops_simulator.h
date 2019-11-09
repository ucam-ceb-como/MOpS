/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Simulator class drives a single reactor simulation and
    takes care of the post-processing.  It requires a Reactor
    object and a Solver object as arguments.

    This is the default simulator for mops.  Simulators for more
    complex systems are under consideration.

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

#ifndef MOPS_SIMULATOR_H
#define MOPS_SIMULATOR_H

#include "mops_params.h"
#include "mops_reactor.h"
#include "mops_psr.h"
#include "mops_timeinterval.h"
#include "mops_solver.h"
#include "mops_mechanism.h"
#include "console_io.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>

namespace Mops
{

// Forward declare the network simulator
class NetworkSimulator;

class Simulator
{
public:
    // Constructors.
    Simulator(void); // Default constructor.

    // Destructors.
    ~Simulator(void); // Default destructor.

    // Copy constructor.
    Simulator(const Mops::Simulator &copy);

    Simulator &operator=(const Simulator &rhs);

    // The NetworkSimulator needs full access to the Simulator
    friend class Mops::NetworkSimulator;

    // SIMULATION SETTINGS.

    // Returns the number of runs to perform
    unsigned int RunCount(void) const;

    // Sets the number of runs to perform.
    void SetRunCount(unsigned int n);

    // Return time vector.
    const timevector &TimeVector() const;

    // Returns the number of iteration to perform per step.
    unsigned int IterCount(void) const;

    // Sets the number of iterations to perform per step.
    void SetIterCount(unsigned int n);

    // Sets the time vector.
    void SetTimeVector(const timevector &times);

    // Return the number of time steps.
    unsigned int TimeStepCount() const;

    // Particle counts should be properties of the Reactor class!

    // Returns the max. stochastic particle count.
    unsigned int MaxPartCount(void) const;

    // Sets the max. stochastic particle count.
    void SetMaxPartCount(unsigned int n);

    //! Returns the max. M0, for initial ensemble scaling.
    double MaxM0(void) const;

    //! Sets max. M0.
    void SetMaxM0(double m0);

    // CONSOLE OUTPUT.

    // Returns the console output interval (in # of steps).
    unsigned int ConsoleInterval() const;

    // Sets the console output interval (# of steps).
    void SetConsoleInterval(unsigned int cint);

    // Returns the vector of console output names.
    const std::vector<std::string> &ConsoleVariables() const;

    // Returns the ith console output name.  Returns "" if i is invalid
    const std::string ConsoleVariable(unsigned int i) const;

    // Adds a variable to the console output.
    void AddConsoleVariable(const std::string &var);

    // Removes the console output variable with the given name.
    void RemoveConsoleVariable(const std::string &var);

    // Removes the ith console output variable.
    void RemoveConsoleVariable(unsigned int i);

    // Returns whether or not messages will be printed to the console.
    bool UseConsoleMsgs() const;

    // Sets whether or not messages are printed to the console.
    void SetUseConsoleMsgs(bool msgs);


    // FILE OUTPUT.

    // Returns the output file name.
    const std::string &OutputFile() const;

    // Sets the output file name.
    void SetOutputFile(const std::string &name);

    // Set simulator to output every iteration.
    void SetOutputEveryIter(bool fout);
    
    //! Set simulator to write the jumps CSV file.
    void SetWriteJumpFile(bool writejumps);

    //! Set the simulator to write the particle binary file
    void SetWriteEnsembleFile(bool writeparticles);

    //! Set simulator to write the PAHs to psl files.
    void SetWritePAH(bool postpocessPAH);

	//! Set simulator to write the PAHs to psl files.
	void SetWritePP(bool postpocessPP);

    // options for Postprocess (only for PAH-PP model)
    //! return the option whehter generate mass spectra
    const bool MassSpectra() const;
    //! set the option whehter generate mass spectra
    void SetMassSpectra(const bool val);

    //! return the option whehter generate mass spectra for the whole ensemble
    const bool MassSpectraEnsemble() const;
    //! set the option whehter generate mass spectra for the whole ensemble
    void SetMassSpectraEnsemble(const bool val);

    //! return the option whehter generate mass spectra for specified gasphase xmer
    const int MassSpectraXmer() const;
    //! set the option whehter generate mass spectra for specified gasphase xmer
    void SetMassSpectraXmer(const int val);

    //! return the option whehter including the xmer in the large soot aggregate
    const bool MassSpectraFrag() const;
    //! set the option whehter including the xmer in the large soot aggregate
    void SetMassSpectraFrag(const bool val);

    // STATISTICAL BOUNDS OUTPUT

    // Set simulator to output data of a given statistical range.
    void SetOutputStatBoundary(
        const Sweep::PropID pid,
        double lower,
        double upper
        );

    // POVRAY OUTPUT.

    // Set number of particle trackings.
    void SetParticleTrackCount(unsigned int ptcount);

	// PARTICLE TRACKING FOR VIDEOS  (alternative to POVRAY)

	//! Set number of track particles 
	void SetTrackBintreeParticleCount(unsigned int val);

	//! Return the (max) number of tracked particles 
	const unsigned int TrackBintreeParticleCount() const;

    // SOLUTION AND POST-PROCESSING.

    // Run the solver for the given reactor and the 
    // given time intervals.
    void RunSimulation(
        Reactor &r,              // Reactor object to solve.
        //const timevector &times, // Vector of time intervals.
        Solver &s,               // Solver to use for simulation.
        size_t seed);

    //! Post-processes binary output files into CSV files.
    void PostProcess(void);

    // Add element for flux analysis postprocessor
    void AddFluxElement(const std::string &elem_str);

    // Clear element list for flux analysis postprocessor
    void ClearFluxElements();

    // READ/WRITE/COPY FUNCTIONS.

    // Writes the simulator to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the simulator data from a binary data stream.
    void Deserialize(std::istream &in);

private:
    // SIMULATION SETTINGS.

    // Number of runs to perform.
    unsigned int m_nruns;

    // Number of internal solver iterations to perform.
    unsigned int m_niter;

    // Simulation time vector.
    timevector m_times;

    // These should be properties of the Reactor class!

    // Max. number of stochastic particles in sweep.
    unsigned int m_pcount;

    //! Max. M0 value, for initial scaling of ensemble.
    double m_maxm0;

    // COMPUTATION TIME.

    clock_t m_cpu_start, m_cpu_mark;
    double m_runtime;


    // CONSOLE OUTPUT PARAMETERS.

    // Interval of console output data (in terms of time steps).
    unsigned int m_console_interval;
    
    // Console column variable names.
    std::vector<std::string> m_console_vars;

    // Console variable mask.
    std::vector<unsigned int> m_console_mask;

    // Set to true if the console is to print messages.
    bool m_console_msgs; 

    // Console output class.
    Console_IO m_console;


    // FILE OUTPUT PARAMETERS.

    // Name of output file.
    std::string m_output_filename;

    // Simulation output file stream.
    mutable std::fstream m_file;

    // Sensitivity output file stream.
    mutable std::fstream m_senfile;

    //! Output stream for the LOI file data
    mutable std::ofstream m_loi_file;

    //! Object to hold LOI data
    std::vector<Mops::fvector> m_loi_data;

    //! The Jacobian for LOI
    double** m_loi_J;

    // Flag controlling iteration output.  If true then output
    // is performed for every iteration at the end of a time step,
    // otherwise output is only performed after all iterations have
    // completed.
    bool m_output_every_iter;

    // Sub-step number at which output will occur.
    unsigned int m_output_step;

    // Sub-step iteration number at which output will occur, 
    // unless m_output_every_iter==true, in which case all
    // iterations are output.
    unsigned int m_output_iter;

    // Flag controlling whether the number of jump events should
    // be written to CSV output. Default false.
    bool m_write_jumps;

    //! Should the ensemble be written to a reusable binary file?
    bool m_write_ensemble_file;

    //! Flag controlling whether post-process the detailed info about every PAH in the particle ensemble. Default false.
    bool m_write_PAH;

	//! Flag controlling whether post-process the detailed info about every primary particle in the particle ensemble. Default false.
	bool m_write_PP;

    //  the method of postprocessing, only for PAH-PP model 
    //!Flag controlling the generation of Mass spectra
    bool m_mass_spectra;

    //! Flag specifying whether the user is interested in the whole ensemble or just the gasphase Xmer.
    bool m_mass_spectra_ensemble;

    //! specifying the largest xmer
    int m_mass_spectra_xmer;

    //! Flag controlling whether considering xmer in the large soot aggregate
    bool m_mass_spectra_frag;

    // STATISTICAL OUTPUT PARAMETERS

    Sweep::Stats::IModelStats::StatBound m_statbound;


    // PARTICLE TRACKING PARAMETERS.

    // Number of particles for which to produce tracking output.  Tracked
    // particles will have their PSL properties output at each time-step in
    // a separate CSV file.  If a primary-particle model is implemented, then
    // TEM-style images will also be generated.
    unsigned int m_ptrack_count;

    // FLUX ANALYSIS ELEMENT LIST.

    // A vector list containing pointers to elements which are wanted to analyse the flux of
    // that element.
    std::vector<std::string> m_flux_elements;

    // LOI THINGS

    //! Set up the LOI calculation
    void setupLOI(const Mops::Reactor &r, Mops::Solver &s);

    //! Solve the LOI Jacobian
    void solveLOIJacobian(
            Mops::Reactor &r,
            Mops::Solver &s,
            const unsigned int istep,
            const double t2);


    // FILE OUTPUT.

    // Opens an output file for the given run number.
    void openOutputFile() const;

    // Closes the output file.
    void closeOutputFile() const;

    // Writes the current reactor state to the output file.  This is
    // of the standard output function form for the Mops::Solver class.
    // It is passed to the Solver classes Solve() function.
    static void fileOutput(
        unsigned int step,  // Internal solver sub-step number for this time-step.
        unsigned int iter,  // Internal solver iteration number for this sub-step.
        const Reactor &r,   // Reactor to output.
        const Solver &s,    // Solver used for calculations.
        void *sim           // Void pointer to cast to this simulator object.
        );

    // Writes the gas-phase conditions of the given reactor to
    // the binary output file.
    void outputGasPhase(const Reactor &r) const;

    // Writes the particle-number count of the given 
    // reactor to the binary output file. 
    void outputParticleNumber(const Reactor &r) const;

    // Writes the particle stats to the binary output file.
    void outputParticleStats(const Reactor &r) const;

    // Writes tracked particles to the binary output file.
    void outputPartTrack(const Reactor &r) const;
    
    // Write sensitivity output to the binary file.
    void outputSensitivity(const Reactor &r) const;

    // Writes the gas-phase reaction rates-of-progress and
    // the species molar production rates due to gas-phase
    // reactions to the binary output file.
    void outputGasRxnRates(const Reactor &r) const;

    // Writes the particle process rates and the molar production
    // rates for each species due to particle processes for the
    // given reactor to the binary output file.
    void outputPartRxnRates(const Reactor &r) const;

	// PARTICLE TRACKING OUTPUT FOR VIDEOS

	// Vector of file names
	std::vector<std::string> m_TrackParticlesName;

	//! (Maximum) number of particles tracked
	unsigned int m_track_bintree_particle_count;

    // CONSOLE OUTPUT.

    // Sets up console output using the given mechanism as a template.
    void setupConsole(const Mechanism &mech);

    // Writes current reactor state to the console.
    void consoleOutput(const Reactor &r) const;

    // POST-PROCESSING ROUTINES.

    // Writes the auxilliary post-processing information using
    // the given file name.  This information is the chemical mechanism
    // and the output time intervals.
    void writeAux(
        const Mops::Mechanism &mech,   // Mechanism which defines the reactor.
        const Mops::timevector &times, // Vector of time intervals.
        const Mops::Solver &solv       // Solver used to perform simulation.
        ) const;

    // Reads auxilliary post-processing information using the
    // given file name.  This information is the chemical mechanism
    // and the output time intervals.
    void readAux(
        Mops::Mechanism &mech,   // Mechanism which defined the reactor.
        Mops::timevector &times, // Vector of time intervals.
        unsigned int &ncput,     // Number of CPU times generated by solver.
        std::vector<std::string> &cput_head // CPU times headings.
        );

    // OUTPUT POINT READING.

    // Reads a gas-phase chemistry data point from the binary file.
    // To allow the averages and confidence intervals to be calculated -
    // the data point is added to a vector of sums, and the squares are
    // added to the vector sumsqr if necessary.
    static void readGasPhaseDataPoint(
        std::istream &in,            // Input stream.
        const Mops::Mechanism &mech, // Chemical mechanism.
        fvector &sum,                // Sums of chemistry data.
        fvector &sumsqr,             // Sums of the squares.
        bool calcsqrs = false        // Set =true to also calculate sums of squares.
        );

    // Reads a particle-number list data point from the binary file.
    // To allow the averages and confidence intervals to be calculated -
    // the data point is added to a vector of sums, and the squares are
    // added to the vector sumsqr if necessary.
    static void readParticleNumberListDataPoint(
	std::istream &in,            // Input stream.
	const Mops::Mechanism &mech, // Chemical mechanism.
	fvector &sum,                // Sums of chemistry data.
	fvector &sumsqr,             // Sums of the squares.
	bool calcsqrs = false        // Set =true to also calculate sums of squares.
	);

    // Reads a CPU timing data from the binary file.
    // To allow the averages and confidence intervals to be calculated
    // the data point is added to a vector of sums, and the squares are
    // added to the vector sumsqr if necessary.
    static void readCTDataPoint(
        std::istream &in,     // Input stream.
        unsigned int N,       // Number of CT points that were written to the file.
        fvector &sum,         // Sums of CPU timing data.
        fvector &sumsqr,      // Sums of the squares.
        bool calcsqrs = false // Set =true to also calculate sums of squares.
        );

    // Reads a particle stats data point from the binary file.
    // To allow the averages and confidence intervals to be calculated
    // the data point is added to a vector of sums, and the squares are
    // added to the vector sumsqr if necessary.
    static void readParticleDataPoint(
        std::istream &in,             // Input stream.
        const Sweep::Mechanism &mech, // Chemical mechanism.
        fvector &sum,                 // Sums of particle stats data.
        fvector &sumsqr,              // Sums of the squares.
        bool calcsqrs = false         // Set =true to also calculate sums of squares.
        );

    // Reads a gas-phase reaction rates stats data point from the binary file.
    // To allow the averages and confidence intervals to be calculated
    // the data point is added to a vector of sums, and the squares are
    // added to the vector sumsqr if necessary.
    static void readGasRxnDataPoint(
        std::istream &in,             // Input stream.
        const Mops::Mechanism &mech,  // Particle mechanism.
        fvector &rates_sum,           // Sums of reaction rates.
        fvector &rates_sumsqr,        // Sums of the squares of reaction rates.
        fvector &fwd_rates_sum,       // Sums of forward reaction rates.
        fvector &fwd_rates_sumsqr,    // Sums of the squares of forward reaction rates.
        fvector &rev_rates_sum,       // Sums of reverse reaction rates.
        fvector &rev_rates_sumsqr,    // Sums of the squares of reverse reaction rates.
        fvector &wdot_sum,            // Sums of species prod. rates.
        fvector &wdot_sumsqr,         // Sums of the squares of species prod. rates.
	fvector &sdot_sum,            // Sums of species prod. rates from surface .
	fvector &sdot_sumsqr,         // Sums of the squares of species prod. rates from surface.
        bool calcsqrs = false         // Set =true to also calculate sums of squares.
        );

    // Reads a particle process rates stats data point from the binary file.
    // To allow the averages and confidence intervals to be calculated
    // the data point is added to a vector of sums, and the squares are
    // added to the vector sumsqr if necessary.
    static void readPartRxnDataPoint(
        std::istream &in,             // Input stream.
        const Sweep::Mechanism &mech, // Particle mechanism.
        fvector &rates_sum,           // Sums of process rates.
        fvector &rates_sumsqr,        // Sums of the squares of process rates.
        fvector &wdot_sum,            // Sums of species prod. rates.
        fvector &wdot_sumsqr,         // Sums of the squares of species prod. rates.
        fvector &jumps_sum,           // Sums of the number of jumps
        fvector &jumps_sumsqr,        // Sums of the squares of the number of jumps
        bool calcsqrs = false         // Set =true to also calculate sums of squares.
        );

    // Reads the tracked particles from the binary file.  The particles are
    // processed so that only a vector of vectors is returned, which contains
    // the PSL data for each tracked particle at that point.
    void readPartTrackPoint(
        std::istream &in,             // Input stream.
        const Sweep::Mechanism &mech, // Particle mechanism.
        std::vector<fvector> &pdata   // Tracked particle output data for this point.
        ) const;

    // OUTPUT VECTOR MANIPULATION.

    // Multiplies all values in a vector by a scaling factor.
    static void multVals(
        fvector &vals, // The values to multiply by the scaling factor.
        double scale     // The scaling factor (numner of runs).
        );

    // Takes vectors of vectors of variable sums and sums of squares, which
    // are converted into the average values and the confidence intervals.
    static void calcAvgConf(
        std::vector<fvector> &avg, // Input=sums, output=averages.
        std::vector<fvector> &err, // Input=sums of squares, output=confidence intervals.
        unsigned int nruns         // Number of runs.
        );

    // Takes a vector of average values and a vector with the confidence
    // bound of each variable and combines them into a single vector:
    // BEFORE:
    // avg = (a1, a2, a3, ..., aN)
    // err = (e1, e2, e3, ..., eN)
    // AFTER:
    // avg = (N, t, a1, e1, a2, e2, a3, e3, ..., aN, eN)
    // The step number and time are insert at the beginning of the avg
    // vector.
    static void buildOutputVector(
        unsigned int step, // Step number.
        double time,         // Step time.
        fvector &avg,      // Averages.
        const fvector &err // Confidence intervals (will be inserted into avg).
        );

    // CSV FILE CREATION.

    // Writes gas-phase conditions profile to a CSV file.
    static void writeGasPhaseCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining chemical species.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of gas-phase time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );

    // Writes the particle-number list count to a CSV file.
    static void writeParticleNumberListCSV(
	const std::string &filename,     // Output file name (incl. extension).
	const Mechanism &mech,           // Mechanism defining chemical species.
	const timevector &times,         // Output time profile.
	std::vector<fvector> &avg,       // Vector of gas-phase time points.
	const std::vector<fvector> &err  // Vector of confidence intervals.
	);

    // Writes computation times profile to a CSV file.
    void writeCT_CSV(
        const std::string &filename,     // Output file name (incl. extension).
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of computation-time time points.
        const std::vector<fvector> &err, // Vector of confidence intervals.
        const std::vector<std::string> &head // CSV column names.
        ) const;

    // Writes particle stats profile to a CSV file.
    static void writeParticleStatsCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining particle ensemble.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of gas-phase time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );

    // Writes particle jump rate profile to a CSV file.
    static void writePartJumpCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Sweep::Mechanism &mech,    // Mechanism defining particle processes.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of process rate time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );

    // Writes gas-phase reaction rates profile to a CSV file.
    static void writeGasRxnCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining kinetics.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of gas-phase reaction time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );

 	// Writes molar prod. rates profile to a CSV file.
    static void writeSurfaceProdRatesCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining kinetics.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of molar prod. rate time points for surface.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );


    // Writes molar prod. rates profile to a CSV file.
    static void writeProdRatesCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining kinetics.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of molar prod. rate time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );

    // Writes particle process rates profile to a CSV file.
    static void writePartProcCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Sweep::Mechanism &mech,    // Mechanism defining particle processes.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of process rate time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );

    // Writes particle tracking for multiple particles to CSV files.
    static void writePartTrackCSV(
        const std::string &filename, // Output file name (excl. extension).
        const Mechanism &mech,       // Mechanism defining particle ensemble.
        const timevector &times,     // Output time profile.
        std::vector<std::vector<fvector> > &track // Vector of tracking data of multiple particles at multiple time points.
        );


    // SAVE POINTS AND PSL POST-PROCESSING.

    // Creates a simulation save point.  The save points can be
    // used to restart an interrupted simulation.
    void createSavePoint(
        const Reactor &r,  // Reactor to output.
        unsigned int step, // Step number.
        unsigned int run   // Run number.
        ) const;

    //! Creates an ensemble binary file.
    void createEnsembleFile(
        const Reactor &r,  // Reactor to output.
        unsigned int step, // Step number.
        unsigned int run   // Run number.
        ) const;


    // Reads a save point file.
    Reactor *const readSavePoint(
        unsigned int step,    // Step number.
        unsigned int run,     // Run number.
        const Mechanism &mech // Mechanism used to define reactor.
        ) const;

    //! Do the grunt work for postprocessing a binary file
    void postProcessSimulation(
        Mechanism &mech,
        timevector &times,
        unsigned int npoints,
        unsigned int ncput,
        std::vector<std::string> cput_head
        ) const;

    // Processes the PSLs at each save point into single files.
    void postProcessPSLs(
        const Mechanism &mech,  // Mechanism use to solve system.
        const timevector &times // Simulation output time intervals.
        ) const;

	// Processes the particle-number list PSLs at each save point into single files.
	void postProcessParticleNumberPSLs(
		const Mechanism &mech,  // Mechanism use to solve system.
		const timevector &times // Simulation output time intervals.
		) const;

    // Processes the PSLs at each save point into single files for PAHs.
    void postProcessPAHPSLs(
        const Mechanism &mech,  // Mechanism use to solve system.
        const timevector &times // Simulation output time intervals.
        ) const;

	// Processes the PSLs at each save point into single files for Primary Particles.
	void postProcessPPPSLs(
		const Mechanism &mech,  // Mechanism use to solve system.
		const timevector &times // Simulation output time intervals.
		) const;

    // post-process the ensemble to find interested information, in this case, mass of Xmer at the end of simulation
    void postProcessXmer(
        const Mechanism &mech,  // Mechanism use to solve system.
        const timevector &times // Simulation output time intervals.
        ) const;

    // post-process the ensemble to find interested information, in this case, info about PAH mass distribution of individual soot aggregate at certain time.
    void postProcessPAHinfo(
        const Mechanism &mech,  // Mechanism use to solve system.
        const timevector &times // Simulation output time intervals.
        ) const;
    // FLUX VIEWER OUTPUT
    // Writes element fluxes to FluxViewer format.
    void writeElementFluxOutput(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining kinetics.
        const timevector &times,         // Output time profile.
        const std::vector<fvector> &agpfwdrates,       // Vector of gas-phase reaction time points.
        const std::vector<fvector> &agprevrates,       // Vector of gas-phase reaction time points.
        const std::vector<fvector> &achem  // Vector of confidence intervals.
        ) const;

    // COMPUTATION TIME CALCULATION.

    // Calculates the time duration from a time mark to the
    // current time.
    double calcDeltaCT(double markt) const;
};
    // calculate M0 for one runs 
    void calculateM0(
        fvector &m_xmer,
        fvector &m_M0,
        double Pcount,
        double PM0
        );
    // calculate M0 for mutiple runs
    void calculateM0(
        fvector &m_mass,
        fvector &m_m0,
        std::vector<std::vector<double> > &m_allmass,
        std::vector<std::vector<double> > &m_allm0
        );

};

#endif
