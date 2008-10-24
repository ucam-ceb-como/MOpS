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
class Simulator
{
public:
    // Constructors.
    Simulator(void); // Default constructor.

    // Destructors.
    ~Simulator(void); // Default destructor.


    // SIMULATION SETTINGS.

    // Returns the number of runs to perform
    unsigned int RunCount(void) const;

    // Sets the number of runs to perform.
    void SetRunCount(unsigned int n);

    // Returns the number of iteration to perform per step.
    unsigned int IterCount(void) const;

    // Sets the number of iterations to perform per step.
    void SetIterCount(unsigned int n);


    // Particle counts should be properties of the Reactor class!

    // Returns the max. stochastic particle count.
    unsigned int MaxPartCount(void) const;

    // Sets the max. stochastic particle count.
    void SetMaxPartCount(unsigned int n);

    // Returns the max. M0, for initial ensemble scaling.
    real MaxM0(void) const;

    // Sets max. M0.
    void SetMaxM0(real m0);

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
    
    // STATISTICAL BOUNDS OUTPUT

    // Set simulator to output data of a given statistical range.
    void SetOutputStatBoundary(
        const Sweep::ParticleCache::PropID pid,
        real lower,
        real upper
        );

    // POVRAY OUTPUT.

    // Set number of particle trackings.
    void SetParticleTrackCount(unsigned int ptcount);

    // SOLUTION AND POST-PROCESSING.

    // Run the solver for the given reactor and the 
    // given time intervals.
    void RunSimulation(
        Reactor &r,              // Reactor object to solve.
        const timevector &times, // Vector of time intervals.
        Solver &s                // Solver to use for simulation.
        );

    // Post-processes binary output files with the given file name
    // into CSV files.
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

    // These should be properties of the Reactor class!

    // Max. number of stochastic particles in sweep.
    unsigned int m_pcount;

    // Max. M0 value, for initial scaling of ensemble.
    real m_maxm0;


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

    // Output file stream.
    mutable std::fstream m_file;

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

    // Writes the particle stats to the binary output file.
    void outputParticleStats(const Reactor &r) const;

    // Writes tracked particles to the binary output file.
    void outputPartTrack(const Reactor &r) const;
    
    // Writes the gas-phase reaction rates-of-progress and
    // the species molar production rates due to gas-phase
    // reactions to the binary output file.
    void outputGasRxnRates(const Reactor &r) const;

    // Writes the particle process rates and the molar production
    // rates for each species due to particle processes for the
    // given reactor to the binary output file.
    void outputPartRxnRates(const Reactor &r) const;


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
        real scale     // The scaling factor (numner of runs).
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
        real time,         // Step time.
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

    // Writes gas-phase reaction rates profile to a CSV file.
    static void writeGasRxnCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining kinetics.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of gas-phase reaction time points.
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

    // Reads a save point file.
    Reactor *const readSavePoint(
        unsigned int step,    // Step number.
        unsigned int run,     // Run number.
        const Mechanism &mech // Mechanism used to define reactor.
        ) const;

    // Processes the PSLs at each save point into single files.
    void postProcessPSLs(
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
        );


    // COMPUTATION TIME CALCULATION.

    // Calculates the time duration from a time mark to the
    // current time.
    double calcDeltaCT(double markt) const;
};
};

#endif
