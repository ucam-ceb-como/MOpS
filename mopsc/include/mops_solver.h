/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The Solver class holds simulation settings for mops and solves reactors.
    This basic solver only solves gas-phase chemistry equations, there is no
    operator splitting to solve gas-phase chemistry coupled to the particle
    population balance using sweep.
*/

#ifndef MOPS_SOLVER_H
#define MOPS_SOLVER_H

#include "mops_params.h"
#include "mops_reactor.h"
#include "mops_timeinterval.h"
#include "mops_mechanism.h"
#include "mops_ode_solver.h"
#include "console_io.h"
#include <vector>
#include <string>
#include <fstream>
#include <time.h>

namespace Mops
{
class Solver
{
public:
    // Constructors.
    Solver(void); // Default constructor.

    // Destructors.
    virtual ~Solver(void); // Default destructor.


    // ERROR TOLERANCES FOR ODE SOLVER.

    // Returns the absolute error tolerance used for ODE
    // calculations.
    real ATOL() const;

    // Sets the absolute error tolerance used for ODE
    // calculations.
    void SetATOL(real atol);

    // Returns the relative error tolerance used for ODE
    // calculations.
    real RTOL() const;

    // Sets the relative error tolerance used for ODE
    // calculations.
    void SetRTOL(real rtol);


    // SWEEP SETTINGS.

    // Returns the number of runs to perform
    unsigned int RunCount(void) const;

    // Sets the number of runs to peform.
    void SetRunCount(unsigned int n);

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


    // SOLUTION AND POST-PROCESSING.

    // Run the solver for the given reactor and the 
    // given time intervals.
    virtual void SolveReactor(
        Reactor &r,              // Reactor object to solve.
        const timevector &times, // Vector of time intervals.
        unsigned int nruns = 1   // Number of runs to perform.
        );

    // Post-processes binary output files with the given file name
    // into CSV files.
    virtual void PostProcess(
        const std::string &filename, // Filename to post-process.
        unsigned int nruns = 1       // Number of runs.
        ) const;

protected:
    // ODE SOLVER.

    ODE_Solver m_ode; // The ODE solver used by the mops solver.


    // SOLVER SETTINGS.

    // Default error tolerances for the ODE solver.
    real m_atol, m_rtol;

    // Number of runs to perform.
    unsigned int m_nruns;

    // Max. number of stochastic particles in sweep.
    unsigned int m_pcount;

    // Max. M0 value, for initial scaling of ensemble.
    real m_maxm0;


    // COMPUTATION TIME.

    clock_t m_cpu_start, m_cpu_mark;
    double m_chemtime;


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


    // FILE OUTPUT.

    // Opens an output file for the given run number.
    void openOutputFile(unsigned int run) const;

    // Closes the output file.
    void closeOutputFile() const;

    // Writes the gas-phase conditions of the given reactor to
    // the binary output file.
    void outputGasPhase(const Reactor &r) const;

    
    // CONSOLE OUTPUT.

    // Sets up console output using the given mechanism as a template.
    void setupConsole(const Mechanism &mech);

    // Writes current reactor state to the console.
    void consoleOutput(const Reactor &r) const;


    // POST-PROCESSING ROUTINES.

    // Writes the auxilliary post-processing information using
    // the given file name.  This information is the chemical mechanism
    // and the output time intervals.
    static void writeAux(
        const std::string &filename,  // Root filename.
        const Mops::Mechanism &mech,  // Mechanism which defines the reactor.
        const Mops::timevector &times // Vector of time intervals.
        );

    // Reads auxilliary post-processing information using the
    // given file name.  This information is the chemical mechanism
    // and the output time intervals.
    static void readAux(
        const std::string &filename, // Root filename.
        Mops::Mechanism &mech,       // Mechanism which defined the reactor.
        Mops::timevector &times      // Vector of time intervals.
        );

    // Reads a gas-phase chemistry data point from the binary file.
    // To allow the averages and confidence intervals to be calculated
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
    // avg = (a1, e1, a2, e2, a3, e3, ..., aN, eN)
    // The step number and time are insert at the beginning of the avg
    // vector.
    static void buildOutputVector(
        unsigned int step, // Step number.
        real time,         // Step time.
        fvector &avg,      // Averages.
        const fvector &err // Confidence intervals (will be inserted into avg).
        );

    // Writes gas-phase conditions profile to a CSV file.
    static void writeGasPhaseCSV(
        const std::string &filename,     // Output file name (incl. extension).
        const Mechanism &mech,           // Mechanism defining chemical species.
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of gas-phase time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        );

    // Writes computation times profile to a CSV file.
    virtual void writeCT_CSV(
        const std::string &filename,     // Output file name (incl. extension).
        const timevector &times,         // Output time profile.
        std::vector<fvector> &avg,       // Vector of computation-time time points.
        const std::vector<fvector> &err  // Vector of confidence intervals.
        ) const;


    // SAVE POINTS.

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


    // COMPUTATION TIME CALCULATION.

    // Calculates the time duration from a time mark to the
    // current time.
    double calcDeltaCT(double markt) const;

private:
    // FILE OUTPUT.

    // Writes the current reactor state to the output file.
    void fileOutput(const Reactor &r) const;
};
};

#endif
