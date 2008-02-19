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
    ~Solver(void); // Default destructor.


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
    std::fstream m_file;

private:

    // CONSOLE OUTPUT.

    // Sets up console output using the given mechanism as a template.
    void setupConsole(const Mechanism &mech);

    // Writes current reactor state to the console.
    void consoleOutput(const Reactor &r) const;


    // FILE OUTPUT.

    // Sets up the file output by outputting an auxilliary file
    // which stores all the information required to post-process the
    // simulation and by opening the output file.
    void beginFileOutput(
        const Mechanism &mech,  // Mechanism which defines the output.
        const timevector &times // Vector of time intervals.
        );

    // Writes the current reactor state to the output file.
    void fileOutput(const Reactor &r);

    // Ends file output by closing all open files.
    void endFileOutput();


    // POST-PROCESSING ROUTINES.

    void buildOutputVector(
        const Reactor &r, // Reactor object holding data to be output.
        fvector &out      // Vector to build.
        ) const;
};
};

#endif
