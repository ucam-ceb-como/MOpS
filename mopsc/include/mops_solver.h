/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Solver class holds simulation settings for mops and solves reactors.
    This basic solver only solves gas-phase chemistry equations, there is no
    operator splitting to solve gas-phase chemistry coupled to the particle
    population balance using sweep.

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
#include <iostream>
#include <time.h>

namespace Mops
{
class Solver
{
public:
    // A function pointer type definition for the solver 
    // output routine.
    typedef void (*OutFnPtr)(
        unsigned int,   // Current internal step number since last Solve() call.
        unsigned int,   // Iteration number for current step.
        const Reactor&, // Reactor being solved.
        const Solver&,  // Reference to current solver.
        void*           // User data object.
        );

    // Constructors.
    Solver(void); // Default constructor.

    // Destructors.
    virtual ~Solver(void); // Default destructor.


    // SOLVER INITIALISATION AND RESET.

    // Initialises the solver to solve the given reactor.
    virtual void Initialise(Reactor &r);

    // Resets the solver to solve the given reactor.
    virtual void Reset(Reactor &r);


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

    // LOI STATUS FOR ODE SOLVER.

    //! Enables LOI status to true.
    void SetLOIStatusTrue();

    //! Sets LOI status to false
    void SetLOIStatusFalse();

    //! Retrieves the LOI status
    bool GetLOIStatus() const;

    //! Sets the cutoff value for LOI comparison
    void SetLOICompValue(double CompValue);

    //! Returns the LOI comparison value
    double ReturnCompValue() const;
    
    // UNDER-RELAXATION.

    // Returns the under-relaxation coefficient.
    real UnderRelaxCoeff(void) const;

    // Sets the under-relaxation coefficient.
    void SetUnderRelaxCoeff(real relax);

    // SOLUTION.

    // Runs the solver for the given reactor, advancing it
    // to the given stop time.  The numerical parameters given
    // are the number of internal steps to take, and the number
    // of internal iterations.  Default values of <=0 will use
    // an adaptive method (NOT YET IMPLEMENTED).  Internal solver
    // output is provided after each step/iteration by passing
    // a function pointer.
    virtual void Solve(
            Reactor &r,   // The reactor to solve.
            real tstop,   // The end time for the step.
            int nsteps,   // Number of internal steps to take.
            int niter,    // Number of internal iterations to take.
            OutFnPtr out, // Output function pointer.
            void *data    // Custom data object which will be passed as argument to out().
        );
    
    //SENSITIVITY.

    //!Returns the double** m_sensitivity matrix with values calculated by CVODES.
    virtual double** GetSensSolution(int n_sensi, int n_species);

    //! Initialises the double** m_sensitivity matrix and sets to identity matrix
    virtual void InitialiseSensMatrix(int n_sensi, int n_species);

    //! Retrieves the number of sensitivity parameters in the problem.
    virtual unsigned int GetNumSens() const;

    //! Destroys the double** m_sensitivity matrix
    virtual void DestroySensMatrix(int n_species);

    //! Retrieves the number of equations solved by the ode solver
    unsigned int GetNEquations() const;


    // COMPUTATION TIME.
    
    // Returns the number of CT time variables tracked by this
    // solver type.
    virtual unsigned int CT_Count(void) const;

    // Outputs internal computation time data to the given
    // binary stream.
    virtual void OutputCT(std::ostream &out) const;

    // Adds the CT descriptions to a vector of strings.
    virtual void CT_Names(
        std::vector<std::string> &names, // Vector of CT names.
        unsigned int start=0 // Optional start index in vector.
        ) const;

    // Attach sensitivity to ODE_Solver by making copy.
    void AttachSensitivity(SensitivityAnalyzer &sensi) const;

    // Outputs sensitivity results to given file stream.
    void OutputSensitivity(std::fstream &fout, const Mops::Reactor &r, void *sim) const;
    //const ODE_Solver &GetODE_Solver() const { return m_ode; };

protected:
    // ODE SOLVER.

    ODE_Solver m_ode; // The ODE solver used by the mops solver.


    // SOLVER SETTINGS.

    // Default error tolerances for the ODE solver.
    real m_atol, m_rtol;

    // SENSITIVITY SETTINGS

    //! Boolean variable for setting LOI status
    bool m_LOIEnable;

    //! Minimum possible LOI value to keep a species in a mechanism.
    double m_LOIComp;

    // Under-relaxation coefficient.
    real m_rlx_coeff;

/*
    // Number of runs to perform.
    unsigned int m_nruns;

    // Max. number of stochastic particles in sweep.
    unsigned int m_pcount;

    // Max. M0 value, for initial scaling of ensemble.
    real m_maxm0;
*/

    // COMPUTATION TIME.

    clock_t m_cpu_start, m_cpu_mark;
    double m_tottime, m_chemtime;

/*
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
*/

    // COMPUTATION TIME CALCULATION.

    // Calculates the time duration from a time mark to the
    // current time.
    double calcDeltaCT(double markt) const;

private:
/*
    // FILE OUTPUT.

    // Writes the current reactor state to the output file.
    void fileOutput(const Reactor &r) const;
*/
};
};

#endif
