// main.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <stdio.h>
#include <exception>
#include <time.h>
#include "sweep.h"
#include "console_io.h"

using namespace std;
using namespace Sweep;

// GLOBAL VARIABLES.
std::ofstream ofile; // Output file stream.
Console_IO console;  // Console output class.


// INPUT FUNCTION DECLARATIONS.

// Reads the gas-phase chemistry profile from a TAB
// separated ASCII file.
void ReadProfile(
    const std::string &file,         // File name.
    Sweep::Solver::GasProfile &prof, // Profile chemical conditions.
    Sprog::SpeciesPtrVector &species // Species definitions read from profile.
    );


// OUTPUT FUNCTION DECLARATIONS.

// Sets up binary file output.
void BeginFileOutput(
    const std::string &file, // File name.
    unsigned int run         // Run number.
    );

// Writes basic system information to the binary output file.
void FileOutput(
    Sweep::real t,                // Time (s).
    const Sweep::Cell &sys,       // System to output.
    const Sweep::Mechanism &mech, // Mechanism which defines the system.
    );

// Prints the console header row.
void PrintConsoleHeader(const Sweep::Mechanism &mech);

// Prints some system information to the console.
void PrintToConsole(
    Sweep::real t,                // Time (s).
    const Sweep::Cell &sys,       // System to output.
    const Sweep::Mechanism &mech, // Mechanism which defines the system.
    );

// Writes the entire particle system to file for the purpose
// of generate particle size distributions.
void WritePSL(
    const std::string &file,      // Root file name for output.
    Sweep::real t,                // Time (s).
    const Sweep::Cell &sys,       // System to output.
    const Sweep::Mechanism &mech, // Mechanism which defines the system.
    );


// MAIN FUNCTION.

int main(int argc, char* argv[])
{
    // Program timer.
    clock_t ct1, ct2;

    // Variables for the solver and the solver object.
    int i=0, n=200, r=0, nruns=10;
    real t = 0.0, dt = 0.00025;
    Solver solver;
    string outfile = "sweep3-all.dat";

    // Read the flame chemistry profile.
    Solver::GasProfile flame;
    Sprog::SpeciesPtrVector species;
    ReadProfile("chem-debug.dat", flame, species);

    // Create a mechanism.
    Mechanism mech;
    mech.SetSpecies(species);

    // Read mechanism from file.
    try {
        MechanismParser::ReadXML("sweep.xml", mech);
    } catch (exception &e) {
        // Failed to read mechanism.
        return -1;
    }

    // Create a particle system.   
    Cell sys(species);
    
    // Run the simulation.
    for (r=1; r<=nruns; r++) {
        ct1 = clock();
        printf("Run number: %d\n", r);
        sys.Reset(3.5e11);
        t = 0.0;

        // Begin output to file and console.
        BeginFileOutput(outfile, r);
        PrinteConsoleHeader(mech);
        PrintToConsole(t, sys, mech);

        for (i=1; i<=n; i++) {
            // Solve time step.
            solver.Run(t, t+dt, flame, sys, mech);

            // Provide some output.
            FileOutput(t, sys, mech);
            PrintToConsole(t, sys, mech);
        }

        // Write the entire particle system to file.
        WritePSL(outfile, t, sys, mech);

        // Stop file output.
        EndFileOutput();

        ct2 = clock();
        cout << (real)(ct2-ct1) / (real)CLOCKS_PER_SEC << " seconds for run." << endl;
    }


    // Post-process output.
    PostProcess(outfile, nruns, mech);
    PostProcessPSLs(outfile, nruns, mech);

	return 0;
}
