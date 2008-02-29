/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    This is the main file for the mops solver.  Includes the driver
    program.
*/

#include "mops.h"
#include "sprog.h"
#include "sweep.h"
#include <vector>

using namespace Mops;
using namespace std;

int main(int argc, char *argv[])
{
    // Command line arguments with default values.
    string chemfile("chem.inp");
    string thermfile("therm.dat");
    string settfile("mops.inx");
    string swpfile("sweep.xml");
    bool fsolve        = true;  // Default is to solve ..
    bool fpostprocess  = false; // .. but not post-process.
    bool foldfmt       = true;  // Settings file format, new format not yet implemented.
    SolverType soltype = GPC;

    // Read command line arguments.
    for (int i=1; i!=argc; ++i) {
        if (strcmp(argv[i], "-c") == 0) {
            // Chemical mechanism file (CK format).
            chemfile = argv[++i];
        } else if (strcmp(argv[i], "-t") == 0) {
            // Thermodynamic properties file (CK format).
            thermfile = argv[++i];
        } else if (strcmp(argv[i], "-r") == 0) {
            // Settings file (F90 mops format).
            settfile = argv[++i];
            foldfmt  = true;
        } else if (strcmp(argv[i], "-s") == 0) {
            // Sweep mechanism file.
            swpfile = argv[++i];
        } else if (strcmp(argv[i], "-p") == 0) {
            // Post-processing switch.  Used to turn PP on.
            fpostprocess = true;
        } else if (strcmp(argv[i], "-po") == 0) {
            // "Post-process only" switch.  Post-processes but doesn't solve.
            // Use this if you have previously run a simulation and want
            // human-readable CSV formatted files with the results.
            fsolve       = false;
            fpostprocess = true;

        // The next statements select the type of solver to use.  The
        // default is to solve gas-phase only, with no particle system.

        } else if (strcmp(argv[i], "-strang") == 0) {
            // Use Strang splitting to couple gas-phase and particle
            // system.
            soltype = Strang;
        } else if (strcmp(argv[i], "-flamepp") == 0) {
            // Post-process a flame gas-phase chemistry
            // profile, just like sweep1.
            soltype = FlamePP;
//        } else if (strcmp(argv[i], "-momic") == 0) {
//            // Use method-of-moments to solve particles (1D only!).
//            soltype = MoMIC;

        } else {
            // Currently nothing else is implemented here.  However, in future
            // the settings file name will not be set with the -r switch and
            // shall be read in this section.
        }
    }

    // Define all the objects required to run the simulation.
    Solver *solver   = NULL; // The solver.
    Reactor *reactor = NULL; // Reactor to solve.
    Mechanism mech;          // Chemical and particle mechanism.
    timevector times;        // A list of output times and step counts.

    // Create the solver.
    switch (soltype) {
        case GPC:
            solver = new Solver();
            break;
        case OpSplit:
            // Not implemented yet.
            return -1;
        case Strang:
            // Not implemented yet.
            solver = new StrangSolver();
            break;
        case MoMIC:
            // Not implemented yet.
            return -1;
        case FlamePP:
            // Post-process a gas-phase profile.
            solver = new Sweep::FlameSolver();
            break;
    }

    // Read the chemical mechanism / profile.
    if (soltype != FlamePP) {
        Sprog::IO::MechanismParser::ReadChemkin(chemfile, mech, thermfile);
    } else {
        dynamic_cast<Sweep::FlameSolver*>(solver)->LoadGasProfile(chemfile, mech);
    }
    
    // Read the particle mechanism.
    mech.ParticleMech().SetSpecies(mech.Species());
    Sweep::MechParser::Read(swpfile, mech.ParticleMech());
    mech.ParticleMech().AddCoagulation();

    // Read the settings file.
    if (foldfmt) {
        // Old file format.
        reactor = Settings_IO::LoadFromXML_V1(settfile, reactor, times, *solver, mech);
        reactor->Initialise(0.0);
    } else {
        // New format not yet implemented.
        delete solver;
        delete reactor;
        return -1;
    }

    // Solve reactor.
    if (fsolve) solver->SolveReactor(*reactor, times, solver->RunCount());

    // Post-process.
    if (fpostprocess) solver->PostProcess(solver->OutputFile(), solver->RunCount());

    // Clear up memory.
    delete solver;
    delete reactor;

    return 0;
}
