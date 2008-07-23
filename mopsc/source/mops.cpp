/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This is the main file for the mops solver.  Includes the driver
    program.

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

#include "mops.h"
#include "sprog.h"
#include "sweep.h"
#include <vector>
#include <stdexcept>

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
    bool foldfmt       = false; // Settings file format, new format not yet implemented.
    SolverType soltype = GPC;
    int diag = 0; // Diagnostics level.

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

        } else if (strcmp(argv[i], "-gpc") == 0) {
            // Solver gas-phase chemistry only.
            soltype = GPC;
        } else if (strcmp(argv[i], "-strang") == 0) {
            // Use Strang splitting to couple gas-phase and particle
            // system.
            soltype = Strang;
        } else if (strcmp(argv[i], "-predcor") == 0) {
            // Use Split-Predictor---Split-Corrector algorithm to 
            // couple gas-phase and particle system.
            soltype = PredCor;
        } else if (strcmp(argv[i], "-flamepp") == 0) {
            // Post-process a flame gas-phase chemistry
            // profile, just like sweep1.
            soltype = FlamePP;
//        } else if (strcmp(argv[i], "-momic") == 0) {
//            // Use method-of-moments to solve particles (1D only!).
//            soltype = MoMIC;

        // the next statements determine diagnostics level (if any).
        } else if (strcmp(argv[i], "-diag1") == 0) {
            diag = 1; // Minimal diagnostics printed to console.
        } else if (strcmp(argv[i], "-diag2") == 0) {
            diag = 2;
        } else if (strcmp(argv[i], "-diag3") == 0) {
            diag = 3;
        } else if (strcmp(argv[i], "-diag4") == 0) {
            diag = 4; // Full diagnostics.

        } else {
            // Currently nothing else is implemented here.  However, in future
            // the settings file name will not be set with the -r switch and
            // shall be read in this section.
            settfile = argv[i];
            foldfmt  = false;
        }
    }

    // Define all the objects required to run the simulation.
    Solver *solver   = NULL; // The solver.
    Reactor *reactor = NULL; // Reactor to solve.
    Mechanism mech;          // Chemical and particle mechanism.
    timevector times;        // A list of output times and step counts.
    Simulator sim;           // The simulator.

    // Create the solver.
    try {
        switch (soltype) {
            case OpSplit:
                // Not implemented yet.
                printf("Attempted to use simple splitting solver, which is not"
                       " yet implemented (mops, main).\n\n");
                return -1;
            case Strang:
                solver = new StrangSolver();
                break;
            case PredCor:
                solver = new PredCorSolver();
                sim.SetOutputEveryIter(true);
                break;
            case MoMIC:
                // Not implemented yet.
                printf("Attempted to use MoMIC solver, which is not"
                       " yet implemented (mops, main).\n\n");
                return -1;
            case FlamePP:
                // Post-process a gas-phase profile.
                solver = new Sweep::FlameSolver();
                break;
            case GPC:
            default:
                solver = new Solver();
                break;
        }
    } catch (std::logic_error &le) {
        printf("mops: Failed to initialise solver due to bad inputs.  Message:\n  ");
        printf(le.what());
        printf("\n\n");
        return -1;
    } catch (std::runtime_error &re) {
        printf("mops: Failed to initialise solver due to a program error.  Message:\n  ");
        printf(re.what());
        printf("\n\n");
        return -1;
    }

    // Read the chemical mechanism / profile.
    try {
        if (soltype != FlamePP) {
            Sprog::IO::MechanismParser::ReadChemkin(chemfile, mech, thermfile, diag);
            if (diag>0) mech.WriteDiagnostics("ckmech.diag");
        } else {
            dynamic_cast<Sweep::FlameSolver*>(solver)->LoadGasProfile(chemfile, mech);
        }
    } catch (std::logic_error &le) {
        printf("mops: Failed to read chemical mechanism/profile due to bad inputs.  Message:\n\n");
        printf(le.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        return -1;
    } catch (std::runtime_error &re) {
        printf("mops: Failed to read chemical mechanism/profile due to a program error.  Message:\n\n");
        printf(re.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        return -1;
    }

    // Read the particle mechanism.
    try {
        if (soltype != GPC) {
            mech.ParticleMech().SetSpecies(mech.Species());
            Sweep::MechParser::Read(swpfile, mech.ParticleMech());
            mech.ParticleMech().AddCoagulation();
        }
    } catch (std::logic_error &le) {
        printf("mops: Failed to read particle mechanism due to bad inputs.  Message:\n  ");
        printf(le.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        return -1;
    } catch (std::runtime_error &re) {
        printf("mops: Failed to read particle mechanism due to a program error.  Message:\n  ");
        printf(re.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        return -1;
    }

    // Read the settings file.
    try {
        if (foldfmt) {
            // Old file format.
            reactor = Settings_IO::LoadFromXML_V1(settfile, reactor, times, sim, *solver, mech);
        } else {
            // New format.
            reactor = Settings_IO::LoadFromXML(settfile, reactor, times, sim, *solver, mech);
        }
    } catch (std::logic_error &le) {
        printf("mops: Failed to load settings file due to bad inputs.  Message:\n  ");
        printf(le.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        return -1;
    } catch (std::runtime_error &re) {
        printf("mops: Failed to load settings file due to a program error.  Message:\n  ");
        printf(re.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        return -1;
    }

    // Solve reactor.
    try {
        if (fsolve) sim.RunSimulation(*reactor, times, *solver);
    } catch (std::logic_error &le) {
        printf("mops: Failed to solve reactor due to bad inputs.  Message:\n  ");
        printf(le.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        delete reactor;
        return -1;
    } catch (std::runtime_error &re) {
        printf("mops: Failed to solve reactor due to a program error.  Message:\n  ");
        printf(re.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        delete reactor;
        return -1;
    }

    // Post-process.
    try {
        if (fpostprocess) sim.PostProcess();
    } catch (std::logic_error &le) {
        printf("mops: Failed to post-process due to bad inputs.  Message:\n  ");
        printf(le.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        delete reactor;
        return -1;
    } catch (std::runtime_error &re) {
        printf("mops: Failed to post-process due to a program error.  Message:\n  ");
        printf(re.what());
        printf("\n\n");
        delete solver; // Must clear memory now.
        delete reactor;
        return -1;
    }

    // Clear up memory.
    delete solver;
    delete reactor;

    printf("mops: Simulation completed successfully!\n");
    printf("mops: Thank you for choosing mops for your particle modelling!\n");
    return 0;
}
