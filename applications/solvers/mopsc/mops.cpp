 /*!
  * @file   mops.cpp
  * @author Matthew Celnik, William Menz
  * @brief  Main source file for MOPS
  *
  *   Licence:
  *
  *      mops is free software; you can redistribute it and/or
  *      modify it under the terms of the GNU Lesser General Public License
  *      as published by the Free Software Foundation; either version 2
  *      of the License, or (at your option) any later version.
  *
  *      This program is distributed in the hope that it will be useful,
  *      but WITHOUT ANY WARRANTY; without even the implied warranty of
  *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *      GNU Lesser General Public License for more details.
  *
  *      You should have received a copy of the GNU Lesser General Public
  *      License along with this program; if not, write to the Free Software
  *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  *      02111-1307, USA.
  *
  *   Contact:
  *      Prof Markus Kraft
  *      Dept of Chemical Engineering
  *      University of Cambridge
  *      New Museums Site
  *      Pembroke Street
  *      Cambridge
  *      CB2 3RA, UK
  *
  *      Email:       mk306@cam.ac.uk
  *      Website:     http://como.cheng.cam.ac.uk
*/



#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <string>

#include "mops.h"
using namespace std;

int main(int argc, char* argv[])
{
    // Declare empty filenames
    string ifile("");       // Main input file
    string cfile("");       // Chemical mechanism
    string tfile("");       // Thermochemical data
    string tranfile("");    // Transport data
    string sfile("");       // Particle mechanism
    string senfile("");     // Sensitivity
    string gpfile("");      // Gas-phase profile
    string surfcfile("");   // Surface chemistry
    string surftfile("");   // Surface thermocehmical data

    // Declare solver options
    size_t rand(0);         // Random seed
    Mops::SolverType soltype = Mops::GPC;
    bool fsurf(false);      // Surface capability on?
    bool fsen(false);       // Sensitivity analysis on?

    // Output options
    int diag(0);            // Diagnostics
    bool fpostproc(false);  // Should the binary files be postprocessed?
    bool fsolve(true);      // Should the system be solved?
    bool fjumps(false);     // Should a jumps file be written?
    bool fpah(false);       // Should full PAHPP data be postprocessed?
	bool fpp(false);       // Should full primary particle data be postprocessed?
    bool fensembles(false); // Should an *.ens file be written?
    bool fnew(false);       // Should the new network interface be used?
    bool fwdotA4(false);    //!< Should postprocess based on the molar rate of
                            //!< production by chemical reaction of the
                            //!< inception species?

    try {

        // Generic options for the program
        po::options_description generic("Generic options");

        generic.add_options()
        ("help,h", "print usage message")
        ("version,v", "print version number")
        ("new,w", "use new network interface")
        ;

        // Paths to input files
        po::options_description inpfiles("Input file options");

        inpfiles.add_options()
        ("mops,r", po::value(&ifile)->default_value("mops.inx"), "path to main input file")
        ("chem,c", po::value(&cfile)->default_value("chem.inp"), "path to chemical mechanism")
        ("therm,t", po::value(&tfile)->default_value("therm.dat"), "path to thermochemical data")
        ("trans,n", po::value(&sfile)->default_value("tran.dat"), "path to transport data")
        ("gasphase,g", po::value(&gpfile)->default_value("gasphase.inp"), "path to gas phase profile")
        ("sweep,s", po::value(&sfile)->default_value("sweep.xml"), "path to particle mechanism")
        ("sensi,q", po::value(&senfile)->default_value("sensi.xml"), "path to sensitivity analysis")
        ("schem", po::value(&surfcfile)->default_value("surfchem.inp"), "path to surface chemical mechanism")
        ("stherm", po::value(&surftfile)->default_value("surftherm.dat"), "path to surface thermochemical data")
        ;

        // Solver options
        po::options_description opt_solver("Solver options");

        opt_solver.add_options()
        ("rand,e", po::value(&rand)->default_value(456), "adjust random seed value")
        ("surf", "turn-on surface chemistry")
        ("opsplit", "use (simple) opsplit solver")
        ("strang", "use strang solver")
        ("predcor", "use predcor solver")
        ("flamepp", "use flamepp solver")
        ;

        // Output options
        po::options_description opt_out("Output options");

        // TODO: Move output options to .inx file specifications
        opt_out.add_options()
        ("postproc,p", "postprocess files")
        ("only,o", "postprocess files only (don't solve)")
        ("diag", po::value(&diag)->default_value(0), "set diagnostics level (0-4)")
        ("ensemble", "write full ensembles to binary files")
        ("ppah", "write full PAHPP data")
		("ppri", "write full primary particle data")
        ("jumps", "write stochastic jumps data")
        ("wdotA4", "postprocess based on the molar rate of production by chemical reaction of the inception species")
        ;

        // Combine sets of program options
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(inpfiles).add(opt_solver).add(opt_out);

        // Parse the command line
        po::variables_map vm;
        po::parsed_options parsed
            = po::command_line_parser(argc, argv).options(cmdline_options).run();
        store(parsed, vm);

        // Print program help and exit
        if (vm.count("help")) {
            cout << cmdline_options << "\n";
            return 0;
        }

        po::notify(vm);

        // Print version number and exit
        if (vm.count("version")) {
            cout << "2.0\n";
            return 0;
        }

        if (vm.count("new")) fnew = true;

        // Assign command-line variables to filenames
        ifile = vm["mops"].as< string >();
        cfile = vm["chem"].as< string >();
        tfile = vm["therm"].as< string >();
        sfile = vm["sweep"].as< string >();
        gpfile = vm["gasphase"].as< string >();
        if (!(vm["sensi"].defaulted())) {
            senfile = vm["sensi"].as< string >();
            fsen = true;
        }
        if (vm["trans"].defaulted()) tranfile = "NOT READ";

        // Check for surface chemistry
        if (vm.count("surf")) fsurf = true;
        if (!fsurf) {
            surfcfile = "NOT READ";
            surftfile = "NOT READ";
        } else {
            surfcfile = vm["schem"].as< std::string >();
            surftfile = vm["stherm"].as< std::string >();
        }

        // Get the solver type
        if (vm.count("opsplit")) soltype = Mops::OpSplit;
        if (vm.count("strang")) soltype = Mops::Strang;
        if (vm.count("predcor")) soltype = Mops::PredCor;
        if (vm.count("flamepp")) soltype = Mops::FlamePP;

        // Get the seed
        rand = vm["rand"].as< size_t >();

        // Get the output options
        if (vm.count("postproc")) fpostproc = true;
        if (vm.count("only")) {fsolve = false; fpostproc = true;}
        diag = vm["diag"].as< int >();
        if (vm.count("ppah")) fpah = true;
		if (vm.count("ppri")) fpp = true;
        if (vm.count("jumps")) fjumps = true;
        if (vm.count("ensemble")) fensembles = true;
        if (vm.count("wdotA4")) fwdotA4 = true;
    }

    // Display any error messages from incorrect command-line flags
    catch(exception& e) {
        std::cerr << "mops: Error getting options. Message:\n"
            << e.what() << "\n";
        return 1;
    }

    // Start main program
    cout << "Using the following files for input:\n" <<
        "  main: " << ifile << "\n" <<
        "  chem: " << cfile << "\n" <<
        "  therm: " << tfile << "\n";
    if (tranfile != "NOT READ") cout << "  trans: " << tranfile << "\n";
    if (soltype != Mops::GPC) cout << "  sweep: " << sfile << "\n";
    if (soltype == Mops::FlamePP) cout << "  gasphase: " << gpfile << "\n";
    if (fsen) cout << "  sensitivity: " << senfile << "\n";
    if (fsurf) {
        cout << "  schem: " << surfcfile << "\n" <<
                "  stherm: " << surftfile << "\n";
    }

    // Define all the objects required to run the simulation.
    Mops::Solver *solver   = NULL; // The solver.
    Mops::Reactor *reactor = NULL; // Reactor to solve.
    Mops::Mechanism mech;          // Chemical and particle mechanism.
    Mops::timevector times;        // A list of output times and step counts.
    Mops::Simulator sim;           // The simulator.
    Mops::ReactorNetwork *net;     // The network.
    Mops::NetworkSimulator *nsim;  // The network simulator.

    // Activate output options
    sim.SetWriteJumpFile(fjumps);
    sim.SetWriteEnsembleFile(fensembles);
    sim.SetWritePAH(fpah);
	sim.SetWritePP(fpp);

    // Create the solver
    solver = Mops::SolverFactory::Create(soltype);

    if (soltype == Mops::PredCor) sim.SetOutputEveryIter(true);

    // Load the (chemical) mechanism
    try {
        if (fsurf) {
            Sprog::IO::MechanismParser::ReadChemkin(
                    cfile,
                    surfcfile,
                    mech.GasMech(),
                    tfile,
                    surftfile,
                    diag,
                    tranfile);
        } else {
            Sprog::IO::MechanismParser::ReadChemkin(
                    cfile,
                    mech.GasMech(),
                    tfile,
                    diag,
                    tranfile);
        }
    } catch (std::logic_error &le) {
        std::cerr << "mops: Failed to read chemical mechanism due to bad inputs. Message:\n"
            << le.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    } catch (std::runtime_error &re) {
        std::cerr << "mops: Failed to read chemical mechanism due to program error. Message:\n"
            << re.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    }

    // Set the sepcies in the particle mechanism
    // (needed even if only gas-phase solver called)
    mech.ParticleMech().SetSpecies(mech.GasMech().Species());

    if (fwdotA4) {
        mech.ParticleMech().setPostprocessingType(Sweep::ParticleModel::wdotA4);
    }

    try {
        // Load the gas profile for flamepp calculations
        if (soltype == Mops::FlamePP)
            dynamic_cast<Sweep::FlameSolver*>(solver)->LoadGasProfile(gpfile, mech);

    } catch (std::logic_error &le) {
        std::cerr << "mops: Failed to read chemical profile due to bad inputs. Message:\n"
            << le.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    } catch (std::runtime_error &re) {
        std::cerr << "mops: Failed to read chemical profile due to program error. Message:\n"
            << re.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    }

    // Write some diagnostics
    if (diag > 0) mech.GasMech().WriteDiagnostics("ckmech.diag");

    try {
        if (soltype != GPC) {
            Sweep::MechParser::Read(sfile, mech.ParticleMech());
        }
    } catch (std::logic_error &le) {
        std::cerr << "mops: Failed to read particle mechanism due to bad inputs. Message:\n  "
            << le.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    } catch (std::runtime_error &re) {
        std::cerr << "mops: Failed to read particle mechanism due to a program error. Message:\n  "
            << re.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    }
	
	//Initialise the rng
	Sweep::rng_type rng = sim.SetRNG(rand, 0);

    // Read the settings file.
    try {
        if (fnew) {
            net = Mops::Settings_IO::LoadNetwork(ifile, times, sim, *solver, mech);
            nsim = new Mops::NetworkSimulator(sim, times);
        } else
            reactor = Mops::Settings_IO::LoadFromXML(ifile, reactor, times, sim, *solver, mech, rng);
    } catch (std::logic_error &le) {
        std::cerr << "mops: Failed to load MOPS settings file due to bad inputs. Message:\n  "
            << le.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    } catch (std::runtime_error &re) {
        std::cerr << "mops: Failed to load MOPS settings file due to a program error. Message:\n  "
            << re.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    }

    // Load the sensitivity analyser if applicable
    try {
        if (fsen) {
            SensitivityAnalyzer *sensi = new SensitivityAnalyzer();
            sensi->SetupProblem(mech, *reactor, senfile);

            if (sensi->isEnable()) solver->AttachSensitivity(*sensi);
            delete sensi;
        }
    } catch (std::logic_error &le) {
        std::cerr << "mops: Failed to load sensitivity file due to bad inputs. Message:\n  "
            << le.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    } catch (std::runtime_error &re) {
        std::cerr << "mops: Failed to load sensitivity file due to a program error. Message:\n  "
            << re.what() << "\n";
        delete solver; // Must clear memory now.
        return -1;
    }

    // Solve the reactor
    try {
        if (fsolve) {
            // Solve the network if enabled..
            if (fnew) {
                nsim->Run(*net, *solver, rand);
            } else {
                sim.SetTimeVector(times);
                sim.RunSimulation(*reactor, *solver, rand, rng);
            }
        }
    } catch (std::logic_error &le) {
        std::cerr << "mops: Failed to solve reactor due to bad inputs.  Message:\n  "
                << le.what() << "\n";
        delete solver; // Must clear memory now.
        delete reactor;
        return -1;
    } catch (std::runtime_error &re) {
        std::cerr << "mops: Failed to solve reactor due to a program error.  Message:\n  "
            << re.what() << "\n";
        delete solver; // Must clear memory now.
        delete reactor;
        return -1;
    }

    // Post-process.
    try {
        if (fpostproc) {
            if (fnew)
                nsim->PostProcess();
            else
                sim.PostProcess();
        }
    } catch (std::logic_error &le) {
        std::cerr << "mops: Failed to post-process due to bad inputs.  Message:\n  "
            << le.what() << "\n";
        delete solver; // Must clear memory now.
        delete reactor;
        return -1;
    } catch (std::runtime_error &re) {
        std::cerr << "mops: Failed to post-process due to a program error.  Message:\n  "
            << re.what() << "\n";
        delete solver; // Must clear memory now.
        delete reactor;
        return -1;
    }

    // Clear up memory.
    delete solver;
    delete reactor;
    if (fnew) {delete net; delete nsim;}

    printf("mops: Simulation completed successfully!\n");
    printf("mops: Thank you for choosing mops for your particle modelling!\n");

    return 0;
}
