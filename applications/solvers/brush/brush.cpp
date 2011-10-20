/*!
 *\file brush.cpp
 *\author Robert I A Patterson
 *
 *\brief Main program for brush
 *
 *  Copyright (C) 2009 Robert I A Patterson.
 *

 Licence:
    This file is part of "brush".

    brush is free software; you can redistribute it and/or
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
    Prof Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
 *
 *
 *\main Brush is a program for the coupled simulation of laminar reacting flow
 * and particle population development in 1d systems.
 */

#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <memory>

#include "camxml.h"
#include "gpc_mech_io.h"
#include "mops_mechanism.h"
#include "swp_mech_parser.h"
#include "mops_timeinterval.h"
#include "mops_settings_io.h"

#include "geometry1d.h"
#include "reactor1d.h"
#include "pred_corr_solver.h"
#include "reset_chemistry.h"
#include "simulator.h"
#include "settings_io.h"

#include "linear_interpolator.hpp"

//! Write summary usage information to stdout
void printUsage();

using namespace Brush;

//! Run brush
int main(int argc, char* argv[])
{
    std::cout << "Starting brush\n";
    //======== Command line arguments ============================

    // Input file names with default values.
    std::string chemfile("chem.inp");
    std::string thermfile("therm.dat");
    std::string settfile("brush.xml");
    std::string swpfile("sweep.xml");
    std::string geomfile("geometry.xml");
    std::string chemsolnfile("chemsoln.dat");
    std::string partsolnfile("partsoln.xml");
    std::string tranfile("");

    // Offset for random number sequence so that independent realisations
    // can be computed in separated program instances.
    size_t randomSeedOffset = 0;

    // Diagnostic level
    int diag = 0;

    // Format for initial chemistry solution data
    ResetChemistry::InputFileType chemFileType = ResetChemistry::Camflow;

    for (int i=1; i!=argc; ++i) {
        if (std::strcmp(argv[i], "-c") == 0) {
            // Chemical mechanism file (CK format).
            chemfile = argv[++i];
        }
        else if (std::strcmp(argv[i], "-d") == 0) {
            // Initial estimate of chemistry solution
            chemsolnfile = argv[++i];
        }
        else if (std::strcmp(argv[i], "-a") == 0) {
            // Initial estimate of particle solution
            partsolnfile = argv[++i];
        }
        else if (strcmp(argv[i], "-t") == 0) {
            // Thermodynamic properties file (CK format).
            thermfile = argv[++i];
        }
        else if (strcmp(argv[i], "-r") == 0) {
            // Species transport properties file (CK format).
            tranfile = argv[++i];
        }
        else if (strcmp(argv[i], "-s") == 0) {
            // Sweep mechanism file.
            swpfile = argv[++i];
        }
        else if (strcmp(argv[i], "-g") == 0) {
            // Initial geometry file.
            geomfile = argv[++i];
        }
        else if (strcmp(argv[i], "-b") == 0) {
            // Initial geometry file.
            settfile = argv[++i];
        }
        else if (strcmp(argv[i], "-v") == 0) {
            // Verbosity
            diag = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "--premix-chem") == 0) {
            // Verbosity
            chemFileType = ResetChemistry::Premix;
        }
        else if (strcmp(argv[i], "-h") == 0) {
            // Print help message and exit
            printUsage();
            return 0;
        }
        else if (strcmp(argv[i], "-e") == 0) {
            // Verbosity
            randomSeedOffset = atoi(argv[++i]);
        }
        else {
            // Options do not make sense so stop give up and quit.
            std::cerr << "Unrecognised option " << argv[i] << " used for brush\n";
            return 1;
        }
    }


    //========= Read numerical settings ===========================
    std::cout << "Reading numerical settings\n";
    // Initial geometry
    std::auto_ptr<Geometry::Geometry1d> pGeom;
    try {
        // Load the XML
        CamXML::Document geomXML;
        geomXML.Load(geomfile);

        // Build the geometry object and hold a smart pointer to it
        pGeom.reset(new Geometry::Geometry1d(*geomXML.Root()));
    }
    catch (std::exception &e) {
        std::cerr << "Failed to load geometry from " << geomfile << ", because\n"
                  << e.what() << '\n';
        return 1;
    }
    // Log the geometry for testing purposes (remove prior to release)
    std::cout << pGeom->printMesh() << '\n';

    // Variables to hold data read from the input file
    Mops::timevector timeIntervals;
    int runs, iterations;
    std::string outputFileBaseName;
    std::vector<std::pair<real, real> > maxPCounts, maxM0s;
    bool splitDiffusion = false;
    bool splitAdvection = false;
    bool weightTransport = false;
    Sweep::Stats::IModelStats::StatBound statBound;

    try {
        // Load the XML
        CamXML::Document settingsXML;
        settingsXML.Load(settfile);
        const CamXML::Element * const root = settingsXML.Root();

        // Time intervals for stepping through the solution process with output
        const CamXML::Element * node = root->GetFirstChild("timeintervals");
        if (node != NULL) {
            Mops::Settings_IO::readTimeIntervals(*node, timeIntervals);
        } else {
            throw std::runtime_error("No time intervals found");
        }

        // Number of paths
        node = root->GetFirstChild("runs");
        if (node != NULL) {
            runs = atoi(node->Data().c_str());
        } else {
            throw std::runtime_error("A <runs> element must be supplied to indicate number of paths");
        }

        // Number of corrector iterations
        node = root->GetFirstChild("iter");
        if (node != NULL) {
            iterations = atoi(node->Data().c_str());
        } else {
            throw std::runtime_error("An <iter> element must be supplied to specify number of corrector iterations");
        }

        // Numerical method for diffusion
        node = root->GetFirstChild("diffusion");
        if ((node != NULL) && ("split" == node->Data())) {
            // simulate diffusion using spltting
            splitDiffusion = true;
        }

        // Numerical method for advection
        node = root->GetFirstChild("advection");
        if ((node != NULL) && ("split" == node->Data())) {
            // simulate advection using spltting
            splitAdvection = true;
        }

        // Numerical method for moving weighted particles between cells
        // should not be used with DSA
        node = root->GetFirstChild("weighttransport");
        if ((node != NULL) && ("weights" == node->Data())) {
            // adjust statistical weights when moving between cells to avoid cloning/killing
            weightTransport = true;
        }

        // Maximum number of computational particles per cell
        maxPCounts = Brush::Settings_IO::readProfile(root, "pcount");

        // Maximum particle concentration
        maxM0s = Brush::Settings_IO::readProfile(root, "maxm0");

        // Output details
        node = root->GetFirstChild("output");
        if(node != NULL) {
            const CamXML::Element * const filenameNode = node->GetFirstChild("filename");
            if(filenameNode != NULL) {
                outputFileBaseName = filenameNode->Data();
            }
            else {
                throw std::runtime_error("A <filename> element must be supplied in the output section");
            }

            // See if the stats bounds were specified in the input file
            const CamXML::Element * const statsBoundNode = node->GetFirstChild("statsbound");
            if(statsBoundNode != NULL) {
                // Read the data from the file
                Sweep::PropID statsboundPropertyID;
                real statsLowerBound;
                real statsUpperBound;
                Mops::Settings_IO::ReadStatsBound(*statsBoundNode, statsboundPropertyID, statsLowerBound, statsUpperBound);

                // Adjust statBound with the newly read data
                statBound.PID = statsboundPropertyID;
                statBound.Lower = statsLowerBound;
                statBound.Upper = statsUpperBound;
            }
            // If no statsbound was specified in the input file, the default constructed value of
            // statBound will be used.
        }
        else {
            throw std::runtime_error("A <output> element must be supplied with details for the required output");
        }

    }
    catch (std::exception &e) {
        std::cerr << "Failed to read settings from " << settfile << ", because\n"
                  << e.what() << '\n';
        return 1;
    }
    std::cout << "Read " << timeIntervals.size() << " time intervals\n";

    //========= Load chemical mechanism ==========================
    // For fixed chemistry simulations the mech does not have to
    // contain reactions, only the species names are required.
    Mops::Mechanism mech;
    try {
        if(diag > 0) {
            std::cout << "Reading chemical species...\n";
        }
        // See if there is a transport data file (not needed for advection only)
        if(tranfile.length() > 0)
            // Read species transport data along with rest of mechanism
            Sprog::IO::MechanismParser::ReadChemkin(chemfile, mech.GasMech(), thermfile, diag, tranfile);
        else
            // Skip species transport data - properties such as thermal conductivity will not be available
            Sprog::IO::MechanismParser::ReadChemkin(chemfile, mech.GasMech(), thermfile, diag);

        if (diag>0)
            mech.GasMech().WriteDiagnostics("ckmech.diag");
    }
    catch (std::exception &e) {
        std::cerr << "Failed to read chemical species from " << chemfile << ", because\n"
                  << e.what() << '\n';
        return 1;
    }

    //========= Load particle mechanism ==========================
    try {
        if(diag > 0) {
            std::cout << "Setting species on particle mechanism...\n";
        }
        mech.ParticleMech().SetSpecies(mech.GasMech().Species());
        if(diag > 0) {
            std::cout << "Reading particle mechanism...\n";
        }
        Sweep::MechParser::Read(swpfile, mech.ParticleMech());
        if(diag > 0) {
            std::cout << "Read particle mechanism with " << mech.ParticleMech().ProcessCount()
                      << " processes\n";
            if(diag > 1) {
                // Get names of the processes
                std::vector<std::string> procNames;
                mech.ParticleMech().GetProcessNames(procNames);

                // Print the process names out one line at a time
                std::vector<std::string>::const_iterator it = procNames.begin();
                const std::vector<std::string>::const_iterator itEnd = procNames.end();
                while(it != itEnd) {
                    std::cout << ' ' << *it++ << '\n';
                }
            }
        }
    }
    catch (std::exception &e) {
        std::cerr << "Failed to read particle mechanism from " << swpfile << ", because\n"
                  << e.what() << '\n';
        return 1;
    }


    //========= Read initial guess for solution ==================
    std::list<ParticlePopulationPoint1d> initialPopulationPoints;
    try {
        // Load the XML
        CamXML::Document initialParticlesXML;
        initialParticlesXML.Load(partsolnfile);
        const CamXML::Element * const root = initialParticlesXML.Root();

        // Now get the XML for each separate population
        typedef std::vector<CamXML::Element*> xml_element_list;
        xml_element_list populations;
        root->GetChildren("population", populations);

        // Read the position and details of each population
        for(xml_element_list::const_iterator it = populations.begin();
            it != populations.end(); ++it) {

            ParticlePopulationPoint1d populationDetails;

            // Get the position
            const CamXML::Attribute *attr = (*it)->GetAttribute("x");
            if(attr) {
                populationDetails.position = atof(attr->GetValue().c_str());
            }
            else {
                throw std::runtime_error("No x (position) attribute for a <population> element.");
            }

            // Get the number density of the particles
            const CamXML::Element * const m0XML = (*it)->GetFirstChild("m0");
            if(m0XML) {
                populationDetails.m0 = atof(m0XML->Data().c_str());
            }
            else {
                throw std::runtime_error("No <m0> (number density) child element for a <population> element.");
            }

            // Read in the population
            populationDetails.particleList =
                    Mops::Settings_IO::ReadInitialParticles(**it, mech.ParticleMech());

            initialPopulationPoints.push_back(populationDetails);
        }
    }
    catch(std::exception &e) {
        std::cerr << "Failed to read initial particle solution from " << partsolnfile
                  << ", because " << e.what() << '\n';
        return 1;
    }


    std::auto_ptr<ResetChemistry> pInitialChem;
    // Initial chemistry is the fixed chemistry
    try {
        if(diag > 0) {
            std::cout << "Reading initial chemistry solution...\n";
        }
        pInitialChem.reset(new ResetChemistry(chemsolnfile, chemFileType, mech.GasMech(), diag));
        if(diag > 0) {
            std::cout << "Read initial chemistry solution\n";
        }
    }
    catch(std::exception &e) {
        std::cerr << "Failed to read initial chemistry solution from " << chemsolnfile
                  << ", because " << e.what() << '\n';
        return 1;
    }

    //========= Build the initial reactor ========================
    Reactor1d initialReactor(*pGeom, mech.GasMech(), mech.ParticleMech(), maxPCounts, maxM0s);

    // Put the initial species concentrations into the reactor.
    // Second argument indicates chemical conditions are fixed and not updated
    // as a result of particle events that affect the gaseous species concentrations.
    // This will need changing when we have proper coupling.
    initialReactor.ReplaceChemistry(*pInitialChem, true);

    // Put the initial particles into the reactor
    initialReactor.ReplaceParticles(initialPopulationPoints.begin(), initialPopulationPoints.end());

    //========= Now run the simulation ===========================
    Simulator sim(runs, iterations, timeIntervals, initialReactor, *pInitialChem,
                  outputFileBaseName, statBound, splitDiffusion, splitAdvection,
                  weightTransport);
    sim.runSimulation(randomSeedOffset);

    //========= Output ===========================================


    //========= Clean up =========================================

    return 0;
}

/*!
 * Method needs to be extended as new options are added.
 */
void printUsage() {
    std::cout << "Usage: brush [OPTION] ...\n";
    std::cout << "Simulate a particle population in an inhomogeneous reactor\n";
    std::cout << "-a INITIAL-PARTICLE-SOLUTION-FILE\n";
    std::cout << "-b BRUSH-SETTINGS-FILE\n";
    std::cout << "-c CHEMICAL-MECHANISM-FILE\n";
    std::cout << "-d INITIAL-CHEMISTRY-SOLUTION-FILE\n";
    std::cout << "-e offset for random number generator seed\n";
    std::cout << "-g GEOMETRY-FILE\n";
    std::cout << "-h output summary usage information\n";
    std::cout << "-r TRANSPORT-DATA-FILE\n";
    std::cout << "-s SWEEP-SETTINGS-FILE\n";
    std::cout << "-t THERMODYNAMICAL-DATA-FILE\n";
    std::cout << "-v verbosity (Integer indicating level of debug info, higher integers mean more output)\n";
    std::cout << "--premix-chem initial chemistry solution file is in Chemkin format\n";
}

