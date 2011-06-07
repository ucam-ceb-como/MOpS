/*!
  \author      Robert I A Patterson
  \author      Sashikumaar Ganesan

  Copyright (C) 2011 Robert I A Patterson & Sashikumaar Ganesan

  \file moonmd_interface.cpp
  \brief Interface between continuum flow solver and stochastic particle
         population balance solver

  Licence:

    This file is free software; you can redistribute it and/or
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
    Dr. Robert Patterson
    Weierstrass Institute
    Mohrenstrasse 39
    10117 Berlin
    Germany
    riap@cantab.net
*/

#include "moonmd_interface.h"

#include "geometry1d.h"
#include "settings_io.h"
#include "reactor1d.h"
#include "reset_chemistry.h"
#include "pred_corr_solver.h"
#include "simulator.h"

#include "mops_settings_io.h"
#include "mops_mechanism.h"

#include "swp_model_stats.h"
#include "mt19937.h"

#include "linear_interpolator.hpp"

#include <memory>
#include <iostream>
#include <stdexcept>
#include <limits>

/*!
 * Function private to this file
 *
 * Work out how long a particle should take to travel through each cell and return the shortest residence
 * time.  A mid point velocity integration is currently used.  This function is useful in determining the number of
 * time splitting steps which the solver must make between particle processes such as coagulation and inception
 * and particle transport.
 *
 *\param[in]         geom                Geometry defining the cells
 *\param[in]         solution_length     Communication is by vectors which must all be at least this long
 *\param[in]         solution_nodes      Distances from start of reactor (in m) at which velocity values are specified
 *\param[in]         velocity            Flow velocity in \f$\mathrm{ms}^{-1}\f$
 */
Brush::real minimumResidenceTime(const Geometry::Geometry1d& geom, const size_t solution_length,
                                 const double solution_nodes[], const double velocity[])
{
    const Utils::LinearInterpolator<double, double> velocityProfile(std::vector<double>(solution_nodes, solution_nodes + solution_length),
                                                                    std::vector<double>(velocity,       velocity + solution_length ));

    double minTime = std::numeric_limits<double>::max();
    for(size_t i = 0; i < geom.numCells(); ++i) {
        // In the 1d case cellVolume == length of cell
        double residenceTime = std::abs(geom.cellVolume(i) / velocityProfile.interpolate(geom.cellCentre(i)));
        minTime = residenceTime < minTime ? residenceTime : minTime;
    }
    return minTime;
}

/*!
 *\param[in]    chemfile     Path to file listing chemical species to consider
 *\param[in]    thermfile    Path to thermodynamic data file
 *\param[in]    settfile     Path to file specifying control parameters for Brush
 *\param[in]    swpfile      Path to XML file specifying the particle mechanism
 *\param[in]    partsolnfile Path to XML file specifying the initial particle population
 *\param[in]    chemsolnfile Path to file specifying the initial continuum phase composition
 *\param[in]    num_grid_nodes   Number of cell boundaries specified in grid_nodes
 *\param[in]    grid_nodes   Boundaries of cells to be used in Brush specified in meters.
 *
 *\pre [grid_nodes, grid_nodes + num_grid_nodes) must be a valid range
 *
 *\return    Pointer to new brush reactor object.   Client must call delete on this pointer.
 */
Brush::MooNMDInterface::particle_reactor_pointer 
 Brush::MooNMDInterface::InitialiseBrush(
    const std::string& chemfile, const std::string& thermfile, const std::string& settfile,
    const std::string& swpfile, const std::string& partsolnfile, const std::string& chemsolnfile,
    const size_t num_grid_nodes, const double grid_nodes[])
{
    std::cout << "Setting up Brush\n";

    // Look at working directory (to help debug tests)
    char *pString = NULL;
    pString = get_current_dir_name();
    std::cout << "Working directory is " << std::string(pString) << std::endl;
    delete[] pString;

    // diagnostic level (could be a user input)
   const int diag = 1;

    // Geometry
    std::auto_ptr<Geometry::Geometry1d> pGeom;
    // Build the geometry object and hold a smart pointer to it
    pGeom.reset(new Geometry::Geometry1d(std::vector<double>(grid_nodes, grid_nodes + num_grid_nodes),
                                          Geometry::dirichlet, Geometry::neumann));

    // Variables to hold data read from the input file
    Mops::timevector timeIntervals;
    int runs, iterations;
    std::string outputFileBaseName;
    std::vector<std::pair<real, real> > maxPCounts, maxM0s;
    bool splitDiffusion = false;
    bool splitAdvection = false;
    Sweep::Stats::IModelStats::StatBound statBound;

    {
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

    //========= Load chemical mechanism ==========================
    // For fixed chemistry simulations the mech does not have to
    // contain reactions, only the species names are required.
    Mops::Mechanism mech;
    {
        if(diag > 0) {
            std::cout << "Reading chemical mechanism...\n";
        }

        Sprog::IO::MechanismParser::ReadChemkin(chemfile, mech, thermfile, diag);

        if (diag>0)
            mech.WriteDiagnostics("ckmech.diag");
    }

    //========= Load particle mechanism ==========================
    {
        if(diag > 0) {
            std::cout << "Setting species on particle mechanism...\n";
        }
        mech.ParticleMech().SetSpecies(mech.Species());
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

    //======== Read the initial state of the continuum phase ================
    if(diag > 0) {
        std::cout << "Reading initial continuum solution...\n";
    }
    Brush::ResetChemistry initChem(chemsolnfile, Brush::ResetChemistry::Camflow, mech, diag);
    if(diag > 0) {
        std::cout << "Read initial continuum solution\n";
    }

    //========= Read initial guess for particle solution ==================
    std::list<ParticlePopulationPoint1d> initialPopulationPoints;
    {
        if(diag > 0) {
            std::cout << "Reading initial particle population...\n";
        }

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
                    Mops::Settings_IO::ReadInitialParticles(**it, mech.ParticleMech(), Sweep::genrand_int);

            initialPopulationPoints.push_back(populationDetails);
        }
    }

    //========= Build the initial reactor ========================
    Reactor1d *pReactor = new Reactor1d(*pGeom, mech, maxPCounts, maxM0s);

    std::cout << "Setting the initial continuum solution\n";
    pReactor->ReplaceChemistry(initChem, true);

    // Put the initial particles into the reactor
    pReactor->ReplaceParticles(initialPopulationPoints.begin(), initialPopulationPoints.end());

    return pReactor;
}

/*!
 * Simulate the particle processes from the current reactor time until t_stop
 * and calculate mean energy and mass source terms.
 *
 *\param[in,out]    reac                The reactor specifying the initial condition and particle dynamics
 *\param[in]         solution_length     Communication is by vectors which must all be at least this long
 *\param[in]         num_species         Number of species for which mass concentrations are supplied
 *\param[in]         solution_nodes      Distances from start of reactor (in m) at which solution values are specified
 *\param[in]         temperature         Reactor temperature in K
 *\param[in]         velocity            Flow velocity in \f$\mathrm{ms}^{-1}\f$
 *\param[in]         mass_concs          Concentrations of species in \f$\mathrm{kg m}^{-3}\f$
 *\param[out]        energy_souce        Energy release by particle processes in \f$\mathrm{J m}^{-3} \mathrm{s}^{-1}\f$
 *\param[out]        mass_conc_souces    Release of species into solution by particle processes in \f$\mathrm{kg m}^{-3}\mathrm{s}^{-1}\f$
 *
 *\pre  All the 1d arrays supplied by the caller must have length at least solution_length, this includes
 *       the output arrays.
 *
 *\pre  The arrays mass_concs and mass_conc_sources are effectively 2d arrays and must have length at least 
 *       solution_length * num_species
 *
 *\return    Pointer to updated brush reactor object.
 *
 * Sign convention for source terms: Positive source terms indicate release of energy or materical by the solid
 * phase.  Thus a positive energy source would indicate exothermic crystallisation and a positive mass source
 * would indicate that crystals were dissolving.
 *
 * The data in mass_concs shall be ordered such that data for species i begins at mass_concs + i * num_species,
 * that is the fastest varying index shall be the index for spatial position.  This is fairly easy to change,
 * but it is important to keep all comments and code consistent!  The data in mass_conc_sources shall be ordered
 * in the same way as that in mass_concs.
 */
Brush::MooNMDInterface::particle_reactor_pointer Brush::MooNMDInterface::RunParticlePhase(
    particle_reactor& reac, const double t_stop,
    const size_t solution_length,
    const size_t num_species,
    const double solution_nodes[],
    const double temperature[],
    const double velocity[],
    const double mass_concs[],
    double energy_source[],
    double mass_conc_sources[],
    std::ostream &moment_output)
{
    //========= Set the chemistry on the reactors ==========
    // First put the data in vectors, because this is what the interface to ResetChemistry expects
    const fvector solutionNodes(solution_nodes, solution_nodes + solution_length);
    const fvector temperatureVector(temperature, temperature + solution_length);
    const fvector velocityVector(velocity, velocity + solution_length);

    // Work out the total density in kg m^-3
    fvector density(solution_length, 0.0);
    // Loop over the spatial positions
    for (unsigned i = 0u; i < solution_length; ++i) {
        // Add up the contribution of all the species at current position
        for (unsigned j = 0u; j < num_species; ++j) {
            density[i] += mass_concs[i + j * solution_length];
        }
    }

    // Now calculate mass fractions (because this is what the ResetChemistry constructor requires)
    std::vector<fvector> massFracs(num_species, fvector(solution_length));
    // Loop over the species
    for (unsigned i = 0u; i < num_species; ++i) {
        // Scale by the total mass density at each point
        for (unsigned j = 0u; j < solution_length; ++j) {
            massFracs[i][j] = mass_concs[j + i * solution_length] / density[j];
        }
    }

    // Dummy vector to fill unused argument positions
    const fvector vec(solution_length);

    ResetChemistry chemInterpolator(solutionNodes, temperatureVector, density, velocityVector,
                                    vec, vec, vec, vec, vec, massFracs);

    // Put the chemistry into the reactor
    reac.ReplaceChemistry(chemInterpolator, true);

    //======== Simulate to update the particle population =======
    Brush::PredCorrSolver solver(chemInterpolator, 0, 0.0, 0.0, false, true);

    // Need to make sure that transport is simulated on intervals that do not
    // always take particles over cell boundaries.
    double minResidence = minimumResidenceTime(reac.getGeometry(), solution_length, solution_nodes, velocity);
    unsigned numSplittings = std::ceil(10 * t_stop / minResidence);
    std::cout << "Num splittings " << numSplittings << ", min residence " << minResidence << std::endl;
    solver.solve(reac, t_stop, numSplittings, 0);
    Brush::Simulator::saveParticleStats(reac, Sweep::Stats::IModelStats::StatBound(), moment_output);

    //======== Estimate the source terms to return to the client =====

    return &reac;
}

