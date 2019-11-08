/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Settings_IO class declared in the
    mops_settings_io.h header file.

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


#include "mops_settings_io.h"
#include "mops_timeinterval.h"
#include "mops_psr.h"
#include "mops_flow_stream.h"
#include "mops_reactor_factory.h"
#include "mops_simulator.h"
#include "mops_solver.h"

#include "camxml.h"
#include "string_functions.h"
#include <iostream>
#include <sstream>
#include <memory>
#include <cstdlib>

using namespace Mops;

//! Anonymous namespace for function only used within this file
namespace {

struct flowNode {
public:
    //! Is the stream an inflow (true) or outflow (false)?
    bool mIn;
    //! Mixture key
    std::string mMixture;
    //! Reactor key
    std::string mReactor;
    //! Fraction of the total in/outflow
    double mFrac;

    //! Constructor
    flowNode():
        mIn(true), mMixture(""), mReactor(""), mFrac(1.0) {};
};

// Returns the temperature in K by reading the value from the given
// XML node and checking the units.
double readTemperature(const CamXML::Element &node)
{
    // Check the temperature units.
    const CamXML::Attribute *attr;
    Sprog::TempUnits units = Sprog::Kelvin;
    attr = node.GetAttribute("units");
    if (attr != NULL) {
        if (attr->GetValue() == "K") {
            units = Sprog::Kelvin;
        } else if (attr->GetValue() == "oC") {
            units = Sprog::Celcius;
        } else
            throw std::runtime_error("Unknown temperature unit "
                    + std::string(attr->GetValue()) + " specified"
                    + " (::readTemperature).");
    }

    // Read the temperature (and convert if necessary).
    double T = Strings::cdble(node.Data());
    switch (units) {
        case Sprog::Celcius:
            T += 273.15;
            break;
        default:
            break;
    }

    return T;
}

// Returns the pressure in Pa by reading the value from the given
// XML node and checking the units.
double readPressure(const CamXML::Element &node)
{
    // Check the pressure units.
    const CamXML::Attribute *attr;
    Sprog::PressureUnits units = Sprog::Bar;
    attr = node.GetAttribute("units");
    if (attr != NULL) {
        if (attr->GetValue() == "Pa") {
            units = Sprog::Pascal;
        } else if (attr->GetValue() == "bar") {
            units = Sprog::Bar;
        } else if (attr->GetValue() == "atm") {
            units = Sprog::Atm;
        } else throw std::runtime_error("Unknown pressure unit "
                + std::string(attr->GetValue()) + " specified"
                + " (::readPressure).");
    }

    // Read the pressure (and convert if necessary).
    double P = Strings::cdble(node.Data());
    switch (units) {
        case Sprog::Bar:
            P *= 1.0e5;
            break;
        case Sprog::Atm:
            P *= 1.01325e5;
            break;
        default:
            break;
    }

    return P;
}

// Reads global simulation settings from the given XML node.
void readGlobalSettings(const CamXML::Element &node,
                                     Simulator &sim, Solver &solver)
{
    const CamXML::Element *subnode;

    // Read the relative error tolerance.
    subnode = node.GetFirstChild("rtol");
    if (subnode != NULL) {
        solver.SetRTOL(Strings::cdble(subnode->Data()));
    }

    // Read the absolute error tolerance.
    subnode = node.GetFirstChild("atol");
    if (subnode != NULL) {
        solver.SetATOL(Strings::cdble(subnode->Data()));
    }

    // Read the number of runs.
    subnode = node.GetFirstChild("runs");
    if (subnode != NULL) {
        sim.SetRunCount((int)Strings::cdble(subnode->Data()));
    }

    // Read the number of iterations.
    subnode = node.GetFirstChild("iter");
    if (subnode != NULL) {
        sim.SetIterCount((int)Strings::cdble(subnode->Data()));
    }

    // Read stochastic particle count.
    subnode = node.GetFirstChild("pcount");
    if (subnode != NULL) {
        sim.SetMaxPartCount((int)Strings::cdble(subnode->Data()));
    }

    // Read predicted maximum M0 value.
    subnode = node.GetFirstChild("maxm0");
    if (subnode != NULL) {
        sim.SetMaxM0(Strings::cdble(subnode->Data())*1.0e6); // Convert from #/cm3 to #/m3.
    }

    // Read predictor-corrector relaxation parameter.
    subnode = node.GetFirstChild("relax");
    if (subnode != NULL) {
        solver.SetUnderRelaxCoeff(Strings::cdble(subnode->Data()));
    }
}


/*!
 * @brief           Helper function to determine if coag kernels are compatible
 *
 * This is here to check if, upon loading an ensemble, that the old and new
 * coagulation kernels are compatible. There are four cases:
 * 1: WPM -> WPM
 * 2: WPM -> DSA
 * 3: DSA -> DSA
 * 4: DSA -> WPM
 * Of these, only number 2 doesn't work. This is because all of the DSA coag
 * kernel calculations assume particles have equal weight.
 *
 * @param old_id    ProcessType ID of binary file's coagulation mechanism
 * @param this_id   ProcessType ID of this simulation
 * @return          Boolean indicating if kernels are compatible
 */
bool checkCoagulationKernel(int old_id, int this_id) {
    bool ans(false);
    if (old_id == Sweep::Processes::Weighted_Additive_Coagulation_ID
            || old_id == Sweep::Processes::Weighted_Constant_Coagulation_ID
            || old_id == Sweep::Processes::Weighted_Transition_Coagulation_ID
            || old_id == Sweep::Processes::Transition_Coagulation_ID) {
        if (old_id == this_id) {
            // This is fine, we've got the same kernel
            ans = true;
        } else {
            if (this_id == Sweep::Processes::Additive_Coagulation_ID
                    || this_id == Sweep::Processes::Constant_Coagulation_ID
                    || this_id == Sweep::Processes::Transition_Coagulation_ID) {
                // Uh oh.. we've gone WPM -> DSA
                ans = false;
            } else {
                // This is fine, we've gone WPMi -> WPMj
                ans = true;
            }
        }
    }
    return ans;
}

/*!
 * @brief           Reads an ensemble .ens file.
 *
 * This function is to be used when loading a *.ens file to initialise
 * particles in the system at t=0. Files are specified in the <reactor> block
 * of the mops.inx file in the following manner:
 *     <population>
 *          <file>silica-fm(0)-SP(100).ens</file>
 *          <m0>3.0e15</m0>
 *     </population>
 *
 */
Sweep::PartPtrList readEnsembleFile(
        const string fname,
        const Sweep::Mechanism &smech, 
	std::vector<unsigned int> &particle_numbers_list) {
    // Open the save point file.
    ifstream fin;
    fin.open(fname.c_str(), ios_base::in | ios_base::binary);

    // Create empty particle pointer list
    Sweep::PartPtrList particles;

    if (fin.good()) {

        // First read-in and check the particle model
        Sweep::Mechanism::ParticleModel filemodel;      // model to be loaded
        filemodel.Deserialize(fin);
        if (filemodel.AggModel() == smech.AggModel()) {
            std::cout << "parser: correct particle model found!\n";
        } else {
            throw runtime_error("Wrong particle model specified in sweep.xml!");
        }

        // For binary-tree based particles, check that the tree status is the same.
        if (filemodel.WriteBinaryTrees() == smech.WriteBinaryTrees()) {
            std::cout << "parser: binary tree status okay!\n";
        } else {
            throw runtime_error("Wrong binary tree output flag in mops.inx.");
        }

        // Now, check the coagulation kernel, assuming only one coagulation process
        int id(0);
        fin.read(reinterpret_cast<char*>(&id), sizeof(id));
        if (checkCoagulationKernel(id, smech.Coagulations()[0]->ID())) {
            std::cout << "parser: coagulation process okay!\n";
        } else {
            throw runtime_error("Conflicting coagulation kernel specification!");
        }

        // Now it's time to load the particle ensemble.
        Sweep::Ensemble fileensemble;
        fileensemble.Deserialize(fin, smech);

        // Begin adding the particles of fileensemble to the simulation's ensemble
        Sweep::Particle *sp;
        for (unsigned int i(0); i != fileensemble.Count(); ++i) {
            // note we need to use Clone, because fileensemble will be deleted once
            // we leave the scope of this function.
            sp = fileensemble.At(i)->Clone();
            particles.push_back(sp);
        }

	// Read hybrid model particle-number info into list 
	if (fileensemble.GetHybridThreshold() > 0)
	{
	    particle_numbers_list.resize(fileensemble.GetHybridThreshold(), 0);
	    for (unsigned int i = 0; i < fileensemble.GetHybridThreshold(); ++i)
	    {
	        particle_numbers_list[i] = fileensemble.NumberAtIndex(i);
	    }
	}

        // Close the input file.
        fin.close();

    } else {
        // Throw error if the output file failed to open.
        throw runtime_error("Failed to open ensemble file "
                            "input (Mops, Simulator::readEnsembleFile).");
    }

    return particles;
}





void readInitialPopulation(
        const CamXML::Element &subnode,
        const Mops::Mechanism &mech,
        Mops::Mixture &mix) {
    // Initialise some storage
    Sweep::PartPtrList fileParticleList, inxParticleList, allParticleList;
    CamXML::Element *m0node, *fnode;
    double initialM0(0.0);

    // Find the overall number density represented by the population
    m0node = subnode.GetFirstChild("m0");
    if (!m0node) {
        throw std::runtime_error("Must specify M0 for initial population"
                                 "(Mops::Settings_IO::readInitialPopulation)");
    }
    initialM0 = std::atof((m0node->Data()).c_str());
    if (initialM0 <= 0.) {
        throw std::runtime_error("m0 for initial particle population may not be negative"
                                 "(Mops::Settings_IO::readInitialPopulation)");
    }

    // Check if a binary file has been specified for particles..
    fnode = subnode.GetFirstChild("file");
    if (fnode) {
        std::string filename;
        filename = fnode->Data();
        std::cout << "parser: binary file " << filename << " specified for input.\n";
	std::vector<unsigned int> particle_numbers;
	particle_numbers.resize(0, 0);
        fileParticleList = readEnsembleFile(filename, mech.ParticleMech(), particle_numbers);
	// Initialise and store hybrid particle-number info for the cell
	// Note this only works if the primary particle model is univariate
	if (particle_numbers.size() > 0)
	{
            mix.Particles().SetHybridThreshold(particle_numbers.size());
            mix.Particles().InitialiseParticleNumberModel();
            mix.Particles().InitialiseDiameters(mix.ParticleModel()->Components()[0]->MolWt(),
            mix.ParticleModel()->Components()[0]->Density()); 
            for (unsigned int i = 0; i < mix.Particles().GetHybridThreshold(); ++i)
            {
                mix.Particles().UpdateTotalsWithIndex(i, particle_numbers[i]);
                mix.Particles().UpdateNumberAtIndex(i, particle_numbers[i]);
                mix.Particles().UpdateTotalParticleNumber(particle_numbers[i]);
            }
	}
    }

    // Choose which file reader to use depending if particle distributions are required
    vector<CamXML::Element*> subitems;
    subnode.GetChildren("particle", subitems); 
    string dflag;
    for (vector<CamXML::Element*>::iterator j = subitems.begin(); j != subitems.end(); ++j) {
	// Get component ID.
	dflag = (*j)->GetAttributeValue("usedistribution");
    }

    // Now read in the list of particles and sum up their statistical weights
    if (dflag == "true") {
	// Use detailed reader to create particles based on distributions specified in input file
	inxParticleList = Settings_IO::ReadInitialParticlesDetailed(subnode, mech.ParticleMech()); 
    }
    else 
	inxParticleList = Settings_IO::ReadInitialParticles(subnode, mech.ParticleMech());

    // Join the particle lists
    allParticleList.splice(allParticleList.end(), fileParticleList);
    allParticleList.splice(allParticleList.end(), inxParticleList);
    Sweep::PartPtrList::const_iterator it = allParticleList.begin();
    const Sweep::PartPtrList::const_iterator itEnd = allParticleList.end();

    // Get the total weights of *all* initialised particles
    double weightSum = 0;
    while(it != itEnd) {
        weightSum += (*it++)->getStatisticalWeight();
    }
    // Add particle-number particles to the total weight
    weightSum += mix.Particles().GetTotalParticleNumber(); 

    mix.SetParticles(allParticleList.begin(), allParticleList.end(), initialM0 / weightSum);
}


Mops::Mixture* readMixture(
        const CamXML::Element &node,
        const Mops::Mechanism &mech,
        const unsigned int pcount,
        const double maxm0,
        std::string &name
        ) {
    std::vector<CamXML::Element*> nodes;
    const CamXML::Element *subnode;
    std::vector<CamXML::Element*>::const_iterator i;
    const CamXML::Attribute *attr;

    // Get the ID of the mixture
    attr = node.GetAttribute("id");
    if (attr) name = attr->GetValue();

    // Create the mixture
    Mixture *mix = NULL;
    mix = new Mops::Mixture(mech.ParticleMech());

    // Read some values
    double val(0.0);
    subnode = node.GetFirstChild("temperature");
    if (subnode != NULL)
        val = readTemperature(*subnode);
    else
        throw std::runtime_error("No temperature specified for the mixture (::readMixture).");
    mix->GasPhase().SetTemperature(val);

    subnode = NULL;
    subnode = node.GetFirstChild("pressure");
    if (subnode != NULL)
        val = readPressure(*(node.GetFirstChild("pressure")));
    else
        throw std::runtime_error("No temperature specified for the mixture (::readMixture).");
    mix->GasPhase().SetPressure(val);

    // Look for components...
    fvector molefracs(mech.GasMech().SpeciesCount(), 0.0);
    node.GetChildren("component", nodes);
    for (i=nodes.begin(); i!=nodes.end(); ++i) {
        // Check the ID attribute.
        attr = (*i)->GetAttribute("id");
        if (attr != NULL) {
            std::string str = attr->GetValue();
            // Get the species index.
            int j = mech.GasMech().FindSpecies(attr->GetValue());

            // Set the species mole fraction.
            if (j >= 0) {
                molefracs[j] = Strings::cdble((*i)->Data()); // from inx input
            } else {
                throw std::runtime_error("No initial condition for species " + str +
                                    " (::readMixture).");
            }
        } else {
            throw std::runtime_error("Initial condition must have ID "
                                " (::readMixture).");
        }
    }
    mix->GasPhase().SetFracs(molefracs);

    // Get viscosity
    subnode = NULL;
    subnode = node.GetFirstChild("viscosity");
    if(subnode != NULL) {
        std::string m = subnode->GetAttributeValue("model");
        if (m.compare("air") == 0) {
            mix->GasPhase().SetViscosityModel(Sprog::iAir);
        } else if (m.compare("chapman-enskog") == 0) {
            // Stop program if not enough transport data
            mix->GasPhase().checkForTransportData();
            mix->GasPhase().SetViscosityModel(Sprog::iChapmanEnskog);
        } else if (m.compare("argon") == 0) {
            mix->GasPhase().SetViscosityModel(Sprog::iArgon);
        } else if (m.compare("hydrogen") == 0) {
            mix->GasPhase().SetViscosityModel(Sprog::iHydrogen);
        } else {
            throw std::runtime_error("Unrecognised viscosity model"
                    " (::readMixture).");
        }
    }

    // Now load some particles
    mix->Particles().SetHybridThreshold(mech.ParticleMech().GetHybridThreshold());
    mix->Particles().Initialise(pcount);
    mix->Particles().InitialiseParticleNumberModel();
    mix->Reset(maxm0);
    subnode = NULL; subnode = node.GetFirstChild("population");
    if (subnode != NULL) {
        if (mech.ParticleMech().ComponentCount() < 1) {
            throw std::runtime_error("Error loading initial population, no particle mechanism loaded. \n"
                    "Ensure a particle solver and sweep file are used (::readMixture).");
        } else {
            readInitialPopulation(*subnode, mech, *mix);
        }
    }

    return mix;
}

Mops::PSR *const readPSR(
        const CamXML::Element &node,
        const Mechanism &mech,
        std::string &name) {

    Mops::PSR *reac = NULL;
    const CamXML::Element *subnode;
    std::vector<CamXML::Element*> nodes;
    std::vector<CamXML::Element*>::const_iterator i;
    const CamXML::Attribute *attr;
    std::string str;

    // Get the ID of the reactor
    attr = node.GetAttribute("id");
    if (attr) name = attr->GetValue();

    // First read the system type.  This is where we initialise
    // the reactor type.
    attr = node.GetAttribute("type");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.compare("psr") == 0) {
            reac = new Mops::PSR(mech);
        } else
            throw std::runtime_error("Only PSR presently supported for networking."
                    " (::readPSR)");
    }

    // Now check for constant temperature.  If not constant temperature
    // then the reactor is solved using the adiabatic energy equation.
    attr = node.GetAttribute("constt");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.compare("true") == 0) {
            // This is a constant temperature reactor.
            reac->SetEnergyEquation(Reactor::ConstT);
        } else {
            // Use the adiabatic energy equation.
            reac->SetEnergyEquation(Reactor::Adiabatic);
        }
    } else {
        // The constt attribute is undefined, so use adiabatic
        // energy equation.
        reac->SetEnergyEquation(Reactor::Adiabatic);
    }
    // Now check if particle terms should be included in the EB.
    // This is only relevant with adiabatic energy equation.
    if (reac->EnergyEquation() == Reactor::Adiabatic)
    {
        attr = node.GetAttribute("includeParticleTerms");
        if (attr != NULL) {
            str = attr->GetValue();
			if (str.compare("true") == 0) {
				if (reac->Mech()->ParticleMech().Components()[0]->Name().compare("Rutile") == 0)
				{
					reac->Mech()->ParticleMech().SetParticleSpeciesIndex(28);
					reac->SetIncludeParticles();
				}
				else
					throw std::runtime_error("Only Rutile primary particle supported for particle contributions to energy balance."
					" (::readPSR)");
            }
        }
    }

    // Now check for constant volume.
    attr = node.GetAttribute("constv");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.compare("true") == 0) {
            // This is a constant volume reactor.
            reac->SetConstV();
        } else {
            // This is a constant pressure reactor.
            reac->SetConstP();
        }
    } else {
        // The constv attribute is undefined, so use constant
        // pressure reactor.
        reac->SetConstP();
    }


     // Now initialise for surface area (surface kinetics purposes).
    attr = node.GetAttribute("surfaceArea");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.length() > 0) {
            stringstream AREA(str);
            double ar;
            AREA >> ar;
            // This is a surface area.
            reac->SetArea(ar);
        } else {
            // This is a default value for area.
            reac->SetArea(0.0);
        }
    } else {
        // The area attribute is undefined, so set it to default
        reac->SetArea(0.0);
    }

    // Now initialise for volume.
    attr = node.GetAttribute("volume");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.length() > 0) {
            stringstream VOLUME(str);
            double vol;
            VOLUME >> vol;
            // This is a volume.
            reac->SetVolume(vol);
        } else {
            // This is a a default value for volume.
            reac->SetVolume(1.0);
        }
    } else {
        // The volume attribute is undefined, so set it to default
        reac->SetVolume(1.0);
    }

    // Get birth and death types
    attr = NULL;
    attr = node.GetAttribute("btype");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.compare("stochastic") == 0)
            reac->SetInflowType(Sweep::Processes::BirthProcess::iStochastic);
        else if (str.compare("continuous") == 0)
            reac->SetInflowType(Sweep::Processes::BirthProcess::iContinuous);
        else
            throw std::runtime_error("Unrecognised birth process type. (::readPSR)");
    }

    attr = NULL;
    attr = node.GetAttribute("dtype");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.compare("cdelete") == 0)
            reac->SetOutflowType(Sweep::Processes::DeathProcess::iContDelete);
        else if (str.compare("cmove") == 0)
            reac->SetOutflowType(Sweep::Processes::DeathProcess::iContMove);
        else if (str.compare("rescale") == 0)
            reac->SetOutflowType(Sweep::Processes::DeathProcess::iContRescale);
        else if (str.compare("adaptive") == 0)
            reac->SetOutflowType(Sweep::Processes::DeathProcess::iContAdaptive);
        else if (str.compare("sdelete") == 0)
            reac->SetOutflowType(Sweep::Processes::DeathProcess::iStochDelete);
        else if (str.compare("smove") == 0)
            reac->SetOutflowType(Sweep::Processes::DeathProcess::iStochMove);
        else
            throw std::runtime_error("Unrecognised death process type. (::readPSR)");
    }

    // TEMPERATURE GRADIENT PROFILE.

    node.GetChildren("dTdt", nodes);
    for (i=nodes.begin(); i!=nodes.end(); ++i) {
        // Get the start time attribute.
        attr = (*i)->GetAttribute("startt");
        if (attr != NULL) {
            Sweep::Maths::Linear fun;
            fun.SetParam(1, Strings::cdble((*i)->Data()));
            reac->Add_dTdt(Strings::cdble(attr->GetValue()), fun);
        }
    }

    //TODO: read PSR-specific settings
    if (reac->SerialType() == Serial_PSR) {

            // Read the residence time first so that the inflow's rate is
            // correctly set.
            subnode = node.GetFirstChild("restime");
            if (subnode != NULL) {
                double tau = Strings::cdble(subnode->Data());
                reac->SetResidenceTime(tau);
            } else {
                throw std::runtime_error("Residence time not specified for PSR!"
                        " (::readPSR)");
            }
        }

    return reac;
}

// Reads the reactor initial settings from the given XML node.
Reactor *const readReactor(const CamXML::Element &node,
                                        const Mechanism &mech,
                                        const unsigned int max_particle_count,
										const double maxM0,
										const Mops::Simulator &sim)
{
    Reactor *reac = NULL;
    const CamXML::Element *subnode;
    std::vector<CamXML::Element*> nodes;
    std::vector<CamXML::Element*>::const_iterator i;
    const CamXML::Attribute *attr;
    std::string str;

    // SYSTEM TYPE AND ATTRIBUTES.

    // First read the system type.  This is where we initialise
    // the reactor type.
    attr = node.GetAttribute("type");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.compare("batch") == 0) {
            // This is a batch (const. P) reactor.
            reac = ReactorFactory::Create(Serial_Batch, mech);
        } else if (str.compare("psr") == 0) {
            // This is a perfectly-stirred reactor.
            reac = ReactorFactory::Create(Serial_PSR, mech);
        } else if (str.compare("shocktube") == 0) {
            // This is a shocktube (const. V batch).
            throw std::runtime_error("Shocktube not currently implemented"
                                " (Mops::Settings_IO::readReactor).");
        } else {
            // Default reactor is a batch reactor.
            reac = ReactorFactory::Create(Serial_Batch, mech);
        }
    }

    // Now check for constant temperature.  If not constant temperature
    // then the reactor is solved using the adiabatic energy equation.
    attr = node.GetAttribute("constt");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.compare("true") == 0) {
            // This is a constant temperature reactor.
            reac->SetEnergyEquation(Reactor::ConstT);
        } else {
            // Use the adiabatic energy equation.
            reac->SetEnergyEquation(Reactor::Adiabatic);
        }
    } else {
        // The constt attribute is undefined, so use adiabatic
        // energy equation.
        reac->SetEnergyEquation(Reactor::Adiabatic);
    }
    // Now check if particle terms should be included in the EB.
    // This is only relevant with adiabatic energy equation.
    if (reac->EnergyEquation() == Reactor::Adiabatic)
    {
        attr = node.GetAttribute("includeParticleTerms");
        if (attr != NULL) {
            str = attr->GetValue();
            if (str.compare("true") == 0) {
				if (reac->Mech()->ParticleMech().Components()[0]->Name().compare("Rutile") == 0)
				{
					reac->Mech()->ParticleMech().SetParticleSpeciesIndex(28);
					reac->SetIncludeParticles();
				}
				else
					throw std::runtime_error("Only Rutile primary particle supported for particle contributions to energy balance."
					" (Mops::Settings_IO::readReactor)");
			}
        }
    }

    // Now check for constant volume.
    attr = node.GetAttribute("constv");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.compare("true") == 0) {
            // This is a constant volume reactor.
            reac->SetConstV();
        } else {
            // This is a constant pressure reactor.
            reac->SetConstP();
        }
    } else {
        // The constv attribute is undefined, so use constant
        // pressure reactor.
        reac->SetConstP();
    }
    
	
	 // Now initialise for surface area (surface kinetics purposes).
    attr = node.GetAttribute("surfaceArea");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.length() > 0) {
			stringstream AREA(str);
			double ar;
			AREA >> ar;
            // This is a surface area.
            reac->SetArea(ar);
        } else {
            // This is a default value for area. 
            reac->SetArea(0.0);
        }
    } else {
        // The area attribute is undefined, so set it to default
        reac->SetArea(0.0);
    }
    
	// Now initialise for volume.
    attr = node.GetAttribute("volume");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.length() > 0) {
			stringstream VOLUME(str);
			double vol;
			VOLUME >> vol;
            // This is a volume. 
            reac->SetVolume(vol);
        } else {
            // This is a a default value for volume. 
            reac->SetVolume(1.0);
        }
    } else {
        // The volume attribute is undefined, so set it to default
        reac->SetVolume(1.0);
    }

    // REACTOR INITIAL CONDITIONS.

    // Create a new Mixture object.
    std::string id;
    Mixture *mix = readMixture(node, mech, max_particle_count, maxM0, id);
    reac->Fill(*mix);

    // Investigate whether doubling is activated
    subnode = node.GetFirstChild("doubling");
    if (subnode != NULL) {
        std::string str_enable = subnode->GetAttributeValue("enable");
        if (str_enable.compare("false") == 0) {
            std::cout << "sweep: Warning! doubling is turned off by user.\n";
            reac->Mixture()->Particles().SetDoubling(false);
        } else reac->Mixture()->Particles().SetDoubling(true);
    }

    // TEMPERATURE GRADIENT PROFILE.

    node.GetChildren("dTdt", nodes);
    for (i=nodes.begin(); i!=nodes.end(); ++i) {
        // Get the start time attribute.
        attr = (*i)->GetAttribute("startt");
        if (attr != NULL) {
            Sweep::Maths::Linear fun;
            fun.SetParam(1, Strings::cdble((*i)->Data()));
            reac->Add_dTdt(Strings::cdble(attr->GetValue()), fun);
        }
    }

    // PSR SPECIFIC SETTINGS.
    fvector molefracs;
    if (reac->SerialType() == Serial_PSR) {

        Mops::PSR* psr = dynamic_cast<PSR*>(reac);

        // Read the residence time first so that the inflow's rate is
        // correctly set.
        subnode = node.GetFirstChild("restime");
        if (subnode != NULL) {
            double tau = Strings::cdble(subnode->Data());
            psr->SetResidenceTime(tau);
        } else
            throw std::runtime_error("Res. time must be defined for a PSR "
                    "(Mops, Settings_IO::readReactor)");

        // Initialise an outflow
        FlowStream *outf = new FlowStream(mech);
        outf->ConnectInflow(*psr);
        psr->SetOutflow(*outf);

        // Read inflow
        subnode = node.GetFirstChild("inflow");
        if (subnode != NULL) {
            FlowStream *inf = new FlowStream(mech);

            // Create a new Mixture object for inflow. Note that the inflow takes
            // a copy of the one loaded from XML.
            mix = readMixture(*subnode, mech, max_particle_count, maxM0, id);
            inf->SetConditions(*mix);

            psr->SetInflow(*inf);

        } else {
            throw std::runtime_error("Inflow conditions must be defined for a PSR "
                                "(Mops, Settings_IO::readReactor)");
        }
    }

    // Check if particle terms are to be included in the energy balance
    if (reac->IncludeParticles()) {
        reac->Mixture()->SetIsAdiabaticFlag(true);
    }
    else {
        reac->Mixture()->SetIsAdiabaticFlag(false);
    }

    return reac;
}


void readNetwork(
        const CamXML::Element &root,
        const Mops::Simulator &sim,
        const Mops::Mechanism &mech,
        Mops::ReactorNetwork &network) {

    // Set up some storage
    CamXML::Document doc;
    CamXML::Element* node;
    std::vector<CamXML::Element*> nodes, inodes;
    std::vector<CamXML::Element*>::const_iterator i, k;
    std::vector<flowNode> flows;
    std::map<unsigned int, std::string> reac_order;
    std::string id, mid, iid;
    const CamXML::Attribute *attr;

    // Look for mixtures first
    root.GetChildren("mixture", nodes);
    unsigned int j(0u);
    if (nodes.size() > 0) {
        for (i=nodes.begin(); i!=nodes.end(); ++i) {
            Mops::Mixture* mix = readMixture(*(*i), mech, sim.MaxPartCount(), sim.MaxM0(), id);

            // Generate a name if one wasn't loaded.
            if (id.size() < 1) {
                std::stringstream sstm;
                sstm << "mix" << j;
                id = sstm.str();
            }
            network.AddMixture(id, *mix);
            j++;
        }
    }


    // Now look for reactors
    nodes.clear();
    root.GetChildren("reactor", nodes);
    if (nodes.size() > 0) {
        for (i=nodes.begin(); i!=nodes.end(); ++i) {
            id = "";
            Mops::PSR* reac = readPSR(*(*i), mech, id);

            // Generate a name if one wasn't loaded.
            if (id.size() < 1) {
                std::stringstream sstm;
                sstm << "reac" << j;
                id = sstm.str();

            }

            network.AddReactor(id, *reac);

            // Check if it has an order attribute
           attr = (*i)->GetAttribute("order");
           if (attr != NULL)
               reac_order[(unsigned int) atoi(attr->GetValue().c_str())] = id;


            // Fill the reactor with its intial state
            node = (*i)->GetFirstChild("initial");
            if (node != NULL) {
                mid = node->Data();
                network.FillReactor(id, mid);
            }

            // Check if doubling should be turned-on
            node = NULL;
            node = (*i)->GetFirstChild("doubling");
            if (node != NULL) {
                std::string str_enable = node->GetAttributeValue("enable");
                if (str_enable.compare("false") == 0) {
                    std::cout << "sweep: Warning! doubling is turned off by user.\n";
                    reac->Mixture()->Particles().SetDoubling(false);
                } else reac->Mixture()->Particles().SetDoubling(true);
            }
            j++;

            // Store the names of the flow
            (*i)->GetChildren("flow", inodes);
            if (inodes.size() < 1) throw std::runtime_error("No flows defined for reactor! (::readNetwork)");
            for (k=inodes.begin(); k!=inodes.end(); ++k) {
                flowNode fn = flowNode();
                std::string val;
                double fval(0.0);

                // Get the mixture id
                fn.mMixture = (*k)->Data();
                fn.mReactor = id;

                // Get the direction
                attr = NULL;
                attr = (*k)->GetAttribute("dir");
                if (attr != NULL) {
                    val = attr->GetValue().c_str();
                    if (val.compare("in")==0) fn.mIn = true;
                    else if (val.compare("out")==0) fn.mIn = false;
                    else {
                        throw std::runtime_error(std::string("Unknown flow direction: ")+val);
                    }
                }

                // Get the flowrate
                attr = NULL;
                attr = (*k)->GetAttribute("frac");
                if (attr != NULL) {
                    fval = (double) std::atof(attr->GetValue().c_str());
                    if (fval > 0.0) fn.mFrac = fval;
                    else {
                        throw std::runtime_error(std::string("Flow fraction must be positive!"));
                    }
                }

                flows.push_back(fn);
            }
            // Check if particle terms are to be included in the energy balance
            if (reac->IncludeParticles()) {
                reac->Mixture()->SetIsAdiabaticFlag(true);
            }
            else {
                reac->Mixture()->SetIsAdiabaticFlag(false);
            }
        }
        // Now that the reactors have been loaded, make node connections.
        for (std::vector<flowNode>::iterator it=flows.begin(); it!=flows.end(); ++it) {
            network.Connect(it->mReactor, it->mMixture, it->mIn, it->mFrac);
        }
    } else {
        throw std::runtime_error("No reactor blocks found in the settings file."
                " (::readNetwork)");
    }

    // Now that we've loaded the reactors, we need to initialise the path
    // in which the reactors are solved.
    if (reac_order.size() == network.ReactorCount()) {
        // First try using the manually specified order
        network.SpecifyPaths(reac_order);
    } else {
        // Incomplete data... try guessing.
        network.FindPaths();
    }

}


/*!
@param[in]     node    XML node for input streams
@param[in]     solver  the solver
*/
void ReadLOIStatus(const CamXML::Element &node, Solver &solver)
{
    //const CamXML::Element *subnode;
    const CamXML::Attribute *attr;
    std::vector<CamXML::Element*> nodes;
    std::vector<CamXML::Element*>::const_iterator i;
    std::string str;
    double m_comp;


    // See if LOI is enabled or not.
    attr = node.GetAttribute("enable");
    if (attr != NULL) {
        str = attr->GetValue();
        if (str.compare("true") == 0) {
            // LOI reduction should take place
            solver.SetLOIStatusTrue();
            }
        else
            solver.SetLOIStatusFalse();
    }

    //Read in the LOI cutoff value
    attr = node.GetAttribute("LOIComp");
    if (attr!= NULL){
        str = attr->GetValue();
        //Set the LOI cutoff value
        m_comp = Strings::cdble(str);
        solver.SetLOICompValue(m_comp);    
    }

    //Need to read in the kept species here...
    if (solver.GetLOIStatus() == true){
        node.GetChildren("Kept_Spec", nodes);
        for (i = nodes.begin(); i != nodes.end(); i++){
            attr = (*i)->GetAttribute("id");
            if (attr != NULL) {
                str = attr->GetValue();
                solver.AddKeptSpecies(str);
            }
        }
    }

}

// Reads simulation postprocess parameters from given XML node.
void readPostprocess(const CamXML::Element &node, Simulator &sim, Mechanism &mech)
{
    const CamXML::Element *subnode;
    // Mass spectra
    subnode = node.GetFirstChild("Mass_spectra");
    if (subnode != NULL) {

        std::string str_enable = subnode->GetAttributeValue("enable");
        if (str_enable.compare("true") == 0) {
            sim.SetMassSpectra(true);
        } else {
            sim.SetMassSpectra(false);
        }

        // Read tabular console output parameters.
        subnode = subnode->GetFirstChild("ensemble");
        if (subnode != NULL) {

            str_enable = subnode->GetAttributeValue("enable");
            if (str_enable.compare("true") == 0) {
                sim.SetMassSpectraEnsemble(true);
            } else {
                sim.SetMassSpectraEnsemble(false);
                const CamXML::Element *xmer = subnode->GetFirstChild("xmer");
                const CamXML::Element *frag = subnode->GetFirstChild("fragment");
                //set the xmer
                sim.SetMassSpectraXmer((int)Strings::cdble(xmer->Data()));
                //set wheter considering fragmentation
                str_enable=(frag->Data());
                if (str_enable.compare("true") == 0) {
                    // By default, frag is turned off
                    sim.SetMassSpectraFrag(true);
                }
            }
        }
        else {
                std::cout << "sweep: Warning! no option for mass spectra has been specified, the default one will be used (Mass spectra for whole ensamble.\n";
                sim.SetMassSpectraEnsemble(true);
        }
    } else {
        throw std::runtime_error("please specify the options about Mass spectra"
                                " (Mops, Settings_IO::readPostprocess()).");
    }
}

// Reads simulation output parameters from given XML node.
void readOutput(const CamXML::Element &node, Simulator &sim, Mechanism &mech)
{
    const CamXML::Element *subnode;
    std::vector<CamXML::Element*> nodes;
    std::vector<CamXML::Element*>::const_iterator i;
    const CamXML::Attribute *attr;
    std::string str;

    // CONSOLE OUTPUT.

    subnode = node.GetFirstChild("console");
    if (subnode != NULL) {
        // Check the console interval.
        attr = subnode->GetAttribute("interval");
        if (attr != NULL) {
            sim.SetConsoleInterval((unsigned int)Strings::cdble(attr->GetValue()));
        }

        // Read tabular console output parameters.
        subnode = subnode->GetFirstChild("tabular");
        if (subnode != NULL) {
            subnode->GetChildren("column", nodes);
            for (i=nodes.begin(); i!=nodes.end(); ++i) {
                sim.AddConsoleVariable((*i)->Data());
            }
        } else {
            throw std::runtime_error("No tabular data defined for console output "
                                "(Mops, Settings_IO::readOutput)");
        }
    } else {
        throw std::runtime_error("No console output defined "
                            "(Mops, Settings_IO::readOutput)");
    }


    // FILE OUTPUT.

    // Read the output file name.
    subnode = node.GetFirstChild("filename");
    if (subnode != NULL) {
        sim.SetOutputFile(subnode->Data());
    } else {
        throw std::runtime_error("No output file name defined "
                            "(Mops, Settings_IO::readOutput)");
    }


    // Check if entire particle binary trees should be dumped (default off)
    // Relevant for Silica and PAHPP models.
    subnode = node.GetFirstChild("writebintree");
    if (subnode != NULL) {
        std::string str_enable = subnode->GetAttributeValue("enable");
        if (str_enable.compare("true") == 0) {
            std::cout << "sweep: Warning! Writing full particle binary trees.\n";
            mech.ParticleMech().SetWriteBinaryTrees(true);
        } else {
            mech.ParticleMech().SetWriteBinaryTrees(false);
        }
    } else {
        mech.ParticleMech().SetWriteBinaryTrees(false);
    }

    // STATISTICAL BOUNDARIES OF OUTPUT

    // Bound parameters
    Sweep::PropID pid = Sweep::iDcol;
    double lower = 0.0; // Default value of the lower bound
    double upper = 1e30; // Default value of the upper bound
    pid = Sweep::iDcol; // Default value of Property ID is iDcol
    // Read the statistical boundaries
    subnode = node.GetFirstChild("statsbound");
    if (subnode != NULL) {
        Settings_IO::ReadStatsBound(*subnode, pid, lower, upper);
    }
    // Set statistical bounds to simulator
    sim.SetOutputStatBoundary(pid, lower, upper);

    // POVRAY OUTPUT.

    // Read POVRAY particle tracking output.
    subnode = node.GetFirstChild("ptrack");
    if (subnode != NULL) {
        std::string str_enable = subnode->GetAttributeValue("enable");
        if (str_enable.compare("true") == 0) {
            std::string str_ptcount = subnode->GetAttributeValue("ptcount");
            sim.SetParticleTrackCount((unsigned int)Strings::cdble(str_ptcount));

            // Also need to ensure full binary trees are written for certain
            // particle models.
            mech.ParticleMech().SetWriteBinaryTrees(true);

        } else if (str_enable.compare("false") == 0) {
            sim.SetParticleTrackCount(0);
        } else {
            throw std::runtime_error("Unknown ptrack enabling keyword in MOPS input "
                                "(Mops, Settings_IO::readOutput)");
        }
    } else {
        // Set number of particle tracking to zero if ptrack is not found.
        sim.SetParticleTrackCount(0);
    }

    // FLUX ANALYSIS OUTPUT

    // Read the flux analysis settings for postprocessing.
    subnode = node.GetFirstChild("fluxanalysis");
    if (subnode != NULL) {
        std::string str_enable = subnode->GetAttributeValue("enable");
        if (str_enable.compare("true") == 0) {
            std::vector<CamXML::Element*> elem_nodes;
            subnode->GetChildren("element", elem_nodes);
            for (unsigned int i = 0; i < elem_nodes.size(); i++) {
                sim.AddFluxElement(elem_nodes.at(i)->GetAttributeValue("id"));
            }
        } else if (str_enable.compare("false") == 0) {
            sim.ClearFluxElements();
        }
    } else {
        sim.ClearFluxElements();
    }

    // PAH KMC POSTPROCESSING

    if (mech.ParticleMech().AggModel() == Sweep::AggModels::PAH_KMC_ID){
        // currently only the data for mass spectra is provided,
        // but the capacity will be extended according to the purpose of the project
        subnode = node.GetFirstChild("postprocess");
        if (subnode != NULL) {
            readPostprocess(*subnode, sim, mech);
        } else {
            // mass spectra data will not be generated, this should be only for the old test files
            sim.SetMassSpectra(false);
        }
    }
}

} // anonymous namespace



// VERSION 2 XML SETTINGS FILES.

// Reads a new-format XML settings file.
Reactor *const Settings_IO::LoadFromXML(const std::string &filename,
                                        Mops::Reactor *reac,
                                        std::vector<TimeInterval> &times,
                                        Simulator &sim, Solver &solver,
                                        Mechanism &mech)
{
    CamXML::Document doc;
    const CamXML::Element *root, *node;

    // We need a blank reactor object to create.  If the pointer
    // passed to this function is valid then it needs to be
    // deleted.
    if (reac != NULL) delete reac;
    reac = NULL;

    // Attempt to load XML document.
    if (doc.Load(filename) == 0) {
        // Get the root element.
        root = doc.Root();

        // Check the mops file version.
        const CamXML::Attribute *attr = root->GetAttribute("version");

        if ((attr!=NULL) && (attr->GetValue() != "2")) {
            throw std::runtime_error("Settings file has wrong version number (expecting 2)"
                                " (Mops::Settings_IO::LoadFromXML).");
        }

        // GLOBAL SETTINGS.

        readGlobalSettings(*root, sim, solver);

        // OUTPUT SETTINGS.
        // wjm34: read output settings before reactor, so we can check if ensemble/g.p.
        // files are consistent with some more simulation settings.

        node = root->GetFirstChild("output");
        if (node != NULL) {
            readOutput(*node, sim, mech);
        } else {
            throw std::runtime_error("Settings file does not contain output"
                                " information (Mops::Settings_IO::LoadFromXML).");
        }

        // REACTOR.

        node = root->GetFirstChild("reactor");
        if (node != NULL) {
			reac = readReactor(*node, mech, sim.MaxPartCount(), sim.MaxM0(), sim);
        } else {
            throw std::runtime_error("Settings file does not contain a reactor definition"
                                " (Mops::Settings_IO::LoadFromXML).");
        }

	
        // MECHANISM REDUCTION.

        node = root->GetFirstChild("LOI");
        if (node!= NULL) {
            ReadLOIStatus(*node, solver);
        }

        // TIME INTERVALS.

        node = root->GetFirstChild("timeintervals");
        if (node != NULL) {
            readTimeIntervals(*node, times);
        } else {
            throw std::runtime_error("Settings file does not contain a time intervals"
                                " (Mops::Settings_IO::LoadFromXML).");
        }

    }

    return reac;
}


/*
 *\param[in]        node        XML node containing time interval data
 *\param[out]       times       vector of times intervals for simulation
 *
 *\exception    std::runtime_error  Missing data
 */
void Settings_IO::readTimeIntervals(const CamXML::Element &node,
                                           std::vector<TimeInterval> &times)
{
    const CamXML::Element *subnode;
    std::vector<CamXML::Element*> nodes;
    std::vector<CamXML::Element*>::const_iterator i;
    const CamXML::Attribute *attr;
    std::string str;

    // Read default splits per step.
    unsigned int splits = 1;
    attr = node.GetAttribute("splits");
    if (attr != NULL) {
        splits = (unsigned int)Strings::cdble(attr->GetValue());
    }


    // Read start time (first interval).
    TimeInterval *ti = NULL;
    subnode = node.GetFirstChild("start");
    if (subnode != NULL) {
        ti = new TimeInterval();
        ti->SetStartTime(Strings::cdble(subnode->Data()));
    } else {
        throw std::runtime_error("Time intervals must have start time defined "
                            "(Mops, Settings_IO::readTimeIntervals)");
    }

    node.GetChildren("time", nodes);
    for (i=nodes.begin(); i!=nodes.end(); ++i) {
        // Need to set the end time of the previous time
        // interval before creating a new one.
        ti->SetEndTime(Strings::cdble((*i)->Data()));

        // Check for steps attribute.
        attr = (*i)->GetAttribute("steps");
        if (attr != NULL) {
            ti->SetStepCount((unsigned int)Strings::cdble(attr->GetValue()));
        } else {
            if (ti != NULL) delete ti;
            throw std::runtime_error("Time interval must define number of steps "
                                "(Mops, Settings_IO::readTimeIntervals)");
        }

        // Check for splits attribute.
        attr = (*i)->GetAttribute("splits");
        if (attr != NULL) {
            ti->SetSplittingStepCount((unsigned int)Strings::cdble(attr->GetValue()));
        } else {
            ti->SetSplittingStepCount(splits);
        }

        // Check for splits attribute.
        attr = (*i)->GetAttribute("subsplits");
        if (attr != NULL) {
            ti->SetSubSplittingStepCount((unsigned int)Strings::cdble(attr->GetValue()));
        }

        // Add the completed time interval to the vector.
        times.push_back(*ti);

        // Create a new time interval and set the start time.
        delete ti; ti = new TimeInterval();
        ti->SetStartTime(Strings::cdble((*i)->Data()));
    }

    // Clear memory used to load time intervals.
    if (ti != NULL) delete ti; ti = NULL;
}


/**
 * Read a list of initial particles specified in an XML file
 *
 *@param[in]    population_xml      xml node containing initial particles as children
 *@param[in]    particle_mech       Mechanism defining the meaning of the particles
 *
 *@return       Container of the particles that form the initial population
 *
 *@exception    std::runtime_error  If particles are specified they must have strictly positive counts
 */
Sweep::PartPtrList Settings_IO::ReadInitialParticles(const CamXML::Element& population_xml,
                                                     const Sweep::Mechanism & particle_mech)
{
    // Accumulate in this container a collection of particles to be inserted into the ensemble
    Sweep::PartPtrList particleList;

    // Now get the XML defining the particles
    std::vector<CamXML::Element*> particleXML;
    population_xml.GetChildren("particle", particleXML);

    // An empty population means nothing to do
    if(particleXML.empty())
    {
        return particleList;
    }

    // See if there is a spatial position for the particles
    const CamXML::Attribute *attr = population_xml.GetAttribute("x");
    double position = 0.0;
    if(attr) {
        position = atof(attr->GetValue().c_str());
    }

    // Now loop through the particles one by one
    std::vector<CamXML::Element*>::const_iterator it(particleXML.begin());
    const std::vector<CamXML::Element*>::const_iterator itEnd(particleXML.end());
    while(it != itEnd)
    {
        // Build a particle from the xml and use an auto_ptr so the particle
        // instance will be deleted if an exception is thrown while processing
        // the rest of the input.
        const std::auto_ptr<Sweep::Particle>
            pParticle(Sweep::Particle::createFromXMLNode(**it, particle_mech));

        // Assume particle population applies at time 0.
        pParticle->setPositionAndTime(position, 0.0);

        // See how many times this particle appears (default is 1)
        int repeatCount = 1;

        // Will initialise to NULL if no count attribute present and if block
        // will be skipped.
        const CamXML::Attribute* const  countAttrib = (*it)->GetAttribute("count");
        if(countAttrib)
        {
            // Read the count specified by the user
            const std::string countString = countAttrib->GetValue();
            repeatCount = atoi(countString.c_str());

            // The count must be a strictly positive integer
            if(repeatCount <= 0)
            {
                // Delete any particles already read into the local list
                // The list itself will be destroyed when the throw std::takes
                // place so there is no need to nullify the pointers
                Sweep::PartPtrList::iterator it = particleList.begin();
                const Sweep::PartPtrList::iterator itEnd = particleList.end();
                while(it != itEnd)
                {
                    delete *it++;
                }
                throw std::runtime_error("Particle counts must be > 0. Mops Reactor::ReadInitialParticles");
            }
        }

        // Add the specified number of copies of the particle to the initial
        // particle container
        while(repeatCount--)
            particleList.push_back(new Sweep::Particle(*pParticle));

        // Move on to the next particle in the xml
        ++it;
    }
    return particleList;
}


// Version of the above reader that will construct (bintree) particles given 
// (mean,std) number of primaries per particle and (mean,std) number of component
// units per particle (assuming lognormal distributions).
// This has not been tested much. 
/**
* Read a list of initial particles specified in an XML file
*
*@param[in]    population_xml      xml node containing initial particles as children
*@param[in]    particle_mech       Mechanism defining the meaning of the particles
*
*@return       Container of the particles that form the initial population
*
*@exception    std::runtime_error  If particles are specified they must have strictly positive counts
*/
Sweep::PartPtrList Settings_IO::ReadInitialParticlesDetailed(const CamXML::Element& population_xml,
	const Sweep::Mechanism & particle_mech)
{
	// Generate a random number generator for use inside this once-off function
	boost::mt19937 rng_temp(0);

	// Accumulate in this container a collection of particles to be inserted into the ensemble
	Sweep::PartPtrList particleList;

	// Now get the XML defining the particles
	std::vector<CamXML::Element*> particleXML;
	population_xml.GetChildren("particle", particleXML);

	// An empty population means nothing to do
	if (particleXML.empty())
	{
		return particleList;
	}

	// See if there is a spatial position for the particles
	const CamXML::Attribute *attr = population_xml.GetAttribute("x");
	double position = 0.0;
	if (attr) {
		position = atof(attr->GetValue().c_str());
	}

	// Now loop through the particles one by one
	std::vector<CamXML::Element*>::const_iterator it(particleXML.begin());
	const std::vector<CamXML::Element*>::const_iterator itEnd(particleXML.end());
	while (it != itEnd)
	{
		// Build a particle from the xml and use an auto_ptr so the particle
		// instance will be deleted if an exception is thrown while processing
		// the rest of the input.
		const std::auto_ptr<Sweep::Particle>
			pParticle(Sweep::Particle::createFromXMLNodeDetailed(**it, particle_mech, rng_temp));

		// Assume particle population applies at time 0.
		pParticle->setPositionAndTime(position, 0.0);

		// See how many times this particle appears (default is 1)
		int repeatCount = 1;

		// Will initialise to NULL if no count attribute present and if block
		// will be skipped.
		const CamXML::Attribute* const  countAttrib = (*it)->GetAttribute("count");
		if (countAttrib)
		{
			// Read the count specified by the user
			const std::string countString = countAttrib->GetValue();
			repeatCount = atoi(countString.c_str());

			// The count must be a strictly positive integer
			if (repeatCount <= 0)
			{
				// Delete any particles already read into the local list
				// The list itself will be destroyed when the throw std::takes
				// place so there is no need to nullify the pointers
				Sweep::PartPtrList::iterator it = particleList.begin();
				const Sweep::PartPtrList::iterator itEnd = particleList.end();
				while (it != itEnd)
				{
					delete *it++;
				}
				throw std::runtime_error("Particle counts must be > 0. Mops Reactor::ReadInitialParticles");
			}
		}

		// Add the specified number of copies of the particle to the initial
		// particle container
		Sweep::Particle *pParticle2;
		while (repeatCount--){
			pParticle2 = Sweep::Particle::createFromXMLNodeDetailed(**it, particle_mech, rng_temp);
			// Assume particle population applies at time 0.
			pParticle2->setPositionAndTime(position, 0.0);
			particleList.push_back(new Sweep::Particle(*pParticle2));
		}

		// Move on to the next particle in the xml
		++it;
	}
	return particleList;
}


/*!
 * Parse XML and set values specifying the boundaries that defined extreme particles to be excluded
 * from particle population statistics.
 *
 * The node must include a "property" attribute and child nodes <lower> and <upper> to contain the
 * numerical values of the bounds.
 *
 * Code extracted from the ReadOutput function also in this file.
 *
 *@param[in]        node            XML element of type <statsbound>
 *@param[in]        property_id     Index of particle property to which limit values apply
 *@param[in]        lower_bound     Lowest  value that a particle may have and still contribute to population statistics
 *@param[in]        upper_bound     Highest value that a particle may have and still contribute to population statistics
 *
 *@exception    std::runtime_error      No lower bound in XML
 *@exception    std::runtime_error      No upper bound in XML
 *@exception    std::runtime_error      No particle property specified in XML
 *@exception    std::runtime_error      Upper bound less than or equal to lower bound
 */
void Settings_IO::ReadStatsBound(const CamXML::Element &node, Sweep::PropID &property_id,
                                 double &lower_bound, double &upper_bound) {

    const CamXML::Element *lobonode = node.GetFirstChild("lower");
    const CamXML::Element *upbonode = node.GetFirstChild("upper");
    // Set the upper and lower bounds if found the entry in MOPS input file
    if (lobonode != NULL) {
        lower_bound = Strings::cdble(lobonode->Data());
    }
    else
        throw std::runtime_error("No lower bound in XML passed to (MOPS, Settings_IO::ReadStatsBound).");

    if (upbonode != NULL) {
        upper_bound = Strings::cdble(upbonode->Data());
        // Check if upper bound is greater than lower bound
        if (upper_bound <= lower_bound) {
            throw std::runtime_error("STATSBOUND Error: Upper bound is less than or equal to lower bound. "
                                     "(MOPS, Settings_IO::ReadStatsBound).");
        }
    }
    else
        throw std::runtime_error("No upper bound in XML passed to (MOPS, Settings_IO::ReadStatsBound).");

    // Set the property id of the particle to enforce the statistical bounds on
    // Currently, only 2 properties are implemented
    std::string prop_str = node.GetAttributeValue("property");
    if ((prop_str.compare("dcol") == 0) ||
        (prop_str.compare("Dcol") == 0)) {
        property_id = Sweep::iDcol;
    } else if ((prop_str.compare("dmob") == 0) ||
               (prop_str.compare("Dmob") == 0)) {
        property_id = Sweep::iDmob;
    } else if ((prop_str.compare("m") == 0) ||
               (prop_str.compare("M") == 0)) {
        property_id = Sweep::iM;
    } else {
        throw std::runtime_error("No particle property specified in XML passed to (MOPS, Settings_IO::ReadStatsBound).");
    }
}



Mops::ReactorNetwork* Settings_IO::LoadNetwork(
        const std::string &filename,
        std::vector<TimeInterval> &times,
        Simulator &sim,
        Solver &solver,
        Mechanism &mech) {

    // Set up some storage
    CamXML::Document doc;
    const CamXML::Element *root, *node;
    Mops::ReactorNetwork* net = new Mops::ReactorNetwork(mech);

    // Attempt to load XML document.
    if (doc.Load(filename) == 0) {
        // Get the root element.
        root = doc.Root();

        // Check the mops file version.
        const CamXML::Attribute *attr = root->GetAttribute("version");

        if ((attr!=NULL) && (attr->GetValue() != "3")) {
            throw std::runtime_error("Settings file has wrong version number (expecting 3)"
                                " (Mops::Settings_IO::LoadNetwork).");
        }

        // GLOBAL SETTINGS.

        readGlobalSettings(*root, sim, solver);

        // OUTPUT SETTINGS.
        // wjm34: read output settings before reactor, so we can check if ensemble/g.p.
        // files are consistent with some more simulation settings.

        node = root->GetFirstChild("output");
        if (node != NULL) {
            readOutput(*node, sim, mech);
        } else {
            throw std::runtime_error("Settings file does not contain output"
                                " information (Mops::Settings_IO::LoadFromXML).");
        }

        // NETWORK.
        readNetwork(*root, sim, mech, *net);
        std::cout << *net;

        // MECHANISM REDUCTION.

        node = root->GetFirstChild("LOI");
        if (node!= NULL) {
            ReadLOIStatus(*node, solver);
        }

        // TIME INTERVALS.

        node = root->GetFirstChild("timeintervals");
        if (node != NULL) {
            readTimeIntervals(*node, times);
        } else {
            throw std::runtime_error("Settings file does not contain a time intervals"
                                " (Mops::Settings_IO::LoadFromXML).");
        }
    }

    return net;
}

