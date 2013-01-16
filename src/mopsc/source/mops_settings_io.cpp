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
        }
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
        }
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

// Reads the reactor initial settings from the given XML node.
Reactor *const readReactor(const CamXML::Element &node,
                                        const Mechanism &mech,
                                        const unsigned int max_particle_count,
										const double maxM0,
										Mops::Simulator &sim)
{
    Reactor *reac = NULL;
    const CamXML::Element *subnode, *subsubnode;
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
    Mixture *mix = new Mixture(mech.ParticleMech());

    // VISCOSITY MODEL

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
            throw std::runtime_error("Unrecognised viscosity model. (Mops, Settings_IO::readReactor).");
        }
    }
    else {
        // Use air by default
        mix->GasPhase().SetViscosityModel(Sprog::iAir);
    }


    fvector molefracs(mech.GasMech().SpeciesCount(), 0.0);

    // Get the temperature.
    subnode = node.GetFirstChild("temperature");
    if(subnode != NULL) {
        mix->GasPhase().SetTemperature(readTemperature(*subnode));
    }
    else {
        throw std::runtime_error("No initial condition for temperature (Mops, Settings_IO::readReactor).");
    }

    // Get the pressure.
    subnode = node.GetFirstChild("pressure");
    if(subnode != NULL) {
        mix->GasPhase().SetPressure(readPressure(*subnode));
    }
    else {
        throw std::runtime_error("No initial condition for pressure (Mops, Settings_IO::readReactor).");
    }

    // Fill the mixture object.
    node.GetChildren("component", nodes);
    for (i=nodes.begin(); i!=nodes.end(); ++i) {
        // Check the ID attribute.
        attr = (*i)->GetAttribute("id");
        if (attr != NULL) {
            str = attr->GetValue();
            // This should be a species initial concentration.

            // Get the species index.
            int j = mech.GasMech().FindSpecies(str);

            // Set the species mole fraction.
            if (j >= 0) {
                molefracs[j] = Strings::cdble((*i)->Data()); // from inx input
            } else {
                throw std::runtime_error("No initial condition for species " + str +
                                    " (Mops, Settings_IO::readReactor).");
            }
        } else {
            throw std::runtime_error("Initial condition must have ID "
                                "(Mops, Settings_IO::readReactor)");
        }

    }

    bool doubling_activated;
    subnode = node.GetFirstChild("DoublingAlgorithm");
    if (subnode != NULL) {
        std::string str_enable = subnode->GetAttributeValue("enable");
        if (str_enable.compare("false") == 0) {
            std::cout << "sweep: Warning! doubling is turned off by user.\n";
            doubling_activated = false;
        } 
        else doubling_activated = true;

    } else doubling_activated = true;

    // Assign the species mole fraction vector to the reactor mixture.
    mix->GasPhase().SetFracs(molefracs);
    std::cout << "ReadReactor: MoleFrac assigned"<< endl;
    mix->Particles().Initialise(max_particle_count, doubling_activated);
	mix->Reset(maxM0);
    reac->Fill(*mix);

    

    // Particles
    subnode = node.GetFirstChild("population");

    // List will be empty unless a population node is present
    if(subnode) {
        Sweep::PartPtrList fileParticleList;
        Sweep::PartPtrList inxParticleList;
        Sweep::PartPtrList allParticleList;
        double initialM0 = 0;

        // Find the overall number density represented by the population
        CamXML::Element* m0Node = subnode->GetFirstChild("m0");
        if(!m0Node) {
            throw std::runtime_error("m0 (number density) must be specified for initial particle population \
                                     (Mops, Settings_IO::readReactor)");
        }
        initialM0 = std::atof((m0Node->Data()).c_str());
        if(initialM0 < 0) {
            throw std::runtime_error("m0 for initial particle population may not be negative \
                                     (Mops, Settings_IO::readReactor)");
        }

        // Check if a binary file has been specified for particles..
        CamXML::Element* fileNode = subnode->GetFirstChild("file");
        if (fileNode) {
            std::string filename;
            filename = fileNode->Data();
            std::cout << "parser: binary file " << filename << " specified for input.\n";
            fileParticleList = sim.ReadEnsembleFile(*reac, filename);
        }

        // Now read in the list of particles and sum up their statistical weights
        inxParticleList = Settings_IO::ReadInitialParticles(*subnode, mech.ParticleMech());

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

        mix->SetParticles(allParticleList.begin(), allParticleList.end(), initialM0 / weightSum);
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

    if (reac->SerialType() == Serial_PSR) {
        subnode = node.GetFirstChild("inflow");

        if (subnode != NULL) {
            // Create a new Mixture object for inflow.
            FlowStream *inf = new FlowStream(mech);
            molefracs.assign(mech.GasMech().SpeciesCount(), 0.0);

            // Get the temperature.
            subsubnode = subnode->GetFirstChild("temperature");
            inf->Mixture()->GasPhase().SetTemperature(readTemperature(*subsubnode));

            // Get the pressure.
            subsubnode = subnode->GetFirstChild("pressure");
            inf->Mixture()->GasPhase().SetPressure(readPressure(*subsubnode));

            // Fill the inflow Mixture object.
            subnode->GetChildren("component", nodes);
            for (i=nodes.begin(); i!=nodes.end(); ++i) {
                // Check the ID attribute.
                attr = (*i)->GetAttribute("id");
                if (attr != NULL) {
                    str = attr->GetValue();

                    // This should be a species initial concentration.

                    // Get the species index.
                    int j = mech.GasMech().FindSpecies(str);

                    // Set the species mole fraction.
                    if (j >= 0) {
                        molefracs[j] = Strings::cdble((*i)->Data());
                    } else {
                        throw std::runtime_error("Unknown species inflow condition "
                                            "(Mops, Settings_IO::readReactor).");
                    }
                } else {
                    throw std::runtime_error("Inflow condition must have ID! "
                                        "(Mops, Settings_IO::readReactor)");
                }
            }

            // Assign the species mole fraction vector
            // to the reactor inflow mixture.
            inf->Mixture()->GasPhase().SetFracs(molefracs);
            dynamic_cast<PSR*>(reac)->SetInflow(*inf);
        } else {
            throw std::runtime_error("Inflow conditions must be defined for a PSR "
                                "(Mops, Settings_IO::readReactor)");
        }

        // Read the residence time.
        subnode = node.GetFirstChild("residencetime");
        if (subnode != NULL) {
        	double tau = Strings::cdble(subnode->Data());
            dynamic_cast<PSR*>(reac)->SetResidenceTime(tau);
        }
    }

    return reac;
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
}

} // anonymous namespace


Reactor *const Settings_IO::LoadFromXML_V1(const std::string &filename,
                                           Mops::Reactor *reac,
                                           std::vector<TimeInterval> &times,
                                           Simulator &sim, Solver &solver,
                                           const Mechanism &mech)
{
    CamXML::Document doc;
    const CamXML::Element *root, *node;
    std::vector<CamXML::Element*> nodes;
    std::vector<CamXML::Element*>::const_iterator i;
    const CamXML::Attribute *attr;
    std::string str;

    // We need a blank reactor object to create.  If the pointer
    // passed to this function is valid then it needs to be
    // deleted.
    if (reac != NULL) delete reac;

    // Attempt to load XML document.
    if (doc.Load(filename) == 0) {
        // Get the root element.
        root = doc.Root();

        // SYSTEM TYPE AND ATTRIBUTES.

        // First read the system type.  This is where we initialise
        // the reactor type.
        attr = root->GetAttribute("system");
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
                throw std::invalid_argument("Shocktube not currently implemented");
            } else {
                // Default reactor is a batch reactor.
                reac = ReactorFactory::Create(Serial_Batch, mech);
            }
        }
	
	std::cout << "Load from old xml Reactor Type Passed"<< endl;

        // Now check for constant temperature.  If not constant temperature
        // then the reactor is solved using the adiabatic energy equation.
        attr = root->GetAttribute("constt");
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


        // Now check for constant volume.
        attr = root->GetAttribute("constv");
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


        // REACTOR INITIAL CONDITIONS.

        // Create a new Mixture object.
        Mixture *mix = new Mixture(mech.ParticleMech());

        fvector molefracs(mech.GasMech().SpeciesCount(), 0.0);

        // Fill the mixture object.
        root->GetChildren("init", nodes);
        for (i=nodes.begin(); i!=nodes.end(); i++) {
            // Check the ID attribute.
            attr = (*i)->GetAttribute("id");
            if (attr != NULL) {
                str = attr->GetValue();

                if (str.compare("T")==0) {
                    // This is the initial temperature.
                    mix->GasPhase().SetTemperature(Strings::cdble((*i)->Data()));
                } else if (str.compare("P")==0) {
                    // This is the initial pressure (Remember to convert to Pa).
                    mix->GasPhase().SetPressure(1.0e5 * Strings::cdble((*i)->Data()));
                } else {
                    // This should be a species initial concentration.

                    // Get the species index.
                    int j = mech.GasMech().FindSpecies(str);

                    // Set the species mole fraction.
                    if (j >= 0) {
                        molefracs[j] = Strings::cdble((*i)->Data());
                    } else {
                        throw std::runtime_error("Unknown species initial condition "
                                            "(Mops, Settings_IO::LoadFromXML_V1).");
                    }
                }
            } else {
                throw std::invalid_argument("Initial condition must have ID!");
            }
        }

        // Assign the species mole fraction vector to the reactor mixture.
        mix->GasPhase().SetFracs(molefracs);
        reac->Fill(*mix);


        // PSR SPECIFIC SETTINGS.

        if (reac->SerialType() == Serial_PSR) {
            // Create a new Mixture object for inflow.
            FlowStream *inf = new FlowStream(mech);
            molefracs.assign(mech.GasMech().SpeciesCount(), 0.0);

            // Fill the inflow Mixture object.
            root->GetChildren("inflow", nodes);
            for (i=nodes.begin(); i!=nodes.end(); i++) {
                // Check the ID attribute.
                attr = (*i)->GetAttribute("id");
                if (attr != NULL) {
                    str = attr->GetValue();

                    if (str.compare("T")==0) {
                        // This is the initial temperature.
                        inf->Mixture()->GasPhase().SetTemperature(Strings::cdble((*i)->Data()));
                    } else if (str.compare("P")==0) {
                        // This is the initial pressure (Remember to convert to Pa).
                        inf->Mixture()->GasPhase().SetPressure(1.0e5 * Strings::cdble((*i)->Data()));
                    } else {
                        // This should be a species initial concentration.

                        // Get the species index.
                        int j = mech.GasMech().FindSpecies(str);

                        // Set the species mole fraction.
                        molefracs[j] = Strings::cdble((*i)->Data());
                    }
                } else {
                    throw std::invalid_argument("Inflow condition must have ID!");
                }
            }

            // Assign the species mole fraction vector
            // to the reactor inflow mixture.
            inf->Mixture()->GasPhase().SetFracs(molefracs);
            dynamic_cast<PSR*>(reac)->SetInflow(*inf);

            // Read the residence time.
            node = root->GetFirstChild("residencetime");
            if (node != NULL) {
                dynamic_cast<PSR*>(reac)->SetResidenceTime(Strings::cdble(node->Data()));
            }
        }


        // GLOBAL SETTINGS.

        // Read the relative error tolerance.
        node = root->GetFirstChild("rtol");
        if (node != NULL) {
            solver.SetRTOL(Strings::cdble(node->Data()));
        }

        // Read the absolute error tolerance.
        node = root->GetFirstChild("atol");
        if (node != NULL) {
            solver.SetATOL(Strings::cdble(node->Data()));
        }

        // Read the number of runs.
        node = root->GetFirstChild("runs");
        if (node != NULL) {
            sim.SetRunCount((int)Strings::cdble(node->Data()));
        }

        // Read stochastic particle count.
        node = root->GetFirstChild("pcount");
        if (node != NULL) {
            sim.SetMaxPartCount((int)Strings::cdble(node->Data()));
        }

        // Read max M0, for scaling.
        node = root->GetFirstChild("maxm0");
        if (node != NULL) {
            sim.SetMaxM0(Strings::cdble(node->Data())*1.0e6); // Convert from #/cm3 to #/m3.
        }

        // TIME INTERVALS.

        // Read splits per step.
        unsigned int splits = 1;
        node = root->GetFirstChild("splitsperstep");
        if (node != NULL) {
            splits = (unsigned int)Strings::cdble(node->Data());
        }

        TimeInterval *ti = NULL;

        root->GetChildren("time", nodes);
        for (i=nodes.begin(); i!=nodes.end(); i++) {
            attr = (*i)->GetAttribute("start");

            if (attr != NULL) {
                // This is the start time.  Create a new time inteval
                // and set the start time.
                ti = new TimeInterval();
                ti->SetStartTime(Strings::cdble((*i)->Data()));
                ti->SetSplittingStepCount(splits);
            } else {
                // This is not the start time.  Need to set the end time of
                // the previous time interval before creating a new one.
                ti->SetEndTime(Strings::cdble((*i)->Data()));

                // Check for steps attribute.
                attr = (*i)->GetAttribute("steps");
                if (attr != NULL) {
                    ti->SetStepCount((unsigned int)Strings::cdble(attr->GetValue()));
                } else {
                    if (ti != NULL) delete ti;
                    throw std::invalid_argument("Time interval must define number of steps.");
                }

                // Add the completed time interval to the vector.
                times.push_back(*ti);

                // Create a new time interval and set the start time.
                delete ti; ti = new TimeInterval();
                ti->SetStartTime(Strings::cdble((*i)->Data()));
            }
        }

        // Clear memory used to load time intervals.
        if (ti != NULL) delete ti; ti = NULL;


        // CONSOLE OUTPUT.

        node = root->GetFirstChild("console");

        if (node != NULL) {
            // Check the console interval.
            attr = node->GetAttribute("interval");
            if (attr != NULL) {
                sim.SetConsoleInterval((unsigned int)Strings::cdble(attr->GetValue()));
            }

            // Read the column variables.
            str = node->Data();
            std::vector<std::string> cvars;
            Strings::split(str, cvars, " ");
            for (unsigned int j=0; j!=cvars.size(); ++j) {
                sim.AddConsoleVariable(cvars[j]);
            }
        }


        // FILE OUTPUT.

        // Read the output file name.
        node = root->GetFirstChild("output");
        if (node != NULL) {
            sim.SetOutputFile(node->Data());
        }
    }

    return reac;
}




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

	std::cout << "Load from NEW xml Global Setting Read"<< endl;

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

        if (mech.ParticleMech().AggModel() == Sweep::AggModels::PAH_KMC_ID){
            //currently only the data for mass spectra is provided, but the capacity will be extrend according to the purpose of the project
            node = root->GetFirstChild("postprocess");
            if (node != NULL) {
                readPostprocess(*node, sim, mech);
            } else {
                // mass spectra data will not be generated, this should be only for the old test files
                sim.SetMassSpectra(false);
            }
        }
        // REACTOR.

        node = root->GetFirstChild("reactor");
        if (node != NULL) {
			reac = readReactor(*node, mech, sim.MaxPartCount(), sim.MaxM0(), sim);
			std::cout << "Load from NEW xml Reactor Read"<< endl;
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

