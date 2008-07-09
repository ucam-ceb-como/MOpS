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

#include "mops_params.h"
#include "mops_settings_io.h"
#include "mops_timeinterval.h"
#include "mops_reactor.h"
#include "mops_psr.h"
#include "mops_mixture.h"
#include "mops_reactor_factory.h"
#include "mops_predcor_solver.h"

#include "camxml.h"
#include "string_functions.h"

#include <vector>
#include <string>
#include <stdexcept>

using namespace Mops;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Settings_IO::Settings_IO(void)
{
}

// Default destructor.
Settings_IO::~Settings_IO(void)
{
}


// VERSION 1 XML SETTINGS FILES.

Reactor *const Settings_IO::LoadFromXML_V1(const std::string &filename, 
                                           Mops::Reactor *reac, 
                                           std::vector<TimeInterval> &times,
                                           Simulator &sim, Solver &solver,
                                           const Mechanism &mech)
{
    CamXML::Document doc;
    const CamXML::Element *root, *node;
    vector<CamXML::Element*> nodes;
    vector<CamXML::Element*>::const_iterator i;
    const CamXML::Attribute *attr;
    string str;

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
                throw invalid_argument("Shocktube not currently implemented");
            } else {
                // Default reactor is a batch reactor.
                reac = ReactorFactory::Create(Serial_Batch, mech);
            }
        }

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

        fvector molefracs(mech.SpeciesCount(), 0.0);

        // Fill the mixture object.
        root->GetChildren("init", nodes);
        for (i=nodes.begin(); i!=nodes.end(); i++) {
            // Check the ID attribute.
            attr = (*i)->GetAttribute("id");
            if (attr != NULL) {
                str = attr->GetValue();

                if (str.compare("T")==0) {
                    // This is the initial temperature.
                    mix->SetTemperature(cdble((*i)->Data()));
                } else if (str.compare("P")==0) {
                    // This is the initial pressure (Remember to convert to Pa).
                    mix->SetPressure(1.0e5 * cdble((*i)->Data()));
                } else {
                    // This should be a species initial concentration.

                    // Get the species index.
                    int j = mech.FindSpecies(str);
                    
                    // Set the species mole fraction.
                    if (j >= 0) {
                        molefracs[j] = cdble((*i)->Data());
                    } else {
                        throw runtime_error("Unknown species initial condition "
                                            "(Mops, Settings_IO::LoadFromXML_V1).");
                    }
                }
            } else {
                throw invalid_argument("Initial condition must have ID!");
            }
        }

        // Assign the species mole fraction vector to the reactor mixture.
        mix->SetFracs(molefracs);
        reac->Fill(*mix);


        // PSR SPECIFIC SETTINGS.

        if (reac->SerialType() == Serial_PSR) {
            // Create a new Mixture object for inflow.
            Mixture *inf = new Mixture(mech.ParticleMech());
            molefracs.assign(mech.SpeciesCount(), 0.0);

            // Fill the inflow Mixture object.
            root->GetChildren("inflow", nodes);
            for (i=nodes.begin(); i!=nodes.end(); i++) {
                // Check the ID attribute.
                attr = (*i)->GetAttribute("id");
                if (attr != NULL) {
                    str = attr->GetValue();

                    if (str.compare("T")==0) {
                        // This is the initial temperature.
                        inf->SetTemperature(cdble((*i)->Data()));
                    } else if (str.compare("P")==0) {
                        // This is the initial pressure (Remember to convert to Pa).
                        inf->SetPressure(1.0e5 * cdble((*i)->Data()));
                    } else {
                        // This should be a species initial concentration.

                        // Get the species index.
                        int j = mech.FindSpecies(str);
                        
                        // Set the species mole fraction.
                        molefracs[j] = cdble((*i)->Data());
                    }
                } else {
                    throw invalid_argument("Inflow condition must have ID!");
                }
            }

            // Assign the species mole fraction vector 
            // to the reactor inflow mixture.
            inf->SetFracs(molefracs);
            dynamic_cast<PSR*>(reac)->SetInflow(*inf);

            // Read the residence time.
            node = root->GetFirstChild("residencetime");
            if (node != NULL) {
                dynamic_cast<PSR*>(reac)->SetResidenceTime(cdble(node->Data()));
            }
        }


        // GLOBAL SETTINGS.

        // Read the relative error tolerance.
        node = root->GetFirstChild("rtol");
        if (node != NULL) {
            solver.SetRTOL(cdble(node->Data()));
        }

        // Read the absolute error tolerance.
        node = root->GetFirstChild("atol");
        if (node != NULL) {
            solver.SetATOL(cdble(node->Data()));
        }

        // Read the number of runs.
        node = root->GetFirstChild("runs");
        if (node != NULL) {
            sim.SetRunCount((int)cdble(node->Data()));
        }

        // Read stochastic particle count.
        node = root->GetFirstChild("pcount");
        if (node != NULL) {
            sim.SetMaxPartCount((int)cdble(node->Data()));
        }

        // Read max M0, for scaling.
        node = root->GetFirstChild("maxm0");
        if (node != NULL) {
            sim.SetMaxM0(cdble(node->Data())*1.0e6); // Convert from #/cm3 to #/m3.
        }

        // TIME INTERVALS.

        // Read splits per step.
        unsigned int splits = 1;
        node = root->GetFirstChild("splitsperstep");
        if (node != NULL) {
            splits = (unsigned int)cdble(node->Data());
        }

        TimeInterval *ti = NULL;

        root->GetChildren("time", nodes);
        for (i=nodes.begin(); i!=nodes.end(); i++) {
            attr = (*i)->GetAttribute("start");

            if (attr != NULL) {
                // This is the start time.  Create a new time inteval
                // and set the start time.
                ti = new TimeInterval();
                ti->SetStartTime(cdble((*i)->Data()));
                ti->SetSplittingStepCount(splits);
            } else {
                // This is not the start time.  Need to set the end time of
                // the previous time interval before creating a new one.
                ti->SetEndTime(cdble((*i)->Data()));

                // Check for steps attribute.
                attr = (*i)->GetAttribute("steps");
                if (attr != NULL) {
                    ti->SetStepCount((unsigned int)cdble(attr->GetValue()));
                } else {
                    if (ti != NULL) delete ti;
                    throw invalid_argument("Time interval must define number of steps.");
                }
                
                // Add the completed time interval to the vector.
                times.push_back(*ti);

                // Create a new time interval and set the start time.
                delete ti; ti = new TimeInterval();
                ti->SetStartTime(cdble((*i)->Data()));
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
                sim.SetConsoleInterval((unsigned int)cdble(attr->GetValue()));
            }

            // Read the column variables.
            str = node->Data();
            vector<string> cvars;
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
                                        const Mechanism &mech)
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
            throw runtime_error("Settings file has wrong version number (expecting 2)"
                                " (Mops::Settings_IO::LoadFromXML).");
        }

        // GLOBAL SETTINGS.

        readGlobalSettings(*root, sim, solver);

        // REACTOR.

        node = root->GetFirstChild("reactor");
        if (node != NULL) {
            reac = readReactor(*node, mech);
        } else {
            throw runtime_error("Settings file does not contain a reactor definition"
                                " (Mops::Settings_IO::LoadFromXML).");
        }

        // TIME INTERVALS.

        node = root->GetFirstChild("timeintervals");
        if (node != NULL) {
            readTimeIntervals(*node, times);
        } else {
            throw runtime_error("Settings file does not contain a time intervals"
                                " (Mops::Settings_IO::LoadFromXML).");
        }

        // OUTPUT SETTINGS.

        node = root->GetFirstChild("output");
        if (node != NULL) {
            readOutput(*node, sim);
        } else {
            throw runtime_error("Settings file does not contain output"
                                " information (Mops::Settings_IO::LoadFromXML).");
        }
    }

    return reac;
}


// V2 SETTINGS FILE SECTIONS.

// Reads global simulation settings from the given XML node.
void Settings_IO::readGlobalSettings(const CamXML::Element &node, 
                                     Simulator &sim, Solver &solver)
{
    const CamXML::Element *subnode;

    // Read the relative error tolerance.
    subnode = node.GetFirstChild("rtol");
    if (subnode != NULL) {
        solver.SetRTOL(cdble(subnode->Data()));
    }

    // Read the absolute error tolerance.
    subnode = node.GetFirstChild("atol");
    if (subnode != NULL) {
        solver.SetATOL(cdble(subnode->Data()));
    }

    // Read the number of runs.
    subnode = node.GetFirstChild("runs");
    if (subnode != NULL) {
        sim.SetRunCount((int)cdble(subnode->Data()));
    }

    // Read the number of iterations.
    subnode = node.GetFirstChild("iter");
    if (subnode != NULL) {
        sim.SetIterCount((int)cdble(subnode->Data()));
    }

    // Read stochastic particle count.
    subnode = node.GetFirstChild("pcount");
    if (subnode != NULL) {
        sim.SetMaxPartCount((int)cdble(subnode->Data()));
    }

    // Read predicted maximum M0 value.
    subnode = node.GetFirstChild("maxm0");
    if (subnode != NULL) {
        sim.SetMaxM0(cdble(subnode->Data())*1.0e6); // Convert from #/cm3 to #/m3.
    }

    // Read predictor-corrector relaxation parameter.
    subnode = node.GetFirstChild("relax");
    if (subnode != NULL) {
        solver.SetUnderRelaxCoeff(cdble(subnode->Data()));
    }
}

// Reads the reactor initial settings from the given XML node.
Reactor *const Settings_IO::readReactor(const CamXML::Element &node,
                                               const Mechanism &mech)
{
    Reactor *reac = NULL;
    const CamXML::Element *subnode, *subsubnode;
    vector<CamXML::Element*> nodes;
    vector<CamXML::Element*>::const_iterator i;
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
            throw runtime_error("Shocktube not currently implemented"
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


    // REACTOR INITIAL CONDITIONS.

    // Create a new Mixture object.
    Mixture *mix = new Mixture(mech.ParticleMech());

    fvector molefracs(mech.SpeciesCount(), 0.0);

    // Get the temperature.
    subnode = node.GetFirstChild("temperature");
    mix->SetTemperature(readTemperature(*subnode));

    // Get the pressure.
    subnode = node.GetFirstChild("pressure");
    mix->SetPressure(readPressure(*subnode));

    // Fill the mixture object.
    node.GetChildren("component", nodes);
    for (i=nodes.begin(); i!=nodes.end(); i++) {
        // Check the ID attribute.
        attr = (*i)->GetAttribute("id");
        if (attr != NULL) {
            str = attr->GetValue();
            // This should be a species initial concentration.

            // Get the species index.
            int j = mech.FindSpecies(str);
            
            // Set the species mole fraction.
            if (j >= 0) {
                molefracs[j] = cdble((*i)->Data());
            } else {
                throw runtime_error("Unknown species initial condition "
                                    "(Mops, Settings_IO::readReactor).");
            }
        } else {
            throw runtime_error("Initial condition must have ID "
                                "(Mops, Settings_IO::readReactor)");
        }
    }

    // Assign the species mole fraction vector to the reactor mixture.
    mix->SetFracs(molefracs);
    reac->Fill(*mix);

    // TEMPERATURE GRADIENT PROFILE.

    node.GetChildren("dTdt", nodes);
    for (i=nodes.begin(); i!=nodes.end(); ++i) {
        // Get the start time attribute.
        attr = (*i)->GetAttribute("startt");
        if (attr != NULL) {
            Sweep::Maths::Linear fun;
            fun.SetParam(1, cdble((*i)->Data()));
            reac->Add_dTdt(cdble(attr->GetValue()), fun);
        }
    }

    // PSR SPECIFIC SETTINGS.

    if (reac->SerialType() == Serial_PSR) {
        subnode = node.GetFirstChild("inflow");

        if (subnode != NULL) {
            // Create a new Mixture object for inflow.
            Mixture *inf = new Mixture(mech.ParticleMech());
            molefracs.assign(mech.SpeciesCount(), 0.0);

            // Get the temperature.
            subsubnode = subnode->GetFirstChild("temperature");
            inf->SetTemperature(readTemperature(*subsubnode));

            // Get the pressure.
            subsubnode = subnode->GetFirstChild("pressure");
            inf->SetPressure(readPressure(*subsubnode));

            // Fill the inflow Mixture object.
            subnode->GetChildren("component", nodes);
            for (i=nodes.begin(); i!=nodes.end(); i++) {
                // Check the ID attribute.
                attr = (*i)->GetAttribute("id");
                if (attr != NULL) {
                    str = attr->GetValue();

                    // This should be a species initial concentration.

                    // Get the species index.
                    int j = mech.FindSpecies(str);
                    
                    // Set the species mole fraction.
                    if (j >= 0) {
                        molefracs[j] = cdble((*i)->Data());
                    } else {
                        throw runtime_error("Unknown species inflow condition "
                                            "(Mops, Settings_IO::readReactor).");
                    }
                } else {
                    throw runtime_error("Inflow condition must have ID! "
                                        "(Mops, Settings_IO::readReactor)");
                }
            }

            // Assign the species mole fraction vector 
            // to the reactor inflow mixture.
            inf->SetFracs(molefracs);
            dynamic_cast<PSR*>(reac)->SetInflow(*inf);
        } else {
            throw runtime_error("Inflow conditions must be defined for a PSR "
                                "(Mops, Settings_IO::readReactor)");
        }

        // Read the residence time.
        subnode = node.GetFirstChild("residencetime");
        if (subnode != NULL) {
            dynamic_cast<PSR*>(reac)->SetResidenceTime(cdble(subnode->Data()));
        }
    }

    return reac;
}

// Reads time intervals from given XML node.
void Settings_IO::readTimeIntervals(const CamXML::Element &node, 
                                           std::vector<TimeInterval> &times)
{
    const CamXML::Element *subnode;
    vector<CamXML::Element*> nodes;
    vector<CamXML::Element*>::const_iterator i;
    const CamXML::Attribute *attr;
    string str;

    // Read default splits per step.
    unsigned int splits = 1;
    attr = node.GetAttribute("splits");
    if (attr != NULL) {
        splits = (unsigned int)cdble(attr->GetValue());
    }


    // Read start time (first interval).
    TimeInterval *ti = NULL;
    subnode = node.GetFirstChild("start");
    if (subnode != NULL) {
        ti = new TimeInterval();
        ti->SetStartTime(cdble(subnode->Data()));
    } else {
        throw runtime_error("Time intervals must have start time defined "
                            "(Mops, Settings_IO::readTimeIntervals)");
    }

    node.GetChildren("time", nodes);
    for (i=nodes.begin(); i!=nodes.end(); ++i) {
        // Need to set the end time of the previous time 
        // interval before creating a new one.
        ti->SetEndTime(cdble((*i)->Data()));

        // Check for steps attribute.
        attr = (*i)->GetAttribute("steps");
        if (attr != NULL) {
            ti->SetStepCount((unsigned int)cdble(attr->GetValue()));
        } else {
            if (ti != NULL) delete ti;
            throw runtime_error("Time interval must define number of steps "
                                "(Mops, Settings_IO::readTimeIntervals)");
        }
        
        // Check for splits attribute.
        attr = (*i)->GetAttribute("splits");
        if (attr != NULL) {
            ti->SetSplittingStepCount((unsigned int)cdble(attr->GetValue()));
        } else {
            ti->SetSplittingStepCount(splits);
        }

        // Add the completed time interval to the vector.
        times.push_back(*ti);

        // Create a new time interval and set the start time.
        delete ti; ti = new TimeInterval();
        ti->SetStartTime(cdble((*i)->Data()));
    }

    // Clear memory used to load time intervals.
    if (ti != NULL) delete ti; ti = NULL;
}

// Reads simulation output parameters from given XML node.
void Settings_IO::readOutput(const CamXML::Element &node, Simulator &sim)
{
    const CamXML::Element *subnode;
    vector<CamXML::Element*> nodes;
    vector<CamXML::Element*>::const_iterator i;
    const CamXML::Attribute *attr;
    string str;

    // CONSOLE OUTPUT.

    subnode = node.GetFirstChild("console");
    if (subnode != NULL) {
        // Check the console interval.
        attr = subnode->GetAttribute("interval");
        if (attr != NULL) {
            sim.SetConsoleInterval((unsigned int)cdble(attr->GetValue()));
        }

        // Read tabular console output parameters.
        subnode = subnode->GetFirstChild("tabular");
        if (subnode != NULL) {
            subnode->GetChildren("column", nodes);
            for (i=nodes.begin(); i!=nodes.end(); ++i) {
                sim.AddConsoleVariable((*i)->Data());
            }
        } else {
            throw runtime_error("No tabular data defined for console output "
                                "(Mops, Settings_IO::readOutput)");
        }
    } else {
        throw runtime_error("No console output defined "
                            "(Mops, Settings_IO::readOutput)");
    }


    // FILE OUTPUT.

    // Read the output file name.
    subnode = node.GetFirstChild("filename");
    if (subnode != NULL) {
        sim.SetOutputFile(subnode->Data());
    } else {
        throw runtime_error("No output file name defined "
                            "(Mops, Settings_IO::readOutput)");
    }
}

// Returns the temperature in K by reading the value from the given
// XML node and checking the units.
real Settings_IO::readTemperature(const CamXML::Element &node)
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
    real T = cdble(node.Data());
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
real Settings_IO::readPressure(const CamXML::Element &node)
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
    real P = cdble(node.Data());
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
