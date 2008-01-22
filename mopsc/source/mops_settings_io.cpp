#include "mops_params.h"
#include "mops_settings_io.h"
#include "mops_timeinterval.h"
#include "mops_reactor.h"
#include "mops_mixture.h"
#include "camxml.h"
#include "string_functions.h"

#include <vector>
#include <string>
#include <exception>

using namespace Mops;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Settings_IO::Settings_IO(void)
{
}

// Default destructor.
Settings_IO::~Settings_IO(void)
{
}


// XML SETTINGS FILES.

Reactor * Settings_IO::LoadFromXML_V1(const std::string &filename, 
                                 Mops::Reactor *reac, 
                                 std::vector<TimeInterval> &times,
                                 Settings &settings,
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
                reac = new Reactor();
            } else if (str.compare("psr") == 0) {
                // This is a perfectly-stirred reactor.
                throw invalid_argument("PSR not currently implemented");
            } else if (str.compare("shocktube") == 0) {
                // This is a shocktube (const. V batch).
                throw invalid_argument("Shocktube not currently implemented");
            } else {
                // Default reactor is a batch reactor.
                reac = new Reactor();
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

        // Assign the mechanism to the reactor.
        reac->SetMechanism(&mech);


        // REACTOR INITIAL CONDITIONS.

        // Create a new Mixture object.
        Mixture *mix = new Mixture();
        mix->SetSpecies(&mech.Species());

        vector<real> molefracs(reac->Mechanism()->Species().size(), 0.0);

        // Fill the mixture object.
        root->GetChildren("init", nodes);
        for (i=nodes.begin(); i!=nodes.end(); i++) {
            // Check the ID attribute.
            attr = (*i)->GetAttribute("id");
            if (attr != NULL) {
                str = attr->GetValue();

                if (str.compare("T")==0) {
                    // This is the initial temperature.
                    mix->SetTemperature(atof((*i)->Data().c_str()));
                } else if (str.compare("P")==0) {
                    // This is the initial pressure (Remember to convert to Pa).
                    mix->SetPressure(1.0e5 * atof((*i)->Data().c_str()));
                } else {
                    // This should be a species initial concentration.

                    // Get the species index.
                    int j = mech.FindSpecies(str);
                    
                    // Set the species mole fraction.
                    molefracs[j] = atof((*i)->Data().c_str());
                }
            } else {
                throw invalid_argument("Initial condition must have ID!");
            }
        }

        // Assign the species mole fraction vector to the reactor mixture.
        mix->SetFracs(molefracs);
        reac->Fill(mix);


        // GLOBAL SETTINGS.

        // Read the relative error tolerance.
        node = root->GetFirstChild("rtol");
        if (node != NULL) {
            settings.SetRTOL(atof(node->Data().c_str()));
        }

        // Read the absolute error tolerance.
        node = root->GetFirstChild("atol");
        if (node != NULL) {
            settings.SetATOL(atof(node->Data().c_str()));
        }


        // TIME INTERVALS.

        TimeInterval *ti;

        root->GetChildren("time", nodes);
        for (i=nodes.begin(); i!=nodes.end(); i++) {
            attr = (*i)->GetAttribute("start");

            if (attr != NULL) {
                // This is the start time.  Create a new time inteval
                // and set the start time.
                ti = new TimeInterval();
                ti->SetStartTime(atof((*i)->Data().c_str()));
            } else {
                // This is not the start time.  Need to set the end time of
                // the previous time interval before creating a new one.
                ti->SetEndTime(atof((*i)->Data().c_str()));

                // Check for steps attribute.
                attr = (*i)->GetAttribute("steps");
                if (attr != NULL) {
                    ti->SetStepCount((unsigned int)atof(attr->GetValue().c_str()));
                } else {
                    if (ti != NULL) delete ti;
                    throw invalid_argument("Time interval must define number of steps.");
                }
                
                // Add the completed time interval to the vector.
                times.push_back(*ti);

                // Create a new time interval and set the start time.
                delete ti; ti = new TimeInterval();
                ti->SetStartTime(atof((*i)->Data().c_str()));
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
                settings.SetConsoleInterval(atof(attr->GetValue().c_str()));
            }

            // Read the column variables.
            str = node->Data();
            vector<string> cvars;
            Strings::split(str, cvars, " ");
            for (int j=0; j<cvars.size(); j++) {
                settings.AddConsoleVariable(cvars[j]);
            }
        }


        // FILE OUTPUT.

        // Read the output file name.
        node = root->GetFirstChild("output");
        if (node != NULL) {
            settings.SetOutputFile(node->Data());
        }
    }

    return reac;
}