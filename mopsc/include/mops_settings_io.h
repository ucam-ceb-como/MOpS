/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The Settings_IO class facillitates the reading of mops simulation
    settings from input files.  Currently only XML formatted input
    files are supported.
*/

#ifndef MOPS_SETTINGS_IO_H
#define MOPS_SETTINGS_IO_H

#include "mops_params.h"
#include "mops_reactor.h"
#include "mops_settings.h"
#include "mops_timeinterval.h"
#include "camxml.h"
#include <vector>
#include <string>

namespace Mops
{
class Settings_IO
{
public:
    // Constructors.
    Settings_IO(void); // Default constructor.

    //Destructors.
    ~Settings_IO(void); // Default destructor.

    // An enumeration of different possible reactor types that
    // may be specified in the settings file.
    enum ReactorType {Batch, PSR, ShockTube};

    // SETTINGS FILE READING.

    // Loads an XML document into the class.  This operation needs
    // to be performed before settings can be acquired.
    static Reactor * LoadFromXML_V1(
        const std::string &filename,      // Input file name.
        Reactor *reac,                    // The reactor to be simulated.
        std::vector<TimeInterval> &times, // Vector of output time intervals.
        Settings &settings,               // General settings incl. output settings.
        const Mechanism &mech             // Mechanism used to define reactor.
        );
};
};

#endif
