/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Settings_IO class facillitates the reading of mops simulation
    settings from input files.  Currently only XML formatted input
    files are supported.

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

#ifndef MOPS_SETTINGS_IO_H
#define MOPS_SETTINGS_IO_H

#include "mops_params.h"
#include "mops_reactor.h"
#include "mops_solver.h"
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
    // enum ReactorType {Batch, PSR, ShockTube};

    // SETTINGS FILE READING.

    // Loads an XML document into the class.  This operation needs
    // to be performed before settings can be acquired.
    static Reactor *const LoadFromXML_V1(
        const std::string &filename,      // Input file name.
        Reactor *reac,                    // The reactor to be simulated.
        std::vector<TimeInterval> &times, // Vector of output time intervals.
        Solver &solver,                   // General settings incl. output settings.
        const Mechanism &mech             // Mechanism used to define reactor.
        );

    // Read a new-format XMl file settings file.
    static Reactor *const LoadFromXML(
        const std::string &filename,      // Input file name.
        Reactor *reac,                    // The reactor to be simulated.
        std::vector<TimeInterval> &times, // Vector of output time intervals.
        Solver &solver,                   // General settings incl. output settings.
        const Mechanism &mech             // Mechanism used to define reactor.
        );

private:

    // V2 SETTINGS FILE SECTIONS.

    // Reads global simulation settings from the given XML node.
    static void readGlobalSettings(
        const CamXML::Element &node, // Root XML node containing simulation settings.
        Solver &solver               // Solver object into which to read global settings.
        );

    // Reads the reactor initial settings from the given XML node.
    static Reactor *const readReactor(
        const CamXML::Element &node, // XML node containing reactor.
        const Mechanism &mech        // Mechanism to define reactor mixture.
        );

    // Reads time intervals from given XML node.
    static void readTimeIntervals(
        const CamXML::Element &node,     // XML node containing time intervals.
        std::vector<TimeInterval> &times // Vector of output time intervals.
        );

    // Reads simulation output parameters from given XML node.
    static void readOutput(
        const CamXML::Element &node, // XML node containing output parameters.
        Solver &solver               // General settings incl. output settings.
        );

    // Returns the temperature in K by reading the value from the given
    // XML node and checking the units.
    static real readTemperature(const CamXML::Element &node);

    // Returns the pressure in Pa by reading the value from the given
    // XML node and checking the units.
    static real readPressure(const CamXML::Element &node);
};
};

#endif
