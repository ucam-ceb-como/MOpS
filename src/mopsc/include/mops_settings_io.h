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

#include "mops_reactor_network.h"
#include "mops_mixture.h"
#include "mops_timeinterval.h"
#include "swp_particle.h"
#include "camxml.h"
#include <vector>
#include <list>
#include <string>

namespace Sweep {
    class Mechanism;
}

namespace Mops
{
    // forward declarations
    class Simulator;
    class Solver;
    class Reactor;
    class Mechanism;
    
namespace Settings_IO
{
    // An enumeration of different possible reactor types that
    // may be specified in the settings file.
    // enum ReactorType {Batch, PSR, ShockTube};

    // SETTINGS FILE READING.

    // Loads an XML document into the class.  This operation needs
    // to be performed before settings can be acquired.
    Reactor *const LoadFromXML_V1(
        const std::string &filename,      // Input file name.
        Reactor *reac,                    // The reactor to be simulated.
        std::vector<TimeInterval> &times, // Vector of output time intervals.
        Simulator &sim,                   // General settings incl. output settings. 
        Solver &solver,                   // The reactor solver (to set numerical params).
        const Mechanism &mech             // Mechanism used to define reactor.
        );

    // Read a new-format XMl file settings file.
    Reactor *const LoadFromXML(
        const std::string &filename,      // Input file name.
        Reactor *reac,                    // The reactor to be simulated.
        std::vector<TimeInterval> &times, // Vector of output time intervals.
        Simulator &sim,                   // General settings incl. output settings. 
        Solver &solver,                   // The reactor solver (to set numerical params).
        Mechanism &mech                   // Mechanism used to define reactor.
        );

    // Reads time intervals from given XML node.
    void readTimeIntervals(
        const CamXML::Element &node,     // XML node containing time intervals.
        std::vector<TimeInterval> &times // Vector of output time intervals.
        );

    //! Read initial particles from a file into a list
    Sweep::PartPtrList ReadInitialParticles(const CamXML::Element &node,
                                            const Sweep::Mechanism & particle_mech);

    //! Read initial particles from a file into a list
    //! This function is used to create detailed particles by sampling from a distribution
    Sweep::PartPtrList ReadInitialParticlesDetailed(const CamXML::Element &node,
		const Sweep::Mechanism & particle_mech);

    //! Read limits that define extreme particles to be excluded from particle population statistics
    void ReadStatsBound(const CamXML::Element &node, Sweep::PropID &property_id,
                        double &lower_bound, double &upper_bound);

    //! Read a reactor network from an XML file.
    Mops::ReactorNetwork* LoadNetwork(
            const std::string &filename,
            std::vector<TimeInterval> &times,
            Simulator &sim,
            Solver &solver,
            Mechanism &mech);

} //namespace Settings_IO
} //namespace Mops

#endif
