/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    A class which reads mechanisms from XML files into Mechanism
    objects.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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

#ifndef SWEEP_MECH_PARSER_H
#define SWEEP_MECH_PARSER_H

#include "swp_mechanism.h"
#include "swp_process.h"
#include "camxml.h"
#include <string>

namespace Sweep
{
class MechParser
{
public:
    // Reads the given XML file into the given mechanism object.
    static void Read(
        const std::string &filename, // Name of XML file.
        Mechanism &mech              // Mechanism to construct.
        );

private:
    // Reads a version 1 sweep mechanism.  This is also the default if no
    // version is specified.
    static void readV1(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // Reads components from a sweep mechanism XML file.
    static void readComponents(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // Reads tracker variables from a sweep mechanism XML file.
    static void readTrackers(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // Reads inception processes from a sweep mechanism XML file.
    static void readInceptions(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // Reads surface reactions from a sweep mechanism XML file.
    static void readSurfRxns(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech              // Mechanism to construct from XML.
        );

    // Reads condensation processes from a sweep mechanism XML file.
    static void readCondensations(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // Reads reactants into a process.
    static void readReactants(
        CamXML::Element &xml,    // CamXML element containing the process definition.
        Processes::Process &proc // Process into which to read the reactants.
        );

    // Reads products into a process.
    static void readProducts(
        CamXML::Element &xml,    // CamXML element containing the process definition.
        Processes::Process &proc // Process into which to read the products.
        );

    // Reads reactant masses and diameters.  This is required by inceptions
    // and condensation to calculate collision rates.
    static void readReactantMDs(
        CamXML::Element &xml, // CamXML element containing the process definition.
        fvector &mass, // Vector of reactant masses.
        fvector &diam  // Vector of reactant diameters.
        );

    // Reads composition changes into a particle process.
    static void readCompChanges(
        CamXML::Element &xml,            // CamXML element containing the process definition.
        Processes::ParticleProcess &proc // Process into which to read the composition changes.
        );

    // Reads tracker variable changes into a particle process.
    static void readTrackChanges(
        CamXML::Element &xml,            // CamXML element containing the process definition.
        Processes::ParticleProcess &proc // Process into which to read the composition changes.
        );
};
};

#endif
