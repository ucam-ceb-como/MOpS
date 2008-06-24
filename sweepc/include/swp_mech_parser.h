/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    A class which reads mechanisms from XML files into Mechanism
    objects.
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
