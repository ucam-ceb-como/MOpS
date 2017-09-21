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
#include "swp_pah_inception.h"
#include "swp_dimer_inception.h"
#include "swp_silicon_inception.h"
#include "swp_surface_reaction.h"
#include "swp_condensation.h"
#include "swp_maths_functional.h"
#include "swp_silica_interparticle.h"
#include "camxml.h"
#include <string>

namespace Sweep
{

//forward declaration
namespace Processes {
    class DiffusionProcess;
    class AdvectionProcess;
    class ConstantInception;
}

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

    // PARTICLE DEFINITIONS.

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

    // INCEPTIONS.

    // Reads inception processes from a sweep mechanism XML file.
    static void readInceptions(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // Reads an inception process from an XML element.
    static void readInception(
        CamXML::Element &xml,     // CamXML document pre-constructed from file.
        Processes::DimerInception &icn // Inception to construct from XML.
        );

    //! Reads a silicon inception process from an XML element.
    static void readSiliconInception(
        CamXML::Element &xml,     // CamXML document pre-constructed from file.
        Processes::SiliconInception &icn // Inception to construct from XML.
        );


    /*!
    * reads a PAH inception process
    @param[in] xml      XML node of type pahinception
    @param[in,out] icn      Inception process to initialise
    */
    static void readPAHInception(CamXML::Element &xml, Processes::PAHInception &icn);

    /*!
    * reads all PAH inception process
    * node that do not contain pahinception will be ignored
    @param[in] xml      XML node with one or more children of type pahinception
    @param[in,out] mech  Mechanism to which the processes will be added
    */
    static void readPAHInceptions(CamXML::Document &xml, Sweep::Mechanism &mech);

    /*!
    * reads all constant inception processes
    * nodes that do not contain constant inception processes will be ignored
    @param[in] xml      XML node with one or more children of type constantinception
    @param[in,out] mech  Mechanism to which the processes will be added
    */
    static void readConstantInceptions(CamXML::Document &xml, Sweep::Mechanism &mech);

    // SURFACE REACTIONS.

    // Reads surface reactions from a sweep mechanism XML file.
    static void readSurfRxns(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // Reads a surface reaction from a sweep mechanism XML file.
    static void readSurfRxn(
        CamXML::Element &xml,           // CamXML document pre-constructed from file.
        Processes::SurfaceReaction &rxn // Reaction to construct from XML.
        );

    // CONDENSATIONS.

    // Reads condensation processes from a sweep mechanism XML file.
    static void readCondensations(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // Reads a condensation process from a sweep mechanism XML file.
    static void readCondensation(
        CamXML::Element &xml,         // CamXML document pre-constructed from file.
        Processes::Condensation &cond // Condensation to construct from XML.
        );

	// INTER PARTICLE REACTIONS.

    // Reads inter particle reaction processes from a sweep mechanism XML file.
    static void readInterParticles(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // Reads an inter particle process from a sweep mechanism XML file.
    static void readInterParticle(
        CamXML::Element &xml,         // CamXML document pre-constructed from file.
        Processes::InterParticle &intpart // Interparticle reaction to construct from XML.
        );


    // COAGULATION
    //! Allow user to choose coagulation kernel via a sweep mechanism XML file.
    static void readCoagulation(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // COAGULATION
    //! Allow user to choose coagulation kernel via a sweep mechanism XML file.
    static void readFragmentation(
        CamXML::Document &xml, // CamXML document pre-constructed from file.
        Mechanism &mech        // Mechanism to construct from XML.
        );

    // TRANSPORT
    //! Read any diffusion processes from the XML document
    static void readDiffusionProcs(
        CamXML::Document &xml,
        Mechanism &mech
        );

    //! Read one diffusion process from its XML description
    static void readDiffusionProc(
        CamXML::Element &xml,
        Processes::DiffusionProcess& tran
        );

    //! Read any advection processes from the XML document
    static void readAdvectionProcs(
        CamXML::Document &xml,
        Mechanism &mech
        );

    //! Read one advection process from its XML description
    static void readAdvectionProc(
        CamXML::Element &xml,
        Processes::AdvectionProcess& tran
        );

    // REACTION SHARED COMPONENTS.

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


    // AUXILLIARY BITS AND BOBS.

    // Reads a maths functional from XML and creates an object to contain it.
    static Sweep::Maths::Functional *const readFunctional(
        CamXML::Element &xml // CamXML element containing the functional def'n.
        );

};
};

#endif
