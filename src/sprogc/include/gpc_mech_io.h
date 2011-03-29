/*
  Author(s):      Matthew Celnik (msc37),
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This class reads/writes Sprog mechanism from file.  Currently only CHEMKIN formatted
    files are supported.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
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

#ifndef GPC_MECH_PARSER_H
#define GPC_MECH_PARSER_H

#include <string>
#include "gpc_mech.h"

namespace Sprog
{
namespace IO
{
//! Read species properties and reaction data from CHEMKIN format files
class MechanismParser
{
public:
    // Reads a CHEMKIN input file.
    static void ReadChemkin(
        const std::string &filename,    // File name of the CHEMKIN input file.
        Sprog::Mechanism &mech,         // Mechanism object to build using data in file.
        const std::string &thermofile,  // File name of thermo data file (optional).
        const int verbose=0             // Set >0 to print parser messages to console.
        );

	//! Read in a mechanism and species properties including transport properties
	static void ReadChemkin(
        const std::string &filename,
        Sprog::Mechanism &mech,
        const std::string &thermofile,
		const std::string &transFile,
		const int verbose = 0);

	//! method to read transport data
	static void ReadTransport(const std::string &, Sprog::Mechanism & );

private:
    // CHEMKIN FILE RELATED THINGS.

    // Enumeration of status flags for reading CHEMKIN files.
    enum CK_PARSE_STATUS {FindKey,
                          BeginParseEl, ParseEl, ParseElWt, ParseElComment,
                          BeginParseSp, ParseSp, ParseSpComment,
                          BeginParseRxn, ParseRxn, ParseRxnComment, BeginParseRxnComment,
                          ParseTherm,
                          End, Fail=-1};

    // Structure to aid reading of CHEMKIN input files.
    struct CK_STATUS {
        CK_PARSE_STATUS Status;
        bool ReadElements;
        bool ReadSpecies;
        bool ReadReactions;
        bool ReadThermo;
        std::string ThermoFile;
        Sprog::Kinetics::ARRHENIUS Scale;  // Arrhenius constant scaling factors.
    };

    // Structure to locate keyworded sections in the CK file.
    struct KEY_POS{
        unsigned int begin;  // First character in the keyword.
        unsigned int end;    // First character in the END keyword.
        unsigned int length; // Keyword length.
        unsigned int line;   // Line number of this keyword.
    };

    //! Reads a CHEMKIN file stream into a std::string
    static void loadCK_File(
        std::ifstream &fin, // File stream to read.
        std::string &out    // Output string to hold file data.
        );

    //! Parses a file stream of a CHEMKIN input file.
    static void parseCK(
        std::ifstream &fin,     // File stream.
        Sprog::Mechanism &mech, // Mechanism into which to read file.
        CK_STATUS &status,      // Status & parsing information.
        int verbose=0           // Set >0 to print parser messages to console.
        );

    // Get positions of a CK keyword and END keyword.
    static KEY_POS getCK_KeyPos(
        const std::string &key,  // Keyword to be found.
        const std::string &ckstr // CK file string to search.
        );

    // Extract element names from CHEMKIN string.
    static void extractCK_Elements(
        const std::string &ckstr, // CK file string to search.
        std::string &elements,    // Return string containing elements.
        unsigned int &lineno      // The line number of the ELEM/ELEMENTS keyword.
        );


    //// extract only reactions string from chemkin string
    //std::string extract_CK_reactions_str(std::string &ckstr);
    //// extract only thermo string from chemkin string
    //std::string extract_CK_thermo_str(std::string &ckstr);

    // Parse the element data in a CHEMKIN input file.
    static unsigned int parseCK_Elements(
        const std::string &elements, // CHEMKIN string for element definitions.
        Sprog::Mechanism &mech,      // Mechanism to receive element information.
        unsigned int lineno,         // Line number of ELEM/ELEMENTS keyword.
        CK_STATUS &status,           // Status and parsing information.
        int verbose=0               // Set >0 to print parser messages to console.
        );

    // Extract species names from CHEMKIN string.
    static void extractCK_Species(
        const std::string &ckstr, // CK file string to search.
        std::string &species,     // Return string containing species.
        unsigned int &lineno      // The line number of the SPEC/SPECIES keyword.
        );

    // Parse the species data in a CHEMKIN input file
    // (elements must already have been read).
    static unsigned int parseCK_Species(
        const std::string &species, // Species names in a string.
        Sprog::Mechanism &mech,     // Mechanism to receive species information.
        unsigned int lineno,        // Line number of SPEC/SPECIES keyword.
        CK_STATUS &status,          // Status and parsing information.
        int verbose=0               // Set >0 to print parser messages to console.
        );

    // Extract thermo data string from CHEMKIN string.
    static void extractCK_Thermo(
        const std::string &ckstr, // CK file string to search.
        std::string &thermo,      // Return string containing thermo data.
        unsigned int &lineno      // The line number of the THER/THERMO keyword.
        );

    //! Parse the thermo data from a separate file (given by filename).
    static unsigned int parseCK_Thermo(
        const std::string &thermo, // File name.
        Sprog::Mechanism &mech,    // Mechanism to receive thermo data.
        unsigned int lineno,       // Line number of THER/THERMO keyword.
        CK_STATUS &status,         // Status and parsing information.
        int verbose=0              // Set >0 to print parser messages to console.
        );

    // Extract reaction data string from CHEMKIN string.
    static void extractCK_Reactions(
        const std::string &ckstr, // CK file string to search.
        std::string &reac,        // Return string containing reaction data.
        unsigned int &lineno      // The line number of the REAC/REACTIONS keyword.
        );

    // Parse the reaction data from a CHEMKIN input file.
    static unsigned int parseCK_Reactions(
        const std::string &reac, // File stream.
        Sprog::Mechanism &mech,  // Mechanism to receive reaction data.
        unsigned int lineno,     // The line number of the REAC/REACTIONS keyword.
        CK_STATUS &status,       // Status and parsing information.
        int verbose=0            // Set >0 to print parser messages to console.
        );

    // Parse a single reaction object from a CHEMKIN formatted string.  Returns pointer to
    // new reaction if successful.
    static Sprog::Kinetics::Reaction *const parseCK_Reaction(
        const std::string &rxndef, // String containing the reaction definition.
        Sprog::Mechanism &mech,    // Mechanism into which the reaction will be inserted (for species definitions).
        unsigned int lineno,       // Line number of reaction definition.
        CK_STATUS &status,         // Status and parsing information.
        int verbose=0              // Set >0 to print parser messages to console.
        );

    // Reads strings of reactant/product species within a reaction string and separates
    // the species from the stoichiometry.  Also checks for third-bodies and fall-off
    // reactions.
    static void parseCK_RxnSpStoich(
        const std::string &sp,            // String containing the reactant/product species.
        const Sprog::Mechanism &mech,     // Mechanism in which the species are defined.
        std::vector<Sprog::Stoich> &mui,  // Vector of integer coefficient stiochiometry.
        std::vector<Sprog::Stoichf> &muf, // Vector of real coefficient stoichiometry.
        bool &isthirdbody,                // Returns true if a third body was detected.
        bool &isfalloff,                  // Returns true if this reaction has fall-off parameters.
        std::string &thirdbody,           // Returns the name of the third body, if detected.
        unsigned int lineno               // Line number on which the reaction is defined.
        );

    // Reads auxilliary reaction data from a CHEMKIN formatted string.  Returns true
    // if auxilliary data was found in the string.
    static bool parseCK_RxnAux(
        const std::string &rxndef,              // String containing the auxilliary information.
        Sprog::Kinetics::Reaction *last_rxn,    // Reaction object into which to read auxilliary information.
        const Sprog::Mechanism &mech,           // Mechanism for which the reaction is defined.
        const Sprog::Kinetics::ARRHENIUS scale, // Scaling parameters for Arrhenius cofficients to convert to correct units.
        CK_STATUS &status,                      // Status and parsing information.
        unsigned int lineno,                    // Line number on which the reaction aux info is defined.
        int verbose=0                           // Set >0 to print parser messages to console.
        );

    // Parse the units data from the REACTION line in a chemkin formatted mechanism file.
    static void parseCK_Units(
        const std::string &rxndef,          // String containing the REACTION statement.
        Sprog::Kinetics::ARRHENIUS &scale); // Scaling factors.

};
};
};

#endif
