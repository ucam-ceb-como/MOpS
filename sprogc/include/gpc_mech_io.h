/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This class reads/writes Sprog mechanism from file.  Currently only CHEMKIN formatted
    files are supported.
*/

#ifndef GPC_MECH_PARSER_H
#define GPC_MECH_PARSER_H

#include <string>
#include "gpc_mech.h"

namespace Sprog
{
namespace IO
{
class MechanismParser
{
public:
    // Reads a CHEMKIN input file.
    static void ReadChemkin(
        const std::string &filename,    // File name of the CHEMKIN input file.
        Sprog::Mechanism &mech,         // Mechanism object to build using data in file.
        const std::string &thermofile); // File name of thermo data file (optional).

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

    // Parses a file stream of a CHEMKIN input file.
    static void parseCK(std::ifstream &fin,     // File stream.
                        Sprog::Mechanism &mech, // Mechanism into which to read file.
                        CK_STATUS &status);     // Status & parsing information.

    // Parse the element data in a CHEMKIN input file.
    static void parseCK_Elements(std::ifstream &fin,     // File stream.
                                 Sprog::Mechanism &mech, // Mechanism to receive element information.
                                 CK_STATUS &status);       // Status and parsing information.

    // Parse the species data in a CHEMKIN input file (elements must already have been read).
    static void parseCK_Species(std::ifstream &fin,     // File stream.
                                Sprog::Mechanism &mech, // Mechanism to receive species information.
                                CK_STATUS &status);     // Status and parsing information.

    // Parse the thermo data from a CHEMKIN input file (elements and species must already have
    // been read).
    static void parseCK_Thermo(std::ifstream &fin,     // File stream.
                               Sprog::Mechanism &mech, // Mechanism to receive thermo data.
                               CK_STATUS &status);     // Status and parsing information.

    // Parse the thermo data from a separate file (given by filename).  Elements and species
    // must already have been read.
    static void parseCK_Thermo(std::string &filename,  // File name.
                               Sprog::Mechanism &mech, // Mechanism to receive thermo data.
                               CK_STATUS &status);     // Status and parsing information.

    // Parse the reaction data from a CHEMKIN input file.
    static void parseCK_Reactions(std::ifstream &fin,     // File stream.
                                  Sprog::Mechanism &mech, // Mechanism to receive reaction data.
                                  CK_STATUS &status);     // Status and parsing information.

    // Parse a single reaction object from a CHEMKIN formatted string.  Returns pointer to
    // new reaction if successful.
    static Sprog::Kinetics::Reaction *const parseCK_Reaction(
        const std::string &rxndef, // String containing the reaction definition.
        Sprog::Mechanism &mech,    // Mechanism into which the reaction will be inserted (for species definitions).
        CK_STATUS &status);        // Status and parsing information.

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
        std::string &thirdbody);          // Returns the name of the third body, if detected.

    // Reads auxilliary reaction data from a CHEMKIN formatted string.  Returns true
    // if auxilliary data was found in the string.
    static bool parseCK_RxnAux(
        const std::string &rxndef,              // String containing the auxilliary information.
        Sprog::Kinetics::Reaction *last_rxn,    // Reaction object into which to read auxilliary information.
        const Sprog::Mechanism &mech,           // Mechanism for which the reaction is defined.
        const Sprog::Kinetics::ARRHENIUS scale, // Scaling parameters for Arrhenius cofficients to convert to correct units.
        CK_STATUS &status);                     // Status and parsing information.

    // Parse the units data from the REACTION line in a chemkin formatted mechanism file.
    static void parseCK_Units(const std::string &rxndef,          // String containing the REACTION statement.
                              Sprog::Kinetics::ARRHENIUS &scale); // Scaling factors.

};
};
};

#endif
