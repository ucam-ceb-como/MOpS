/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This class reads/writes Sprog mechanism from file.  Different file formats
    are supported.
*/

#ifndef GPC_MECH_IO_H
#define GPC_MECH_IO_H

#include <string>
#include "gpc_mech.h"

namespace Sprog
{
namespace IO
{
class Mechanism_IO
{
public:
    // Constructors.
    Mechanism_IO(void); // Default constructor.

    // Destructors.
    virtual ~Mechanism_IO(void); // Default destructor.

    // Reads a CHEMKIN input file.
    static void ReadChemkin(const std::string &filename, Sprog::Mechanism &mech, const std::string &thermofile);

private:
    // Enumeration of status flags for reading CHEMKIN files.
    enum CK_PARSE_STATUS {FindKey, 
                          BeginParseEl, ParseEl, ParseElWt, ParseElComment,
                          BeginParseSp, ParseSp, ParseSpComment,
                          BeginParseRxn, ParseRxn, ParseRxnComment,
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
    };

    // Parses a file stream of a CHEMKIN input file.
    static void parseCK(std::ifstream &fin, Sprog::Mechanism &mech, CK_STATUS &stat);
    static void parseCK_Elements(std::ifstream &fin, Sprog::Mechanism &mech, CK_STATUS &stat);
    static void parseCK_Species(std::ifstream &fin, Sprog::Mechanism &mech, CK_STATUS &stat);
    static void parseCK_Thermo(std::ifstream &fin, Sprog::Mechanism &mech, CK_STATUS &stat);
    static void parseCK_Thermo(std::string &filename, Sprog::Mechanism &mech, CK_STATUS &stat);
    static void parseCK_Reactions(std::ifstream &fin, Sprog::Mechanism &mech, CK_STATUS &stat);
    static Sprog::Kinetics::Reaction *const parseCK_Reaction(const std::string &rxndef, 
        Sprog::Mechanism &mech, CK_STATUS &stat);
    static void parseCK_RxnSpStoich(const std::string &sp, const Sprog::Mechanism &mech,
        std::vector<Sprog::Stoich> &mui, std::vector<Sprog::Stoichf> &muf, 
        bool &isthirdbody, bool &isfalloff, std::string &thirdbody);
    static bool parseCK_RxnAux(const std::string &rxndef, Sprog::Kinetics::Reaction *last_rxn, 
        const Sprog::Mechanism &mech, CK_STATUS &stat);
};
};
};

#endif