#include "morestrings.h"
#include "gpc_mech_io.h"
#include "gpc_mech.h"

#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>

using namespace Sprog::IO;
using namespace std;

Mechanism_IO::Mechanism_IO(void)
{
}

Mechanism_IO::~Mechanism_IO(void)
{
}

// Reads a mechanism from a CHEMKIN input file.
void Mechanism_IO::ReadChemkin(const std::string &filename, Sprog::Mechanism &mech, const std::string &thermofile)
{
    // Open file for reading.
    ifstream fin; //(filename.c_str(), ios::in);
    fin.open(filename.c_str(), ios::in);

    if (fin.good()) {
        CK_STATUS stat;
        stat.Status = FindKey;
        stat.ReadElements = false;
        stat.ReadSpecies = false;
        stat.ReadReactions = false;
        stat.ReadThermo = false;
        stat.ThermoFile = thermofile;

        try {
            parseCK(fin, mech, stat);
        } catch (exception &e) {
            fin.close();
            throw e;
        }
        
        fin.close();
        return;
    } else {
        // Failed to open file.
        throw invalid_argument(string("Could not open CHEMKIN file: ").append(filename));
    }
}


// CHEMKIN PARSING ROUTINES.

// Loads a mechanism from a CHEMKIN formatted file stream.
void Mechanism_IO::parseCK(std::ifstream &fin, Sprog::Mechanism &mech, Sprog::IO::Mechanism_IO::CK_STATUS &stat)
{
    char c;
    string line;
    real val;
    int i=0;

    // Positions in the file stream at the beginning of the elements, species
    // reactions and thermo data.
    streamoff iel=0, isp=0, irxn=0, ithrm=0;

    Element * last_el;
    Species * last_sp;

    // Locate in the file stream the starting point of the elements, species
    // reactions and thermo data.
    while (fin.good()) {
        fin.get(c);
        i++;
        if ((c=='\n') || (c=='\r')) {
            // This is the end of a line.
            line = convertToCaps(line);
            if (line.substr(0,4).compare("ELEM")==0) {
                // This is the start of the elements.
                if (iel == 0) iel = (std::streamoff)i;
            } else if (line.substr(0,4).compare("SPEC")==0) {
                // This is the start of the species.
                if (isp == 0) isp = (std::streamoff)i;
            } else if (line.substr(0,4).compare("REAC")==0) {
                // This is the start of the reactions.
                if (irxn == 0) irxn = (std::streamoff)i;
            } else if (line.substr(0,4).compare("THER")==0) {
                // This is the start of the thermo data.
                if (ithrm == 0) ithrm = (std::streamoff)i;
            }
            line.clear();

            if ((iel!=0) && (isp!=0) && (irxn!=0) && ((ithrm!=0) || (stat.ThermoFile!=""))) {
                break;
            }
        } else {
            // Append character to line.
            line.append(&c, 1);
        }
    }

    // Read the elements.
    fin.seekg(iel);
    parseCK_Elements(fin, mech, stat);
    stat.ReadElements = true;

    // Read the species.
    fin.seekg(isp);
    parseCK_Species(fin, mech, stat);
    stat.ReadSpecies = true;

    // Read the thermo data.
    if (stat.ThermoFile == "") {
        // Read thermo data from this file.
        fin.seekg(ithrm);
        parseCK_Thermo(fin, mech, stat);
    } else {
        // Read thermo data from thermo file.
        parseCK_Thermo(stat.ThermoFile, mech, stat);
    }
    stat.ReadThermo = true;

    // Read the reactions.
    fin.seekg(irxn);
    parseCK_Reactions(fin, mech, stat);
}

void Mechanism_IO::parseCK_Elements(std::ifstream &fin, Sprog::Mechanism &mech, 
                                    Sprog::IO::Mechanism_IO::CK_STATUS &stat)
{
    char c;
    string tag;
    real val;

    Element * last_el;

    stat.Status = BeginParseEl;

    while((stat.Status!=End) && (stat.Status!=Fail) && (fin.good())) {
        // Get the next character from the input file.
        fin.get(c);

        switch(stat.Status) {
            case BeginParseEl:
                // Here we are looking for the start of an element name.
                if (isLetter(c)) {
                    // This is the start of an element name.
                    tag = c;
                    stat.Status = ParseEl;

                    // Clear the last read element.
                    last_el = NULL;
                } else if (c == '/') {
                    // This character means that the last element's Mol. Wt. is defined in
                    // the input file.  If the last element is valid, then we must read
                    // the mol. wt. from the file.
                    if (last_el != NULL) {
                        stat.Status = ParseElWt;
                        tag = "";
                    } else {
                        // Oh dear.  There is an element mol. wt. definition without an
                        // element.  This is an error.
                        stat.Status = Fail;
                        throw range_error("Invalid / character in element definition.");
                    }
                } else if (c == '!') {
                    // This is a comment.
                    stat.Status = ParseElComment;
                } else if (!isWhiteSpace(c)) {
                    // This must be an invalid character.
                    stat.Status = Fail;
                    throw range_error("Invalid character in element definition.");
                }
                break;
            case ParseEl:
                // We are now reading an element name.
                if (isLetterOrNum(c)) {
                    // This is a valid character for the element name.
                    tag.append(&c, 1);
                } else if (isWhiteSpace(c)) {
                    // We have read to the end of the element name.  We need to check
                    // that this element name is not actually a keyword.
                    if (convertToCaps(tag.substr(0,4)) == "ELEM") {
                        stat.Status = BeginParseEl;
                    } else if (convertToCaps(tag.substr(0,4)) == "SPEC") {
                        stat.Status = End;
                    } else if (convertToCaps(tag.substr(0,4)) == "REAC") {
                        stat.Status = End;
                    } else if (convertToCaps(tag.substr(0,4)) == "END") {
                        stat.Status = End;
                    } else {
                        // This is an element name, so add a new
                        // element to the mechanism and set its name.
                        last_el = mech.AddElement();
                        last_el->SetName(tag);
                        // Attempt to find the element weight in the library of
                        // know elements.
                        last_el->SetMolWtFromLibrary();
                        // Start searching for the next element.
                        stat.Status = BeginParseEl;
                    }
                    tag = "";
                } else if (c == '/') {
                    // This character means that the element Mol. Wt. is defined in
                    // the input file.  First task is to initialise a new element.
                    last_el = mech.AddElement();
                    last_el->SetName(tag);
                    // We must now read Mol. Wt from file.
                    stat.Status = ParseElWt;
                    tag = "";
                }
                break;
            case ParseElWt:
                // We are reading an element mol. wt. from the file.
                if (c == '/') {
                    // This is the end of the mol. wt. definition.  Attempt to convert it
                    // to a number.
                    val = atof(tag.c_str());
                    if (val > 0.0) {
                        // Set the mol. wt. of the last read element.
                        last_el->SetMolWt(val);
                    } else {
                        // The mol. wt. is invalid.
                        stat.Status = Fail;
                        throw range_error(string("An invalid mol. wt. was found for element ").append(last_el->Name()));
                    }
                } else {
                    // Add character to number.
                    tag.append(&c, 1);
                }
                break;
            case ParseElComment:
                // We are reading a comment in the element list.  A comment is only
                // terminated by a line ending.
                if ((c == '\n') || (c == '\r')) {
                    stat.Status = BeginParseEl;
                }
        }
    }
}

void Mechanism_IO::parseCK_Species(std::ifstream &fin, Sprog::Mechanism &mech, Sprog::IO::Mechanism_IO::CK_STATUS &stat)
{
    char c;
    string tag;

    Species * last_sp;

    stat.Status = BeginParseSp;

    while((stat.Status!=End) && (stat.Status!=Fail) && (fin.good())) {
        // Get the next character from the input file.
        fin.get(c);

        switch(stat.Status) {
            case BeginParseSp:
                // We are searching for the beginning of a species name.
                if (isLetter(c)) {
                    // This is the start of a species name.
                    tag = c;
                    stat.Status = ParseSp;

                    // Clear the last read species.
                    last_sp = NULL;
                } else if (c == '!') {
                    // This is a comment.
                    stat.Status = ParseSpComment;
                } else if (!isWhiteSpace(c)) {
                    // This must be an invalid character.
                    stat.Status = Fail;
                    throw range_error("Invalid character in species definition.");
                }
                break;
            case ParseSp:
                // We are now reading a species name.
                if (isLetterOrNum(c) || (c=='(') || (c==')') ||
                    (c=='*') || (c=='-') || (c=='+')) {
                    // This is a valid character for the species name.
                    tag.append(&c, 1);
                } else if (isWhiteSpace(c)) {
                    // We have read to the end of the species name.  We need to check
                    // that this species name is not actually a keyword.
                    if (convertToCaps(tag.substr(0,4)) == "ELEM") {
                        // This is an error. The elements must be defined before
                        // the species.
                        stat.Status = Fail;
                        throw range_error("Element definition after species definitions.");
                    } else if (convertToCaps(tag.substr(0,4)) == "SPEC") {
                        stat.Status = BeginParseSp;
                    } else if (convertToCaps(tag.substr(0,4)) == "REAC") {
                        stat.Status = End;
                    } else if (convertToCaps(tag.substr(0,4)) == "END") {
                        stat.Status = End;
                    } else {
                        // This is a species name, so add a new
                        // species to the mechanism and set its name.
                        last_sp = mech.AddSpecies();
                        last_sp->SetName(tag);
                        // Begin searching for the next species name.
                        stat.Status = BeginParseSp;
                    }
                }
                break;
            case ParseSpComment:
                // We are reading a comment in the species list.  A comment is only
                // terminated by a line ending.
                if ((c == '\n') || (c == '\r')) {
                    stat.Status = BeginParseSp;
                }
        }
    }
}

// Reads CHEMKIN formatted thermo data for all species in the given mechanism from
// the supplied file stream.
void Mechanism_IO::parseCK_Thermo(std::ifstream &fin, Sprog::Mechanism &mech, 
                                  Sprog::IO::Mechanism_IO::CK_STATUS &stat)
{
    int i, isp;
    char c, line[200];
    string tag, spname;

    real trange[3], lowT, highT, commT;
    Thermo::THERMO_PARAMS up; // Coeffs for upper T interval.
    Thermo::THERMO_PARAMS lp; // Coeffs for lower T interval.
    string els[4];
    int nels[4];

    Species * sp;

    stat.Status = ParseTherm;
    tag.resize(200,' ');

    // The first line holds the temperature ranges for 2 sets of coefficients.
    fin.getline(&line[0], 200);
    tag.clear();
    tag = line;
    trange[0] = atof(tag.substr(0,10).c_str());
    trange[1] = atof(tag.substr(10,10).c_str());
    trange[2] = atof(tag.substr(20,10).c_str());

    while((stat.Status!=End) && (stat.Status!=Fail) && (fin.good())) {
        // Get species name line.
        fin.getline(&line[0], 200);
        tag.clear();
        tag = line;

        // Get species name from line.
        c = line[0];
        for(i=1; (i<18)&&(!isWhiteSpace(c)); i++) {
            c = line[i];
        }
        spname = tag.substr(0,i-1);

        if (spname.compare("END") != 0) {
            // Find species in the mechanism.
            sp = mech.GetSpecies(spname);

            if (sp != NULL) {
                // Species was found in the mechanism, so we need to
                // read the thermo data.
                
                // Get elemental composition from this (1st) line.
                els[0] = tag.substr(24, 2);
                els[1] = tag.substr(29, 2);
                els[2] = tag.substr(34, 2);
                els[3] = tag.substr(39, 2);
                nels[0] = (int)atof(tag.substr(26,3).c_str());
                nels[1] = (int)atof(tag.substr(31,3).c_str());
                nels[2] = (int)atof(tag.substr(36,3).c_str());
                nels[3] = (int)atof(tag.substr(41,3).c_str());

                // Get the thermo temperature ranges from this (1st) line.
                lowT = atof(tag.substr(45,10).c_str());
                highT = atof(tag.substr(55,10).c_str());
                commT = atof(tag.substr(65,10).c_str());
                if (commT == 0) commT = trange[1];

                // Get 2nd line and read thermo coeffs.
                fin.getline(&line[0], 200);
                tag = line;
                up.Count = 7;
                up.Params[0] = atof(tag.substr(0,15).c_str());
                up.Params[1] = atof(tag.substr(15,15).c_str());
                up.Params[2] = atof(tag.substr(30,15).c_str());
                up.Params[3] = atof(tag.substr(45,15).c_str());
                up.Params[4] = atof(tag.substr(60,15).c_str());

                // Get 3rd line and read thermo coeffs.
                fin.getline(&line[0], 200);
                tag = line;
                up.Params[5] = atof(tag.substr(0,15).c_str());
                up.Params[6] = atof(tag.substr(15,15).c_str());
                lp.Count = 7;
                lp.Params[0] = atof(tag.substr(30,15).c_str());
                lp.Params[1] = atof(tag.substr(45,15).c_str());
                lp.Params[2] = atof(tag.substr(60,15).c_str());

                // Get 4th line and read thermo coeffs.
                fin.getline(&line[0], 200);
                tag = line;
                lp.Params[3] = atof(tag.substr(0,15).c_str());
                lp.Params[4] = atof(tag.substr(15,15).c_str());
                lp.Params[5] = atof(tag.substr(30,15).c_str());
                lp.Params[6] = atof(tag.substr(45,15).c_str());

                // Now save the data to the species:
                
                // Add the elemental composition to the species.
                for (i=0; i<4; i++) {
                    // Need to check the element name is valid.
                    try {
                        if (!isWhiteSpace(*els[i].substr(0,1).c_str())) {
                            if (!isWhiteSpace(*els[i].substr(1,1).c_str())) {
                                sp->AddElement(els[i], nels[i]);
                            } else {
                                sp->AddElement(els[i].substr(0,1), nels[i]);
                            }                    
                        }
                    } catch (std::invalid_argument &ia) {
                        throw range_error(string("Invalid element defined for species ").append(spname));
                        cout << spname;
                    }
                }
                
                // Add the thermo parameters.
                sp->SetThermoStartTemperature(lowT);
                sp->AddThermoParams(commT, lp);
                sp->AddThermoParams(highT, up);
            } else {
                // Species was not found in the mechanism, so we need to
                // skip lines to the next species.
                fin.getline(&line[0], 200);
                fin.getline(&line[0], 200);
                fin.getline(&line[0], 200);
            }
        } else {
            stat.Status = End;
        }
    }
}

// Reads CHEMKIN formatted thermo data for all species in the given mechanism from
// the file specified by the file name.
void Mechanism_IO::parseCK_Thermo(std::string &filename, Sprog::Mechanism &mech, 
                                  Sprog::IO::Mechanism_IO::CK_STATUS &stat)
{
    // Open file for reading.
    ifstream fin; //(filename.c_str(), ios::in);
    fin.open(filename.c_str(), ios::in);

    if (fin.good()) {
        try {
            char c[80]; fin.getline(&c[0],80); // Discard first line.
            parseCK_Thermo(fin, mech, stat);
        } catch (exception &e) {
            fin.close();
            throw e;
        }
        
        fin.close();
        return;
    } else {
        // Failed to open file.
        throw invalid_argument(string("Could not open THERMO file: ").append(filename));
    }
}

// Reads all the chemical reactions from the CHEMKIN formatted file stream.
void Mechanism_IO::parseCK_Reactions(std::ifstream &fin, Sprog::Mechanism &mech, 
                                     Sprog::IO::Mechanism_IO::CK_STATUS &stat)
{
    char c;
    string tag, rxndef;
    bool fcont = false, fdup = false;

    Kinetics::Reaction * last_rxn;

    stat.Status = ParseRxn;
    tag.clear();
    rxndef.clear();

    while((stat.Status!=End) && (stat.Status!=Fail) && (fin.good())) {
        // Get the next character from the input file.
        fin.get(c);

        switch(stat.Status) {
            case BeginParseRxn:
                // Clear white space at the beginning of the line.
                if (!isWhiteSpace(c)) {
                    rxndef = c;
                    stat.Status = ParseRxn;
                }
            case ParseRxn:
                // We are now reading a reaction definition.  Currently we are building
                // a string rxndef to hold to the reaction info.  This string will be
                // passed to the function parseCK_Reaction for translation.

                if (c=='!') {
                    // This is the beginning of a comment.  All the preceeding line must
                    // be the reaction definition unless a continuation & was present.
                    if (!fcont) {
                        // Before parsing the string we should check that it doesn't
                        // hold auxilliary info for the previous reaction.
                        if (!parseCK_RxnAux(rxndef, last_rxn, stat)) {
                            // No aux info.  First add the last reaction to the
                            // mechanism.
                            mech.AddReaction(last_rxn);

                            // No aux info, so parse the reaction definition.
                            last_rxn = parseCK_Reaction(rxndef, mech, stat);
                        }
                    }
                    stat.Status = ParseRxnComment;
                } else if ((c=='\n') || (c=='\r')) {
                    // This is the end of the line.  If there was a continuation character '&'
                    // then we continue reading the reaction string, otherwise it gets passed
                    // to parseCK_Reaction for translation.
                    if (!fcont) {
                        // Before parsing the string we should check that it doesn't
                        // hold auxilliary info for the previous reaction.
                        if (!parseCK_RxnAux(rxndef, last_rxn, stat)) {
                            // No aux info.  First add the last reaction to the
                            // mechanism.
                            mech.AddReaction(last_rxn);

                            // No aux info, so parse the reaction definition.
                            last_rxn = parseCK_Reaction(rxndef, mech, stat);
                        }
                    }
                }
                break;
            case ParseRxnComment:
                // We are reading a comment in the reaction list.  A comment is only
                // terminated by a line ending.
                if ((c == '\n') || (c == '\r')) {
                    stat.Status = ParseRxn;
                }
        }
    }
}

// Parses a string and builds a Reaction object from the data therein.
Sprog::Kinetics::Reaction *const Mechanism_IO::parseCK_Reaction(const std::string &rxndef, 
                                                                Sprog::Mechanism &mech, 
                                                                Sprog::IO::Mechanism_IO::CK_STATUS &stat)
{
    // The reaction string contains four pieces of information:  the reaction formula, A, n & E.
    Kinetics::Reaction * rxn = NULL;
    Kinetics::FallOffReaction * forxn = NULL;
    Kinetics::ThirdBodyReaction * tbrxn = NULL;

    int i, k, ilast_reac, ifirst_prod, ilast_prod, imu; // String & vector indices.
    bool frev; // Reversible reaction?
    bool fsearch, fplus; // Flags for searching string.
    bool fisfo = false, fistb = false; // Fall-off and third-body reaction flags.
    string strtb; // Third body, if not M.

    // Vectors to hold stoichiometry information.
    vector<Sprog::Stoich> rmui, pmui;
    vector<Sprog::Stoichf> rmuf, pmuf;
  
    vector<string> arrstr; // String vector to hold Arrhenius coefficients.
    Kinetics::ARRHENIUS arr; // Arrhenius coefficients.

    // Locate the delimiter between the reactions and products.
    i = rxndef.find_first_of("<=>");
    if (i != rxndef.npos) {
        // Found reversible delimiter.
        ilast_reac = i - 1;
        ifirst_prod = i + 3;
        frev = true;
    } else {
        i = rxndef.find_first_of("=>");
        if (i != rxndef.npos) {
            // Found irreversible delimiter.
            ilast_reac = i - 1;
            ifirst_prod = i + 2;
            frev = false;
        } else {
            i = rxndef.find_first_of("=");
            if (i != rxndef.npos) {
                // Found reversible delimiter.
                ilast_reac = i - 1;
                ifirst_prod = i + 1;
                frev = true;
            }
        }
    }

    // Find the end of the products by locating a text item aftera space that doesn't
    // follow a "+".
    k = ifirst_prod + 1; fsearch = true; fplus = true;
    while (k < rxndef.length()) {
        if (rxndef.substr(k, 1) == "+") {
            fplus = true; fsearch = true;
        } else if (rxndef.substr(k, 1) == " ") {
            fsearch = true;
        } else if (rxndef.substr(k, 1) == "(") {
            if (fsearch && !fplus) {
                if (rxndef.substr(k+1,1) == "+") {
                    // Just a delimiter.
                    fsearch = false;
                } else {
                    // That's it!
                    ilast_prod = k - 1;
                    break;
                }
            }
        } else {
            if (fsearch && !fplus) {
                // That's it, found the end.
                ilast_prod = k - 1;
                break;
            } else if (fsearch && fplus) {
                // Just another species name.
                fsearch = false; fplus = false;
            }
        }
        k++;
    }

    // Get the reactant stoichiometry and set it in the reaction.
    parseCK_RxnSpStoich(rxndef.substr(0,ilast_reac), mech, rmui, rmuf, fistb, fisfo, strtb);

    // Get the product stoichiometry and set it in the reaction.
    parseCK_RxnSpStoich(rxndef.substr(ifirst_prod-1,ilast_prod-ifirst_prod), mech, pmui, pmuf, fistb, fisfo, strtb);

    // Now we have got the reactants and the products, we must now parse the Arrhenius coefficients.
    // This is quite simple as they must be space delimited.
    split(rxndef.substr(ilast_reac), arrstr, " ");
    arr.A = atof(arrstr[0].c_str());
    arr.n = atof(arrstr[1].c_str());
    arr.E = atof(arrstr[2].c_str());
    
    // Now all the information about the reaction has been acquired we need to initialise the reaction
    // using the correct class and add it to the mechanism.
    if (fisfo) {
        // This is a fall-off reaction.
        forxn = new Kinetics::FallOffReaction();

        // Set third-body if not M.  This needs to be done after the reaction is added
        // to the mechanism so that the species can be correctly identified.
        if (strtb.compare("M") != 0) {
            forxn->SetFallOffThirdBody(mech.FindSpecies(strtb));
        }

        rxn = forxn;
    } else if (fistb) {
        // This is a third-body reaction, but without pressure dependence.
        tbrxn = new Kinetics::ThirdBodyReaction();
        rxn = tbrxn;
    } else {
        // This is a basic reaction.
        rxn = new Kinetics::Reaction();
    }

    // Set basic reaction properties:

    // Reactants.
    for (imu=0; imu<rmui.size(); imu++) rxn->AddReactant(rmui[imu]);
    for (imu=0; imu<rmuf.size(); imu++) rxn->AddReactant(rmuf[imu]);
    rmui.clear(); rmuf.clear();
    
    // Products.
    for (imu=0; imu<pmui.size(); imu++) rxn->AddReactant(pmui[imu]);
    for (imu=0; imu<pmuf.size(); imu++) rxn->AddReactant(pmuf[imu]);
    pmui.clear(); pmuf.clear();

    // Arrhenius coefficients.
    rxn->SetArrhenius(arr);

   return NULL;
}

// Parses a string of species names in a reaction, delimited by '+'s.
void Sprog::IO::Mechanism_IO::parseCK_RxnSpStoich(const std::string &sp,
                                                  const Sprog::Mechanism &mech,
                                                  std::vector<Sprog::Stoich> &mui, 
                                                  std::vector<Sprog::Stoichf> &muf, 
                                                  bool &isthirdbody, bool &isfalloff, 
                                                  std::string &thirdbody)
{
    vector<string> species;
    vector<string>::iterator k;
    int i = 0;
    bool search = true, plus = true;

    // Split the string into individual species names.
    while (i < sp.length()) {
        switch (sp.at(i)) {
            case ' ':
                search = true;
                break;
            case '+':
                if (search && !plus) {
                    // Delimiter between species.
                    plus = true;
                } else if (search && plus) {
                    // Oh dear, cannot start a species name with '+'.
                    throw invalid_argument("Invalid species name beginning with '+'");
                } else if ((sp.at(i-1) == '(') || (sp.at(i+1) == '+')) {
                    // Part of the species name.
                    (*k).append(sp.substr(i,1));
                } else {
                    // Delimiter between species.
                    search = true; plus = true;
                }
                break;
            case '(':
                if (!plus) {
                    // Better hope this is followed immediately by a plus sign, indicating
                    // a pressure dependent reaction.
                    if (sp.at(i+1) == '+') {
                        // Good.  Allow "(" to act as species delimiter, but
                        // save it in the specie's name as the interpreter
                        // will need to know about the pressure dependence.
                        search = false;
                        (*k).append(sp.substr(i,1));
                    } else {
                        // Noooo!  So many things wrong here.  Can't start a
                        // species name with a "(", and even if we could we
                        // have no "+" delimiter to denote a new species, so 
                        // this is technically a space in a species.
                        throw invalid_argument("Invalid species name.");
                    }
                } else if (search && plus) {
                    // Can't start a species name with a bracket.
                    throw invalid_argument("Invalid species name beginning with a bracket.");
                } else {
                    // This is just part of the species name.
                    (*k).append(sp.substr(i,1));
                }
                break;
            default:
                if (search && !plus) {
                    // There is a space in a species name.  This is not allowed.
                    throw invalid_argument("Space in a species name.");
                } else if (search && plus) {
                    // Found the beginning of a new species name.
                    species.push_back(sp.substr(i,1));
                    k = species.end()-1;
                    search = false; plus = false;
                } else {
                    // Continuing to add characters to a species name.
                    (*k).append(sp.substr(i,1));
                }
                break;
        }
    }

    // Now we must separate the stoichiometry from the species' names.  Also we can check here
    // for third-bodies or fall-off reactions.
    string sym, mu;
    int intmu;
    real floatmu;

    for (k=species.begin(); k!=species.end(); k++) {
        // Separate stoichiometry from specie's name.
        i = (*k).find_first_not_of("0123456789.");
        if (i != string.npos) {
            if (i > 0) 
                mu = (*k).substr(0,i); 
            else 
                mu = "";
            sym = (*k).substr(i);
        } else {
            mu = ""; sym = "";
        }

        // Check for third body.
        if ((sym.substr(0,2) == "(+") && (sym.substr(sym.length()-1) == ")")) {
            // This is a third-body indicating pressure dependence.  Set the TB and FO flags.  Also
            // remember the third body symbol if not M.
            isthirdbody = true; isfalloff = true;
            thirdbody = sym.substr(2, sym.length()-1);
            if (thirdbody.compare("M")==0) thirdbody = "";
        } else if (sym != "M") {
            // This is not a third body at all.

            // Get the index of the species from the mechanism.
            i = mech.FindSpecies(sym);
            
            if (i > 0) {
                // Before we save the species stoichiometry we need to check if it is integer or real.
                intmu = atoi(mu.c_str()); floatmu = atof(mu.c_str());
                if (floatmu == (real)intmu) {
                    // This is integer stoichiometry.
                    mui.push_back(Stoich(i, intmu));
                } else {
                    // This is real stoichiometry.
                    muf.push_back(Stoichf(i, floatmu));
                }
            } else {
                // Species not found in mechanism.
                throw invalid_argument("Unknown species in reaction.");
            }
        } else {
            //This is a third-body reaction, but not a fall-off reaction.
            isthirdbody = true;
        }
    }
}

// Parses auxilliary reaction information.
bool Sprog::IO::Mechanism_IO::parseCK_RxnAux(const std::string &rxndef, 
                                             Sprog::Kinetics::Reaction *last_rxn, 
                                             const Sprog::Mechanism &mech, 
                                             Sprog::IO::Mechanism_IO::CK_STATUS &stat)
{
      
}