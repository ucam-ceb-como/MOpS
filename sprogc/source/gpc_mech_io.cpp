/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the MechanismParser class declared in the
    gpc_mech_io.h header file.

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

#include "comostrings.h"
#include "gpc_mech_io.h"
#include "gpc_mech.h"

#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>

using namespace Sprog::IO;
using namespace std;
using namespace Sprog::Kinetics;
using namespace Strings;

// Reads a mechanism from a CHEMKIN input file.
void MechanismParser::ReadChemkin(const std::string &filename, 
                                  Sprog::Mechanism &mech, 
                                  const std::string &thermofile)
{
    // Open file for reading.
    ifstream fin; //(filename.c_str(), ios::in);
    fin.open(filename.c_str(), ios::in);

    if (fin.good()) {
        CK_STATUS status;
        status.Status = FindKey;
        status.ReadElements = false;
        status.ReadSpecies = false;
        status.ReadReactions = false;
        status.ReadThermo = false;
        status.ThermoFile = thermofile;
        status.Scale = ARRHENIUS(1.0, 1.0, 1.0);

        // Clear current mechanism.
        mech.Clear();

        // CHEMKIN files are read in CGS units.
        mech.SetUnits(CGS);

        try {
            parseCK(fin, mech, status);
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
void MechanismParser::parseCK(std::ifstream &fin, 
                           Sprog::Mechanism &mech, 
                           Sprog::IO::MechanismParser::CK_STATUS &status)
{
    char c;
    string line;

    // Convert chemkin file stream to a string and store in chemkinstr
    string chemkinstr;
    do {
        // get a charater to see if file is empty, this is a must.
        // get function consider \r\n sequence as \n
        fin.get(c);
        if (fin.eof()) break;
        // skiping comment until new line or return character is found
        if  (c=='!') {
            do {
                fin.get(c);
            } while (!fin.eof() && (c!='\r') && (c!='\n'));
            if (fin.eof()) break;
        }
        if (c=='\r') {
            // replace return character by new line character, if some system doesn't read \r\n as \n
            c = '\n';
        } else if (c=='\t') {
            // replace tab character by a space bar character
            c = ' ';
        }
        chemkinstr.append(1,c);
        //cout << c ;
    } while (!fin.eof());
    // Convert chemkin file string to an istringstream for passing as argument in the other parsers
    istringstream strin(ios_base::in);
    strin.str(chemkinstr);

    // Locate in the string stream the starting point of the elements, species
    // reactions and thermo data. This is hard coded. It is not the best but it does the job
    // for this specific case.
    int el_index = 0, sp_index = 0, rt_index = 0, tm_index = 0;
    string chemkinCapstr = convertToCaps(chemkinstr);
    el_index = chemkinCapstr.find("ELEMENT", 0);
    sp_index = chemkinCapstr.find("SPECIES", 0);
    rt_index = chemkinCapstr.find("REACION", 0);
    tm_index = chemkinCapstr.find("THERMO", 0);
    if (el_index == string::npos) {
        el_index = chemkinCapstr.find("ELEM", 0);
        if (el_index == string::npos) {
            throw invalid_argument("ELEM or ELEMENT keyword not found in chemkin input file (Sprogc, MechanismParser::parseCK).");
        } else {
            el_index += 4;
        }
    } else {
        el_index += 7;
    }
    if (sp_index == string::npos) {
        sp_index = chemkinCapstr.find("SPEC", 0);
        if (sp_index == string::npos) {
            throw invalid_argument("SPEC or SPECIES keyword not found in chemkin input file (Sprogc, MechanismParser::parseCK).");
        } else {
            sp_index += 4;
        }
    } else {
        sp_index += 7;
    }
    if (rt_index == string::npos) {
        rt_index = chemkinCapstr.find("REAC", 0);
        if (rt_index == string::npos) {
            throw invalid_argument("REAC or REACTION keyword not found in chemkin input file (Sprogc, MechanismParser::parseCK).");
        } else {
            rt_index += 4;
        }
    } else {
        rt_index += 7;
    }
    if (tm_index == string::npos) {
        tm_index = chemkinCapstr.find("THER", 0);
        if (tm_index == string::npos) {
            // Thermo keywarod is not needed.
            //throw invalid_argument("THER or THERMO keyword not found in chemkin input file (Sprogc, MechanismParser::parseCK).");
        } else {
            tm_index += 4;
        }
    } else {
        tm_index += 6;
    }

    // Read the elements.
    strin.clear();
    strin.seekg(el_index);
    parseCK_Elements(strin, mech, status);
    status.ReadElements = true;

    // Read the species.
    strin.clear();
    strin.seekg(sp_index);
    parseCK_Species(strin, mech, status);
    status.ReadSpecies = true;

    // Read the thermo data.
    if (status.ThermoFile == "") {
        // Read thermo data from this file.
        strin.clear();
        strin.seekg(tm_index);
        parseCK_Thermo(strin, mech, status);
    } else {
        // Read thermo data from thermo file.
        parseCK_Thermo(status.ThermoFile, mech, status);
    }
    status.ReadThermo = true;

    // Read the reactions.
    strin.clear();
    strin.seekg(rt_index);
    parseCK_Reactions(strin, mech, status);
    mech.BuildStoichXRef();
    mech.SetUnits(SI);
}

void MechanismParser::parseCK_Elements(std::istream&fin, Sprog::Mechanism &mech, 
                                    Sprog::IO::MechanismParser::CK_STATUS &status)
{
    char c = '\0';
    string tag;
    real val = 0.0;

    Element * last_el = NULL;

    status.Status = BeginParseEl;

    while((status.Status!=End) && (status.Status!=Fail) && (fin.good())) {
        // Get the next character from the input file.
        fin.get(c);

        switch(status.Status) {
            case BeginParseEl:
                // Here we are looking for the start of an element name.
                if (isLetter(c)) {
                    // This is the start of an element name.
                    tag = c;
                    status.Status = ParseEl;

                    // Clear the last read element.
                    last_el = NULL;
                } else if (c == '/') {
                    // This character means that the last element's Mol. Wt. is defined in
                    // the input file.  If the last element is valid, then we must read
                    // the mol. wt. from the file.
                    if (last_el != NULL) {
                        status.Status = ParseElWt;
                        tag = "";
                    } else {
                        // Oh dear.  There is an element mol. wt. definition without an
                        // element.  This is an error.
                        status.Status = Fail;
                        throw range_error("Invalid / character in element definition.");
                    }
                } else if (c == '!') {
                    // This is a comment.
                    status.Status = ParseElComment;
                } else if (!isWhiteSpace(c)) {
                    // This must be an invalid character.
                    status.Status = Fail;
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
                        status.Status = BeginParseEl;
                    } else if (convertToCaps(tag.substr(0,4)) == "SPEC") {
                        status.Status = End;
                    } else if (convertToCaps(tag.substr(0,4)) == "REAC") {
                        status.Status = End;
                    } else if (convertToCaps(tag.substr(0,4)) == "END") {
                        status.Status = End;
                    } else {
                        // This is an element name, so add a new
                        // element to the mechanism and set its name.
                        last_el = mech.AddElement();
                        last_el->SetName(tag);
                        // Attempt to find the element weight in the library of
                        // know elements.
                        last_el->SetMolWtFromLibrary();
                        // Start searching for the next element.
                        status.Status = BeginParseEl;
                    }
                    tag = "";
                } else if (c == '/') {
                    // This character means that the element Mol. Wt. is defined in
                    // the input file.  First task is to initialise a new element.
                    last_el = mech.AddElement();
                    last_el->SetName(tag);
                    // We must now read Mol. Wt from file.
                    status.Status = ParseElWt;
                    tag = "";
                }
                break;
            case ParseElWt:
                // We are reading an element mol. wt. from the file.
                if (c == '/') {
                    // This is the end of the mol. wt. definition.  Attempt to convert it
                    // to a number.
                    val = cdble(tag);
                    if (val > 0.0) {
                        // Set the mol. wt. of the last read element.
                        last_el->SetMolWt(val);
                    } else {
                        // The mol. wt. is invalid.
                        status.Status = Fail;
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
                    status.Status = BeginParseEl;
                }
                break;
            default:
                throw runtime_error("Illegal parsing flag encountered "
                                    "(Sprog, MechanismParser::parseCK_Elements).");
        }
    }
}

void MechanismParser::parseCK_Species(std::istream&fin, 
                                   Sprog::Mechanism &mech, 
                                   Sprog::IO::MechanismParser::CK_STATUS &status)
{
    char c;
    string tag;

    Species * last_sp;

    status.Status = BeginParseSp;

    while((status.Status!=End) && (status.Status!=Fail) && (fin.good())) {
        // Get the next character from the input file.
        fin.get(c);

        switch(status.Status) {
            case BeginParseSp:
                // We are searching for the beginning of a species name.
                if (isLetter(c)) {
                    // This is the start of a species name.
                    tag = c;
                    status.Status = ParseSp;

                    // Clear the last read species.
                    last_sp = NULL;
                } else if (c == '!') {
                    // This is a comment.
                    status.Status = ParseSpComment;
                } else if (!isWhiteSpace(c)) {
                    // This must be an invalid character.
                    status.Status = Fail;
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
                        status.Status = Fail;
                        throw range_error("Element definition after species definitions.");
                    } else if (convertToCaps(tag.substr(0,4)) == "SPEC") {
                        status.Status = BeginParseSp;
                    } else if (convertToCaps(tag.substr(0,4)) == "REAC") {
                        status.Status = End;
                    } else if (convertToCaps(tag.substr(0,4)) == "END") {
                        status.Status = End;
                    } else {
                        // This is a species name, so add a new
                        // species to the mechanism and set its name.
                        last_sp = mech.AddSpecies();
                        last_sp->SetName(tag);
                        // Begin searching for the next species name.
                        status.Status = BeginParseSp;
                    }
                }
                break;
            case ParseSpComment:
                // We are reading a comment in the species list.  A comment is only
                // terminated by a line ending.
                if ((c == '\n') || (c == '\r')) {
                    status.Status = BeginParseSp;
                }
                break;
            default:
                throw runtime_error("Illegal parsing flag encountered "
                                    "(Sprog, MechanismParser::parseCK_Species).");
        }
    }
}

// Reads CHEMKIN formatted thermo data for all species in the given mechanism from
// the supplied file stream.
void MechanismParser::parseCK_Thermo(std::istream&fin, Sprog::Mechanism &mech, 
                                  Sprog::IO::MechanismParser::CK_STATUS &status)
{
    int i;
    char c, line[200];
    string tag, spname;

    real trange[3], lowT, highT, commT;
    Thermo::THERMO_PARAMS up; // Coeffs for upper T interval.
    Thermo::THERMO_PARAMS lp; // Coeffs for lower T interval.
    vector<string> els(4);
    int nels[4];

    Species * sp;

    status.Status = ParseTherm;
    tag.resize(200,' ');

    // The first line holds the temperature ranges for 2 sets of coefficients.
    fin.getline(&line[0], 200);
    while (line[0]=='\0') fin.getline(&line[0], 200);
    tag.clear();
    tag = line;
    trange[0] = cdble(tag.substr(0,10));
    trange[1] = cdble(tag.substr(10,10));
    trange[2] = cdble(tag.substr(20,10));

    while((status.Status!=End) && (status.Status!=Fail) && (fin.good())) {
        // Get species name line.
        fin.getline(&line[0], 200);
        tag.clear();
        tag = line;

        if ((tag != "") && (line[0] != '!')) {
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
                    nels[0] = (int)cdble(tag.substr(26,3));
                    nels[1] = (int)cdble(tag.substr(31,3));
                    nels[2] = (int)cdble(tag.substr(36,3));
                    nels[3] = (int)cdble(tag.substr(41,3));

                    // Get the thermo temperature ranges from this (1st) line.
                    lowT = cdble(tag.substr(45,10));
                    highT = cdble(tag.substr(55,10));
                    commT = cdble(tag.substr(65,10));
                    if (commT == 0) commT = trange[1];

                    // Get 2nd line and read thermo coeffs.
                    fin.getline(&line[0], 200);
                    tag = line;
                    up.Count = 7;
                    up.Params[0] = cdble(tag.substr(0,15));
                    up.Params[1] = cdble(tag.substr(15,15));
                    up.Params[2] = cdble(tag.substr(30,15));
                    up.Params[3] = cdble(tag.substr(45,15));
                    up.Params[4] = cdble(tag.substr(60,15));

                    // Get 3rd line and read thermo coeffs.
                    fin.getline(&line[0], 200);
                    tag = line;
                    up.Params[5] = cdble(tag.substr(0,15));
                    up.Params[6] = cdble(tag.substr(15,15));
                    lp.Count = 7;
                    lp.Params[0] = cdble(tag.substr(30,15));
                    lp.Params[1] = cdble(tag.substr(45,15));
                    lp.Params[2] = cdble(tag.substr(60,15));

                    // Get 4th line and read thermo coeffs.
                    fin.getline(&line[0], 200);
                    tag = line;
                    lp.Params[3] = cdble(tag.substr(0,15));
                    lp.Params[4] = cdble(tag.substr(15,15));
                    lp.Params[5] = cdble(tag.substr(30,15));
                    lp.Params[6] = cdble(tag.substr(45,15));

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
                            throw ia;
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
                status.Status = End;
            }
        }
    }
}

// Reads CHEMKIN formatted thermo data for all species in the given mechanism from
// the file specified by the file name.
void MechanismParser::parseCK_Thermo(std::string &filename, Sprog::Mechanism &mech, 
                                     Sprog::IO::MechanismParser::CK_STATUS &status)
{
    // Open file for reading.
    ifstream fin; //(filename.c_str(), ios::in);
    fin.open(filename.c_str(), ios::in);

    if (fin.good()) {
        try {
            char c[81]; fin.getline(&c[0],80); // Discard first line.
            parseCK_Thermo(fin, mech, status);
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
void MechanismParser::parseCK_Reactions(std::istream&fin, Sprog::Mechanism &mech, 
                                     Sprog::IO::MechanismParser::CK_STATUS &status)
{
    char c;
    string tag, rxndef;
    bool fcont = false; //, fdup = false;

    Kinetics::Reaction * last_rxn = NULL;

    status.Status = ParseRxn;
    tag.clear();
    rxndef.clear();

    while((status.Status!=End) && (status.Status!=Fail) && (fin.good())) {
        // Get the next character from the input file.
        fin.get(c);

        switch(status.Status) {
            case BeginParseRxn:
                // Clear white space at the beginning of the line.
                if (c=='!') {
                    status.Status = BeginParseRxnComment;
                    rxndef = "";
                } else if (!isWhiteSpace(c)) {
                    rxndef = c;
                    status.Status = ParseRxn;
                }
                break;
            case ParseRxn:
                // We are now reading a reaction definition.  Currently we are building
                // a string rxndef to hold to the reaction info.  This string will be
                // passed to the function parseCK_Reaction for translation.

                if (c=='!') {
                    // Break if this is the end of the reaction definitions.
                    if (rxndef.substr(0,3).compare("END")==0) {                    
                        // Before parsing the string we should check that it doesn't
                        // hold auxilliary info for the previous reaction.
                        if (!parseCK_RxnAux(rxndef, last_rxn, mech, status.Scale, status)) {
                            // No aux info.  First add the last reaction to the
                            // mechanism.
                            if (last_rxn != NULL) mech.AddReaction(last_rxn);
                        }
                        status.Status = End;
                        break;
                    }

                    // This is the beginning of a comment.  All the preceeding line must
                    // be the reaction definition unless a continuation & was present.
                    if (!fcont) {
                        // Before parsing the string we should check that it doesn't
                        // hold auxilliary info for the previous reaction.
                        if (!parseCK_RxnAux(rxndef, last_rxn, mech, status.Scale, status)) {
                            // No aux info.  First add the last reaction to the
                            // mechanism.
                            if (last_rxn != NULL) mech.AddReaction(last_rxn);

                            // No aux info, so parse the reaction definition.
                            if (rxndef.length() > 0) last_rxn = parseCK_Reaction(rxndef, mech, status);
                        }

                        rxndef = "";
                    }
                    status.Status = BeginParseRxnComment;
                } else if ((c=='\n') || (c=='\r')) {
                    // Break if this is the end of the reaction definitions.
                    if (rxndef.substr(0,3).compare("END")==0) {                    
                        // Before parsing the string we should check that it doesn't
                        // hold auxilliary info for the previous reaction.
                        if (!parseCK_RxnAux(rxndef, last_rxn, mech, status.Scale, status)) {
                            // No aux info.  First add the last reaction to the
                            // mechanism.
                            if (last_rxn != NULL) mech.AddReaction(last_rxn);
                        }
                        status.Status = End;
                        break;
                    }

                    // This is the end of the line.  If there was a continuation character '&'
                    // then we continue reading the reaction string, otherwise it gets passed
                    // to parseCK_Reaction for translation.
                    if (!fcont) {
                        // Before parsing the string we should check that it doesn't
                        // hold auxilliary info for the previous reaction.
                        if (!parseCK_RxnAux(rxndef, last_rxn, mech, status.Scale, status)) {
                            // No aux info.  First add the last reaction to the
                            // mechanism.
                            if (last_rxn != NULL) {
                                mech.AddReaction(last_rxn);
                                delete last_rxn;
                                last_rxn = NULL;
                            }

                            // No aux info, so parse the reaction definition.
                            if (rxndef.length() > 0) last_rxn = parseCK_Reaction(rxndef, mech, status);
                        }
                        status.Status = BeginParseRxn;
                    }
                } else {
                    rxndef.append(&c, 1);
                }
                break;
            case ParseRxnComment:
                // We are reading a comment in the reaction list.  A comment is only
                // terminated by a line ending.
                if ((c == '\n') || (c == '\r')) {
                    status.Status = ParseRxn;
                }
                break;
            case BeginParseRxnComment:
                // We are reading a comment in the reaction list.  A comment is only
                // terminated by a line ending.
                if ((c == '\n') || (c == '\r')) {
                    status.Status = BeginParseRxn;
                }
                break;
            default:
                throw runtime_error("Illegal parsing flag encountered "
                                    "(Sprog, MechanismParser::parseCK_Reactions).");
        }
    }
}

// Parses a string and builds a Reaction object from the data therein.
Sprog::Kinetics::Reaction *const MechanismParser::parseCK_Reaction(const std::string &rxndef, 
                                                                Sprog::Mechanism &mech, 
                                                                Sprog::IO::MechanismParser::CK_STATUS &status)
{
    // The reaction string contains four pieces of information:  the reaction formula, A, n & E.
    Kinetics::Reaction * rxn = NULL;

    string::size_type i=0, k=0, ilast_reac=0; // String & vector indices.
    string::size_type ifirst_prod=0, ilast_prod=0;
    unsigned int imu=0;
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
    i = rxndef.find_first_of("<");
    if (i != rxndef.npos) {
        // Found reversible delimiter.
        ilast_reac = i - 1;
        ifirst_prod = i + 3;
        frev = true;
    } else {
        i = rxndef.find_first_of(">");
        if (i != rxndef.npos) {
            // Found irreversible delimiter.
            ilast_reac = i - 2;
            ifirst_prod = i + 1;
            frev = false;
        } else {
            i = rxndef.find_first_of("=");
            if (i != rxndef.npos) {
                // Found reversible delimiter.
                ilast_reac = i - 1;
                ifirst_prod = i + 1;
                frev = true;
            } else {
                // No delimiter found.
                throw invalid_argument("No reactant/product delimiter found in reaction definition.");
            }
        }
    }

    // Find the end of the products by locating a text item after a space that doesn't
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
                ilast_prod = rxndef.find_last_not_of(" ", k-1)+1;
                break;
            } else if (fsearch && fplus) {
                // Just another species name.
                fsearch = false; fplus = false;
            }
        }
        k++;
    }

    // Get the reactant stoichiometry and set it in the reaction.
    parseCK_RxnSpStoich(rxndef.substr(0,ilast_reac+1), mech, rmui, rmuf, fistb, fisfo, strtb);

    // Get the product stoichiometry and set it in the reaction.
    parseCK_RxnSpStoich(rxndef.substr(ifirst_prod,ilast_prod-ifirst_prod+1), 
                        mech, pmui, pmuf, fistb, fisfo, strtb);

    // Now we have got the reactants and the products, we must now parse the Arrhenius coefficients.
    // This is quite simple as they must be space delimited.
    split(rxndef.substr(ilast_prod), arrstr, " ");
    arr.A = cdble(arrstr[0]) * status.Scale.A;
    arr.n = cdble(arrstr[1]);
    arr.E = cdble(arrstr[2]) * status.Scale.E;
    
    // Now all the information about the reaction has been acquired so we need to 
    // initialise the reaction object.
    rxn = new Kinetics::Reaction();
    
    if (fisfo) {
        // Set fall-off third body if not M.
        if (strtb.compare("M") != 0) {
            rxn->SetFallOffThirdBody(mech.FindSpecies(strtb));
        }
    } else if (fistb) {
        // This is a third-body reaction, but without pressure dependence.
        rxn->SetUseThirdBody(true);
    }

    // Set basic reaction properties:

    // Reactants.
    for (imu=0; imu<rmui.size(); imu++) rxn->AddReactant(rmui[imu]);
    for (imu=0; imu<rmuf.size(); imu++) rxn->AddReactant(rmuf[imu]);
    rmui.clear(); rmuf.clear();
    
    // Products.
    for (imu=0; imu<pmui.size(); imu++) rxn->AddProduct(pmui[imu]);
    for (imu=0; imu<pmuf.size(); imu++) rxn->AddProduct(pmuf[imu]);
    pmui.clear(); pmuf.clear();

    // Arrhenius coefficients.
    rxn->SetArrhenius(arr);

    // Is reaction reversible?
    rxn->SetReversible(frev);

   return rxn;
}

// Parses a string of species names in a reaction, delimited by '+'s.
void Sprog::IO::MechanismParser::parseCK_RxnSpStoich(const std::string &sp,
                                                  const Sprog::Mechanism &mech,
                                                  std::vector<Sprog::Stoich> &mui, 
                                                  std::vector<Sprog::Stoichf> &muf, 
                                                  bool &isthirdbody, bool &isfalloff, 
                                                  std::string &thirdbody)
{
    vector<string> species;
    vector<string>::iterator k;
    string str;
    string::size_type i = 0, j = 0;
    // bool search = true, plus = true, fosym = false;

    // Locate the first species delimiter.
    i = sp.find_first_not_of(" ");
    if (i != sp.npos) {
        j = sp.find_first_of("+",i+1);
    } else {
        j = sp.npos;
    }

    // Loop over all delimiters until we reach the end of the string.
    while (j != sp.npos) {
        // There is a special case for fall-off reactions of the form (+M).
        // We must check for that.
        if (sp.substr(j-1,1) == "(") {
            // This is a fall-off reaction.
            isfalloff = true;
            isthirdbody = true;

            // Must remove ( bracket before saving previous species symbol, hence j-i-1.
            str = sp.substr(i, j-i-1);
            species.push_back(str.substr(0, str.find_last_not_of(" ")+1));

            // Save fall-off species name and reset iterators.
            i = j + 1; j = sp.find_first_of(")");
            thirdbody.assign(sp.begin()+i, sp.begin()+j);
            if (thirdbody.compare("M")==0) thirdbody = "";

        } else {
            // There is no bracket to worry about as this is not a fall-off symbol.
            str = sp.substr(i, j-i);
            species.push_back(str.substr(0, str.find_last_not_of(" ")+1));
        }

        // Find next delimiter.
        i = sp.find_first_not_of(" ", j+1);
        if (i != sp.npos) {
            j = sp.find_first_of("+",i+1);
        } else {
            j = sp.npos;
        }
    }

    // Save the last species.
    if ((i != sp.npos) && (i < sp.length())) {
        str = sp.substr(i);
        species.push_back(str.substr(0, str.find_last_not_of(" ")+1));
    }

    // Now we must separate the stoichiometry from the species' names.  Also we can check here
    // for third-bodies or fall-off reactions.
    string sym, mu;
    int intmu;
    real floatmu;

    for (k=species.begin(); k!=species.end(); k++) {
        // Separate stoichiometry from specie's name.
        i = (*k).find_first_not_of("0123456789.");
        if (i != sp.npos) {
            if (i > 0) 
                mu = (*k).substr(0,i); 
            else 
                mu = "1";
            sym = (*k).substr(i);
        } else {
            mu = "1"; sym = "";
        }

        if (sym.compare("M")==0) {
            // This is a third-body.
            isthirdbody = true;
        } else if (sym != "M") {
            // This is not a third body at all.

            // Get the index of the species from the mechanism.
            i = mech.FindSpecies(sym);
            
            if (i >= 0) {
                // Before we save the species stoichiometry we need to check if it is integer or real.
                intmu = atoi(mu.c_str()); floatmu = cdble(mu);
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
        }
    }
}

// Parses auxilliary reaction information.
bool Sprog::IO::MechanismParser::parseCK_RxnAux(const std::string &rxndef, 
                                             Sprog::Kinetics::Reaction *last_rxn, 
                                             const Sprog::Mechanism &mech,
                                             Sprog::Kinetics::ARRHENIUS scale,
                                             Sprog::IO::MechanismParser::CK_STATUS &status)
{
    string str, key, vals;
    vector<string> params;
    string::size_type i, j, k;
    int sp;
    real val, foparams[Kinetics::FALLOFF_PARAMS::MAX_FALLOFF_PARAMS];
    bool foundaux = false;

    // Find first non-blank character in aux information.
    i = rxndef.find_first_not_of(' ');
    
    // Find delimiters /.
    j = rxndef.find_first_of('/', i+1);
    k = rxndef.find_first_of('/', j+1);

    while (i != rxndef.npos) {
        // Now i-j is the keyword (or species name), j-k are the parameters.

        // Get the keyword.
        key = rxndef.substr(i, j-i);
        key = key.substr(0, key.find_first_of(' '));

        if (j != rxndef.npos) {
            // Check for unterminated set of parameters.
            if (k == rxndef.npos) {
                // An unterminated set of parameters.  This is an error.
                status.Status = Fail;
                throw invalid_argument("Unterminated set of parameters in auxilliary data.");
            }

            // Get the parameters.
            vals = rxndef.substr(j+1, k-j-1);
            split(vals, params, " ");
            foundaux = true;
        } else {
            // There are no parameters with this keyword.
            vals = "";
            params.clear();

            if (key.compare("DUPLICATE")==0) foundaux = true;
            return foundaux;
        }

        // Decide next action based on the keyword.
        if (key.compare("LOW") == 0) {
            // Low-pressure limit for fall-off reaction (3 parameters).
            if (params.size() > 2) {
                last_rxn->SetLowPressureLimit(ARRHENIUS(cdble(params[0]) * scale.A,
                                                        cdble(params[1]) * scale.n,
                                                        cdble(params[2]) * scale.E));
            } else {
                // Insufficient parameters.
                throw invalid_argument("Insufficient parameters for LOW keyword in auxilliary data.");
            }
        } else if (key.compare("HIGH") == 0) {
            // High-pressure limit for chemically activated reaction (3 parameters).
            // NOT IN USE.
        } else if (key.compare("TROE") == 0) {
            // TROE form pressure-dependent reaction (3 or 4 parameters).
            if (params.size() > 2) {
                foparams[0] = cdble(params[0]);
                foparams[1] = cdble(params[1]);
                foparams[2] = cdble(params[2]);
                foparams[4] = 0.0;

                if (params.size() > 3) {
                    // 4-parameter TROE form.
                    foparams[3] = cdble(params[3]);
                    last_rxn->SetFallOffParams(Troe4, foparams);
                } else {
                    // 3-parameter TROE form.
                    foparams[3] = 0.0;
                    last_rxn->SetFallOffParams(Troe3, foparams);
                }
            } else {
                // Insufficient parameters.
                throw invalid_argument("Insufficient parameters for TROE keyword in auxilliary data.");
            }
        } else if (key.compare("SRI") == 0) {
            // SRI pressure-dependent reaction (3 or 5 parameters).
            if (params.size() > 4) {
                // 5-parameter form.
                foparams[0] = cdble(params[0]);
                foparams[1] = cdble(params[1]);
                foparams[2] = cdble(params[2]);
                foparams[3] = cdble(params[3]);
                foparams[4] = cdble(params[4]);
                last_rxn->SetFallOffParams(SRI, foparams);
            } else if (params.size() > 2) {
                //3-parameter form.
                foparams[0] = cdble(params[0]);
                foparams[1] = cdble(params[1]);
                foparams[2] = cdble(params[2]);
                foparams[3] = 1.0;
                foparams[4] = 0.0;
                last_rxn->SetFallOffParams(SRI, foparams);
            } else {
                // Insufficient parameters.
                throw invalid_argument("Insufficient parameters for SRI keyword in auxilliary data.");
            }
        } else if (key.compare("USER") == 0) {
            // User-defined pressure fall-off reaction (5 parameters).
            // NOT IN USE.
        } else if (key.compare("REV") == 0) {
            // Reverse rate Arrhenius parameters (3 required).
            if (params.size() > 2) {
                last_rxn->SetRevArrhenius(ARRHENIUS(cdble(params[0]) * scale.A,
                                                    cdble(params[1]) * scale.n,
                                                    cdble(params[2]) * scale.E));
            } else {
                // Insufficient parameters.
                throw invalid_argument("Insufficient parameters for REV keyword in auxilliary data.");
            }
        } else if (key.compare("LT") == 0) {
            // Landau-Teller reaction parameters (2 required).
            if (params.size() > 1) {
                last_rxn->SetLTCoeffs(LTCOEFFS(cdble(params[0]), cdble(params[1])));
            } else {
                // Insufficient parameters.
                throw invalid_argument("Insufficient parameters for LT keyword in auxilliary data.");
            }
        } else if (key.compare("RLT") == 0) {
            // Reverse Landau-Teller parameters (2 required).
             if (params.size() > 1) {
                last_rxn->SetRevLTCoeffs(LTCOEFFS(cdble(params[0]), cdble(params[1])));
            } else {
                // Insufficient parameters.
                throw invalid_argument("Insufficient parameters for RLT keyword in auxilliary data.");
            }
       } else if (key.compare("TDEP") == 0) {
            // Species temperature dependence (1 parameter).
            // NOT IN USE.
        } else if (key.compare("EXCI") == 0) {
            // Energy loss (1 parameter).
            // NOT IN USE.
        } else if (key.compare("JAN") == 0) {
            // Rate constant fit expression (9 parameters).
            // NOT IN USE.
        } else if (key.compare("FIT1") == 0) {
            // Rate constant fit expression (4 parameters).
            // NOT IN USE.
        } else if (key.compare("HV") == 0) {
            // Radiation wavelength specification (1 parameter).
            // NOT IN USE.
        } else if (key.compare("MOME") == 0) {
            // Momentum transfer collision frequency for electrons (no params).
            // NOT IN USE.
        } else if (key.compare("FORD") == 0) {
            // Forward reaction order of species.
            // NOT IN USE.
        } else if (key.compare("RORD") == 0) {
            // Reverse reaction order of species.
            // NOT IN USE.
        } else if (key.compare("UNITS") == 0) {
            // Specify different units for a particular reaction.
            // NOT IN USE.
        } else {
            // Not a keyword at all, but a species name (hopefully!) for third-body efficiencies.
            sp = mech.FindSpecies(key);
            if ((sp >= 0) && (params.size() > 0)) {
                val = cdble(params[0]);
                last_rxn->AddThirdBody(sp, val);
            } else {
                // Unrecognised species or missing efficiency.
                throw invalid_argument("Error parsing auxilliary data.  Unknown keyword.");
            }
        }

        // Find first non-blank character in aux information, after current position.
        i = rxndef.find_first_not_of(' ', k+1);
        
        // Find delimiters /.
        j = rxndef.find_first_of('/', i+1);
        k = rxndef.find_first_of('/', j+1);
    }

    return foundaux;
}

// Parses the units from the REACTION line of a CHEMKIN formatted file.
void Sprog::IO::MechanismParser::parseCK_Units(const std::string &rxndef, Sprog::Kinetics::ARRHENIUS &scale)
{
    // Split the string into parts.
    vector<string> parts;
    split(rxndef, parts, " ");

    // Default scaling.
    scale.A = 1.0;
    scale.n = 1.0;
    scale.E = 4.184e7; // Convert cal/mol -> J/mol.

    // The first part will be the REACTION keyword, so we ignore that.
    vector<string>::iterator i;
    for (i=parts.begin()+1; i!=parts.end(); i++) {
        if ((*i).compare("CAL/MOLE")==0) {
            // Default scaling.
        } else if ((*i).compare("KCAL/MOLE")==0) {
            scale.E = 4184.0e7;
        } else if ((*i).compare("JOULES/MOLE")==0) {
            scale.E = 1.0e7;
        } else if ((*i).compare("KJOULES/MOLE")==0) {
            scale.E = 1.0e10;
        } else if ((*i).compare("KJOU/MOLE")==0) {
            scale.E = 1.0e10;
        } else if ((*i).compare("KJOU/MOL")==0) {
            scale.E = 1.0e10;
        } else if ((*i).compare("KELVINS")==0) {
            scale.E = R;
        } else if ((*i).compare("EVOLTS")==0) {
            scale.E = 1.60217646e-12;
        } else if ((*i).compare("MOLES")==0) {
            // Default scaling.
        } else if ((*i).compare("MOLECULES")==0) {
            scale.A = 1.0 / NA;
        }
    }

}
