/*
  Author(s):      Matthew Celnik (msc37), 
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
#include "gpc_string.h"

#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>

using namespace Sprog::IO;
using namespace std;
using namespace Sprog::Kinetics;
using namespace Strings;

// Reads a mechanism from a CHEMKIN input file.
void MechanismParser::ReadChemkin(const std::string &filename, 
                                  Sprog::Mechanism &mech, 
                                  const std::string &thermofile,
                                  int verbose)
{
    string msg;

    // Open file for reading.
    ifstream fin; //(filename.c_str(), ios::in);
    fin.open(filename.c_str(), ios::in);

    if (fin.good()) {
        if (verbose>0) {
            msg = "sprog: Successfully opened CHEMKIN file: " + filename + ".\n";
            printf(msg.c_str());
        }

        CK_STATUS status;
        status.Status = FindKey;
        status.ReadElements = false;
        status.ReadSpecies = false;
        status.ReadReactions = false;
        status.ReadThermo = false;
        status.ThermoFile = thermofile;
        status.Scale = ARRHENIUS(1.0, 1.0, 4.184E7);

        // Clear current mechanism.
        mech.Clear();

        // CHEMKIN files are read in CGS units.
        mech.SetUnits(CGS);

        try {
                parseCK(fin, mech, status, verbose);
        } catch (logic_error &le) {
            fin.close();
            throw;
        } catch (runtime_error &re) {
            fin.close();
            throw;
        }
        
        fin.close();
        return;
    } else {
        // Failed to open file.
        throw std::invalid_argument("Could not open CHEMKIN file: " + filename);
    }
}



void MechanismParser::ReadChemkin(const std::string &filename, 
                                  Sprog::Mechanism &mech, 
                                  const std::string &thermofile,
                                  const std::string &transFile)
{
	// member function added by vniod to enable the reading of transport data.
	// This function utilizes the ReadChemkin function which was provided earlier.

	//std::string mchFile = filename;
	//Sprog::Mechanism mch = mech;
	//std::string thermo = thermofile;
	//int verb = verbose;
	//std::string trans = transFile;

	ReadChemkin(filename,mech,thermofile);

	ReadTransport(transFile,mech);
}

void MechanismParser::ReadTransport(const string& transFile, Sprog::Mechanism &mech){
	// function added by vinod to read the transport data
	ifstream inf;	
	int i;
	inf.open(transFile.c_str());
	if(inf.good()){
		std::string transString;
		vector<std::string> transpVector;
		map< string,vector<string> > transpMap;
		while(!inf.eof()){
			getline(inf,transString);
			if (!isEmpty(transString)){								
				Strings::split(transString,transpVector," \t");						
				i = mech.FindSpecies(convertToCaps(transpVector[0]));
				if( i>= 0){
					//cout << transpVector[0] <<  endl;
					transpMap.insert(pair< std::string,vector<string> >(convertToCaps(transpVector[0]),transpVector));
				}							
				
			}
			
		}
		//cout << "Transport data read successfully\n";		
		mech.setSpeciesTransport(transpMap,mech); 
		// pass the map and mech object forf setting the transport object for species
		
	}else{
		cout << "Transport file not specified \n";
		//exit(1);
	}

}



// CHEMKIN PARSING ROUTINES.

// Reads a CHEMKIN file stream into a std::string.  Removes
// comments and fixes line endings as well.
void MechanismParser::loadCK_File(std::ifstream &fin, std::string &out)
{
    // Current position in the file stream (for later reset).
    int current_pos = fin.tellg();

    // Seek to beginning of stream
    fin.clear();
    fin.seekg (0, std::ios::beg);

    // Convert CK file stream to a string.
    char c;
    // Must try to get the first character in the stream first 
    // because the stream might be empty.  This sets good=false.
    fin.get(c);
    while (fin.good()) {
        // Skipping comment until new line or return character is found.
        if  (c=='!') {
            do {
                fin.get(c);
            } while (fin.good() && (c!='\r') && (c!='\n'));
            if (!fin.good()) break;
        }

        // Fix line endings and tabs.
        if (c=='\r') {
            // replace return character by new line character, if some 
            // systems don't read \r\n as \n
            c = '\n';
        } else if (c=='\t') {
            // replace tab character by a space bar character
            c = ' ';
        }

        // Store this character.
        out.append(1,c);
        
        // Get next character.
        fin.get(c);
    };

    // Reset position of of stream to where it was before
    // we read it.
    fin.clear();
    fin.seekg(current_pos);

    // Convert the string to capital letters.
    out = Strings::convertToCaps(out);

    // Now correct 0.0D+0 exponentials to 0.0E+0.
    std::string::size_type pos = out.find("D+", 0);
    while(pos != out.npos) {
        out.replace(pos, 1, "E");
        pos = out.find("D+", pos+2);
    }
    pos = out.find("D-", 0);
    while(pos != out.npos) {
        out.replace(pos, 1, "E");
        pos = out.find("D-", pos+2);
    }
}

// Loads a mechanism from a CHEMKIN formatted file stream.
void MechanismParser::parseCK(std::ifstream &fin, 
                           Sprog::Mechanism &mech, 
                           Sprog::IO::MechanismParser::CK_STATUS &status,
                           int verbose)
{
    string msg;

    // READ FILE INTO A STRING.

    string ck;
    loadCK_File(fin, ck);

    // PARSE THE STRING.

    string substr = "";
    unsigned int lineno=0, n=0;

    // Read the elements.
    if (verbose>0) printf("sprog: Parsing CHEMKIN file for elements...\n");
    extractCK_Elements(ck, substr, lineno);
    n = parseCK_Elements(substr, mech, lineno, status, verbose);
    status.ReadElements = true;
    if (verbose>0) {
        msg = "sprog: Read ";
        msg = msg + cstr(n) + " elements.\n";
        printf(msg.c_str());
    }

    // Read the species.
    if (verbose>0) printf("sprog: Parsing CHEMKIN file for species...\n");
    extractCK_Species(ck, substr, lineno);
    n = parseCK_Species(substr, mech, lineno, status, verbose);
    status.ReadSpecies = true;
    if (verbose>0) {
        msg = "sprog: Read ";
        msg = msg + cstr(n) + " species.\n";
        printf(msg.c_str());
    }

    // Read the thermo data.
    if (status.ThermoFile == "") {
        // Read thermo data from this file. // This could be bugs
        if (verbose>0) printf("sprog: Parsing CHEMKIN file for thermo data...\n");
        extractCK_Thermo(ck, substr, lineno);
    } else {
        // Read thermo data from thermo file.
        if (verbose>0) printf("sprog: Parsing THERMO file for thermo data...\n");
        // Read the file into a string.
        string thermo = "";
        ifstream tfin;
        tfin.open(status.ThermoFile.c_str(), ios::in);
        loadCK_File(tfin, thermo);
        tfin.close();
        // Extract the thermo data.
        extractCK_Thermo(thermo, substr, lineno);
    }
    n = parseCK_Thermo(substr, mech, lineno, status, verbose);
    status.ReadThermo = true;
    if (n < mech.SpeciesCount()) {
        // Insufficient species thermodynamic properties
        // read from input file.

        // Determine which species are still undefined by counting
        // the element counts of each species.  If a species is 
        // undefined, then it will have no elements.
        for (Mechanism::const_sp_iterator i=mech.SpBegin(); i!=mech.SpEnd(); ++i) {
            if ((*i)->Composition().size() == 0) {
                msg = "sprog: Thermo properties not defined for species ";
                msg = msg + (*i)->Name() + ".\n";
                printf(msg.c_str());
            }
        }

        // Throw an exception now.
        throw invalid_argument("Insufficient species thermo defined in CHEMKIN files.\n"
            "Function: (Sprog, MechanismParser::parseCK)");
    }
    if (verbose>0) {
        msg = "sprog: Read ";
        msg = msg + cstr(n) + " species thermo.\n";
        printf(msg.c_str());
    }

    // Read the reactions.
    if (verbose>0) printf("sprog: Parsing CHEMKIN file for reactions...\n");
    extractCK_Reactions(ck, substr, lineno);
    n = parseCK_Reactions(substr, mech, lineno, status, verbose);
    mech.BuildStoichXRef();
    mech.SetUnits(SI);
    if (verbose>0) {
        msg = "sprog: Read ";
        msg = msg + cstr(n) + " reactions.\n";
        printf(msg.c_str());
    }
}

// Get positions of a CK keyword and END keyword.
MechanismParser::KEY_POS MechanismParser::getCK_KeyPos(const std::string &key, 
                                                       const std::string &ckstr)
{
    KEY_POS pos;

    // Try to locate the full keyword.
    pos.begin = ckstr.find(key, 0);

    // Try to locate smaller, 4-character keyword.
    if ((pos.begin == std::string::npos) || (pos.begin >= ckstr.length())) {
        pos.begin  = ckstr.find(key.substr(0,4), 0);
        pos.length = 4;
    } else {
        pos.length = key.length();
    }

    // Find the END keyword.
    pos.end = ckstr.find("END", pos.begin);

    // Now try to find the other keywords before the END.  If
    // keywords are present before END then this is an error.
    /*if (ckstr.find("ELEM", pos.begin+pos.length) < pos.end){
        throw std::invalid_argument("Erroneous ELEM/ELEMENTS keyword found within " +
                                    key + " definition "
                                    "(Sprog, MechanismParser::getCK_KeyPos).");
    }
    if (ckstr.find("SPEC", pos.begin+pos.length) < pos.end){
        throw std::invalid_argument("Erroneous SPEC/SPECIES keyword found within " +
                                    key + " definition "
                                    "(Sprog, MechanismParser::getCK_KeyPos).");
    }
    if (ckstr.find("REAC", pos.begin+pos.length) < pos.end){
        throw std::invalid_argument("Erroneous REAC/REACTIONS keyword found within " +
                                    key + " definition "
                                    "(Sprog, MechanismParser::getCK_KeyPos).");
    }
    if (ckstr.find("THER", pos.begin+pos.length) < pos.end){
        throw std::invalid_argument("Erroneous THER/THERMO keyword found within " +
                                    key + " definition "
                                    "(Sprog, MechanismParser::getCK_KeyPos).");
    }*/

    // Determine the keyword line number.
    pos.line = 1;
    string::size_type linepos = ckstr.find_last_of("\n", pos.begin);
    while (linepos < ckstr.length()) {
        ++pos.line;
        linepos = ckstr.find_last_of("\n", linepos-1);
    }

    // No erroneous keywords were found.
    return pos;
}

// Extract element names from ck string.
void MechanismParser::extractCK_Elements(const std::string &ckstr, 
                                         std::string &elements,
                                         unsigned int &lineno)
{
    KEY_POS pos;

    // Find the beginning and end of the elements.
    pos = getCK_KeyPos("ELEMENTS", ckstr);

    // Check for valid keyword positions.
    if ((pos.begin == std::string::npos) || (pos.begin >= ckstr.length())) {
        throw std::invalid_argument("ELEM/ELEMENTS keyword not found in "
                                    "CHEMKIN input file "
                                    "(Sprog, MechanismParser::extractCK_Elements).");
    }
    if ((pos.end == std::string::npos) || (pos.end >= ckstr.length())) {
        throw std::invalid_argument("END keyword missing after ELEM/ELEMENTS "
                                    "keywords in CHEMKIN input file "
                                    "(Sprog, MechanismParser::extractCK_Elements).");
    }
    
    // Copy only the elements part from chemkin string, ignoring
    // the keywords.
    elements = ckstr.substr(pos.begin + pos.length, pos.end - pos.begin - pos.length);
    lineno   = pos.line;
}

// Parse the element data in a CHEMKIN input file.
unsigned int MechanismParser::parseCK_Elements(const std::string &elements, 
                                       Sprog::Mechanism &mech,
                                       unsigned int lineno,
                                       Sprog::IO::MechanismParser::CK_STATUS &status,
                                       int verbose)
{
    string tag;
    real val = 0.0;
    Element *last_el = NULL;
    status.Status = BeginParseEl;
    unsigned int n = 0; // Element count.

    // Get iterator to beginning of elements declaration.
    string::const_iterator c = elements.begin();

    while((status.Status!=End) && (status.Status!=Fail) && (c!=elements.end())) {
        switch(status.Status) {
            case BeginParseEl:
                // Here we are looking for the start of an element name.
                if (isLetter(*c)) {
                    // This is the start of an element name.
                    tag = *c;
                    status.Status = ParseEl;

                    // Clear the last read element.
                    last_el = NULL;
                } else if (*c == '/') {
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
                        throw invalid_argument("Invalid / character in element definition "
                                               "on line number " + cstr(lineno) +
                                               " (Sprog, MechanismParser::parseCK_Elements).");
                    }
                } else if (*c == '!') {
                    // This is a comment.
                    status.Status = ParseElComment;
                } else if (*c == '\n') {
                    // This is a new line.
                    ++lineno;
                } else if (!isWhiteSpace(*c)) {
                    // This must be an invalid character.
                    status.Status = Fail;
                    string msg = "Invalid character in element definition on line number "+
                        cstr(lineno)+" of CHEMKIN input file.\n"
                        "Function: (Sprog, MechanismParser::parseCK_Elements).";
                    throw std::invalid_argument(msg);

                }
                break;
            case ParseEl:
                // We are now reading an element name.
                if (isLetterOrNum(*c)) {
                    // This is a valid character for the element name.
                    tag.append(&(*c), 1);
                } else if (isWhiteSpace(*c)) {
                    if (*c == '\n') ++lineno;

                    // This tag contains a complete element name, so add a new
                    // element to the mechanism and set its name.
                    last_el = mech.AddElement();
                    last_el->SetName(tag);
                    // Attempt to find the element weight in the library of
                    // know elements.
                    last_el->SetMolWtFromLibrary();
                    // Start searching for the next element.
                    status.Status = BeginParseEl;
                    tag = "";
                    // Increment element count.
                    ++n;
                    if (verbose>1) printf((string("sprog: Added element ") + last_el->Name() + ".\n").c_str());
                } else if (*c == '/') {
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
                if (*c == '/') {
                    // This is the end of the mol. wt. definition.  Attempt to convert it
                    // to a number.
                    val = cdble(tag);
                    if (val > 0.0) {
                        // Set the mol. wt. of the last read element.
                        last_el->SetMolWt(val);
                    } else {
                        // The mol. wt. is invalid.
                        status.Status = Fail;
                        throw invalid_argument("An invalid mol. wt. was found for element " + 
                                               last_el->Name() + " on line number " + 
                                               cstr(lineno) +
                                               "(Sprog, MechanismParser::parseCK_Elements).");
                    }
                } else {
                    // Add character to number.
                    tag.append(&(*c), 1);
                }
                break;
            case ParseElComment:
                // We are reading a comment in the element list.  A comment is only
                // terminated by a line ending.
                if (*c == '\n') {
                    ++lineno;
                    status.Status = BeginParseEl;
                }
                break;
            default:
                throw runtime_error("Illegal parsing flag encountered "
                                    "(Sprog, MechanismParser::parseCK_Elements).");
        }

        // Increment string location.
        ++c;
    }

    return n;
}

// Extract species names from CHEMKIN string.
void MechanismParser::extractCK_Species(const std::string &ckstr, 
                                        std::string &species,
                                        unsigned int &lineno)
{
    KEY_POS pos;

    // Find the beginning and end of the elements.
    pos = getCK_KeyPos("SPECIES", ckstr);

    // Check for valid keyword positions.
    if ((pos.begin == std::string::npos) || (pos.begin >= ckstr.length())) {
        throw std::invalid_argument("SPEC/SPECIES keyword not found in "
                                    "CHEMKIN input file.\n"
                                    "Function: (Sprog, MechanismParser::extractCK_Species).");
    }
    if ((pos.end == std::string::npos) || (pos.end >= ckstr.length())) {
        throw std::invalid_argument("END keyword missing after SPEC/SPECIES "
                                    "keywords in CHEMKIN input file.\n"
                                    "Function: (Sprog, MechanismParser::extractCK_Species).");
    }
    
    // Copy only the elements part from chemkin string, ignoring
    // the keywords.
    species = ckstr.substr(pos.begin + pos.length, pos.end - pos.begin - pos.length);
    lineno   = pos.line;
}

// Parse the species data in a CHEMKIN input file 
// (elements must already have been read).
unsigned int MechanismParser::parseCK_Species(const std::string &species, 
                                      Sprog::Mechanism &mech,
                                      unsigned int lineno,
                                      Sprog::IO::MechanismParser::CK_STATUS &status,
                                      int verbose)
{
    string tag;
    Species *last_sp;
    status.Status = BeginParseSp;
    unsigned int n = 0; // Species count.

    // Get iterator to first character in the species defintion.
    string::const_iterator c = species.begin();

    // Loop over all characters in the species definition.
    while((status.Status!=End) && (status.Status!=Fail) && (c!=species.end())) {
        switch(status.Status) {
            case BeginParseSp:
                // We are searching for the beginning of a species name.
                if (isLetter(*c)) {
                    // This is the start of a species name.
                    tag = *c;
                    status.Status = ParseSp;

                    // Clear the last read species.
                    last_sp = NULL;
                } else if (*c == '!') {
                    // This is a comment.
                    status.Status = ParseSpComment;
                } else if (*c == '\n') {
                    // A new line.
                    ++lineno;
                } else if (!isWhiteSpace(*c)) {
                    // This must be an invalid character.
                    status.Status = Fail;
                    throw invalid_argument("Invalid character in species definition on"
                        "line number " + cstr(lineno) + "of CHEMKIN input file.\n"
                        "Function: (Sprog, MechanismParser::parseCK_Species).");
                }
                break;
            case ParseSp:
                // We are now reading a species name.
                if (isLetterOrNum(*c) || (*c=='(') || (*c==')') ||
                    (*c=='*') || (*c=='-') || (*c=='+')) {
                    // This is a valid character for the species name.
                    tag.append(&(*c), 1);
                } else if (isWhiteSpace(*c)) {
                    // This tag contains a full species name, so add a new
                    // species to the mechanism and set its name.
                    ++n;
                    last_sp = mech.AddSpecies();
                    last_sp->SetName(tag);
                    if (verbose>1) printf((string("sprog: Added species ") + tag + ".\n").c_str());
                    // Begin searching for the next species name.
                    status.Status = BeginParseSp;
                }
                break;
            case ParseSpComment:
                // We are reading a comment in the species list.  A comment is only
                // terminated by a line ending.
                if (*c == '\n') {
                    ++lineno;
                    status.Status = BeginParseSp;
                }
                break;
            default:
                throw runtime_error("Illegal parsing flag encountered "
                                    "(Sprog, MechanismParser::parseCK_Species).");
        }

        // Increment iterator.
        ++c;
    }

    return n;
}

// Extract thermo data from CHEMKIN string.
void MechanismParser::extractCK_Thermo(const std::string &ckstr, 
                                       std::string &thermo,
                                       unsigned int &lineno)
{
    KEY_POS pos;

    // Find the beginning and end of the elements.
    pos = getCK_KeyPos("THERMO", ckstr);

    // Check for valid keyword positions.
    if ((pos.begin == std::string::npos) || (pos.begin >= ckstr.length())) {
        throw std::invalid_argument("THER/THERMO keyword not found in "
                                    "CHEMKIN input file.\n"
                                    "Function: (Sprog, MechanismParser::extractCK_Species).");
    }
    if ((pos.end == std::string::npos) || (pos.end >= ckstr.length())) {
        throw std::invalid_argument("END keyword missing after THER/THERMO "
                                    "keywords in CHEMKIN input file.\n"
                                    "Function: (Sprog, MechanismParser::extractCK_Species).");
    }
    
    // Unlike other CHEMKIN data, the thermo has a strict
    // 80 character format.  We must locate the first character
    // on the line following the THERMO keyword.  This involves searching
    // for the line ending \n.
    int i = -1;
    while (ckstr[pos.begin+pos.length+(++i)] != '\n');
    ++i;

    // Copy only the thermo part from chemkin string, ignoring
    // the keywords.
    thermo = ckstr.substr(i + pos.length, pos.end - i - pos.length);
    lineno = pos.line;
}

// Reads CHEMKIN formatted thermo data for all species in the given mechanism from
// the supplied file stream.
unsigned int MechanismParser::parseCK_Thermo(const std::string &thermo, Sprog::Mechanism &mech, 
                                    unsigned int lineno,
                                    Sprog::IO::MechanismParser::CK_STATUS &status,
                                    int verbose)
{
    unsigned int i, j1, j2, n=0;
    char c;
    string tag, spname;

    real trange[3], lowT, highT, commT;
    Thermo::THERMO_PARAMS up; // Coeffs for upper T interval.
    Thermo::THERMO_PARAMS lp; // Coeffs for lower T interval.
    vector<string> els(4);
    int nels[4];

    Species *sp;

    status.Status = ParseTherm;
    tag.resize(81,' ');

    // The first line holds the temperature ranges for 
    // 2 sets of coefficients.
    j1 = j2 = 0;
    tag = "\0";
    while ((tag == "\0") && (j2<thermo.length())){
        j1 = j2;
        j2 = thermo.find("\n", j1) + 1;
        tag = thermo.substr(j1, j2-j1); ++lineno;
    }
    trange[0] = cdble(tag.substr(0,10));
    trange[1] = cdble(tag.substr(10,10));
    trange[2] = cdble(tag.substr(20,10));

    while ((status.Status!=End) && (status.Status!=Fail) && (j2 < thermo.length())) {
        // Get the next line, ignoring blank lines.
        tag = "\0";
        while ((tag=="\0") && (j2<thermo.length())) {
            j1 = j2;
            j2 = thermo.find("\n", j1) + 1;
            tag = thermo.substr(j1, j2-j1-1); ++lineno;
        }

        // Check for end of data.
        if (j2>=thermo.length()) break;

        if ((tag != "") && (tag[0] != '!')) {
            // Get species name from line.
            c = tag[0];
            for(i=1; (i<18)&&(!isWhiteSpace(c)); ++i) {
                c = tag[i];
            }
            spname = tag.substr(0,i-1);

            // Find species in the mechanism.
            sp = mech.GetSpecies(spname);

            if ((sp != NULL) && (sp->Composition().size() == 0)) {
                // Species was found in the mechanism, so we need to
                // read the thermo data.
                if (verbose>2) printf((string("sprog: Reading ") + spname + " thermo.\n").c_str());
                ++n;
                
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
                tag = "\0";
                while ((tag=="\0") && (j2<thermo.length() )) {
                    j1 = j2;
                    j2 = thermo.find("\n", j1) + 1;
                    tag = thermo.substr(j1, j2-j1-1); ++lineno;
                }
                up.Count = 7;
                up.Params[0] = cdble(tag.substr(0,15));
                up.Params[1] = cdble(tag.substr(15,15));
                up.Params[2] = cdble(tag.substr(30,15));
                up.Params[3] = cdble(tag.substr(45,15));
                up.Params[4] = cdble(tag.substr(60,15));

                // Get 3rd line and read thermo coeffs.
                tag = "\0";
                while ((tag=="\0") && (j2<thermo.length())) {
                    j1 = j2;
                    j2 = thermo.find("\n", j1) + 1;
                    tag = thermo.substr(j1, j2-j1-1); ++lineno;
                }
                up.Params[5] = cdble(tag.substr(0,15));
                up.Params[6] = cdble(tag.substr(15,15));
                lp.Count = 7;
                lp.Params[0] = cdble(tag.substr(30,15));
                lp.Params[1] = cdble(tag.substr(45,15));
                lp.Params[2] = cdble(tag.substr(60,15));

                // Get 4th line and read thermo coeffs.
                tag = "\0";
                while ((tag=="\0") && (j2<thermo.length())) {
                    j1 = j2;
                    j2 = thermo.find("\n", j1) + 1;
                    tag = thermo.substr(j1, j2-j1-1); ++lineno;
                }
                lp.Params[3] = cdble(tag.substr(0,15));
                lp.Params[4] = cdble(tag.substr(15,15));
                lp.Params[5] = cdble(tag.substr(30,15));
                lp.Params[6] = cdble(tag.substr(45,15));

                // Now save the data to the species:
                
                // Add the elemental composition to the species.
                if (verbose>3) printf("sprog: Elemental composition:\n");
                for (i=0; i!=4; ++i) {
                    // Need to check the element name is valid.
                    try {
                        if (!isWhiteSpace(*els[i].substr(0,1).c_str())) {
                            if (!isWhiteSpace(*els[i].substr(1,1).c_str())) {
                                sp->AddElement(els[i], nels[i]);
                                if (verbose>3) printf((string("sprog: ") + cstr(nels[i]) + "x " + els[i] + ".\n").c_str());
                            } else {
                                sp->AddElement(els[i].substr(0,1), nels[i]);
                                if (verbose>3) printf((string("sprog: ") + cstr(nels[i]) + "x " + els[i].substr(0,1) + ".\n").c_str());
                            }                    
                        }
                    } catch (std::invalid_argument &ia) {
                        throw ia;
                    }
                }
                
                if (verbose>3) {
                    printf("sprog: Lower range parameters:\n");
                    printf((string("sprog: a1 = ") + cstr(lp.Params[0]) + ".\n").c_str());
                    printf((string("sprog: a2 = ") + cstr(lp.Params[1]) + ".\n").c_str());
                    printf((string("sprog: a3 = ") + cstr(lp.Params[2]) + ".\n").c_str());
                    printf((string("sprog: a4 = ") + cstr(lp.Params[3]) + ".\n").c_str());
                    printf((string("sprog: a5 = ") + cstr(lp.Params[4]) + ".\n").c_str());
                    printf((string("sprog: a6 = ") + cstr(lp.Params[5]) + ".\n").c_str());
                    printf((string("sprog: a7 = ") + cstr(lp.Params[6]) + ".\n").c_str());

                    printf("sprog: Upper range parameters:\n");
                    printf((string("sprog: a1 = ") + cstr(up.Params[0]) + ".\n").c_str());
                    printf((string("sprog: a2 = ") + cstr(up.Params[1]) + ".\n").c_str());
                    printf((string("sprog: a3 = ") + cstr(up.Params[2]) + ".\n").c_str());
                    printf((string("sprog: a4 = ") + cstr(up.Params[3]) + ".\n").c_str());
                    printf((string("sprog: a5 = ") + cstr(up.Params[4]) + ".\n").c_str());
                    printf((string("sprog: a6 = ") + cstr(up.Params[5]) + ".\n").c_str());
                    printf((string("sprog: a7 = ") + cstr(up.Params[6]) + ".\n").c_str());
                }

                // Add the thermo parameters.
                sp->SetThermoStartTemperature(lowT);
                sp->AddThermoParams(commT, lp);
                sp->AddThermoParams(highT, up);
            } else {
                // Species was not found in the mechanism, so we need to
                // skip lines to the next species.
                if (verbose>2) printf((string("sprog: Skipping ") + spname + " thermo.\n").c_str());
                j2 = thermo.find("\n", j1=j2) + 1;
                j2 = thermo.find("\n", j1=j2) + 1;
                j2 = thermo.find("\n", j1=j2) + 1;
                lineno += 3;
            }
        }
    }

    return n;
}

// Extract element names from ck string.
void MechanismParser::extractCK_Reactions(const std::string &ckstr, 
                                          std::string &reac,
                                          unsigned int &lineno)
{
    KEY_POS pos;

    // Find the beginning and end of the elements.
    pos = getCK_KeyPos("REACTIONS", ckstr);

    // Check for valid keyword positions.
    if ((pos.begin == std::string::npos) || (pos.begin >= ckstr.length())) {
        throw std::invalid_argument("REAC/REACTIONS keyword not found in "
                                    "CHEMKIN input file "
                                    "(Sprog, MechanismParser::extractCK_Reactions).");
    }
    if ((pos.end == std::string::npos) || (pos.end >= ckstr.length())) {
        throw std::invalid_argument("END keyword missing after REAC/REACTIONS "
                                    "keywords in CHEMKIN input file "
                                    "(Sprog, MechanismParser::extractCK_Reactions).");
    }
    
    // Copy only the elements part from chemkin string, ignoring
    // the keywords.
    reac   = ckstr.substr(pos.begin + pos.length, pos.end - pos.begin - pos.length);
    lineno = pos.line;
}

// Reads all the chemical reactions from the CHEMKIN formatted file stream.
unsigned int MechanismParser::parseCK_Reactions(const std::string &reac, Sprog::Mechanism &mech, 
                                        unsigned int lineno,
                                        Sprog::IO::MechanismParser::CK_STATUS &status,
                                        int verbose)
{
    string tag="", rxndef="";
    bool fcont = false; //, fdup = false;
    unsigned int n = 0; // Reaction count.

    Kinetics::Reaction *last_rxn = NULL;

    status.Status = ParseRxn;

    // First try to read the reaction units.
    unsigned int i = reac.find("\n");
    tag = reac.substr(0, i);
    parseCK_Units(tag, status.Scale);

    // Get iterator to beginning of reactions declaration.
    string::const_iterator c = reac.begin()+i;
//    ++lineno;

    while((status.Status!=End) && (status.Status!=Fail) && (c!=reac.end())) {
        switch(status.Status) {
            case BeginParseRxn:
                fcont = false;
                // Clear white space at the beginning of the line.
                if (!isWhiteSpace(*c)) {
                    rxndef = *c;
                    status.Status = ParseRxn;
                } else if (*c=='\n') {
                    ++lineno;
                }
                break;
            case ParseRxn:
                // We are now reading a reaction definition.  Currently we are building
                // a string rxndef to hold to the reaction info.  This string will be
                // passed to the function parseCK_Reaction for translation.

                if (*c=='\n') {
                    // This is the end of the line.  If there was a continuation character '&'
                    // then we continue reading the reaction string, otherwise it gets passed
                    // to parseCK_Reaction for translation.
                    if (!fcont) {
                        // Before parsing the string we should check that it doesn't
                        // hold auxilliary info for the previous reaction.
                        if (!parseCK_RxnAux(rxndef, last_rxn, mech, status.Scale, status, lineno, verbose)) {
                            // No aux info.  First add the last reaction to the
                            // mechanism.
                            if (last_rxn != NULL) {
                                mech.AddReaction(last_rxn);
                                delete last_rxn;
                                last_rxn = NULL;
                                ++n;
                            }

                            // No aux info, so parse the reaction definition.
                            if (rxndef.length() > 0) {
                                last_rxn = parseCK_Reaction(rxndef, mech, lineno, status, verbose);
                            }
                        }
                        status.Status = BeginParseRxn;
                    } else {
                        // Reset continuation flag to avoid reading all
                        // lines at once.
                        fcont = false;
                    }
                    // Increment line number.
                    ++lineno;
                } else if (*c == '&') {
                    // Continuation character located.
                    fcont = true;
                } else {
                    if (fcont && !isWhiteSpace(*c)) {
                        // There should not be any non-white space
                        // characters after a continuation & character.
                        throw invalid_argument("Illegal character after continuation '&' "
                                "on line number " + cstr(lineno) + ".\n"
                                "Function: (Sprog, MechanismParser::parseCK_Reactions).");
                    } else {
                        rxndef.append(&(*c), 1);
                    }
                }
                break;
            default:
                throw runtime_error("Illegal parsing flag encountered "
                                    "(Sprog, MechanismParser::parseCK_Reactions).");
        }

        // Increment string iterator.
        ++c;
    }

    // Add the last read reaction, if exists.
    if (last_rxn != NULL) {
        ++n;
        mech.AddReaction(last_rxn);
        delete last_rxn;
        last_rxn = NULL;
    }

    return n;
}

// Parses a string and builds a Reaction object from the data therein.
Sprog::Kinetics::Reaction *const MechanismParser::parseCK_Reaction(
                                    const std::string &rxndef, 
                                    Sprog::Mechanism &mech, unsigned int lineno,
                                    Sprog::IO::MechanismParser::CK_STATUS &status,
                                    int verbose)
{
    // The reaction string contains four pieces of information:  
    // the reaction formula, A, n & E.
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
    if (i < rxndef.length()) {
        // Found reversible delimiter.
        ilast_reac = i - 1;
        ifirst_prod = i + 3;
        frev = true;
    } else {
        i = rxndef.find_first_of(">");
        if (i < rxndef.length()) {
            // Found irreversible delimiter.
            ilast_reac = i - 2;
            ifirst_prod = i + 1;
            frev = false;
        } else {
            i = rxndef.find_first_of("=");
            if ((i!=rxndef.npos) && (i<rxndef.length())) {
                // Found reversible delimiter.
                ilast_reac = i - 1;
                ifirst_prod = i + 1;
                frev = true;
            } else {
                // No delimiter found.
                throw invalid_argument("No reactant/product delimiter found "
                        "in reaction definition on line " + cstr(lineno) + 
                        " of the CHEMKIN input file.\n"
                        "Function: (Sprog, MechanismParser::parseCK_Reaction)");
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
    parseCK_RxnSpStoich(rxndef.substr(0,ilast_reac+1), mech, rmui, rmuf, 
                        fistb, fisfo, strtb, lineno);

    // Get the product stoichiometry and set it in the reaction.
    parseCK_RxnSpStoich(rxndef.substr(ifirst_prod,ilast_prod-ifirst_prod+1), 
                        mech, pmui, pmuf, fistb, fisfo, strtb, lineno);

    // Now we have got the reactants and the products, we must now parse the Arrhenius coefficients.
    // This is quite simple as they must be space delimited.
    split(rxndef.substr(ilast_prod), arrstr, " ");
    arr.A = cdble(arrstr[0]) * status.Scale.A;
    arr.n = cdble(arrstr[1]);
    arr.E = cdble(arrstr[2]) * status.Scale.E;

    if (verbose>2) {
        printf((string("sprog: Parsing (line ") + cstr(lineno) + 
               ") " + rxndef.substr(0,ilast_prod) + 
               " (A=" + arrstr[0] + ", "
                 "n=" + arrstr[1] + ", "
                 "E=" + arrstr[2] + ").\n").c_str());
    }

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
                                                  std::string &thirdbody, 
                                                  unsigned int lineno)
{

    vector<string> species;
    vector<string>::iterator k;
    string str;
    string::size_type i = 0, j = 0;
    // bool search = true, plus = true, fosym = false;

    // Locate the first species delimiter.
    i = sp.find_first_not_of(" ");
    if ((i!=sp.npos) && (i<sp.length())) {
        j = sp.find_first_of("+",i+1);
    } else {
        j = sp.npos;
    }

    // Loop over all delimiters until we reach the end of the string.
    while ((j!=sp.npos) && (j<sp.length())) {
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
            i = j + 1; j = sp.find_first_of(")",i);
            thirdbody.assign(sp.begin()+i, sp.begin()+j);
            if (thirdbody.compare("M")==0) thirdbody = "";

        } else {
            // There is no bracket to worry about as this is not a fall-off symbol.
            str = sp.substr(i, j-i);
            species.push_back(str.substr(0, str.find_last_not_of(" ")+1));
        }

        // Find next delimiter.
        i = sp.find_first_not_of(" ", j+1);
        if ((i!=sp.npos) && (i<sp.length())) {
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

    // Now we must separate the stoichiometry from the species' names.  
    // Also we can check here for third-bodies or fall-off reactions.
    string sym, mu;
    int intmu;
    real floatmu;

    for (k=species.begin(); k!=species.end(); k++) {
        // Separate stoichiometry from specie's name.
        i = (*k).find_first_not_of("0123456789.");
        if ((i!=sp.npos) && (i<sp.length())) {
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
                // Before we save the species stoichiometry we need to check
                // if it is integer or real.
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
                throw invalid_argument(
                    "Unknown species in reaction on line " + cstr(lineno) + 
                    "of the CHEMKIN input file.\n"
                    "Function: (Sprog, MechanismParser::parseCK_RxnSpStoich)");
            }
        }
    }
}

// Parses auxilliary reaction information.
bool Sprog::IO::MechanismParser::parseCK_RxnAux(const std::string &rxndef, 
                                             Sprog::Kinetics::Reaction *last_rxn, 
                                             const Sprog::Mechanism &mech,
                                             Sprog::Kinetics::ARRHENIUS scale,
                                             Sprog::IO::MechanismParser::CK_STATUS &status,
                                             unsigned int lineno, int verbose)
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

    while ((i != rxndef.npos) && (i < rxndef.length())) {
        // Now i-j is the keyword (or species name), j-k are the parameters.

        // Get the keyword.
        key = rxndef.substr(i, j-i);
        key = key.substr(0, key.find_first_of(' '));

        if ((j != rxndef.npos) && (j < rxndef.length())) {
            // Check for unterminated set of parameters.
            if ((k == rxndef.npos) && (k < rxndef.length())) {
                // An unterminated set of parameters.  This is an error.
                status.Status = Fail;
                throw invalid_argument(
                    "Unterminated set of parameters in auxilliary data on "
                    "line " + cstr(lineno) + " of the CHEMKIN input file.\n"
                    "Function: (Sprog, MechanismParser::parseCK_RxnAux)");
            }

            // Get the parameters.
            vals = rxndef.substr(j+1, k-j-1);
            split(vals, params, " ");
            foundaux = true;
        } else {
            // There are no parameters with this keyword.
            vals = "";
            params.clear();

            if ((key.compare("DUPLICATE")==0) || (key.compare("DUP")==0)) {
                if (verbose>2) {
                    printf((string("sprog: DUPLICATE keyword found on line ") + 
                           cstr(lineno) + ".\n").c_str());
                }
                foundaux = true;
            }
            return foundaux;
        }

        // Decide next action based on the keyword.
        if (key.compare("LOW") == 0) {
            // Low-pressure limit for fall-off reaction (3 parameters).
            if (params.size() > 2) {
                if (verbose>3) {
                    printf((string("sprog: Parsing (line ") + cstr(lineno) + 
                        ") low pressure limit parameters ("
                        "A=" + params[0] + ", "
                        "n=" + params[1] + ", "
                        "E=" + params[2] + ").\n").c_str());
                }
                last_rxn->SetLowPressureLimit(ARRHENIUS(cdble(params[0]) * scale.A,
                                                        cdble(params[1]) * scale.n,
                                                        cdble(params[2]) * scale.E));
            } else {
                // Insufficient parameters.
                throw invalid_argument(
                    "Insufficient parameters for LOW keyword in auxilliary data on "
                    "line " + cstr(lineno) + " of the CHEMKIN input file.\n"
                    "Function: (Sprog, MechanismParser::parseCK_RxnAux)");
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
                    if (verbose>3) {
                        printf((string("sprog: Parsing (line ") + cstr(lineno) + 
                            ") TROE-4 fall-off parameters ("
                            "a=" + params[0] + ", "
                            "T***=" + params[1] + ", "
                            "T*=" + params[2] + ", "
                            "T**=" + params[3] + ").\n").c_str());
                    }
                } else {
                    // 3-parameter TROE form.
                    if (verbose>3) {
                        printf((string("sprog: Parsing (line ") + cstr(lineno) + 
                            ") TROE-3 fall-off parameters ("
                            "a=" + params[0] + ", "
                            "T***=" + params[1] + ", "
                            "T*=" + params[2] + ").\n").c_str());
                    }
                    foparams[3] = 0.0;
                    last_rxn->SetFallOffParams(Troe3, foparams);
                }
            } else {
                // Insufficient parameters.
                throw invalid_argument(
                    "Insufficient parameters for TROE keyword in auxilliary data on "
                    "line " + cstr(lineno) + " of the CHEMKIN input file.\n"
                    "Function: (Sprog, MechanismParser::parseCK_RxnAux)");
            }
        } else if (key.compare("SRI") == 0) {
            // SRI pressure-dependent reaction (3 or 5 parameters).
            if (params.size() > 4) {
                if (verbose>3) {
                    printf((string("sprog: Parsing 5-param form SRI parameters on line ") + 
                           cstr(lineno) + ".\n").c_str());
                }
                // 5-parameter form.
                foparams[0] = cdble(params[0]);
                foparams[1] = cdble(params[1]);
                foparams[2] = cdble(params[2]);
                foparams[3] = cdble(params[3]);
                foparams[4] = cdble(params[4]);
                last_rxn->SetFallOffParams(SRI, foparams);
            } else if (params.size() > 2) {
                if (verbose>3) {
                    printf((string("sprog: Parsing 3-param form SRI parameters on line ") + 
                           cstr(lineno) + ".\n").c_str());
                }
                //3-parameter form.
                foparams[0] = cdble(params[0]);
                foparams[1] = cdble(params[1]);
                foparams[2] = cdble(params[2]);
                foparams[3] = 1.0;
                foparams[4] = 0.0;
                last_rxn->SetFallOffParams(SRI, foparams);
            } else {
                // Insufficient parameters.
                throw invalid_argument(
                    "Insufficient parameters for SRI keyword in auxilliary data on "
                    "line " + cstr(lineno) + " of the CHEMKIN input file.\n"
                    "Function: (Sprog, MechanismParser::parseCK_RxnAux)");
            }
        } else if (key.compare("USER") == 0) {
            // User-defined pressure fall-off reaction (5 parameters).
            // NOT IN USE.
        } else if (key.compare("REV") == 0) {
            // Reverse rate Arrhenius parameters (3 required).
            if (params.size() > 2) {
                if (verbose>3) {
                    printf((string("sprog: Parsing reversible Arrhenius parameters on line ") + 
                           cstr(lineno) + ".\n").c_str());
                }
                last_rxn->SetRevArrhenius(ARRHENIUS(cdble(params[0]) * scale.A,
                                                    cdble(params[1]) * scale.n,
                                                    cdble(params[2]) * scale.E));
            } else {
                // Insufficient parameters.
                throw invalid_argument(
                    "Insufficient parameters for REV keyword in auxilliary data on "
                    "line " + cstr(lineno) + " of the CHEMKIN input file.\n"
                    "Function: (Sprog, MechanismParser::parseCK_RxnAux)");
            }
        } else if (key.compare("LT") == 0) {
            // Landau-Teller reaction parameters (2 required).
            if (params.size() > 1) {
                if (verbose>3) {
                    printf((string("sprog: Parsing Landau-Teller parameters on line ") + 
                           cstr(lineno) + ".\n").c_str());
                }
                last_rxn->SetLTCoeffs(LTCOEFFS(cdble(params[0]), cdble(params[1])));
            } else {
                // Insufficient parameters.
                throw invalid_argument(
                    "Insufficient parameters for LT keyword in auxilliary data on "
                    "line " + cstr(lineno) + " of the CHEMKIN input file.\n"
                    "Function: (Sprog, MechanismParser::parseCK_RxnAux)");
            }
        } else if (key.compare("RLT") == 0) {
            // Reverse Landau-Teller parameters (2 required).
             if (params.size() > 1) {
                if (verbose>3) {
                    printf((string("sprog: Parsing reverse Landau-Teller parameters on line ") + 
                           cstr(lineno) + ".\n").c_str());
                }
                last_rxn->SetRevLTCoeffs(LTCOEFFS(cdble(params[0]), cdble(params[1])));
            } else {
                // Insufficient parameters.
                throw invalid_argument(
                    "Insufficient parameters for RLT keyword in auxilliary data on "
                    "line " + cstr(lineno) + " of the CHEMKIN input file.\n"
                    "Function: (Sprog, MechanismParser::parseCK_RxnAux)");
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
            // Not a keyword at all, but a species name (hopefully!) 
            // for third-body efficiencies.
            sp = mech.FindSpecies(key);
            if ((sp >= 0) && (params.size() > 0)) {
                if (verbose>3) {
                    printf((string("sprog: Parsing (line ") + cstr(lineno) + 
                        ") third-body enhancement of " + key + 
                        " (" + params[0] + ").\n").c_str());
                }
                val = cdble(params[0]);
                last_rxn->AddThirdBody(sp, val);
            } else {
                // Unrecognised species or missing efficiency.
                throw invalid_argument(
                    "Error parsing auxilliary reaction data.  Unknown keyword on "
                    "line " + cstr(lineno) + " of the CHEMKIN input file.\n"
                    "Function: (Sprog, MechanismParser::parseCK_RxnAux)");
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
void Sprog::IO::MechanismParser::parseCK_Units(const std::string &rxndef, 
                                               Sprog::Kinetics::ARRHENIUS &scale)
{
    // Split the string into parts.
    vector<string> parts;
    split(rxndef, parts, " ");

    // Default scaling.
    scale.A = 1.0;
    scale.n = 1.0;
    scale.E = 4.184e7; // Convert cal/mol -> J/mol.

    vector<string>::iterator i;
    for (i=parts.begin(); i!=parts.end(); ++i) {
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
