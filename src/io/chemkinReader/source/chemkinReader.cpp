/*
 * chemkinReader.cpp
 *
 *  Created on: Jun 22, 2011
 *      Author: lrm29
 */

#include "chemkinReader.h"
#include "stringFunctions.h"
#include "thermoParser.h"
#include "transportParser.h"
#include "reactionParser.h"
#include <iostream>
#include <sstream>
#include <map>

using namespace std;
using namespace boost;

const regex IO::ChemkinReader::elementListRegex
("ELEM(?:|ENT|ENTS)\\s+(.*?)\\s+END");

const regex IO::ChemkinReader::elementSingleRegex
("(\\w+)");

const regex IO::ChemkinReader::speciesListRegex
("SPEC(?:|IE|IES)\\s+(.*?)\\s+END");

const regex IO::ChemkinReader::mainSurfaceListRegex
("MATERIAL\\s+(.*?)\\s+(SITE(?:|S)\\s*\\/\\s*(.*?)\\s+END)\\s+REACTION");

const regex IO::ChemkinReader::surfaceSpeciesListRegex
("SITE(?:|S)\\s*\\/\\s*(.*?)\\s*\\/\\s+SDEN\\s*\\/\\s*(.*?)\\s*\\/\\s+(.*?)\\s+END");

/*
const regex IO::ChemkinReader::surfaceSpeciesListRegex
("MATERIAL\\s+(.*?)\\s+SITE(?:|S)\\s*\\/\\s*(.*?)\\s*\\/\\s+SDEN\\s*\\/\\s*(.*?)\\s*\\/\\s+(.*?)\\s+END");
*/

const regex IO::ChemkinReader::speciesSingleRegex
("([^ \\n]+)"); // correct
//("(\\w+(?:|\\(.*?\\)))");

const regex IO::ChemkinReader::surfaceSpeciesSingleRegex
("([^ \\n]+)");
//("([^ \\n]+)""(?:\\/(\\d+)\\/)");
//("([^ \\n\\/]+)(?:|\\/(\\w+)\\/)");
//("([^ \\n\\/]+(?:|\\/(\\w+)\\/))"); 
 
const regex IO::ChemkinReader::surfaceSpeciesSingleSubRegex
("([^ \\n]+)""(?:\\/(\\d+)\\/)");
//("(.*?)\\/(.*?)\\/|(.*?)");
//("(.*?)(?:|\\/(.*?)\\/)"); 

const regex IO::ChemkinReader::reactionListRegex
("REAC(?:|TION|TIONS)\\s+(.*?)\\s+END");

const regex IO::ChemkinReader::unitsRegex
("REAC(?:|TION|TIONS)\\s+\\b(CAL/MOLE|KCAL/MOLE|JOULES/MOL|KJOULES/MOLE|KJOU/MOL|KJOU/MOLE|KELVINS|EVOLTS|MOLES|MOLECULES)\\b");

const regex IO::ChemkinReader::surfaceUnitsRegex
("REAC(?:|TION|TIONS)\\s+(.*?)\\s+(CAL/MOLE|KCAL/MOL|JOULES/MOL|KJOULES/MOL|KJOU/MOL|KJOU/MOL|KELVINS|EVOLTS|MOLES|MOLECULES)");


IO::ChemkinReader::ChemkinReader
(
    const string chemfile,
    const string thermfile,
    const string transfile
)
:
    chemfile_(chemfile),
    chemSurfFile_("NOT READ"), // Added by mm864
    thermfile_(thermfile),
    thermSurfFile_("NOT READ"), // Added by mm864
    transfile_(transfile),
    chemfilestring_(convertToCaps(replaceComments(fileToString(chemfile_)))),
    globalUnits_("NO GLOBAL UNITS"), 
    surfUnits_("NO GLOBAL UNITS") // Added by mm864	
{
    checkChemFile();
}


IO::ChemkinReader::ChemkinReader
(
    const string chemfile,
	const string chemSurfFile, // Added by mm864
    const string thermfile,
	const string thermSurfFile, // Added by mm864
    const string transfile
)
:
    chemfile_(chemfile),
	chemSurfFile_(chemSurfFile), // Added by mm864
    thermfile_(thermfile),
    thermSurfFile_(thermSurfFile), // Added by mm864
	transfile_(transfile),
    chemfilestring_(convertToCaps(replaceComments(fileToString(chemfile_)))),
	chemSurfFilestring_(convertToCaps(replaceComments(fileToString(chemSurfFile_)))), // Added by mm864
    globalUnits_("NO GLOBAL UNITS"), 
    surfUnits_("NO GLOBAL UNITS") // Added by mm864 	
{
    checkChemSurfFile();
}

bool IO::ChemkinReader::checkChemFile()
{
    cout << "Checking the format of the chem.inp file." << endl;

    const regex fileStructure
    (
        "ELEM(?:|ENT|ENTS).*?END.*?"
        "SPEC(?:|IE|IES).*?END.*?"
        "REAC(?:|TION|TIONS).*?END"
    );

    //! \todo NEED TO ADD REGEXES FOR (LT)|(RLT) IN REACTIONPARSER ASAP
    const regex unsupported
    (
        "(TDEP)|(EXCI)|(JAN)|(FIT1)|"
        "(HV)|(MOME)|(FORD)|(RORD)|(UNITS)|(HIGH)|(USER)|"
        "(LT)|(RLT)"
    );

    if(!regex_search(chemfilestring_, fileStructure))
    {
        throw regex_error("chem.inp needs the structure:\n ELEMENTS END\n SPECIES END\n REACTIONS END.");
    }

    smatch result;
    if(regex_search(chemfilestring_, result, unsupported))
    {
        if (result[1] == "TDEP") throw regex_error("TDEP not supported yet.");
        if (result[2] == "EXCI") throw regex_error("EXCI not supported yet.");
        if (result[3] == "JAN") throw regex_error("JAN not supported yet.");
        if (result[4] == "FIT1") throw regex_error("FIT1 not supported yet.");
        if (result[5] == "HV") throw regex_error("HV not supported yet.");
        if (result[6] == "MOME") throw regex_error("MOME not supported yet.");
        if (result[7] == "RORD") throw regex_error("RORD not supported yet.");
        if (result[8] == "UNITS") throw regex_error("UNITS not supported yet.");
        if (result[9] == "HIGH") throw regex_error("HIGH not supported yet.");
        if (result[10] == "USER") throw regex_error("USER not supported yet.");
        if (result[11] == "LT") throw regex_error("LT not supported yet.");
        if (result[12] == "RLT") throw regex_error("RLT not supported yet.");
    }

	cout << "chem.inp file format check PASSED." << endl;
	
    return true;

}

bool IO::ChemkinReader::checkChemSurfFile()
{
    cout << "Checking the format of the chem.inp file." << endl;

    const regex fileStructure
    (
        "ELEM(?:|ENT|ENTS).*?END.*?"
        "SPEC(?:|IE|IES).*?END.*?"
        "REAC(?:|TION|TIONS).*?END"
    );

     const regex surfacefileStructure
    (
        "MATERIAL.*?(?:|END).*?"
        "SITE(?:|S).*?END.*?"
        "REAC(?:|TION|TIONS).*?END"
    );

    //! \todo NEED TO ADD REGEXES FOR (LT)|(RLT) IN REACTIONPARSER ASAP
    const regex unsupported
    (
        "(TDEP)|(EXCI)|(JAN)|(FIT1)|"
        "(HV)|(MOME)|(FORD)|(RORD)|(UNITS)|(HIGH)|(USER)|"
        "(LT)|(RLT)"
    );

    if(!regex_search(chemfilestring_, fileStructure))
    {
        throw regex_error("chem.inp needs the structure:\n ELEMENTS END\n SPECIES END\n REACTIONS END.");
    }

    smatch result;
    if(regex_search(chemfilestring_, result, unsupported))
    {
        if (result[1] == "TDEP") throw regex_error("TDEP not supported yet.");
        if (result[2] == "EXCI") throw regex_error("EXCI not supported yet.");
        if (result[3] == "JAN") throw regex_error("JAN not supported yet.");
        if (result[4] == "FIT1") throw regex_error("FIT1 not supported yet.");
        if (result[5] == "HV") throw regex_error("HV not supported yet.");
        if (result[6] == "MOME") throw regex_error("MOME not supported yet.");
        if (result[7] == "RORD") throw regex_error("RORD not supported yet.");
        if (result[8] == "UNITS") throw regex_error("UNITS not supported yet.");
        if (result[9] == "HIGH") throw regex_error("HIGH not supported yet.");
        if (result[10] == "USER") throw regex_error("USER not supported yet.");
        if (result[11] == "LT") throw regex_error("LT not supported yet.");
        if (result[12] == "RLT") throw regex_error("RLT not supported yet.");
		
    }

	/*
	* Added by mm864 for chemsurf.inp file 
	*/
	if(!regex_search(chemSurfFilestring_, surfacefileStructure))
    {
        throw regex_error("chemsurf.inp needs the structure:\n MATERIAL END\n SITE END\n REACTIONS END.");
    }

    smatch Surfresult;
    if(regex_search(chemSurfFilestring_, Surfresult, unsupported))
    { // You need to modify later (leave it now)
        if (result[1] == "TDEP") throw regex_error("TDEP not supported yet.");
        if (result[2] == "EXCI") throw regex_error("EXCI not supported yet.");
        if (result[3] == "JAN") throw regex_error("JAN not supported yet.");
        if (result[4] == "FIT1") throw regex_error("FIT1 not supported yet.");
        if (result[5] == "HV") throw regex_error("HV not supported yet.");
        if (result[6] == "MOME") throw regex_error("MOME not supported yet.");
        if (result[7] == "RORD") throw regex_error("RORD not supported yet.");
        if (result[8] == "UNITS") throw regex_error("UNITS not supported yet.");
        if (result[9] == "HIGH") throw regex_error("HIGH not supported yet.");
        if (result[10] == "USER") throw regex_error("USER not supported yet.");
        if (result[11] == "LT") throw regex_error("LT not supported yet.");
        if (result[12] == "RLT") throw regex_error("RLT not supported yet.");
		if (result[13] == "DCOL") throw regex_error("DCOL not supported yet.");
		if (result[14] == "LANG") throw regex_error("LANG not supported yet.");
		if (result[15] == "LHDE") throw regex_error("LHDE not supported yet.");
		if (result[16] == "LHNU") throw regex_error("LHNU not supported yet.");
		if (result[17] == "YIELD") throw regex_error("YIELD not supported yet.");
		if (result[18] == "USRPROG") throw regex_error("USRPROG not supported yet.");
		if (result[19] == "NUCL") throw regex_error("NUCL not supported yet.");
		if (result[20] == "LHPR") throw regex_error("LHPR not supported yet.");
		if (result[21] == "BOHM") throw regex_error("BOHM not supported yet.");
		if (result[22] == "ENRGDEP") throw regex_error("ENRGDEP not supported yet.");
    }

    cout << "chem.inp and chemsurf.inp file format check PASSED." << endl;

    return true;

}

void IO::ChemkinReader::check()
{
    cout << "Chemistry file: " << chemfile_ << endl;
	cout << "Surface Chemistry file: " << chemSurfFile_ << endl;
    cout << "Thermo file: " << thermfile_ << endl;
	cout << "Surface Thermo file: " << thermSurfFile_ << endl;
    cout << "Trans file: " << transfile_ << endl;

    if (globalUnits_ != "NO GLOBAL UNITS" && globalUnits_ != "CAL/MOLE")
    {
        throw std::logic_error(globalUnits_+" are not supported yet.");
    }
    cout << "Global Units are " << globalUnits_ << endl;
	
	ofstream outputReactions("reactionsParsed");
	outputReactions << reactions_ << endl;
    ofstream outputSpecies("speciesParsed");
    outputSpecies << elements_ << endl;
	outputSpecies << phase_ << endl;
    outputSpecies << species_ << endl;
    
    

    cout << "Data output to speciesParsed and reactionsParsed." << endl;
}

void IO::ChemkinReader::read()
{

    readGlobalUnits();
    if (chemSurfFile_ != "NOT READ"){
	readPhase();
    }
    readElements();
    readSpecies();
    
	
     if (transfile_ != "NOT READ" )
    {
        TransportParser transportParser(transfile_);
        transportParser.parse(species_);
    }
    
     // Moved by mm864 to check only gas species, seems that the transport data for surface species is absent 
     if (chemSurfFile_ != "NOT READ"){
       readSurfaceSpecies();
     }
    /*
    ThermoParser thermoParser(thermfile_ );
    thermoParser.parse(species_);

    ThermoParser thermoParser(thermSurfFile_);
    thermoSurfParser.parse(species_);
    */

    ThermoParser thermoParser(thermfile_, thermSurfFile_);
    thermoParser.parse(species_);
    readReactions();
    /* Read surface units must be placed after read gas reaction otherwise if the surface reaction unit is not default,
     *  everything in gas reaction will be converted as well. 
     */ 

    if (chemSurfFile_ != "NOT READ"){
    readSurfaceUnits();	
    readSurfReactions();
    }
}

/*
void IO::ChemkinReader::read()
{

    readGlobalUnits();
	readPhase();
    readElements();
    readSpecies();
	readSurfaceSpecies();
	
	
    ThermoParser thermoParser(thermfile_);
    thermoParser.parse(species_);
	ThermoParser thermoSurfParser(thermSurfFile_);
	thermoSurfParser.parse(species_);

    if (transfile_ != "NOT READ")
    {
        TransportParser transportParser(transfile_);
        transportParser.parse(species_);
    }

    readReactions();

}
*/

// Not that the surface element would be defined in chem.inp
void IO::ChemkinReader::readElements() {

    smatch result;
    regex_search(chemfilestring_, result, elementListRegex);
    string elementString = result[1];

    string::const_iterator start = elementString.begin();
    string::const_iterator end = elementString.end();
    match_results<string::const_iterator> what;

    while (regex_search(start, end, what, elementSingleRegex)) {
        elements_.push_back(Element(what[1], 0.0));
        start = what[0].second;
    }
	cout << "elements read" << endl; 
}

void IO::ChemkinReader::readSpecies() {

    smatch result;
    regex_search(chemfilestring_, result, speciesListRegex);
    string speciesString = result[1];
    phase_.push_back(Phase("gas","", 0.0)); // stores the gas phase Phase object

    string::const_iterator start = speciesString.begin();
    string::const_iterator end = speciesString.end();
    match_results<string::const_iterator> what1; 
    int sp_index =0;
    std::map<std::string, int> sp_map; 
    while (regex_search(start, end, what1, speciesSingleRegex)){
      species_.push_back(Species(what1[1]));
      species_.back().setSiteOccupancy(0);
      species_.back().setPhaseName("gas");
      sp_map.insert(pair<string, int> (what1[1], sp_index));
      start = what1[0].second;
      sp_index = sp_index + 1; 

      // cout << what1[1] << endl; // for checking
    }	
    phase_.back().setSpecies(sp_map);
       	
    cout << "species read" << endl; 
}


void IO::ChemkinReader::readSurfaceSpecies() {
  
  smatch result;
  regex_search(chemSurfFilestring_, result, mainSurfaceListRegex);
  string siteList = result[2];
  //cout << siteList << endl; // for checking
    
  string::const_iterator start0 = siteList.begin();
  string::const_iterator end0 = siteList.end();
  match_results<string::const_iterator> what0;

  while (regex_search(start0, end0, what0, surfaceSpeciesListRegex)){;
    string spList = what0[3];
    string phName = what0[1];
    //cout << phName << endl; // for checking
    //cout << spList << endl;

    string::const_iterator start = spList.begin();
    string::const_iterator end = spList.end();
    match_results<string::const_iterator> what;

    
    while (regex_search(start, end, what, surfaceSpeciesSingleRegex)){

      string spList1 = what[1];
      size_t found = spList1.find('/'); 

      if (found!=string::npos){
	spList1 = spList1; 
      }
      else {
	spList1 = spList1.append("/1/");
      }


      string::const_iterator start1 = spList1.begin();
      string::const_iterator end1 = spList1.end();
      match_results<string::const_iterator> what1; 
		
    
      while (regex_search(start1, end1, what1, surfaceSpeciesSingleSubRegex)){
		string surfaceSpeciesString = what1[1];
		string siteO = what1[2];
		species_.push_back(Species(surfaceSpeciesString));
		species_.back().setPhaseName(phName);
      
		if (siteO.length() > 0)
		{
			stringstream ss(siteO);
			int sOccup;
			ss >> sOccup;
			species_.back().setSiteOccupancy(sOccup);
			cout << sOccup << endl; // for checking
		}
		else {species_.back().setSiteOccupancy(1);}
		start1 = what1[0].second;

		cout << surfaceSpeciesString << endl; // for checking
		
		}
		start = what[0].second;

	}	

    start0 = what0[0].second;  
  }
  cout << "surface species read" << endl; 
  
}


void IO::ChemkinReader::readPhase() {

  smatch result;
	
  regex_search(chemSurfFilestring_, result, mainSurfaceListRegex);
  string siteList = result[2];
  //cout << siteList << endl;
    
  string::const_iterator start0 = siteList.begin();
  string::const_iterator end0 = siteList.end();
  match_results<string::const_iterator> what0;
 
  while (regex_search(start0, end0, what0, surfaceSpeciesListRegex)){;
    string spList = what0[3];
    string siteName = what0[1];
    string sitedenString = what0[2];
    stringstream ss(sitedenString);
    double siteD;
    ss >> siteD;
    phase_.push_back(Phase(siteName,"", siteD));
    // for checking
    //cout << spList << endl;
    //cout << siteName << endl;
    //cout << siteD << endl;
   
    string::const_iterator start = spList.begin();
    string::const_iterator end = spList.end();
    match_results<string::const_iterator> what;

    int sp_index =0;
	std::map<std::string, int> sp_map;
    while (regex_search(start, end, what, surfaceSpeciesSingleRegex)){

      string spList1 = what[1];
      size_t found = spList1.find('/'); 

      if (found!=string::npos){

	spList1 = spList1; 

      }

      else {

	spList1 = spList1.append("/1/");
      }
      
		string::const_iterator start1 = spList1.begin();
		string::const_iterator end1 = spList1.end();
		match_results<string::const_iterator> what1; 
		
		
		while (regex_search(start1, end1, what1, surfaceSpeciesSingleSubRegex)){
		  string surfaceSpeciesString = what1[1];
		  //string siteOc = what1[2]; (Not in use in phase)
		  sp_map.insert(pair<string, int> (surfaceSpeciesString, sp_index));
		  start1 = what1[0].second;	
		  sp_index = sp_index + 1;
		  
		  //cout << "*" << surfaceSpeciesString << "*" << endl;	// for checking
		  //cout << "*" << siteOc << "*" << endl; // for checking
		  
    		}
		
		start = what[0].second;
    }		
	phase_.back().setSpecies(sp_map);
    start0 = what0[0].second;
   
  }
 cout << "phase read" << endl; 
}
 



void IO::ChemkinReader::readReactions() {

    smatch result;
    regex_search(chemfilestring_, result, reactionListRegex);
    string reactionString = result[1];

    ReactionParser reactionParser(reactionString);
    reactionParser.parse(reactions_);

	cout << "reaction read, size all: " << reactions_.size() << endl; 
}

void IO::ChemkinReader::readSurfReactions() {

	smatch result;
    regex_search(chemSurfFilestring_, result, reactionListRegex);
    string reactionString = result[1];

    ReactionParser reactionParser(reactionString);
    reactionParser.setSurfaceReactionUnit(surfUnits_);
    reactionParser.parse(reactions_);
	
	
	cout << "surface reaction read, size all: " << reactions_.size() << endl; 
}


void IO::ChemkinReader::readGlobalUnits()
{

    smatch units;
    string::const_iterator start = chemfilestring_.begin();
    string::const_iterator end = chemfilestring_.end();

    while (regex_search(start, end, units, unitsRegex))
    {
        if (globalUnits_ != "NO GLOBAL UNITS")
            throw std::logic_error("Units are already specified as " + globalUnits_);
       
        globalUnits_ = units[1];
        start = units[0].second;
    }

	cout << "Global Units Read" << endl; 
	cout << "Units for gas reaction " << globalUnits_ << endl;
}

void IO::ChemkinReader::readSurfaceUnits()
{	

    	smatch surf_units;
	regex_search(chemSurfFilestring_, surf_units, surfaceUnitsRegex);
      
        surfUnits_ = surf_units[2];

	if (surfUnits_.length() == 0)
	{ surfUnits_ = "NO GLOBAL UNITS"; 
	}
	cout << "Surface Units Read" << endl; 
	cout << "Units for surface reaction " << surfUnits_ << endl;

}
