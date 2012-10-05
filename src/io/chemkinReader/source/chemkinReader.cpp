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

using namespace std;
using namespace boost;

const regex IO::ChemkinReader::elementListRegex
("ELEM(?:|ENT|ENTS)\\s+(.*?)\\s+END");

const regex IO::ChemkinReader::elementSingleRegex
("(\\w+)");

const regex IO::ChemkinReader::speciesListRegex
("SPEC(?:|IE|IES)\\s+(.*?)\\s+END");

const regex IO::ChemkinReader::speciesSingleRegex
("\\s+");

const regex IO::ChemkinReader::reactionListRegex
("REAC(?:|TION|TIONS)\\s+(.*?)\\s+END");

const regex IO::ChemkinReader::unitsRegex
("REAC(?:|TION|TIONS)\\s+\\b(CAL/MOLE|KCAL/MOLE|JOULES/MOLE|KJOULES/MOLE|KJOU/MOL|KJOU/MOLE|KELVINS|EVOLTS|MOLES|MOLECULES)\\b");

IO::ChemkinReader::ChemkinReader
(
    const string chemfile,
    const string thermfile,
    const string transfile
)
:
    chemfile_(chemfile),
    thermfile_(thermfile),
    transfile_(transfile),
    chemfilestring_(convertToCaps(replaceComments(fileToString(chemfile_)))),
    globalUnits_("NO GLOBAL UNITS")
{
    checkChemFile();
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
        if (result[7] == "FORD") throw regex_error("FORD not supported yet.");
        if (result[8] == "RORD") throw regex_error("RORD not supported yet.");
        if (result[9] == "UNITS") throw regex_error("UNITS not supported yet.");
        if (result[10] == "HIGH") throw regex_error("HIGH not supported yet.");
        if (result[11] == "USER") throw regex_error("USER not supported yet.");
        if (result[12] == "LT") throw regex_error("LT not supported yet.");
        if (result[13] == "RLT") throw regex_error("RLT not supported yet.");
    }

    cout << "chem.inp file format check PASSED." << endl;

    return true;

}

void IO::ChemkinReader::check()
{
    cout << "Chemistry file: " << chemfile_ << endl;
    cout << "Thermo file: " << thermfile_ << endl;
    cout << "Trans file: " << transfile_ << endl;

    if (globalUnits_ != "NO GLOBAL UNITS" && globalUnits_ != "CAL/MOLE")
    {
        throw std::logic_error(globalUnits_+" are not supported yet.");
    }
    cout << "Global Units are " << globalUnits_ << endl;

    ofstream outputSpecies("speciesParsed");
    ofstream outputReactions("reactionsParsed");
    outputSpecies << elements_ << endl;
    outputSpecies << species_ << endl;
    outputReactions << reactions_ << endl;

    cout << "Data output to speciesParsed and reactionsParsed." << endl;
}

void IO::ChemkinReader::read()
{

    readGlobalUnits();

    readElements();
    readSpecies();

    ThermoParser thermoParser(thermfile_);
    thermoParser.parse(species_);

    if (transfile_ != "NOT READ")
    {
        TransportParser transportParser(transfile_);
        transportParser.parse(species_);
    }

    readReactions();

}

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

}

void IO::ChemkinReader::readSpecies() {

    smatch result;
    regex_search(chemfilestring_, result, speciesListRegex);
    string speciesString = result[1];

    string::const_iterator start = speciesString.begin();
    string::const_iterator end = speciesString.end();

    sregex_token_iterator i(start, end, speciesSingleRegex, -1);
    sregex_token_iterator j;
    while (i != j) {
        species_.push_back(Species(*i++));
        //species_.back().checkElementsInSpecies(elements_);
    }

}

void IO::ChemkinReader::readReactions() {

    smatch result;
    regex_search(chemfilestring_, result, reactionListRegex);
    string reactionString = result[1];

    ReactionParser reactionParser(reactionString);
    reactionParser.parse(reactions_);

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
        cout << units[1]<<endl;
        globalUnits_ = units[1];
        start = units[0].second;
    }

}
