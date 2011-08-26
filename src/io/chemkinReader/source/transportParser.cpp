/*
 * transportParser.cpp
 *
 *  Created on: Jun 23, 2011
 *      Author: gigadot
 *     License: Apache 2.0
 */

#include "transportParser.h"
#include "stringFunctions.h"
#include "species.h"

using namespace std;
using namespace boost;

const string IO::TransportParser::transportRegex
(
    "\\s+([0-2]+?)\\s+"
    "((?:[0-9]|\\.)+(?:\\.[0-9]*)?)\\s+"
    "((?:[0-9]|\\.)+(?:\\.[0-9]*)?)\\s+"
    "((?:[0-9]|\\.)+(?:\\.[0-9]*)?)\\s+"
    "((?:[0-9]|\\.)+(?:\\.[0-9]*)?)\\s+"
    "((?:[0-9]|\\.)+(?:\\.[0-9]*)?)"
);

// Empty default constructor, can be removed but leave it there just in case.
IO::TransportParser::TransportParser
(
    const string tranfile
)
:
    tranfile_(tranfile),
    transportfilestring_(convertToCaps(fileToString(tranfile)))
{}

void IO::TransportParser::parse(vector<Species>& species)
{

    cout << "Transport file: " << tranfile_ << endl;

    for (size_t i = 0; i != species.size(); ++i)
    {
        smatch specieTransportData = findSpecies(species[i]);
        setSpecieData(species[i], specieTransportData);
    }
    cout << "End of Transport file: " << tranfile_ << endl;
}


const smatch IO::TransportParser::findSpecies(const IO::Species& specie)
{

    smatch what;
    regex reg("\\b"+IO::regex_escape(specie.name())+transportRegex);

    string::const_iterator start = transportfilestring_.begin();
    string::const_iterator end = transportfilestring_.end();

    if(!regex_search(start, end, what, reg))
    {
        throw regex_error
        (
            "Species " + specie.name() + " not found in tran.dat. using "
            "\\b"+IO::regex_escape(specie.name())+transportRegex
        );
    }

    return what;

}

void IO::TransportParser::setSpecieData
(
    IO::Species& specie,
    const smatch& specieTransportData
)
{

    specie.transport().setMoleculeIndex(from_string<int>(specieTransportData[1]));
    specie.transport().setPotentialWellDepth(from_string<double>(specieTransportData[2]));
    specie.transport().setCollisionDiameter(from_string<double>(specieTransportData[3]));
    specie.transport().setDipoleMoment(from_string<double>(specieTransportData[4]));
    specie.transport().setPolarizability(from_string<double>(specieTransportData[5]));
    specie.transport().setRotRelaxationNumber(from_string<double>(specieTransportData[6]));

}

