/*
* reactionParser.cpp
*
*  Created on: Jun 23, 2011
*      Author: gigadot
*     License: Apache 2.0
*/

#include <boost/algorithm/string.hpp>
#include "boost/algorithm/string/trim.hpp"

#include "reactionParser.h"
#include "stringFunctions.h"
#include "reaction.h"

using namespace std;
using namespace boost;

const regex IO::ReactionParser::reactionSingleRegex
    (
    "(.*?)\\s*"
    "(<=>|=>|=)\\s*"
    "(.*?)"
    "\\s+((?:[0-9]+|\\.)\\.*[0-9]*(?:[eEgG][-+]?[0-9]*)*)"
    "\\s+(.*?)"
    "\\s+(.*?)$|\\n"
    );

const regex IO::ReactionParser::blankLine
    (
    "\\s*\\n*$"
    );

const regex IO::ReactionParser::DUPLICATE
    (
    "DUPLICATE|DUP"
    );

const regex IO::ReactionParser::LOW
    (
    "(LOW)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)\\s*\\/"
    );

const regex IO::ReactionParser::TROE
    (
    "(TROE)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)(?:|\\s+(.*?))\\s*\\/"
    );

const regex IO::ReactionParser::SRI
    (
    "(SRI)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)(?:|\\s+(.*?)\\s+(.*?))\\s*\\/"
    );

const regex IO::ReactionParser::REV
    (
    "(REV)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)\\s*\\/"
    );

const regex IO::ReactionParser::pressureDependent
    (
    "\\(\\+(.*?)\\)"
    );

const regex IO::ReactionParser::EaUnitsRegex
    ("REAC(?:|TION|TIONS)\\s+.*\\b(CAL/MOLE|KCAL/MOLE|JOULES/MOLE|KJOULES/MOLE|KJOU/MOL|KJOU/MOLE|KELVINS|EVOLTS)\\b");

const regex IO::ReactionParser::AUnitsRegex
    ("REAC(?:|TION|TIONS)\\s+.*\\b(MOLES|MOLECULES)\\b");


// Empty default constructor, can be removed but leave it there just in case.
IO::ReactionParser::ReactionParser
    (
    const string reactionString
    )
    :
reactionString_(reactionString),
    globalEaUnits_("NO GLOBAL UNITS"),
    globalAUnits_("NO GLOBAL UNITS"),
    scale_Ea(1.0),
    scale_A(1.0)
{
    split(reactionStringLines_, reactionString, boost::is_any_of("\n|\r\n"));

    // Check for a blank line and erase.
    for (size_t i=0; i<reactionStringLines_.size(); ++i)
    {
        if(isBlankLine(reactionStringLines_[i]))
        {
            reactionStringLines_.erase(reactionStringLines_.begin()+i);
            --i;
        }
    }
}

void IO::ReactionParser::parse(vector<IO::Reaction>& reactions)
{

    readGlobalUnits();

    for (size_t i=0; i<reactionStringLines_.size(); ++i)
    {

        cout << reactionStringLines_[i] << endl;

        Reaction reaction;

        // Check for pressure dependency now as it screws up reactionSingleRegex.
        if (checkForPressureDependentReaction(reactionStringLines_[i]))
        {
            reactionStringLines_[i] = regex_replace(reactionStringLines_[i], pressureDependent, "");
            reaction.setPressureDependent();
        }

        smatch what;
        string::const_iterator start = reactionStringLines_[i].begin();
        string::const_iterator end = reactionStringLines_[i].end();

        if (regex_search(start, end, what, reactionSingleRegex))
        {
            reaction.setReactants(parseReactionSpecies(what[1]));

            if (what[2] == "=>") reaction.setReversible(false);

            reaction.setProducts(parseReactionSpecies(what[3]));

            // Do conversion to SI units here. Assume for the moment that it
            // is not a 3rd body reaction, we will deail with that case later.
            reaction.setArrhenius
                (
                from_string<double>(what[4]) * pow(scale_A, reaction.getReactantStoich() - 1.0),
                from_string<double>(what[5]),
                from_string<double>(what[6]) * scale_Ea
                );

            while (i < reactionStringLines_.size()-1)
            {

                cout << "Next = " << reactionStringLines_[i+1] << endl;

                start = reactionStringLines_[i+1].begin();
                end = reactionStringLines_[i+1].end();

                if (regex_search(start, end, reactionSingleRegex))
                {
                    break;
                }
                else if (regex_search(start, end, DUPLICATE))
                {
                    reaction.setDuplicate();
                    // Skip one line when looking for the next reaction.
                    ++i;
                    //break;
                }
                else if (regex_search(start, end, REV))
                {
                    vector<double> reverseArrhenius = parseLOWTROEREV(reactionStringLines_[i+1], REV);
                    // Do conversion to SI units here. Assume for the moment that it
                    // is not a 3rd body reaction, we will deail with that case later.
                    reaction.setArrhenius(reverseArrhenius[0] * pow(scale_A, reaction.getReactantStoich() - 1.0),reverseArrhenius[1],reverseArrhenius[2] * scale_Ea,true);
                    // Skip one line when looking for the next reaction.
                    ++i;
                    //break;
                }
                else if (reaction.hasThirdBody() || reaction.isPressureDependent())
                {
                    // It could be a third body or Lindemann reaction so need to scale the forward
                    IO::Arrhenius arr( reaction.getArrhenius() );
                    reaction.setArrhenius( arr.A * scale_A, arr.n, arr.E );
                    // and reverse A's, if present, a bit more.
                    if( reaction.hasREV() )
                    {
                        arr = reaction.getArrhenius(true);
                        reaction.setArrhenius( arr.A * scale_A, arr.n, arr.E, true );
                    }

                    string lineType = findLineType(reactionStringLines_[i+1]);
                    if (lineType == "THIRDBODY")
                    {
                        reaction.setThirdBodies(parseThirdBodySpecies(reactionStringLines_[i+1]));
                        ++i;
                    } else if (lineType == "LOW")
                    {
                        vector<double> lowParams = parseLOWTROEREV(reactionStringLines_[i+1], LOW);
                        // Convert to SI units.
                        lowParams[0] *= pow(scale_A, reaction.getReactantStoich());
                        lowParams[2] *= scale_Ea;
                        reaction.setLOW(lowParams);

                        ++i;

                    }else if (lineType == "TROE")
                    {
                        reaction.setTROE(parseLOWTROEREV(reactionStringLines_[i+1], TROE));

                        // It's not a Lindemann reaction so need to scale the forward
                        IO::Arrhenius arr( reaction.getArrhenius() );
                        reaction.setArrhenius( arr.A / scale_A, arr.n, arr.E );
                        // and reverse A's, if present, back again.
                        if( reaction.hasREV() )
                        {
                            arr = reaction.getArrhenius(true);
                            reaction.setArrhenius( arr.A / scale_A, arr.n, arr.E, true );
                        }
                        ++i;
                    } else if (lineType == "SRI")
                    {
                        reaction.setSRI(parseLOWTROEREV(reactionStringLines_[i+1], SRI));

                        // It's not a Lindemann reaction so need to scale the forward
                        IO::Arrhenius arr( reaction.getArrhenius() );
                        reaction.setArrhenius( arr.A / scale_A, arr.n, arr.E );
                        // and reverse A's, if present, back again.
                        if( reaction.hasREV() )
                        {
                            arr = reaction.getArrhenius(true);
                            reaction.setArrhenius( arr.A / scale_A, arr.n, arr.E, true );
                        }
                        ++i;
                    }
                }
                else
                {
                    throw std::logic_error("Reaction "+reactionStringLines_[i+1]+" is not supported.");
                }

            }

            reactions.push_back(reaction);

        }

    }

}

multimap<string, double>
    IO::ReactionParser::parseReactionSpecies(string reactionSpecies)
{

    std::multimap<std::string, double> reactionSpeciesMap;

    regex splitSpecies("\\+");
    regex splitStoichiometry("([0-9]*)([A-Z].*)");

    sregex_token_iterator i(reactionSpecies.begin(), reactionSpecies.end(), splitSpecies,-1);
    sregex_token_iterator j;

    while (i != j)
    {
        // *i Gives a reactant species. Now get its stoichiometry.
        smatch splitStoic;
        // ++ iterates to the next reactant!
        string reactionSpecie = *i++;

        string::const_iterator start = reactionSpecie.begin();
        string::const_iterator end = reactionSpecie.end();

        regex_search(start, end, splitStoic, splitStoichiometry);

        string speciesName = splitStoic[2];

        if (splitStoic[1] == "")
        {
            reactionSpeciesMap.insert
                (
                pair<string,double>
                (
                trim_copy(speciesName),
                1.0
                )
                );
        }
        else
        {
            reactionSpeciesMap.insert
                (
                pair<string,double>
                (
                trim_copy(speciesName),
                from_string<double>(splitStoic[1])
                )
                );
        }

    }

    return reactionSpeciesMap;

}

multimap<string, double>
    IO::ReactionParser::parseThirdBodySpecies(const string& thirdBodies)
{

    string trim_copymed = trim_copy(thirdBodies);
    multimap<string, double> thirdBodyMap;
    // Split the next line using / as delimiter.
    regex splitThirdBodies("\\/");

    sregex_token_iterator j
        (
        trim_copymed.begin(),
        trim_copymed.end(),
        splitThirdBodies,
        -1
        );
    sregex_token_iterator k;

    while (j != k)
    {
        string name = *j++;
        string trim_copyName = trim_copy(name);
        double efficiencyFactor = from_string<double>(*j++);
        thirdBodyMap.insert(pair<string,double>(trim_copyName,efficiencyFactor));
    }

    return thirdBodyMap;

}

bool
    IO::ReactionParser::isBlankLine(const string& line)
{

    string::const_iterator start = line.begin();
    string::const_iterator end = line.end();

    return regex_match(start, end, blankLine);

}

string
    IO::ReactionParser::findLineType(const string& line)
{

    string::const_iterator start = line.begin();
    string::const_iterator end = line.end();

    if (regex_search(start, end, LOW))
        return "LOW";
    if (regex_search(start, end, TROE))
        return "TROE";
    if (regex_search(start, end, SRI))
        return "SRI";

    return "THIRDBODY";

}

vector<double>
    IO::ReactionParser::parseLOWTROEREV(const string& line, const regex& reg)
{

    vector<double> vec;
    smatch result;

    regex_search(line.begin(), line.end(), result, reg);

    for (size_t i=2; i<result.size(); ++i)
    {
        if (result[i] != "") vec.push_back(from_string<double>(result[i]));
    }

    return vec;

}

bool
    IO::ReactionParser::checkForPressureDependentReaction(const string& line)
{
    if (!regex_search(line.begin(), line.end(), pressureDependent))
    {
        return false;
    }
    return true;
}

void IO::ReactionParser::readGlobalUnits()
{
    smatch units;
    string::const_iterator start = reactionString_.begin();
    string::const_iterator end = reactionString_.end();

    while (regex_search(start, end, units, EaUnitsRegex))
    {
        if (globalEaUnits_ != "NO GLOBAL UNITS")
            throw std::logic_error("Units are already specified as " + globalEaUnits_);
        cout << units[1]<<endl;
        globalEaUnits_ = units[1];
        start = units[0].second;

        if ("CAL/MOLE" == globalEaUnits_) {
            scale_Ea = 4.184e7;
        } else if ("KCAL/MOLE" == globalEaUnits_) {
            scale_Ea = 4.184e10;
        } else if ("JOULES/MOLE" == globalEaUnits_) {
            scale_Ea = 1.0e7;
        } else if ("KJOULES/MOLE" == globalEaUnits_ || "KJOU/MOLE" == globalEaUnits_ || "KJOU/MOL" == globalEaUnits_) {
            scale_Ea = 1.0e10;   
        } else if ("KELVINS" == globalEaUnits_) {
            scale_Ea = 8.3144621e7; // R in cgs
        } else if ("EVOLTS" == globalEaUnits_) {
            scale_Ea = 1.60217646e-12;
        } else {
            throw std::logic_error(globalEaUnits_+" are not supported.");
        }
    }

    if ("NO GLOBAL UNITS" == globalEaUnits_)
    {
        globalEaUnits_ = "CAL/MOLE";
        scale_Ea = 4.184e7;
    }

    scale_Ea *= 1e-7;

    cout << "Global Units for Ea are " << globalEaUnits_ << endl;

    while (regex_search(start, end, units, AUnitsRegex))
    {
        if (globalAUnits_ != "NO GLOBAL UNITS")
            throw std::logic_error("Units are already specified as " + globalAUnits_);
        cout << units[1]<<endl;
        globalAUnits_ = units[1];
        start = units[0].second;
        if ("MOLES" == globalAUnits_) {
            scale_A = 1.0e-6;
        } else if ("MOLECULES" == globalAUnits_) {
            scale_A = 1.0e-6 / 6.0221415e23;
        } else {
            throw std::logic_error(globalAUnits_+" are not supported.");
        }
    }

    if ("NO GLOBAL UNITS" == globalAUnits_)
    {
        globalAUnits_ = "MOLES";
        scale_A = 1.0e-6;
    }

    cout << "Global Units for A are " << globalAUnits_ << endl;

}