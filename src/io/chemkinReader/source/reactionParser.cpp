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

const regex IO::ReactionParser::STICK // Checked by mm864!
(
    "(STICK)"
);

const regex IO::ReactionParser::COV // Checked by mm864!
(
    "(COV)\\s*\\/\\s*(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s*\\/"
);

const regex IO::ReactionParser::FORD // Checked by mm864!
(
    "(FORD)\\s*\\/\\s*(.*?)\\s+(.*?)\\s*\\/"
);

const regex IO::ReactionParser::MWON// Checked by mm864!
(
    "(MWON)"
);

const regex IO::ReactionParser::MWOFF // Checked by mm864!
(
    "(MWOFF)"
);


// Empty default constructor, can be removed but leave it there just in case.
IO::ReactionParser::ReactionParser
(
    const string reactionString
)
:
    reactionString_(reactionString)
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

void IO::ReactionParser::setSurfaceReactionUnit(const std::string &units){

surfaceUnits = units; 
}

void IO::ReactionParser::parse(vector<IO::Reaction>& reactions)
{

    for (size_t i=0; i<reactionStringLines_.size(); ++i)
    {

        // comment out printing reaction mechanism 
        // cout << reactionStringLines_[i] << endl;

        Reaction reaction;

        // Check for pressure dependency now as it screws up reactionSingleRegex. It is only for gas phase reaction (added comment by mm864)
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
		
	    	if ((surfaceUnits == "KJOULES/MOL") || (surfaceUnits == "KJOULES/MOLE")){	
            	reaction.setArrhenius
           	 (
                from_string<double>(what[4]),
                from_string<double>(what[5]),
                from_string<double>(what[6])*239.005736 // KJOULES/mol to cal/mol
            	);
			}
			
			else if ((surfaceUnits == "JOULES/MOL") || (surfaceUnits == "JOULES/MOLE")){	
            	reaction.setArrhenius
           	 (
                from_string<double>(what[4]),
                from_string<double>(what[5]),
                from_string<double>(what[6])*239.005736/1000 // JOULES/mol to cal/mol
            	);
			}

		else {
	    	reaction.setArrhenius
            	(
                from_string<double>(what[4]),
                from_string<double>(what[5]),
                from_string<double>(what[6])
            	);


		}


            while (i < reactionStringLines_.size()-1)
            {
                // Comment out hard coded print to terminal of reaction mechanism:
                // cout << "Next = " << reactionStringLines_[i+1] << endl;

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
                    reaction.setArrhenius(reverseArrhenius[0],reverseArrhenius[1],reverseArrhenius[2],true);
                    // Skip one line when looking for the next reaction.
                    ++i;
                    //break;
                }
                else if (reaction.hasThirdBody() || reaction.isPressureDependent())
                {
                    string lineType = findLineType(reactionStringLines_[i+1]);
                    if (lineType == "THIRDBODY")
                    {
                        reaction.setThirdBodies(parseThirdBodySpecies(reactionStringLines_[i+1]));
                        ++i;
                    }
                    if (lineType == "LOW")
                    {
                        reaction.setLOW(parseLOWTROEREV(reactionStringLines_[i+1], LOW));
                        ++i;
                    }
                    if (lineType == "TROE")
                    {
                        reaction.setTROE(parseLOWTROEREV(reactionStringLines_[i+1], TROE));
                        ++i;
                    }
                    if (lineType == "SRI")
                    {
                        reaction.setSRI(parseLOWTROEREV(reactionStringLines_[i+1], SRI));
                        ++i;
                    }
                }
		else if (regex_search(start, end, COV))
		{
		   smatch result;
		   string::const_iterator start = reactionStringLines_[i+1].begin();
		   string::const_iterator end = reactionStringLines_[i+1].end();
		   regex_search(start, end, result, COV);

		   std::string speciesName =  result[2];
		   double eta = from_string<double>(result[3]);
		   double miu = from_string<double>(result[4]);
		   double epsilon = from_string<double>(result[5]);
		   
			if ((surfaceUnits == "KJOULES/MOL") || (surfaceUnits == "KJOULES/MOLE")) {	
			epsilon = from_string<double>(result[5])*1000; // Convert to Joules 
			}
		   
		   reaction.setCOV(eta, miu, epsilon, speciesName);  
		   // Skip one line when looking for the next reaction.
		   ++i;
		}
		else if (regex_search(start, end, FORD))
		  {
		    smatch result;
		    string::const_iterator start = reactionStringLines_[i+1].begin();
			string::const_iterator end = reactionStringLines_[i+1].end(); 
		    regex_search(start, end, result, FORD);

		    std::string speciesName =  result[2];
		    double coeff = from_string<double>(result[3]);

		    reaction.setFORD(coeff, speciesName);  

		    // Skip one line when looking for the next reaction.
		    ++i;
		  }
				
		else if (regex_search(start, end, STICK))
		  {
		  reaction.setSTICK();

                 // Skip one line when looking for the next reaction.
		  ++i;
                 //break;
		  }
				
		else if (regex_search(start, end, MWON))
		  {
		    reaction.setMWON();
                    // Skip one line when looking for the next reaction.
                    ++i;
                    //break;
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
                    from_string<double>(splitStoic[1]) // this is the stoichiometry number
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
