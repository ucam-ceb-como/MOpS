/*
 * reactionParser.h
 *
 *  Created on: Jun 23, 2011
 *      Author: Laurence R. McGlashan
 *     License: GPL
 */

#ifndef REACTIONPARSER_H_
#define REACTIONPARSER_H_

#include "boost/regex.hpp"
#include <string>
#include <map>

namespace IO {

    class Reaction;

    class ReactionParser
    {

            static const boost::regex reactionSingleRegex;
            static const boost::regex blankLine;
            static const boost::regex DUPLICATE;
            static const boost::regex LOW;
            static const boost::regex TROE;
            static const boost::regex SRI;
            static const boost::regex REV;
            static const boost::regex pressureDependent;

            const std::string reactionString_;
            std::vector<std::string> reactionStringLines_;

        public:

            ReactionParser(const std::string reactionString);

            ~ReactionParser(){}

            void parse(std::vector<IO::Reaction>& reactions);

            std::multimap<std::string, double> parseReactionSpecies(std::string reactants);

            std::multimap<std::string, double> parseThirdBodySpecies(const std::string& thirdBodies);

            bool isBlankLine(const std::string& line);

            std::string findLineType(const std::string& line);

            std::vector<double> parseLOWTROEREV(const std::string& line, const boost::regex& reg);

            bool checkForPressureDependentReaction(const std::string& line);
    };

}

#endif /* REACTIONPARSER_H_ */
