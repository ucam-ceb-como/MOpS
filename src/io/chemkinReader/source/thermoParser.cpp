/*
 * thermoParser.cpp
 *
 *  Created on: Jun 23, 2011
 *      Author: gigadot
 *     License: Apache 2.0
 */

#include "thermoParser.h"
#include "stringFunctions.h"
#include <istream>
#include "boost/algorithm/string/erase.hpp"
#include "boost/algorithm/string/replace.hpp"
#include "boost/algorithm/string/trim.hpp"

using namespace std;
using namespace boost;

IO::ThermoParser::ThermoParser(const string thermo_file)
:
thermo_file_(convertToCaps(thermo_file)),
thermo_file_string_(convertToCaps(fileToString(thermo_file))),
lines_(fileToStrings(thermo_file)),
globalLowT_(-1),
globalCommonT_(-1),
globalHighT_(-1)
{}

void IO::ThermoParser::parse(vector<Species>& species) {

    getGlobalTemperatures();

    cout << "Parsing NASA thermo file: " << thermo_file_ << endl;
    parseAllThermoData();

    ensureSpeciesNamesAreValid();
    ensureNoDuplicates();

    for (size_t i = 0; i != species.size(); ++i) {
        if (!setThermoDataFor(species[i])) {
            throw runtime_error("Thermo data for Species " + species[i].name() + " not found in " + thermo_file_);
        }
    }
    cout << "End of Parsing NASA thermo file: " << thermo_file_ << endl;

}

void IO::ThermoParser::ensureNoDuplicates() {
    for (size_t i = 0; i < thermos_.size(); ++i) {
        string aname = thermos_[i].getSpeciesName();
        for (size_t j = i; j < thermos_.size(); ++j) {
            string bname = thermos_[j].getSpeciesName();
            if ((i != j) && (aname.compare(bname) == 0)) {
                //throw runtime_error("Species name '" + aname + "' is found more than once in thermo file.");
            }
        }
    }
}
void IO::ThermoParser::ensureSpeciesNamesAreValid() {
    const regex space("\\s");
    for (size_t i = 0; i < thermos_.size(); ++i) {
        string name = thermos_[i].getSpeciesName();
        if (regex_search(name, space)) {
            //throw runtime_error("Invalid species name '" + name + "'. Species name must contain no spaces.");
        }
    }
}

bool IO::ThermoParser::setThermoDataFor(IO::Species& species) {
    for (size_t i = 0; i < thermos_.size(); ++i) {
        if (convertToCaps(thermos_[i].getSpeciesName()).compare(species.name()) == 0) {
            species.thermo() = thermos_[i];
            return true;
        }
    }
    return false;
}

void IO::ThermoParser::parseAllThermoData() {

    vector<string> thermo_lines = getThermoSection(lines_);

    const regex empty("\\s*");

    for (unsigned int i = 0; i < thermo_lines.size(); i++) {
        if (isSectionMatchedNASA(thermo_lines, i)) {
            parseNASASection(thermo_lines[i], thermo_lines[i + 1], thermo_lines[i + 2], thermo_lines[i + 3]);
            i += 3;
        } else if (!regex_match(thermo_lines[i], empty, regex_constants::match_single_line)) {
            // Only print unmatched lines if it is not empty.
            cout << "Unmatched : " << thermo_lines[i] << endl;
        }
    }

    //    // this is where it print out info at the moment
    //    for (unsigned int i = 0; i < thermos_.size(); i++) {
    //        cout << thermos_[i] << endl;
    //    }

}

bool IO::ThermoParser::parseNASASection(string l1, string l2, string l3, string l4) {

    cout << l1 << endl;
    // line 1

    string speciesString = l1.substr(0, 18);
    string speciesName = extractSpeciesName(speciesString);

    Thermo thermo(speciesName);
    thermo.setNote(trim_copy(l1.substr(18, 6)));
    thermo.setPhase(l1.substr(44, 1));
    thermo.setTLow(from_string<double>(trim_copy(l1.substr(45, 10))));
    thermo.setTHigh(from_string<double>(trim_copy(l1.substr(55, 10))));

    if (trim_copy(l1.substr(65, 8)) == "")
    {
        thermo.setTCommon(globalCommonT_);
    }
    else
    {
        thermo.setTCommon(from_string<double>(trim_copy(l1.substr(65, 8))));
    }
    string elements_string = convertToCaps(l1.substr(24, 20));
    thermo.setElements(parseElements(elements_string));
    // line 2, 3 4

    double ah1 = from_string<double>(boost::erase_all_copy(l2.substr(0, 15)," "));
    double ah2 = from_string<double>(boost::erase_all_copy(l2.substr(15, 15)," "));
    double ah3 = from_string<double>(boost::erase_all_copy(l2.substr(30, 15)," "));
    double ah4 = from_string<double>(boost::erase_all_copy(l2.substr(45, 15)," "));
    double ah5 = from_string<double>(boost::erase_all_copy(l2.substr(60, 15)," "));
    double ah6 = from_string<double>(boost::erase_all_copy(l3.substr(0, 15)," "));
    double ah7 = from_string<double>(boost::erase_all_copy(l3.substr(15, 15)," "));
    double al1 = from_string<double>(boost::erase_all_copy(l3.substr(30, 15)," "));
    double al2 = from_string<double>(boost::erase_all_copy(l3.substr(45, 15)," "));
    double al3 = from_string<double>(boost::erase_all_copy(l3.substr(60, 15)," "));
    double al4 = from_string<double>(boost::erase_all_copy(l4.substr(0, 15)," "));
    double al5 = from_string<double>(boost::erase_all_copy(l4.substr(15, 15)," "));
    double al6 = from_string<double>(boost::erase_all_copy(l4.substr(30, 15)," "));
    double al7 = from_string<double>(boost::erase_all_copy(l4.substr(45, 15)," "));
    thermo.setUpperTemperatureCoefficients(ah1, ah2, ah3, ah4, ah5, ah6, ah7);
    thermo.setLowerTemperatureCoefficients(al1, al2, al3, al4, al5, al6, al7);
    thermos_.push_back(thermo);
    return false;
}

map<string, int> IO::ThermoParser::parseElements(string elements_string) {

    std::map<std::string, int> elem_count_map;
    if (elements_string.length() % 5 == 0) {
        for (unsigned int i = 0; i < elements_string.length() / 5; i++) {
            string elem = trim_copy(elements_string.substr(i * 5, 3));
            string count_str = trim_copy(elements_string.substr(3 + i * 5, 2));
            if (count_str.length() == 0) continue;
            // check if count string is a number
            int count = from_string<int>(count_str);
            if (count == 0) continue;
            if (elem_count_map.find(elem) != elem_count_map.end()) {
               // throw runtime_error("Are you tried to play trick with me? Duplicated element found: " + elem);
            }
            elem_count_map.insert(pair<string, int> (elem, count));
        }
    } else {
        throw runtime_error("Invalid element string found for value: " + elements_string);
    }
    return elem_count_map;

}

bool IO::ThermoParser::isSectionMatchedNASA(std::vector<std::string> lines, unsigned int offset) {
    const regex nasaLine1("[^!]{79}1.*");
    const regex nasaLine2("[^!]{79}2.*");
    const regex nasaLine3("[^!]{79}3.*");
    const regex nasaLine4("[^!]{79}4.*");
    if (lines.size() >= (offset + 4)) {
        return regex_match(lines.at(offset), nasaLine1, regex_constants::match_single_line) &&
                regex_match(lines.at(offset + 1), nasaLine2, regex_constants::match_single_line) &&
                regex_match(lines.at(offset + 2), nasaLine3, regex_constants::match_single_line) &&
                regex_match(lines.at(offset + 3), nasaLine4, regex_constants::match_single_line);
    } else {
        return false;
    }
}

std::vector<std::string> IO::ThermoParser::getThermoSection(std::vector<std::string> lines) {
    vector<string> thermo_lines;
    int begin = -1;
    int end = -1;
    const regex thermoTag("\\s*THER(?:|MO).*");
    const regex endTag("\\s*END.*");

    // this loop only for finding index of THERMO and END keywords in the lines
    for (unsigned int i = 0; i < lines.size(); i++) {
        if ((begin < 0) && regex_match(lines.at(i), thermoTag, regex_constants::match_single_line)) {
            begin = (int) i;
        } else if ((end < 0) && regex_match(lines.at(i), endTag, regex_constants::match_single_line)) {
            end = (int) i;
            break;
        }
    }

    // check if the indexes are correct
    if (-1 < begin && begin < end) {
        for (int i = ++begin; i < end; i++) {
            thermo_lines.push_back(lines[i]);
        }
    } else {
        cout << "Error: Thermo file does not contain valid starting and ending tags." << endl;
    }

    return thermo_lines;
}

std::string
IO::ThermoParser::extractSpeciesName(const string& speciesString)
const
{
    const regex speciesNameRegex("^\\s*(.*?)\\s+");
    string speciesName;
    smatch what;
    string::const_iterator start = speciesString.begin();
    string::const_iterator end = speciesString.end();

    if(regex_search(start, end, what, speciesNameRegex))
    {
        speciesName = what[1];
    } else
    {
        throw std::runtime_error("Can't find a species name. There is probably no"
                                 "space after column 18.");
    }
    return speciesName;
}

void
IO::ThermoParser::getGlobalTemperatures()
{
    const regex globalTRegex("\\s*THER(?:|MO)\\s+(?:|ALL)\\s*([0-9\\.]+)\\s+([0-9\\.]+)\\s+([0-9\\.]+)\\s*");
    string speciesName;
    smatch what;

    string::const_iterator start = thermo_file_string_.begin();
    string::const_iterator end = thermo_file_string_.end();

    if(regex_search(start, end, what, globalTRegex))
    {
        globalLowT_ = from_string<double>(what[1]);
        globalCommonT_ = from_string<double>(what[2]);
        globalHighT_ = from_string<double>(what[3]);
    } else
    {
        throw std::runtime_error("Could not find list of global temperatures.");
    }
}
