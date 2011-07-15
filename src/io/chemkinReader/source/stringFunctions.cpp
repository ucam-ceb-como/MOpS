/*
 * stringFunctions.cpp
 *
 *  Created on: Jun 24, 2011
 *      Author: Laurence R. McGlashan
 */

#include "stringFunctions.h"
#include "boost/regex.hpp"
#include "boost/algorithm/string.hpp"
#include <string>

using namespace std;
using namespace boost;

const string
IO::fileToString(const string& fileName)
{
    ifstream fin(fileName.c_str(), ios::in);

    string fileInString((istreambuf_iterator<char>(fin)),
                         istreambuf_iterator<char>());

    return fileInString;
}

const vector<string>
IO::fileToStrings(const string fileName)
{
    vector<string> lines;
    ifstream fin(fileName.c_str(), ios::in);
    string line;
    while (getline(fin, line)) {
        lines.push_back(line);
    }
    return lines;
}

string
IO::regex_escape(const string& string_to_escape)
{
    static const regex re_boostRegexEscape("[\\^\\.\\$\\|\\(\\)\\[\\]\\*\\+\\?\\/\\\\]");
    const string rep("\\\\\\1&");
    string result = regex_replace(string_to_escape, re_boostRegexEscape, rep, match_default | format_sed);
    return result;
}


string
IO::replaceComments(string stringToReplace)
{
    regex commentRegex("(!.*?)\\n|(!.*?)$");
    string format_string = " \n";
    stringToReplace = regex_replace(stringToReplace, commentRegex, format_string, match_default | format_sed);
    return stringToReplace;
}

string
IO::convertToCaps(const string& str)
{
    return boost::to_upper_copy(str);
}

string
IO::trim(const string &str)
{
    int b = str.find_first_not_of(" \t");
    if (b < 0)
        return "";
    int e = str.find_last_not_of(" \t");

    return str.substr(b, e - b + 1);
}
