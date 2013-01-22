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

namespace {

/*!
 * New getline function to handle reading of LF or CRLF line endings.
 * Function taken from StackOverflow with permalink:
 * http://stackoverflow.com/a/6089413/1021919
 *
 * @param is    Input filestream
 * @param t     Line to get
 * @return
 */
std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

}


using namespace std;
using namespace boost;

/*!
 * Read a file and convert it to a string, INCLUDING the line
 * ending characters. Also convert all to caps. Use safe getline
 * to permit parsing of CRLF lines endings.
 *
 * @param fileName  Name of file to read
 * @return          File as a string
 */
string IO::fileToString(const string& fileName)
{
    ifstream fin(fileName.c_str(), ios::in);

    string fs;
    string line;
    while (::safeGetline(fin, line)) {
        fs.append(convertToCaps(line));
        fs.append("\n");
    }
    return convertToCaps(fs);
}

/*!
 * Read a file and convert it to a vector of strings, with each line a new
 * entry in the vector. Use safe getline to permit parsing of CRLF line
 * endings.
 *
 * @param fileName  Name of file to read
 * @return          File as a vector of strings
 */
vector<string> IO::fileToStrings(const string fileName)
{
    vector<string> lines;
    ifstream fin(fileName.c_str(), ios::in);
    string line;
    while (::safeGetline(fin, line)) {
        lines.push_back(convertToCaps(line));
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
IO::replaceComments(const string& stringToReplace)
{
    regex commentRegex("(!.*?)\\n|(!.*?)$");
    string format_string = " \n";
    string result = regex_replace(stringToReplace, commentRegex, format_string, match_default | format_sed);
    return result;
}

string
IO::convertToCaps(const string& str)
{
    return boost::to_upper_copy(str);
}

vector<string>
IO::convertToCaps(const vector<string>& str)
{
    for (size_t i=0; i<str.size(); ++i)
    {
        convertToCaps(str[i]);
    }
    return str;
}

//! Check the format of the number.
void
IO::checkNumberFormat(std::string& t)
{
    const boost::regex numberFormatRegex("[A-CI-Z]");
    boost::smatch what;

    std::string::const_iterator start = t.begin();
    std::string::const_iterator end = t.end();

    if (boost::regex_search(start, end, what, numberFormatRegex))
    {
        throw std::runtime_error("from_string<> is complaining about the format"
                                 " of a string you've given it.");
    } else
    {
         // This is a bastard: Some numbers have D or G in them instead of E which
         // from_string can't deal with.
         boost::replace_all(t,"D","E");
         boost::replace_all(t,"F","E");
         boost::replace_all(t,"G","E");
         boost::replace_all(t,"H","E");
    }
}
