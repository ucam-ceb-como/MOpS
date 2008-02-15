/*
  Author(s):      Matthew Celnik (msc37)
  Project:        comostrings

  File purpose:
    Some additional subroutines for manipulating STL strings.
*/

#ifndef COMOSTRINGS_SFUNCTIONS_H
#define COMOSTRINGS_SFUNCTIONS_H

#include <string>
#include <vector>

namespace Strings
{
    // Splits a string into substrings using the given delimiter
    // characters.  The delimiter characters are supplied as a
    // string, all the characters of which are used as delimiter.
    void split(
        const std::string &str,         // The string to split.
        std::vector<std::string> &subs, // Output vector of substrings. 
        const std::string &delims       // List of delimiters.
        );

    // Returns the index of a string in a vector of strings. Returns
    // < 0 if not found.
    int findinlist(
        std::string str,               // String to locate.
        std::vector<std::string> &list // Vector of strings in which 
                                       // to locate target string.
        );

    // Converts a number less than 10 into a single character.  If
    // number is greater than 9 then returns 'X'.
    char &cchar(
        unsigned int i, // Integer to convert.
        char &c         // Output character.
        );

    // Converts any type into a string representation.
    template<class T> std::string cstr(T &t);

    // Converts a string into a double if possible.
    double cdble(const std::string &str);

    // Removes all the white space characters from the string.
    std::string removeWhiteSpace(const std::string &str);

    // Returns true if a character is a whitespace character.
    bool isWhiteSpace(const char c);

    // Returns true if a character is a Roman letter.
    bool isLetter(const char c);

    // Returns true if a character is a Roman letter or a number.
    bool isLetterOrNum(const char c);

    // Converts a string to all capital letters.
    std::string convertToCaps(const std::string &str);
};

#endif
