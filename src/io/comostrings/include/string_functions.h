/*
  Author(s):      Matthew Celnik (msc37)
  Project:        comostrings
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Some additional subroutines for manipulating STL strings.

  Licence:
    This file is part of the "comostrings" library.

    comostrings is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
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
        std::vector<std::string> &string_list // Vector of strings in which 
                                       // to locate target string.
        );

    // Converts a number less than 10 into a single character.  If
    // number is greater than 9 then returns 'X'.
    char &cchar(
        unsigned int i, // Integer to convert.
        char &c         // Output character.
        );

    // Converts any type into a string representation.
    template<class T> std::string cstr(const T &t);

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
	//following functions are added by vinod

	std::string trim(const std::string &str);
	bool isEmpty(const std::string &str);


};

#endif
