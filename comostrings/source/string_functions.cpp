/*

  Author(s):      Matthew Celnik (msc37)
  Project:        comostrings (mopssuite)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the functions declared in the string_functions.h
    header file.

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

#include "string_functions.h"
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>

using namespace Strings;

void Strings::split(const std::string &str,
                    std::vector<std::string> &subs,
                    const std::string &delims)
{
    // Clear the vector of substrings.
    subs.clear();

    // Locate the first substring in the string. Iterator i1 now
    // points to the first character of the substring and i2 points to the
    // delimiter after the substring.
    std::string::size_type i1 = str.find_first_not_of(delims, 0);
    std::string::size_type i2 = str.find_first_of(delims, i1);

    // Loop over the whole string.
    while((i1!=std::string::npos) || (i2!=std::string::npos)) {
        // Add the substring to the list.
        subs.push_back(str.substr(i1,i2-i1));

        // Locate the next substring.
        i1 = str.find_first_not_of(delims, i2);
        i2 = str.find_first_of(delims, i1);
    }
}

int Strings::findinlist(std::string str, std::vector<std::string> &string_list)
{
    // Very straightforward, just loop over the list and locate
    // the first element that matches the given string.
    int i;
    for (i=0; i<(int)string_list.size(); i++) {
        if (str.compare(string_list[i]) == 0) return i;
    }
    return -1;
}

char &Strings::cchar(unsigned int i, char &c)
{
    if (i<10) {
        switch (i) {
            case 0: c = '0';
            case 1: c = '1';
            case 2: c = '2';
            case 3: c = '3';
            case 4: c = '4';
            case 5: c = '5';
            case 6: c = '6';
            case 7: c = '7';
            case 8: c = '8';
            case 9: c = '9';
        }
    } else {
        c = 'X';
    } 
    return c;
};


template<class T> std::string Strings::cstr(const T &t)
{
    // This routine takes advantage of the STL ostringstream class.
    std::ostringstream ss;

    if ((ss << t).fail()) {
        return "X";
    } else {
        return ss.str();
    }
}
template std::string Strings::cstr<int>(const int&);
template std::string Strings::cstr<unsigned int>(const unsigned int&);
template std::string Strings::cstr<float>(const float&);
template std::string Strings::cstr<double>(const double&);
template std::string Strings::cstr<long double>(const long double&);

double Strings::cdble(const std::string &str)
{
    // First remove all the white space from the string.
    std::string nospace = removeWhiteSpace(str);

    // Now convert to a double.
    return atof(nospace.c_str());
}

std::string Strings::removeWhiteSpace(const std::string &str)
{
    std::string nospace = "";
    for (unsigned int i=0; i<str.length(); i++) {
        if (!isWhiteSpace(str[i])) {
            nospace.append(str.substr(i,1));
        }
    }
    return nospace;
}

bool Strings::isWhiteSpace(const char c)
{
    return (c==' ') || (c=='\n') || (c=='\t') || (c=='\r');
}

bool Strings::isLetter(const char c)
{
    return (c=='a') || (c=='b') || (c=='c') || (c=='d') || (c=='e') ||
           (c=='f') || (c=='g') || (c=='h') || (c=='i') || (c=='j') ||
           (c=='k') || (c=='l') || (c=='m') || (c=='n') || (c=='o') ||
           (c=='p') || (c=='q') || (c=='r') || (c=='s') || (c=='t') ||
           (c=='u') || (c=='v') || (c=='w') || (c=='x') || (c=='y') ||
           (c=='z') ||
           (c=='A') || (c=='B') || (c=='C') || (c=='D') || (c=='E') ||
           (c=='F') || (c=='G') || (c=='H') || (c=='I') || (c=='J') ||
           (c=='K') || (c=='L') || (c=='M') || (c=='N') || (c=='O') ||
           (c=='P') || (c=='Q') || (c=='R') || (c=='S') || (c=='T') ||
           (c=='U') || (c=='V') || (c=='W') || (c=='X') || (c=='Y') ||
           (c=='Z');
}

bool Strings::isLetterOrNum(const char c)
{
    return isLetter(c) ||
           (c=='0') || (c=='1') || (c=='2') || (c=='3') || (c=='4') ||
           (c=='5') || (c=='6') || (c=='7') || (c=='8') || (c=='9');
}

std::string Strings::convertToCaps(const std::string &str)
{
    std::string::const_iterator i;
    std::string caps;

    for (i=str.begin(); i!=str.end(); i++) {
        switch (*i) {
            case 'a': caps.append("A"); break;
            case 'b': caps.append("B"); break;
            case 'c': caps.append("C"); break;
            case 'd': caps.append("D"); break;
            case 'e': caps.append("E"); break;
            case 'f': caps.append("F"); break;
            case 'g': caps.append("G"); break;
            case 'h': caps.append("H"); break;
            case 'i': caps.append("I"); break;
            case 'j': caps.append("J"); break;
            case 'k': caps.append("K"); break;
            case 'l': caps.append("L"); break;
            case 'm': caps.append("M"); break;
            case 'n': caps.append("N"); break;
            case 'o': caps.append("O"); break;
            case 'p': caps.append("P"); break;
            case 'q': caps.append("Q"); break;
            case 'r': caps.append("R"); break;
            case 's': caps.append("S"); break;
            case 't': caps.append("T"); break;
            case 'u': caps.append("U"); break;
            case 'v': caps.append("V"); break;
            case 'w': caps.append("W"); break;
            case 'x': caps.append("X"); break;
            case 'y': caps.append("Y"); break;
            case 'z': caps.append("Z"); break;
            default: caps.append(i,i+1); break;
        }
    }
    return caps;

}

std::string Strings::trim(const std::string &str){
	int b = str.find_first_not_of(" \t");
	if( b < 0)
		return "";
	int e = str.find_last_not_of(" \t");
	
	return str.substr(b,e-b+1);
}

bool Strings::isEmpty(const std::string &str){
	return (trim(str)).empty();
}
