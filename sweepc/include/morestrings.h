/*
  Author(s):      Matthew Celnik (msc37)
  Project:        none

  File purpose:
    Some additional subroutines for STL strings that aren't provided.
*/

#ifndef MORE_STRINGS_H
#define MORE_STRINGS_H

#include <string>
#include <vector>
#include <stdlib.h>
#include <sstream>

using namespace std;

/* Splits a string into sub strings using the given delimiter
   characters. */
inline void split(string &str, vector<string> &subs, string &delims)
{
    // Clear the vector of substrings.
    subs.clear();

    // Locate the first substring in the string. Iterator i1 now
    // points to the first character of the substring and i2 points to the
    // delimiter after the substring.
    string::size_type i1 = str.find_first_not_of(delims, 0);
    string::size_type i2 = str.find_first_of(delims, i1);

    // Loop over the whole string.
    while((i1!=string::npos) || (i2!=string::npos)) {
        // Add the substring to the list.
        subs.push_back(str.substr(i1,i2-i1));

        // Locate the next substring.
        i1 = str.find_first_not_of(delims, i2);
        i2 = str.find_first_of(delims, i1);
    }
};

/* Returns the index of a string in a vector of strings. Returns
   < 0 is not found. */
inline int findinlist(string &str, vector<string> &list)
{
    // Very straightforward, just loop over the list and locate
    // the first element that matches the given string.
    int i;
    for (i=0; i<(int)list.size(); i++) {
        if (str.compare(list[i]) == 0) return i;
    }
    return -1;
}

/* Converts a number less than 10 into a single character. */
inline char &cchar(unsigned int i, char &c)
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
}

///* Converts an integer to a string. */
//inline string cstr(int i)
//{
//    string str = "";
//    
//    // Build string in reverse order.
//    int a = abs(i), p = 1, n = 0;
//    div_t d = div(a, 10);
//    char c;
//    while (d.quot != 0) {
//        n = d.rem / p;
//        str.append(&cchar(n, c));
//        p *= 10;
//        d = div(a, p*10);
//    }
//    if (i<0) str.append("-");
//
//    // Reverse string to get correct string order.
//    char tmp;
//    n=0, a=(int)str.length()-1;
//    while(n<a) {
//        tmp    = str[n];
//        str[n] = str[a];
//        str[a] = tmp;
//        n++; a--;
//    }
//    return str;
//}

/* Converts any type into a string representation. */
template<class T> string cstr(T &t)
{
    // This routine takes advantage of the STL ostringstream class.
    ostringstream ss;

    if ((ss << t).fail()) {
        return "X";
    } else {
        return ss.str();
    }
};
#endif