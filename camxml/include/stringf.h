/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        camxml (Cambridge c++ XML library)

  File purpose:
    Additional string manipulation functions for creating 
    XML wide string in xmlelement.h and xmlelement.cpp
*/


#ifndef CAMXML_STRING_FUNCTIONS_H
#define CAMXML_STRING_FUNCTIONS_H

#include <string>
//#include <cassert>
//#include <cstdlib>
//#include <iostream>
//#include <iomanip>
using namespace std;

namespace ComoString {
    // Check if there is any white space to the right of wstring.
    inline bool rtrimable(std::wstring & str);
    // Check if there is any white space to the left of wstring.
    inline bool ltrimable(std::wstring & str);
    // Remove white space from both side of wstring.
    std::wstring trim(std::wstring str);
    // Return a wstring containing 4 x depth space bar charaters.
    // only use in xmlelement.cpp to construct nicely XML wstring.
    std::wstring insertSpace(const int depth);
}



#endif
