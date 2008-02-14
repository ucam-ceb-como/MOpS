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
#include <sstream>
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
    // Convert string to integer
    int string2int (const std::string &str);
    // Convert integer to string
    std::string int2string (int n);
    /* In the future use this to modified string2int to return error if conversion fails
    class BadConversion : public std::runtime_error {
    public:
       BadConversion(const std::string& s)
         : std::runtime_error(s)
         { }
    };
     
    inline std::string stringify(double x)
    {
       std::ostringstream o;
       if (!(o << x))
         throw BadConversion("stringify(double)");
       return o.str();
    }*/

}



#endif
