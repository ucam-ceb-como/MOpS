/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        camxml (Cambridge c++ XML library)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    Additional string manipulation functions for creating 
    XML wide string in xmlelement.h and xmlelement.cpp

  Licence:
    This file is part of the "camxml" library.

    camxml is free software; you can redistribute it and/or
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


#ifndef CAMXML_STRING_FUNCTIONS_H
#define CAMXML_STRING_FUNCTIONS_H

#include <string>
//#include <cassert>
//#include <cstdlib>
#include <sstream>
//#include <iomanip>

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
