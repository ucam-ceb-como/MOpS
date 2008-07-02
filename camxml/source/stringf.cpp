/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        camxml (Cambridge c++ XML library)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    Implementation of functions declared in the strungf.h header file.

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

#include "stringf.h"

namespace ComoString {
    inline bool rtrimable(std::wstring & str) {
	    if (str.length() == 0) {
		    return false;
	    } else if ((str[str.length()-1] == L' ')||(str[str.length()-1] == L'\n')||(str[str.length()-1] == L'\t')||(str[str.length()-1] == L'\r')) {
		    return true;
	    } else {
		    return false;
	    }
    }
    inline bool ltrimable(std::wstring & str) {
	    if (str.length() == 0) {
		    return false;
	    } else if ((str[0] == L' ')||(str[0] == L'\n')||(str[0] == L'\t')||(str[0] == L'\r')) {
		    return true;
	    } else {
		    return false;
	    }
    }
    std::wstring trim(std::wstring str) {

        while (rtrimable(str)) {
            str.erase(str.length()-1,1);
        }
        while (ltrimable(str)) {
            str.erase(0,1);
        }
        return str;
    }

    std::wstring insertSpace(const int depth) {
        std::wstring space = L"";
	    for (int i=0; i < depth; i++) {
		    space.append(L"    ");
	    }
        return space;
    }
    int string2int (const std::string &str) {
        std::stringstream ss(str);
        int n;
        ss >> n;
        return n;
    }

    std::string int2string (int n) {
        std::stringstream ss;
        ss << n;
        return ss.str();
    }
}
