/*
  Author(s):      Matthew Celnik (msc37)
  Project:        camxml (Cambridge c++ XML library)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Attribute class declared in the xmlattribute.h
    header file.

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

#include "xmlattribute.h"

using namespace CamXML;

Attribute::Attribute(void)
{
}

Attribute::Attribute(const std::string &name, const std::string &value)
{
    m_name = ComoUnicode::StringToWString(name);
    m_value = ComoUnicode::StringToWString(value);
}
Attribute::Attribute(const std::wstring &name, const std::wstring &value)
{
    m_name = name;
    m_value = value;
}

Attribute::~Attribute(void)
{
}
