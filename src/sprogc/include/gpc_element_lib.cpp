/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Element library holds a list of known elements.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
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

#ifndef GPC_ELEMENT_LIB_CPP
#define GPC_ELEMENT_LIB_CPP

#include "gpc_element.h"

// ELEMENT LIBRARY.
const Sprog::Element Sprog::Element::m_lib[Element::m_nlib] = {
    Sprog::Element("H",    1.0079), 
    Sprog::Element("HE",   4.0026),
    Sprog::Element("LI",   6.941),
    Sprog::Element("BE",   9.01218),
    Sprog::Element("B",   10.81),
    Sprog::Element("C",   12.011),
    Sprog::Element("N",   14.0067), 
    Sprog::Element("O",   15.9994), 
    Sprog::Element("F",   18.9984), 
    Sprog::Element("NE",  20.179), 
    Sprog::Element("NA",  22.98977), 
    Sprog::Element("MG",  24.305), 
    Sprog::Element("AL",  26.9815), 
    Sprog::Element("SI",  28.0855), 
    Sprog::Element("P",   30.9738), 
    Sprog::Element("S",   32.06), 
    Sprog::Element("CL",  35.453), 
    Sprog::Element("AR",  39.948),
    Sprog::Element("K",   39.0983), 
    Sprog::Element("CA",  40.08), 
    Sprog::Element("SC",  44.9559), 
    Sprog::Element("TI",  47.88), 
    Sprog::Element("CR",  50.9415), 
    Sprog::Element("MN",  51.996), 
    Sprog::Element("FE",  55.847), 
    Sprog::Element("CO",  58.9332), 
    Sprog::Element("NI",  58.69), 
    Sprog::Element("CU",  63.546), 
    Sprog::Element("ZN",  65.38), 
    Sprog::Element("GA",  69.72), 
    Sprog::Element("GE",  72.59), 
    Sprog::Element("AS",  74.9216), 
    Sprog::Element("SE",  78.96), 
    Sprog::Element("BR",  79.904), 
    Sprog::Element("KR",  83.80), 
    Sprog::Element("RB",  85.4678), 
    Sprog::Element("SR",  87.62), 
    Sprog::Element("Y",   88.9059), 
    Sprog::Element("ZR",  91.22), 
    Sprog::Element("NB",  92.9064), 
    Sprog::Element("MO",  95.94), 
    Sprog::Element("TC", (98)), 
    Sprog::Element("RU", 101.07), 
    Sprog::Element("RH", 102.9055), 
    Sprog::Element("PD", 106.42), 
    Sprog::Element("AG", 107.868), 
    Sprog::Element("CD", 112.41), 
    Sprog::Element("IN", 114.82), 
    Sprog::Element("SN", 118.69), 
    Sprog::Element("SB", 121.75), 
    Sprog::Element("TE", 127.60), 
    Sprog::Element("I",  126.9045), 
    Sprog::Element("XE", 131.29), 
    Sprog::Element("CS", 132.9054), 
    Sprog::Element("BA", 137.33), 
    Sprog::Element("LA", 138.9055), 
    Sprog::Element("HF", 178.49), 
    Sprog::Element("TA", 180.9479), 
    Sprog::Element("W",  183.85), 
    Sprog::Element("RE", 186.207), 
    Sprog::Element("OS", 190.2), 
    Sprog::Element("IR", 192.22), 
    Sprog::Element("PT", 195.08), 
    Sprog::Element("AU", 196.9665), 
    Sprog::Element("HG", 200.59), 
    Sprog::Element("TL", 204.383), 
    Sprog::Element("PB", 207.2), 
    Sprog::Element("BI", 208.9804), 
    Sprog::Element("PO", (209)), 
    Sprog::Element("AT", (210)), 
    Sprog::Element("RN", (222)), 
    Sprog::Element("FR", (223)), 
    Sprog::Element("RA", 226.0254), 
    Sprog::Element("AC", 227.0278)
};

#endif // GPC_ELEMENT_LIB_CPP
