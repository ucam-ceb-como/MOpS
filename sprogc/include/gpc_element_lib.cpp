/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    The Element library holds a list of known elements.  This file is included in the gpc_element.h
    header.
*/

#ifndef GPC_ELEMENT_LIB_H
#define GPC_ELEMENT_LIB_H

#include "gpc_params.h"
#include "gpc_element.h"
#include <vector>
#include <string>

// ELEMENT LIBRARY.
const Sprog::Element Sprog::Element::m_lib[Element::m_nlib] = {
    Sprog::Element("O",15.9994), 
    Sprog::Element("H",1.0079), 
    Sprog::Element("N", 14.0067), 
    Sprog::Element("C", 12.011),
    Sprog::Element("AR",39.948)
};

#endif // GPC_ELEMENT_LIB_H
