/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    The Element library holds a list of known elements.
*/

#ifndef GPC_ELEMENT_LIB_CPP
#define GPC_ELEMENT_LIB_CPP

#include "gpc_element.h"

// ELEMENT LIBRARY.
const Sprog::Element Sprog::Element::m_lib[Element::m_nlib] = {
    Sprog::Element("O",15.9994), 
    Sprog::Element("H",1.0079), 
    Sprog::Element("N", 14.0067), 
    Sprog::Element("C", 12.011),
    Sprog::Element("AR",39.948)
};

#endif // GPC_ELEMENT_LIB_CPP
