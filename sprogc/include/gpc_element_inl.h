/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    Inline member functions of Element class.  This file is included
    at the end of the gpc_element.h file, which contains the definition
    of the element.
*/

#ifndef GPC_ELEMENT_INL_H
#define GPC_ELEMENT_INL_H

// ELEMENT NAME.
inline const std::string Element::Name(void) const {return m_name;};

// MOLECULAR WEIGHT.
inline real Element::MolWt(void) const {return m_molwt;};

#endif
