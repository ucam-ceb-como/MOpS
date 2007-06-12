/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    Inline member functions of Element class.
*/

#ifndef GPC_ELEMENT_INL_H
#define GPC_ELEMENT_INL_H

inline const std::string Element::Name(void) const {return m_name;};
inline void Element::SetName(const std::string &name) {m_name=name;};
inline const real Element::MolWt(void) const {return m_molwt;};
inline void Element::SetMolWt(const real molwt) {m_molwt=molwt;};

#endif