/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    Inline function definitions for Species class.
*/


#ifndef GPC_SPECIES_INL_H
#define GPC_SPECIES_INL_H

inline const std::string &Species::Name() const {return m_name;};
inline void Species::SetName(const std::string &name) {m_name=name;};
inline const ElCompVector &Species::Composition() const {return m_elcomp;};
inline const real Species::MolWt() const {return m_molwt;};
inline void Species::SetMolWt(const real molwt) {m_molwt=molwt;};
inline unsigned int Species::ThermoRangeCount(void) const {return m_thermoparams.size();};
inline void Species::SetThermoStartTemperature(const real T) {m_T1 = T;};

#endif