/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    Inline function definitions for Species class.
*/

#ifndef GPC_SPECIES_INL_H
#define GPC_SPECIES_INL_H

//  SPECIES NAME.
inline const std::string &Species::Name() const {return m_name;};

// SPECIES COMPOSITION.
inline const ElCompVector &Species::Composition() const {return m_elcomp;};
inline unsigned int Species::ComponentCount(void) const {return m_elcomp.size();};

// SPECIES MOLECULAR WEIGHT.
inline real Species::MolWt() const {return m_molwt;};

// THERMODYNAMIC FITTING PARAMETERS.
inline unsigned int Species::ThermoRangeCount(void) const {return m_thermoparams.size();};
inline void Species::SetThermoStartTemperature(const real T) {m_T1 = T;};

// PARENT MECHANISM.
inline const Sprog::Mechanism *const Species::Mechanism(void) const {return m_mech;};

// CLONING.
inline Species *const Species::Clone(void) const {return new Species(*this);};

#endif
