/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of inline member functions of the Component class.
*/

#ifndef SWEEP_COMPONENT_INL_H
#define SWEEP_COMPONENT_INL_H

#include "swp_params.h"
#include "swp_component.h"
#include <string>

// MOLECULAR WEIGHT

// Returns component molecular weight (g/mol).
inline Sweep::real Sweep::Component::MolWt() const {return m_molwt;};

// Sets the molecular weight (g/mol).
inline void Sweep::Component::SetMolWt(const Sweep::real molwt) {m_molwt = molwt;};

// DENSITY.

// Returns component density (g/cm3).
inline Sweep::real Sweep::Component::Density() const {return m_density;};

// Sets the density (g/cm3).
inline void Sweep::Component::SetDensity(const Sweep::real dens) {m_density = dens;};

// COMPONENT NAME.

// Returns component symbol or name.
inline const std::string &Sweep::Component::Name() const {return m_name;};

// Sets the symbol or name.
inline void Sweep::Component::SetName(const std::string &name) {m_name = name;};

#endif
