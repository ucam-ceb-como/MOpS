/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of inline member functions of the Component class.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#ifndef SWEEP_COMPONENT_INL_H
#define SWEEP_COMPONENT_INL_H

#include "swp_params.h"
#include "swp_component.h"
#include <string>

// MOLECULAR WEIGHT

// Returns component molecular weight (kg/mol).
inline double Sweep::Component::MolWt() const {return m_molwt;};

// Sets the molecular weight (kg/mol).
inline void Sweep::Component::SetMolWt(const double molwt) {m_molwt = molwt;};

// DENSITY.

//! Returns the coalescence threshold 
inline double Sweep::Component::CoalescThresh() const {return m_coalesc_thresh;};

//! Returns the growhtfact 
inline double Sweep::Component::GrowthFact() const {return m_growthfact;};

//! Returns the minpah 
inline double Sweep::Component::MinPAH() const {return m_minPAH;};

//! Returns component density (kg/m3).
inline double Sweep::Component::Density() const {return m_density;};

//! If a jump process reduces the total number of 6-member rings (excludes 5-member rings) in a PAH below a certain threshold it is removed. Return this threshold in terms of the total number of 6-member rings.
inline double Sweep::Component::ThresholdOxidation() const {return m_thresholdOxidation;};

//! Allow PAHs in soot particles to point to the same memory location after a doubling event.
inline double Sweep::Component::SharedPointers() const {return m_sharedPointers;};

//! Allow particles composed of only single PAHs to be respresented with weighted particles
inline double Sweep::Component::WeightedPAHs() const { return m_weightedPAHs; };

//! Sets the density (g/cm3).
inline void Sweep::Component::SetDensity(const double dens) {m_density = dens;};

//! Sets the coalescence threshold 
inline void Sweep::Component::SetCoalescThresh(const double ct) {m_coalesc_thresh = ct;};

//! Sets the growhtfact 
inline void Sweep::Component::SetGrowthFact(const double gf) {m_growthfact = gf;};

//! Sets the minimum number of PAHs
inline void Sweep::Component::SetMinPAH(const int mp) {m_minPAH = mp;};

//! If a jump process reduces the total number of 6-member rings (excludes 5-member rings) in a PAH below a certain threshold it is removed. Sets this threshold in terms of the total number of 6-member rings.
inline void Sweep::Component::SetThresholdOxidation(const int to) {m_thresholdOxidation = to;};

//! Allow PAHs in soot particles to point to the same memory location after a doubling event.
inline void Sweep::Component::SetSharedPointers(const int sp) {m_sharedPointers = sp;};

//! Allow particles composed of only single PAHs to be respresented with weighted particles
inline void Sweep::Component::SetWeightedPAHs(const int wpah) { m_weightedPAHs = wpah; };

// COMPONENT NAME.

// Returns component symbol or name.
inline const std::string &Sweep::Component::Name() const {return m_name;};

// Sets the symbol or name.
inline void Sweep::Component::SetName(const std::string &name) {m_name = name;};

// PHASE

// Return phase
inline const std::string &Sweep::Component::Phase() const
{
	return m_phase;
};

// Set phase
inline void Sweep::Component::SetPhase(const std::string &phase)
{
	m_phase = phase;
};

// Return element
inline const std::string &Sweep::Component::Element() const
{
	return m_element;
};

// Set element
inline void Sweep::Component::SetElement(const std::string &element)
{
	m_element = element;
};

#endif
