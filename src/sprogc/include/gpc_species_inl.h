/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Inline function definitions for Species class.

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

#ifndef GPC_SPECIES_INL_H
#define GPC_SPECIES_INL_H

/*
 * These function frames have been defined in gpc_species.h  
 */

//  SPECIES NAME.
inline const std::string &Species::Name() const {return m_name;};
inline const int Species::SiteOccupancy(void) const {return site_occupancy;};

// SPECIES PHASE NAME.
inline const std::string &Species::PhaseName() const {return m_phaseName;};

// SPECIES COMPOSITION.
inline const ElCompVector &Species::Composition() const {return m_elcomp;};
inline unsigned int Species::ComponentCount(void) const {return m_elcomp.size();};

// SPECIES MOLECULAR WEIGHT.
inline double Species::MolWt() const {return m_molwt;};

// THERMODYNAMIC FITTING PARAMETERS.
inline unsigned int Species::ThermoRangeCount(void) const {return m_thermoparams.size();};
inline void Species::SetThermoStartTemperature(const double T) {m_T1 = T;};

// PARENT MECHANISM.
inline const Sprog::Mechanism *const Species::Mechanism(void) const {return m_mech;};

// SPECIES LOOKUP.
inline int Species::Find(const std::string &name, const SpeciesPtrVector &list)
{
    // Loop over species to find index.
    unsigned int i;
    for (i=0; i<list.size(); ++i) {
        if (*list[i] == name) {
            // Found species!
            return i;
        }
    }

    // We are here because the species wasn't found.
    return -1;
}

// CLONING.
inline Species *const Species::Clone(void) const {return new Species(*this);};

#endif
