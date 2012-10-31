/*
  Author(s):      Martin Martin (mm864)
  Project:        sprog (gas-phase and surface chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2012 Martin Martin.

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

#ifndef GPC_PHASE_INL_H
#define GPC_PHASE_INL_H
#include <vector>
#include <string>
#include "gpc_phase.h"

//  PHASE NAME.
inline const std::string &Phase::Name() const {return m_name;};

// PHASE ID
inline const std::string &Phase::ID() const {return m_id;};

// PHASE COMPOSITION.
inline const SpCompVector &Phase::Composition() const {return m_spcomp;};
inline unsigned int Phase::ComponentCount(void) const {return m_spcomp.size();};

// PHASE SITE DENSITY.
inline const double Phase::SiteDen() const {return m_siteden;};


// PARENT MECHANISM.
inline const Sprog::Mechanism *const Phase::Mechanism(void) const {return m_mech;};

// PHASE LOOKUP.
inline int Phase::Find(const std::string &name, const PhasePtrVector &list)
{
    // Loop over phase to find index.
    unsigned int i;
    for (i=0; i<list.size(); ++i) {
        if (*list[i] == name) {
            // Found phase!
            return i;
        }
    }

    // We are here because the phase wasn't found.
    return -1;
}

// CLONING.
inline Phase *const Phase::Clone(void) const {return new Phase(*this);};

#endif
