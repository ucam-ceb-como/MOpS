/*
  Author(s):      Martin Martin (mm864)
  Project:        sprog (gas-phase and surface chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2012 Martin Martin.

  File purpose:
    Implementation of the Stoichiometry class declared in the
    gpc_stoich.h header file.

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


#ifndef GPC_DELTA_STOICH_INL_H
#define GPC_DELTA_STOICH_INL_H
#include <vector>
#include <string>


// PHASE NAME 

inline const std::string &DeltaStoich::Name() const {return m_spName;};


// TOTAL STOICHIOMETRY FOR EACH SPECIES

inline const double &DeltaStoich::TotalStoich(void) const {return total_stoich_sp;};

// REACTANT STOICHIOMETRY FOR EACH SPECIES

inline const double &DeltaStoich::ReacStoich() const {return reac_stoich_sp;};

// PRODUCT STOICHIOMETRY FOR EACH SPECIES

inline const double &DeltaStoich::ProdStoich() const {return prod_stoich_sp;};

// CLONING.
inline DeltaStoich *const DeltaStoich::Clone(void) const {return new DeltaStoich(*this);};

#endif
