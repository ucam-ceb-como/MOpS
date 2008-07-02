/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Inline function definitions for the Reaction class.

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

#ifndef GPC_REACTION_INL_H
#define GPC_REACTION_INL_H

#include "gpc_stoich.h"
#include <vector>
#include <string>

class Reaction; // Forward declaration of Reaction class.

// REACTION NAME/DESRIPTION.
inline const std::string &Reaction::Name() const {return m_name;};
inline void Reaction::SetName(const std::string &name) {m_name = name;};

// REVERSIBLE REACTION?.
inline bool Reaction::IsReversible() const {return m_reversible;};
inline void Reaction::SetReversible(const bool isrev) {m_reversible=isrev;};

// REACTANTS.
inline const std::vector<Stoich> &Reaction::Reactants() const {return m_reac;};
inline const std::vector<Stoichf> &Reaction::FReactants() const {return m_freac;};

// PRODUCTS.
inline const std::vector<Stoich> &Reaction::Products() const {return m_prod;};
inline const std::vector<Stoichf> &Reaction::FProducts() const {return m_fprod;};

// STOICHIOMETRY.
inline real Reaction::TotalStoich() const {return m_dstoich;};
inline real Reaction::ReactantStoich() const {return m_dreac;};
inline real Reaction::ProductStoich() const {return m_dprod;};

// ARRHENIUS COEFFICIENTS.
inline const ARRHENIUS &Reaction::Arrhenius() const {return m_arrf;};

// Returns the explicit reverse Arrhenius parameters as a constant
// pointer.  If the reaction doesn't have explicit reverse parameters, then
// this function returns a NULL pointer.  This provides a way of checking
// for explicit reverse parameters.
inline const ARRHENIUS *const Reaction::RevArrhenius(void) const {return m_arrr;};


// LANDAU TELLER COEFFICIENTS.
inline const LTCOEFFS *const Reaction::LTCoeffs() const {return m_lt;};
inline const LTCOEFFS *const Reaction::RevLTCoeffs() const {return m_revlt;};


// THIRD BODIES.
inline bool Reaction::UseThirdBody() const {return m_usetb;};
inline void Reaction::SetUseThirdBody(bool usetb) {m_usetb = usetb;};

#endif
