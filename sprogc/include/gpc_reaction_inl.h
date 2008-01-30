/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    Inline function definitions for the Reaction class.
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
