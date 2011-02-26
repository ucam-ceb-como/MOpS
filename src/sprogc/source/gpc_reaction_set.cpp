/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ReactionSet class declared in the
    gpc_reaction_set.h header file.

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

#include "gpc_reaction_set.h"
#include "gpc_reaction.h"
#include "gpc_mech.h"
#include "gpc_stoich.h"
#include "gpc_idealgas.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <memory.h>
#include <iostream>
using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ReactionSet::ReactionSet()
{
    m_mech = NULL;
}

// Copy constructor.
ReactionSet::ReactionSet(const Sprog::Kinetics::ReactionSet &rxn)
{
    m_mech = NULL;
    *this = rxn;
}

// Stream-reading constructor.
ReactionSet::ReactionSet(std::istream &in)
{
    Deserialize(in);
}

// Destructor.
ReactionSet::~ReactionSet()
{
    releaseMemory();
}


// OPERATOR OVERLOADING.

// Assignment operator.
ReactionSet &ReactionSet::operator=(const ReactionSet &rxns)
{
    // Check for self assignment!
    if (this != &rxns) {
        // Clear current memory.
        releaseMemory();

        // Copy the reaction list.  Use the Clone() member function
        // to ensure reactions of the correct type are added.
        RxnPtrVector::const_iterator i;
        for (i=rxns.m_rxns.begin(); i!=rxns.m_rxns.end(); i++) {
            m_rxns.push_back((*i)->Clone());
        }

        // Build reversible reaction map.  Loop over incoming map to
        // get the reaction indices, but remember to use the pointers
        // to the new reactions!
        unsigned int j;
        RxnMap::const_iterator jrxn;
        for (jrxn=rxns.m_rev_rxns.begin(); jrxn!=rxns.m_rev_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_rev_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build forward Landau Teller reaction map.
        for (jrxn=rxns.m_lt_rxns.begin(); jrxn!=rxns.m_lt_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_lt_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build reverse Landau Teller reaction map.
        for (jrxn=rxns.m_revlt_rxns.begin(); jrxn!=rxns.m_revlt_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_revlt_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build third-body reaction map.
        for (jrxn=rxns.m_tb_rxns.begin(); jrxn!=rxns.m_tb_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_tb_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build fall-off reaction map.
        for (jrxn=rxns.m_fo_rxns.begin(); jrxn!=rxns.m_fo_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_fo_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }
    }

    return *this;
}

// Compound assignment operator:  Adds the contents of one reaction set
// to this one.
ReactionSet &ReactionSet::operator+=(const ReactionSet &rxns)
{
    // It is currently easier to not allow self-compounding here.
    if (this != &rxns) {
        // Save the current number of reactions, we'll need it later.
        int n = m_rxns.size();

        // Loop over all incoming reactions and copy them to the vector.
        RxnPtrVector::const_iterator i;
        for (i=rxns.m_rxns.begin(); i!= rxns.m_rxns.end(); i++) {
            m_rxns.push_back((*i)->Clone());
        }

        // Build reversible reaction map.  Loop over incoming map to
        // get the reaction indices, but remember to use the pointers
        // to the new reactions!
        unsigned int j;
        RxnMap::const_iterator jrxn;
        for (jrxn=rxns.m_rev_rxns.begin(); jrxn!=rxns.m_rev_rxns.end(); jrxn++) {
            j = n + (*jrxn).first; // Note n!
            m_rev_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build forward Landau Teller reaction map.
        for (jrxn=rxns.m_lt_rxns.begin(); jrxn!=rxns.m_lt_rxns.end(); jrxn++) {
            j = n + (*jrxn).first;
            m_lt_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build reverse Landau Teller reaction map.
        for (jrxn=rxns.m_revlt_rxns.begin(); jrxn!=rxns.m_revlt_rxns.end(); jrxn++) {
            j = n + (*jrxn).first;
            m_revlt_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build third-body reaction map.
        for (jrxn=rxns.m_tb_rxns.begin(); jrxn!=rxns.m_tb_rxns.end(); jrxn++) {
            j = n + (*jrxn).first;
            m_tb_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build fall-off reaction map.
        for (jrxn=rxns.m_fo_rxns.begin(); jrxn!=rxns.m_fo_rxns.end(); jrxn++) {
            j = n + (*jrxn).first;
            m_fo_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }
    }

    return *this;
}

// Addition operator:  Adds the contents of two reaction sets together.
const ReactionSet ReactionSet::operator+(const ReactionSet &rxns) const
{
    ReactionSet rs(*this);
    rs += rxns;
    return rs;
}

// Subscripting operator:  Provides a different way to access a particular
// reaction by index in the list.
Reaction *const ReactionSet::operator[](unsigned int i)
{
    return m_rxns.at(i);
}

// Subscripting operator:  Provides a different way to access a particular
// reaction by index in the list.
const Reaction *const ReactionSet::operator[](unsigned int i) const
{
    return m_rxns.at(i);
}


// REACTIONS.

// Returns the number of reactions in the set.
unsigned int ReactionSet::Count(void) const
{
    return m_rxns.size();
}

// Returns the vector of all reactions.
const RxnPtrVector &ReactionSet::Reactions() const
{
    return m_rxns;
}

// Returns a pointer to the ith reaction.  Returns NULL if i is invalid.
const Reaction *const ReactionSet::Reactions(unsigned int i) const
{
    if (i < m_rxns.size()) {
        return m_rxns[i];
    } else {
        return NULL;
    }
}

// Adds a reaction to the set.
Reaction *const ReactionSet::AddReaction(const Sprog::Kinetics::Reaction &rxn)
{
    // Clone the reaction and add it to the vector.
    Reaction *pr = rxn.Clone();
    m_rxns.push_back(pr);


    // Check for reverse parameters.
    if (rxn.RevArrhenius() != NULL) {
        m_rev_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    // Check for forward LT parameters.
    if (rxn.LTCoeffs() != NULL) {
        m_lt_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    // Check for reverse LT parameters.
    if (rxn.RevLTCoeffs() != NULL) {
        m_revlt_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    // Check for third-bodies.
    if (rxn.ThirdBodies().size() > 0) {
        m_tb_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    // Check for fall-off parameters.
    if (rxn.FallOffType() != None) {
        m_fo_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    return pr;
}


// TIDYING UP.

// Clears all reactions from the set.
void ReactionSet::Clear()
{
    releaseMemory();
}


// MOLAR PRODUCTION RATES RATES.

// Calculates the molar production rates of all species given the rate
// of progress of each reaction.
real ReactionSet::GetMolarProdRates(const fvector &rop,
                                    fvector &wdot) const 
{
    unsigned int k;
    const RxnStoichMap *mu;
    RxnStoichMap::const_iterator i;
    real wtot = 0.0;

    // Assign sufficient memory for output.
    wdot.resize(m_mech->SpeciesCount());

    // Loop over all species in mechanism.
    for (k=0; k!=m_mech->SpeciesCount(); ++k) {
        // Reset prod. rate of this species to zero.
        wdot[k] = 0.0;

        // Get the map of reaction stoichiometry for this species.
        mu = &m_mech->GetStoichXRef(k);

        // Sum up the production rate of this species from all reactions.
        //   first  = index of reaction.
        //   second = stoichiometry of species in reaction.
        for (i=mu->begin(); i!= mu->end(); i++) {
            wdot[k] += (*i).second * rop[(*i).first];
        }

        // Sum up total production rate.
        wtot += wdot[k];
    }

    return wtot;
}

// Calculates the molar production rates of all species.
real ReactionSet::GetMolarProdRates(const Sprog::Thermo::GasPhase &gas,
                                    fvector &wdot) const
{
    static fvector rop;
    GetRatesOfProgress(gas, rop);
    return GetMolarProdRates(rop, wdot);
}

// Calculates the molar production rates of all species.
real ReactionSet::GetMolarProdRates(real T, real density, const real *const x, 
                                    unsigned int n, 
                                    const Sprog::Thermo::ThermoInterface &thermo, 
                                    fvector &wdot) const
{
    static fvector rop;
    GetRatesOfProgress(T, density, x, n, thermo, rop);
    return GetMolarProdRates(rop, wdot);
}





// returns the molar production rate given the species mixture
void ReactionSet::GetMolarProdRates(Sprog::Thermo::Mixture &mix, fvector &wdot) const{		

   	fvector kfrwd,krev,rop,Gs;
	Sprog::Thermo::IdealGas ig(*mix.Species());
	ig.CalcGs_RT(mix.Temperature(),Gs);
 	GetRateConstants(mix.Temperature(),mix.Density(),&(mix.MoleFractions()[0]),m_mech->SpeciesCount(),Gs,kfrwd,krev);
	GetRatesOfProgress(mix.Density(),&(mix.MoleFractions()[0]),m_mech->SpeciesCount(),kfrwd,krev,rop);
	GetMolarProdRates(rop,wdot);
 

 
}


// REACTION RATES OF PROGRESS.

// Calculates the rate of progress of each reaction given the
// precalculated rate constants.
void ReactionSet::GetRatesOfProgress(real density,
                                     const real *const x,
                                     unsigned int n,
                                     const fvector &kforward,
                                     const fvector &kreverse,
                                     fvector &rop,
                                     fvector &rfwd,
                                     fvector &rrev) const
{
    unsigned int i;
    int j, k;

    // Resize output vector to sufficient length.
    rop.resize(m_rxns.size(), 0.0);

    if (n >= m_mech->Species().size()) {
        // Assign temp memory.
        rfwd.resize(m_rxns.size(), 0.0);
        rrev.resize(m_rxns.size(), 0.0);

        // Loop over all reactions.
        for (i=0; i!=m_rxns.size(); ++i) {
            // Use rfwd to store forward rates of production, 
            // and rrev to store reverse rates.
            rfwd[i] = kforward[i];
            rrev[i] = kreverse[i];

            // Integer reactants.
            for (k=0; k!=m_rxns[i]->ReactantCount(); ++k) {
                // As the stoichiometry is integer, it is more computationally efficient
                // to multiply the values together than to use the pow() function.
                for (j=0; j!=m_rxns[i]->Reactants()[k].Mu(); ++j) {
                    rfwd[i] *= density * x[m_rxns[i]->Reactants()[k].Index()];
                }
            }

            // Integer products.
            for (k=0; k!=m_rxns[i]->ProductCount(); ++k) {
                // As the stoichiometry is integer, it is more computationally efficient
                // to multiply the values together than to use the pow() function.
                for (j=0; j!=m_rxns[i]->Products()[k].Mu(); ++j) {
                    rrev[i] *= density * x[m_rxns[i]->Products()[k].Index()];
                }
            }

            // Real reactants.
            for (k=0; k!=m_rxns[i]->FReactantCount(); ++k) {
                // Now the stoichiometry is non-integer, we must use the pow() function.
                rfwd[i] *= pow(density * x[m_rxns[i]->FReactants()[k].Index()], 
                              m_rxns[i]->FReactants()[k].Mu()); 

            }

            // Real products.
            for (k=0; k!=m_rxns[i]->FProductCount(); ++k) {
                // Now the stoichiometry is non-integer, we must use the pow() function.
                rrev[i] *= pow(density * x[m_rxns[i]->FProducts()[k].Index()], 
                              m_rxns[i]->FProducts()[k].Mu()); 

            }

            // Calculate the net rates of production.
            rop[i] = rfwd[i] - rrev[i];
        }
    }
}

// Calculates the rate of progress of each reaction given the
// precalculated rate constants.
void ReactionSet::GetRatesOfProgress(real density,
                                     const real *const x,
                                     unsigned int n,
                                     const fvector &kforward,
                                     const fvector &kreverse,
                                     fvector &rop) const
{
    static fvector rfwd, rrev;
    GetRatesOfProgress(density, x, n, kforward, kreverse, rop, rfwd, rrev);
}

// Returns the rates of progress of all reactions given the mixture
// object.
void ReactionSet::GetRatesOfProgress(const Sprog::Thermo::GasPhase &mix, 
                                     const fvector &kforward, 
                                     const fvector &kreverse, 
                                     fvector &rop) const
{
    GetRatesOfProgress(mix.Density(), &(mix.MoleFractions()[0]), 
                       m_mech->Species().size(),
                       kforward, kreverse, rop);
}



// Calculates the rate of progress of each reaction.
void ReactionSet::GetRatesOfProgress(const Sprog::Thermo::GasPhase &gas, fvector &rop) const
{
    static fvector kf, kr;
    GetRateConstants(gas, kf, kr);
    GetRatesOfProgress(gas, kf, kr, rop);
}

// Calculates the rate of progress of each reaction.
void ReactionSet::GetRatesOfProgress(const Sprog::Thermo::GasPhase &gas,
                                     fvector &rop,
                                     fvector &rfwd,
                                     fvector &rrev) const
{
    static fvector kf, kr;
    GetRateConstants(gas, kf, kr);
    GetRatesOfProgress(gas.Density(), &(gas.MoleFractions()[0]),
                       m_mech->Species().size(),
                       kf, kr, rop, rfwd, rrev);
}

// Calculates the rate of progress of each reaction.
void ReactionSet::GetRatesOfProgress(real T, real density, const real *const x, 
                                     unsigned int n, 
                                     const Sprog::Thermo::ThermoInterface &thermo, 
                                     fvector &rop) const
{
    static fvector kf, kr;
    GetRateConstants(T, density, x, n, thermo, kf, kr);
    GetRatesOfProgress(density, x, n, kf, kr, rop);
}


// RATE CONSTANTS.

// Calculates the forward and reverse rate constants of all reactions
// given the mixture temperature, density and species mole fractions.
void ReactionSet::GetRateConstants(real T,
                                   real density,
                                   const real *const x,
                                   unsigned int n,
                                   const fvector &Gs,
                                   fvector &kf,
                                   fvector &kr) const 
{
    static fvector tbconcs;
    static bool fallocated = false;

    // Check that we have been given enough species concentrations.
    if (n < m_mech->Species().size()) {
        return;
    }

    // SETUP WORKSPACE.

    // Check that we have been given enough species concentrations.
    if (n < m_mech->Species().size()) {
        return;
    } else {
        if (!fallocated) {
            // Allocate temporary memory.
            tbconcs.resize(m_rxns.size(), 0.0);
            fallocated = true;
        }
        kf.resize(m_rxns.size(), 0.0);
        kr.resize(m_rxns.size(), 0.0);
    }

    // Calculate concentration-independent rate constants.
    calcRateConstantsT(T, Gs, kf, kr);

    // Calculate third-body concentrations for all reactions.  These
    // values will be multiplied by the rate constants, therefore if
    // a reaction does not have third-bodies then tbconcs is set to 1.0.
    calcTB_Concs(density, x, n, tbconcs);

    // Calculate the pressure-dependent fall-off terms in the rate
    // constant expressions.  This function multiplies the rate constants
    // by the fall-off terms.  This function may also change the values in
    // the tbconcs vector.
    calcFallOffTerms(T, density, x, n, tbconcs, kf, kr);

    // Apply third-body concentrations to rate constants.  It is important
    // to apply the fall-off terms before doing this, as that routine may
    // change the values in the tbconcs vector.
    for (RxnMap::const_iterator im=m_tb_rxns.begin(); im!=m_tb_rxns.end(); ++im) {
        unsigned int j = (*im).first;
        kf[j] *= tbconcs[j];
        kr[j] *= tbconcs[j];
    }
}

// Calculates the forward and reverse rate constants 
// of all reactions given a mixture object.
void ReactionSet::GetRateConstants(const Sprog::Thermo::GasPhase &mix, 
                                   const std::vector<real> &Gs,
                                   std::vector<real> &kforward, 
                                   std::vector<real> &kreverse) const
{
    GetRateConstants(mix.Temperature(), mix.Density(), &(mix.MoleFractions()[0]), 
                     m_mech->Species().size(), Gs, kforward, kreverse);
}

// Calculates the forward and reverse rate constants of all reactions
// given the mixture temperature, density and species mole fractions.
void ReactionSet::GetRateConstants(real T,
                                   real density,
                                   const real *const x,
                                   unsigned int n,
                                   const Sprog::Thermo::ThermoInterface &thermo,
                                   fvector &kforward,
                                   fvector &kreverse) const 
{
    static fvector Gs;
    thermo.CalcGs_RT(T, Gs);
    GetRateConstants(T, density, x, n, Gs, kforward, kreverse);
}

// Calculates the forward and reverse rate constants 
// of all reactions given a mixture object.
void ReactionSet::GetRateConstants(const Sprog::Thermo::GasPhase &mix, 
                                   std::vector<real> &kforward, 
                                   std::vector<real> &kreverse) const
{
    static fvector Gs;
    mix.Gs_RT(Gs);
    GetRateConstants(mix.Temperature(), mix.Density(), &(mix.MoleFractions()[0]), 
                     m_mech->Species().size(), Gs, kforward, kreverse);
}


// Calculates the concentration-independent portions
// of the rates constants.
void ReactionSet::calcRateConstantsT(real T, const fvector &Gs, 
                                     fvector &kf, fvector &kr) const
{
    real lnT=0.0, invRT=0.0, Patm_RT=0.0, T_1_3=0.0, T_2_3=0.0;
    RxnPtrVector::const_iterator i;
    RxnMap::const_iterator im;
    int j, k;

    // Precalculate some temperature parameters.
    lnT = log(T);
    switch (m_mech->Units()) {
        case SI :
            invRT = 1.0 / (R * T);
            Patm_RT = 101325.0 * invRT;
            break;
        case CGS :
            invRT = 1.0 / (R_CGS * T);
            Patm_RT = 1013250.0 * invRT;
            break;
        default:
            // Something has gone wrong to end up here.
            invRT = 0.0;
            Patm_RT = 0.0;
    }
    if (m_lt_rxns.size() > 0) {
        T_1_3 = 1.0 / pow(T, ONE_THIRD);
        T_2_3 = T_1_3 * T_1_3;
    }

    // Calculate classic Arrhenius forward rate expression.
    for (i=m_rxns.begin(),j=0; i!=m_rxns.end(); ++i,++j) {
        kf[j] = (*i)->Arrhenius().A * 
                      exp(((*i)->Arrhenius().n * lnT) - 
                          ((*i)->Arrhenius().E * invRT));
    }

    // Landau-Teller rate expressions.
    for (im=m_lt_rxns.begin(); im!=m_lt_rxns.end(); ++im) {
        j = (*im).first;
        kf[j] *= exp(((*im).second->LTCoeffs()->B / T_1_3) + 
                           ((*im).second->LTCoeffs()->C / T_2_3));
    }

    // Reverse rate constants.
    for (i=m_rxns.begin(),j=0; i!=m_rxns.end(); ++i,++j) {
        if ((*i)->RevArrhenius() != NULL) {
            // This reaction has explicit reverse rate parameters.
            kr[j] = (*i)->RevArrhenius()->A * 
                          exp(((*i)->RevArrhenius()->n * lnT) - 
                              ((*i)->RevArrhenius()->E * invRT));
        } else if ((*i)->IsReversible()) {
            // Must find the reverse rate constants by equilibrium.
            kr[j] = 0.0;

            // Calculate the Gibbs free energy change for reaction i and sum up
            // the stoichiometric coefficients.
            for (k=0; k!=(*i)->ReactantCount(); ++k) {
                // Integer Reactants.
                kr[j] += (*i)->Reactants()[k].Mu() * Gs[(*i)->Reactants()[k].Index()];
            }
            for (k=0; k!=(*i)->ProductCount(); ++k) {
                // Integer Products.
                kr[j] -= (*i)->Products()[k].Mu() * Gs[(*i)->Products()[k].Index()];
            }
            for (k=0; k!=(*i)->FReactantCount(); ++k) {
                // Real Reactants.
                kr[j] += (*i)->FReactants()[k].Mu() * Gs[(*i)->FReactants()[k].Index()];
            }
            for (k=0; k!=(*i)->FProductCount(); ++k) {
                // Real Products.
                kr[j] -= (*i)->Products()[k].Mu() * Gs[(*i)->Products()[k].Index()];
            }

            // Calculate the reverse rate constant.
            kr[j]  = exp(min(kr[j], log(1.0e250)));
            kr[j] *= pow(Patm_RT, (*i)->TotalStoich());
            kr[j]  = kf[j] / max(kr[j], 1.0e-250);
        }
    }

    // Explicit reverse Landau-Teller parameters.
    for (im=m_revlt_rxns.begin(); im!=m_revlt_rxns.end(); ++im) {
        j = (*im).first;
        kr[j] *= exp(((*im).second->RevLTCoeffs()->B / T_1_3) + 
                           ((*im).second->RevLTCoeffs()->C / T_2_3));
    }
}

// Calculates third-body concentrations for all reactions.  These
// values will be multiplied by the rate constants, therefore if
// a reaction does not have third-bodies the tbconcs is set to 1.0.
void ReactionSet::calcTB_Concs(real density, const real *const x, 
                               unsigned int n, fvector &tbconcs) const
{
    RxnPtrVector::const_iterator i;
    RxnMap::const_iterator im;
    int j, k;

    // Third body concentrations for each reaction.
    for (i=m_rxns.begin(),j=0; i!=m_rxns.end(); ++i,++j) {
        if ((*i)->UseThirdBody()) {
            // Calculate enhanced third body concentration using the enhancement
            // factors defined for this reaction.
            tbconcs[j] = 0.0;
            for (k=0; k<(*i)->ThirdBodyCount(); ++k) {
                tbconcs[j] += ((*i)->ThirdBody(k).Mu() - 1.0) * 
                              x[(*i)->ThirdBody(k).Index()];
            }
            tbconcs[j] = density * (1.0 + tbconcs[j]);
        } else {
            // This reaction has no third body requirement.
            tbconcs[j] = 1.0;
        }
    }
}

// Calculates the pressure-dependent fall-off terms in the rate
// constant expressions.  This function multiplies the rate constants
// by the fall-off terms.  This function may also change the values in
// the tbconcs vector.
void ReactionSet::calcFallOffTerms(real T, real density, const real *const x, 
                                   unsigned int n, fvector &tbconcs, 
                                   fvector &kf, fvector &kr) const
{
    real lowk=0.0, pr=0.0, logpr=0.0, lnT=log(T), invRT=0.0;
    RxnPtrVector::const_iterator i;
    RxnMap::const_iterator im;
    int j;
    const Reaction *rxn;

    switch (m_mech->Units()) {
        case SI :
            invRT = 1.0 / (R * T);
            break;
        case CGS :
            invRT = 1.0 / (R_CGS * T);
            break;
        default:
            // Something has gone wrong to end up here.
            invRT = 0.0;
    }

    // Pressure dependent fall-off reactions.
    for (im=m_fo_rxns.begin(); im!=m_fo_rxns.end(); ++im) {
        j   = (*im).first;  // Reaction index of imth fall-off reaction.
        rxn = (*im).second;

        // Calculate low pressure limit.
        lowk = rxn->FallOffParams().LowP_Limit.A * 
               exp((rxn->FallOffParams().LowP_Limit.n * lnT) - 
                    (rxn->FallOffParams().LowP_Limit.E * invRT));

        // Calculate reduced pressure.
        if (rxn->FallOffParams().ThirdBody >= 0) {
            // A particular species is to be used as the third body.
            pr = lowk * density * x[rxn->FallOffParams().ThirdBody] / kf[j];
        } else {
            // Use all species as third bodies.
            pr = lowk * tbconcs[j] / kf[j];
            tbconcs[j] = 1.0;        }

        // Calculate rate constants based on equation form.
        logpr = log10(pr);
        pr    = pr / (1.0 + pr);
        switch (rxn->FallOffType()) {
            case Troe3: // 3-parameter Troe form.
                pr *= rxn->FTROE3(T, logpr);
                kf[j] *= pr;
                kr[j] *= pr;
                break;
            case Troe4: // 4-parameter Troe form.
                pr *= rxn->FTROE4(T, logpr);
                kf[j] *= pr;
                kr[j] *= pr;
                break;
            case SRI: // SRI form.
                pr *= rxn->FSRI(T, logpr);
                kf[j] *= pr;
                kr[j] *= pr;
                break;
            case Custom: // A custom function is defined to calculate the fall-off form.
                rxn->FallOffFn()(*rxn, lowk, tbconcs[j], T, kf[j], kr[j]);
                break;
            case Lindemann:
                // Need do nothing for Lindemann form.
                break;
            default:
                break;
        }
    }
}


// JACOBIAN EVALUATION.

// Calculates the Jacobian matrix for a constant volume, adiabatic
// homogeneous mixture. J[j][i] is the Jacobian entry for variable
// i with respect to i: dFi/dYj.  It is assumed that the Jacobian
// matrix array J has already been allocated for NSP+2 variables
// (all species, temperature and density).
void ReactionSet::CalcJacobian(real T, real density, real *const x,
                               unsigned int n, 
                               const Sprog::Thermo::ThermoInterface &thermo, 
                               real pfac, real **J,
                               bool constV, bool constT) const
{
    bool fallocated=false;
    fvector tbconcs, kfT, krT, kf, kr, Gs, rop0, 
                   rop1, wdot0, wdot1, Hs, xdot0, xdot1;
    real wtot0=0.0, wtot1=0.0, //invrho=1.0 / density, 
         Tdot0=0.0, Tdot1=0.0, Cp=0.0, xsave=0.0, dx=0.0, invdx=0.0;

    // SETUP WORKSPACE.

    // Check that we have been given enough species concentrations.
    if (n < m_mech->Species().size()) {
        return;
    } else {
        if (!fallocated) {
            // Allocate temporary memory.
            tbconcs.resize(m_rxns.size(), 0.0);
            kfT.resize(m_rxns.size(), 0.0);
            krT.resize(m_rxns.size(), 0.0);
            kf.resize(m_rxns.size(), 0.0);
            kr.resize(m_rxns.size(), 0.0);
            Gs.resize(m_mech->SpeciesCount(), 0.0);
            rop0.resize(m_rxns.size(), 0.0);
            rop1.resize(m_rxns.size(), 0.0);
            wdot0.resize(m_mech->SpeciesCount(), 0.0);
            wdot1.resize(m_mech->SpeciesCount(), 0.0);
            xdot0.resize(m_mech->SpeciesCount(), 0.0);
            xdot1.resize(m_mech->SpeciesCount(), 0.0);
            fallocated = true;
        }
    }

    // CALCULATE UNPERTURBED VALUES.

    // Calculate Gibbs free energies.
    thermo.CalcGs_RT(T, Gs);

    // Calculate concentration-independent rate constants.
    calcRateConstantsT(T, Gs, kfT, krT);
    
    // Copy conc-independent rate constants into total rate
    // constant vectors.  As the vectors are equal lengths we
    // can just use the memcpy function.
    memcpy(&kf[0], &kfT[0], sizeof(real)*m_rxns.size());
    memcpy(&kr[0], &krT[0], sizeof(real)*m_rxns.size());

    // Calculate third-body concentrations for all reactions.  These
    // values will be multiplied by the rate constants, therefore if
    // a reaction does not have third-bodies then tbconcs is set to 1.0.
    calcTB_Concs(density, x, n, tbconcs);

    // Calculate the pressure-dependent fall-off terms in the rate
    // constant expressions.  This function multiplies the rate constants
    // by the fall-off terms.  This function may also change the values in
    // the tbconcs vector.
    calcFallOffTerms(T, density, x, n, tbconcs, kf, kr);

    // Apply third-body concentrations to rate constants.  It is important
    // to apply the fall-off terms before doing this, as that routine may
    // change the values in the tbconcs vector.
    for (RxnMap::const_iterator im=m_tb_rxns.begin(); im!=m_tb_rxns.end(); ++im) {
        unsigned int j = (*im).first;
        kf[j] *= tbconcs[j];
        kr[j] *= tbconcs[j];
    }

    // Calculate unperturbed reaction rates-of-progress and molar
    // production rates of species.
    GetRatesOfProgress(density, x, n, kf, kr, rop0);
    wtot0 = GetMolarProdRates(rop0, wdot0);

    // Copy rates-of-progess to working vector using memcpy.
    memcpy(&rop1[0], &rop0[0], sizeof(real)*m_rxns.size());

    // Calculate unperturbed dx/dt.
    for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
        xdot0[i] = ((wdot0[i] - (x[i]*wtot0)) / density);
    }

    // Calculate unperturbed temperature term dT/dt.
    Tdot0 = 0.0;
    if (!constT) {
        if (constV) {
            // Use internal energies for constant volume.
            thermo.CalcUs_RT(T, Hs);
            Cp = thermo.CalcBulkCv_R(T, x, n);
        } else {
            // Use enthalpies for constant pressure.
            thermo.CalcHs_RT(T, Hs);
            Cp = thermo.CalcBulkCp_R(T, x, n);
        }
        for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
            Tdot0 += wdot0[i] * Hs[i];
        }
        Tdot0 *= - T / (density * Cp);
    }

    // FINITE DIFFERENCING W.R.T. SPECIES MOLE FRACTIONS.

    dx = 0.0; invdx=0.0;
    for (unsigned int k=0; k!=m_mech->SpeciesCount(); ++k) {
        // PERTURB VARIABLE.

        // Perturb mole fraction of species k.
        xsave = x[k];
        // dx    = sqrt(1.0e-16 * max(1.0e-8, abs(xsave)));
        dx    = max(sqrt(pfac) * abs(xsave), 3.6390968218251355e-021);
        invdx = 1.0 / dx;
        x[k] += dx;

        // RECALCULATE RATE EXPRESSIONS.

        // Copy back in the concentration-independent parts of
        // the rate constant expressions.
        memcpy(&kf[0], &kfT[0], sizeof(real)*m_rxns.size());
        memcpy(&kr[0], &krT[0], sizeof(real)*m_rxns.size());

        // Calculate concentration-dependent terms.
        // TODO:  Optimise this to only calculate for affected reactions.
        calcTB_Concs(density, x, n, tbconcs);
        calcFallOffTerms(T, density, x, n, tbconcs, kf, kr);
        for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); ++i) {
            unsigned int j = i->first;
            kf[j] *= tbconcs[j];
            kr[j] *= tbconcs[j];
        }

        // RECALCULATE RATES-OF-PROGRESS.

        // Copy unperturbed values.
        memcpy(&rop1[0], &rop0[0], sizeof(real)*m_rxns.size());
        // memcpy(&wdot1[0], &wdot0[0], sizeof(real)*m_mech->SpeciesCount());

        // Recalculate third-body and fall-off reaction rates-of-progress.
        for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); ++i) {
            unsigned int j = i->first;
            rop1[j] = (i->second)->RateOfProgress(density, x, n, kf[j], kr[j]);
        }
        for (RxnMap::const_iterator i=m_fo_rxns.begin(); i!=m_fo_rxns.end(); ++i) {
            unsigned int j = i->first;
            rop1[j] = (i->second)->RateOfProgress(density, x, n, kf[j], kr[j]);
        }

        // Recalculate other reaction rates-of-progress for which
        // the kth species is a reactant.
        const RxnStoichMap &mu = m_mech->GetStoichXRef(k);
        for (RxnStoichMap::const_iterator i=mu.begin(); i!=mu.end(); ++i) {
            unsigned int j = i->first;
            rop1[j] = m_rxns[j]->RateOfProgress(density, x, n, kf[j], kr[j]);
        }

        // RECALCULATE MOLAR PRODUCTION RATES.

        // Calculate perturbed molar production rates.
        // TODO:  Can this be optimised to avoid looping over all species?
        wtot1 = GetMolarProdRates(rop1, wdot1);

        // Calculate dx/dt.
        for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
            xdot1[i] = ((wdot1[i] - (x[i]*wtot1)) / density);
        }

        // Calculate temperature term dT/dt.
        Tdot1 = 0.0;
        if (!constT) {
            for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
                Tdot1 += wdot1[i] * Hs[i];
            }
            Tdot1 *= - T / (density * Cp);
        }

        // CALCULATE JACOBIAN ENTRIES.

        // Calculate Jacobian entries for species mole fractions.
        for (unsigned int j=0; j!=m_mech->SpeciesCount(); ++j) {
            J[k][j] = invdx * (xdot1[j] - xdot0[j]);
        }

        // Calculate temperature Jacobian entry.
        J[k][m_mech->SpeciesCount()] = (Tdot1 - Tdot0) * invdx;

        // Calculate density Jacobian entry.
        if (constV) {
            J[k][m_mech->SpeciesCount()+1] = (wtot1 - wtot0) * invdx;
        } else {
            J[k][m_mech->SpeciesCount()+1] = 0.0;
        }

        // UNPERTURB.

        // Put species mole fraction back to correct value.
        x[k] = xsave;

    } // (Loop over all species).


    // FINITE DIFFERENCING W.R.T. DENSITY.
    // Density affects (almost) all reaction rates.

    // Perturb density.
    real Dpert = density;
    // dx    = sqrt(1.0e-16 * max(1.0e-8, abs(density)));
    dx    = sqrt(pfac) * abs(density);
    invdx = 1.0 / dx;
    Dpert += dx;

    // Copy conc-independent rate constants into total rate
    // constant vectors.  As the vectors are equal lengths we
    // can just use the memcpy function.
    memcpy(&kf[0], &kfT[0], sizeof(real)*m_rxns.size());
    memcpy(&kr[0], &krT[0], sizeof(real)*m_rxns.size());

    // Recalculate all rate constants.
    calcTB_Concs(Dpert, x, n, tbconcs);
    calcFallOffTerms(T, Dpert, x, n, tbconcs, kf, kr);
    for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); ++i) {
        unsigned int j = i->first;
        kf[j] *= tbconcs[j];
        kr[j] *= tbconcs[j];
    }

    // Recalculate reaction rates-of-progress and molar
    // production rates of species.
    GetRatesOfProgress(Dpert, x, n, kf, kr, rop1);
    wtot1 = GetMolarProdRates(rop1, wdot1);

    // Calculate dx/dt.
    for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
        xdot1[i] = ((wdot1[i] - (x[i]*wtot1)) / Dpert);
    }

    // Calculate temperature term dT/dt.
    Tdot1 = 0.0;
    if (!constT) {
        for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
            Tdot1 += wdot1[i] * Hs[i];
        }
        Tdot1 *= - T / (Dpert * Cp);
    }

    // CALCULATE JACOBIAN ENTRIES.

    // Calculate Jacobian entries for species mole fractions.
    for (unsigned int j=0; j!=m_mech->SpeciesCount(); ++j) {
        J[m_mech->SpeciesCount()+1][j] = invdx * (xdot1[j] - xdot0[j]);
    }

    // Calculate temperature Jacobian entry.
    J[m_mech->SpeciesCount()+1][m_mech->SpeciesCount()] = (Tdot1 - Tdot0) * invdx;

    // Calculate density Jacobian entry.
    if (constV) {
        J[m_mech->SpeciesCount()+1][m_mech->SpeciesCount()+1] = (wtot1 - wtot0) * invdx;
    } else {
        J[m_mech->SpeciesCount()+1][m_mech->SpeciesCount()+1] = 0.0;
    }


    // FINITE DIFFERENCING W.R.T. TEMPERATURE.
    // Temperature affects all reaction rates.

    // Perturb temperature.
    real Tpert = T;
    // dx    = sqrt(1.0e-16 * max(1.0e-8, abs(T)));
    dx    = sqrt(pfac) * abs(T);
    invdx = 1.0 / dx;
    Tpert += dx;

    // Recalculate all rate constants.
    thermo.CalcGs_RT(Tpert, Gs);
    calcRateConstantsT(Tpert, Gs, kf, kr);
    calcTB_Concs(density, x, n, tbconcs);
    calcFallOffTerms(Tpert, density, x, n, tbconcs, kf, kr);
    for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); ++i) {
        unsigned int j = i->first;
        kf[j] *= tbconcs[j];
        kr[j] *= tbconcs[j];
    }

    // Recalculate reaction rates-of-progress and molar
    // production rates of species.
    GetRatesOfProgress(density, x, n, kf, kr, rop1);
    wtot1 = GetMolarProdRates(rop1, wdot1);

    // Calculate dx/dt.
    for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
        xdot1[i] = ((wdot1[i] - (x[i]*wtot1)) / density);
    }

    // Calculate temperature term dT/dt.
    Tdot1 = 0.0;
    if (!constT) {
        if (constV) {
            // Use internal energies for constant volume.
            thermo.CalcUs_RT(Tpert, Hs);
            Cp = thermo.CalcBulkCv_R(Tpert, x, n);
        } else {
            // Use enthalpies for constant pressure.
            thermo.CalcHs_RT(Tpert, Hs);
            Cp = thermo.CalcBulkCp_R(Tpert, x, n);
        }
        for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
            Tdot1 += wdot1[i] * Hs[i];
        }
        Tdot1 *= - Tpert / (density * Cp);
    }

    // CALCULATE JACOBIAN ENTRIES.

    // Calculate Jacobian entries for species mole fractions.
    for (unsigned int j=0; j!=m_mech->SpeciesCount(); ++j) {
        J[m_mech->SpeciesCount()][j] = invdx * (xdot1[j] - xdot0[j]);
    }

    // Calculate temperature Jacobian entry.
    J[m_mech->SpeciesCount()][m_mech->SpeciesCount()] = (Tdot1 - Tdot0) * invdx;

    // Calculate density Jacobian entry.
    if (constV) {
        J[m_mech->SpeciesCount()][m_mech->SpeciesCount()+1] = (wtot1 - wtot0) * invdx;
    } else {
        J[m_mech->SpeciesCount()][m_mech->SpeciesCount()+1] = 0.0;
    }
}

/*!
The Jacobian matrix is originally calculated for [dOmegai/dxj]/rho
For LOI it is needed to be [dOmegai/dCj], therefore a change of basis by 
vector multiplication is needed and is implemented below. It is assumed that this
will be deleted when the Jacobian is implemented in the correct manner.

@param[in]      T           The mixture temperature.
@param[in]      density     Mixture molar density.
@param[in]      x           Species mole fractions.
@param[in]      n           Number of values in x array.
@param[in]      thermo      Thermodynamics interface.
@param[in]      pfac        Perturbation factor for calculating J entries.
@param[in, out] J           Jacobian matrix array.
@param[in]      constV      Volume (constant).
@param[in]      constT      Temperature (constant).
*/
void ReactionSet::RateJacobian(
        real T,            
        real density,      
        real *const x,    
        unsigned int n,   
        const Sprog::Thermo::ThermoInterface &thermo, 
        real pfac,         
        real **J,          
        bool constV, bool constT
        ) const
{
    bool fallocated=false;
    fvector tbconcs, kfT, krT, kf, kr, rop0, rop1, wdot0, wdot1, Gs, Hs;
    real wtot0=0.0, wtot1=0.0, 
         Tdot0=0.0, Tdot1=0.0, Cp=0.0, xsave=0.0, dconc=0.0, invdconc=0.0;

    // SETUP WORKSPACE.

    // Check that we have been given enough species concentrations.
    if (n < m_mech->Species().size()) {
        return;
    } else {
        if (!fallocated) {
            // Allocate temporary memory.
            tbconcs.resize(m_rxns.size(), 0.0);
            kfT.resize(m_rxns.size(), 0.0);
            krT.resize(m_rxns.size(), 0.0);
            kf.resize(m_rxns.size(), 0.0);
            kr.resize(m_rxns.size(), 0.0);
            Gs.resize(m_mech->SpeciesCount(), 0.0);
            rop0.resize(m_rxns.size(), 0.0);
            rop1.resize(m_rxns.size(), 0.0);
            wdot0.resize(m_mech->SpeciesCount(), 0.0);
            wdot1.resize(m_mech->SpeciesCount(), 0.0);
            fallocated = true;
        }
    }

    // CALCULATE UNPERTURBED VALUES.

    // Calculate Gibbs free energies.
    thermo.CalcGs_RT(T, Gs);

    // Calculate concentration-independent rate constants.
    calcRateConstantsT(T, Gs, kfT, krT);
    
    // Copy conc-independent rate constants into total rate
    // constant vectors.  As the vectors are equal lengths we
    // can just use the memcpy function.
    memcpy(&kf[0], &kfT[0], sizeof(real)*m_rxns.size());
    memcpy(&kr[0], &krT[0], sizeof(real)*m_rxns.size());

    // Calculate third-body concentrations for all reactions.  These
    // values will be multiplied by the rate constants, therefore if
    // a reaction does not have third-bodies then tbconcs is set to 1.0.
    calcTB_Concs(density, x, n, tbconcs);

    // Calculate the pressure-dependent fall-off terms in the rate
    // constant expressions.  This function multiplies the rate constants
    // by the fall-off terms.  This function may also change the values in
    // the tbconcs vector.
    calcFallOffTerms(T, density, x, n, tbconcs, kf, kr);

    // Apply third-body concentrations to rate constants.  It is important
    // to apply the fall-off terms before doing this, as that routine may
    // change the values in the tbconcs vector.
    for (RxnMap::const_iterator im=m_tb_rxns.begin(); im!=m_tb_rxns.end(); ++im) {
        unsigned int j = (*im).first;
        kf[j] *= tbconcs[j];
        kr[j] *= tbconcs[j];
    }

    // Calculate unperturbed reaction rates-of-progress and molar
    // production rates of species.
    GetRatesOfProgress(density, x, n, kf, kr, rop0);

    GetMolarProdRates(rop0, wdot0);

    // Copy rates-of-progess to working vector using memcpy.
    memcpy(&rop1[0], &rop0[0], sizeof(real)*m_rxns.size());

    // Calculate unperturbed temperature term dT/dt.
    Tdot0 = 0.0;
    if (!constT) {
        if (constV) {
            // Use internal energies for constant volume.
            thermo.CalcUs_RT(T, Hs);
            Cp = thermo.CalcBulkCv_R(T, x, n);
        } else {
            // Use enthalpies for constant pressure.
            thermo.CalcHs_RT(T, Hs);
            Cp = thermo.CalcBulkCp_R(T, x, n);
        }
        for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
            Tdot0 += wdot0[i] * Hs[i];
        }
        Tdot0 *= - T / (density * Cp);
    }

    // FINITE DIFFERENCING W.R.T. SPECIES MOLE CONCENTRATIONS.

    dconc = 0.0; invdconc=0.0;
    for (unsigned int k=0; k!=m_mech->SpeciesCount(); ++k) {
        // PERTURB VARIABLE.

        // Perturb mole fraction of species k.
        xsave = x[k];
        dconc    = max(sqrt(pfac) * abs(xsave), 3.6390968218251355e-021);
        invdconc = 1.0 / dconc;
        x[k] += dconc;

        // RECALCULATE RATE EXPRESSIONS.

        // Copy back in the concentration-independent parts of
        // the rate constant expressions.
        memcpy(&kf[0], &kfT[0], sizeof(real)*m_rxns.size());
        memcpy(&kr[0], &krT[0], sizeof(real)*m_rxns.size());

        // Calculate concentration-dependent terms.
        // TODO:  Optimise this to only calculate for affected reactions.
        calcTB_Concs(density, x, n, tbconcs);
        calcFallOffTerms(T, density, x, n, tbconcs, kf, kr);
        for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); ++i) {
            unsigned int j = i->first;
            kf[j] *= tbconcs[j];
            kr[j] *= tbconcs[j];
        }

        // RECALCULATE RATES-OF-PROGRESS.

        // Copy unperturbed values.
        memcpy(&rop1[0], &rop0[0], sizeof(real)*m_rxns.size());
        // memcpy(&wdot1[0], &wdot0[0], sizeof(real)*m_mech->SpeciesCount());

        // Recalculate third-body and fall-off reaction rates-of-progress.
        for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); ++i) {
            unsigned int j = i->first;
            rop1[j] = (i->second)->RateOfProgress(density, x, n, kf[j], kr[j]);
        }
        for (RxnMap::const_iterator i=m_fo_rxns.begin(); i!=m_fo_rxns.end(); ++i) {
            unsigned int j = i->first;
            rop1[j] = (i->second)->RateOfProgress(density, x, n, kf[j], kr[j]);
        }

        // Recalculate other reaction rates-of-progress for which
        // the kth species is a reactant.
        const RxnStoichMap &mu = m_mech->GetStoichXRef(k);
        for (RxnStoichMap::const_iterator i=mu.begin(); i!=mu.end(); ++i) {
            unsigned int j = i->first;
            rop1[j] = m_rxns[j]->RateOfProgress(density, x, n, kf[j], kr[j]);
        }

        // RECALCULATE MOLAR PRODUCTION RATES.

        // Calculate perturbed molar production rates.
        // TODO:  Can this be optimised to avoid looping over all species?
        GetMolarProdRates(rop1, wdot1);

        // Calculate temperature term dT/dt.
        Tdot1 = 0.0;
        if (!constT) {
            for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
                Tdot1 += wdot1[i] * Hs[i];
            }
            Tdot1 *= - T / (density * Cp);
        }

        // CALCULATE JACOBIAN ENTRIES.

        // Calculate Jacobian entries for species molar concentrations.
        for (unsigned int j=0; j!=m_mech->SpeciesCount(); ++j) {
            J[k][j] = invdconc * (wdot1[j] - wdot0[j]);
        }

        // Calculate temperature Jacobian entry.
        J[k][m_mech->SpeciesCount()] = (Tdot1 - Tdot0) * invdconc;

        // Calculate density Jacobian entry.
        if (constV) {
            J[k][m_mech->SpeciesCount()+1] = (wtot1 - wtot0) * invdconc;
        } else {
            J[k][m_mech->SpeciesCount()+1] = 0.0;
        }

        // UNPERTURB.

        // Put species mole fraction back to correct value.
        x[k] = xsave;

    } // (Loop over all species).


    // FINITE DIFFERENCING W.R.T. DENSITY.
    // Density affects (almost) all reaction rates.

    // Perturb density.
    real Dpert = density;
    // dconc    = sqrt(1.0e-16 * max(1.0e-8, abs(density)));
    dconc    = sqrt(pfac) * abs(density);
    invdconc = 1.0 / dconc;
    Dpert += dconc;

    // Copy conc-independent rate constants into total rate
    // constant vectors.  As the vectors are equal lengths we
    // can just use the memcpy function.
    memcpy(&kf[0], &kfT[0], sizeof(real)*m_rxns.size());
    memcpy(&kr[0], &krT[0], sizeof(real)*m_rxns.size());

    // Recalculate all rate constants.
    calcTB_Concs(Dpert, x, n, tbconcs);
    calcFallOffTerms(T, Dpert, x, n, tbconcs, kf, kr);
    for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); ++i) {
        unsigned int j = i->first;
        kf[j] *= tbconcs[j];
        kr[j] *= tbconcs[j];
    }

    // Recalculate reaction rates-of-progress and molar
    // production rates of species.
    GetRatesOfProgress(Dpert, x, n, kf, kr, rop1);
    wtot1 = GetMolarProdRates(rop1, wdot1);

    // Calculate temperature term dT/dt.
    Tdot1 = 0.0;
    if (!constT) {
        for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
            Tdot1 += wdot1[i] * Hs[i];
        }
        Tdot1 *= - T / (Dpert * Cp);
    }

    // CALCULATE JACOBIAN ENTRIES.

    // Calculate temperature Jacobian entry.
    J[m_mech->SpeciesCount()+1][m_mech->SpeciesCount()] = (Tdot1 - Tdot0) * invdconc;

    // Calculate density Jacobian entry.
    if (constV) {
        J[m_mech->SpeciesCount()+1][m_mech->SpeciesCount()+1] = (wtot1 - wtot0) * invdconc;
    } else {
        J[m_mech->SpeciesCount()+1][m_mech->SpeciesCount()+1] = 0.0;
    }


    // FINITE DIFFERENCING W.R.T. TEMPERATURE.
    // Temperature affects all reaction rates.

    // Perturb temperature.
    real Tpert = T;
    // dconc    = sqrt(1.0e-16 * max(1.0e-8, abs(T)));
    dconc    = sqrt(pfac) * abs(T);
    invdconc = 1.0 / dconc;
    Tpert += dconc;

    // Recalculate all rate constants.
    thermo.CalcGs_RT(Tpert, Gs);
    calcRateConstantsT(Tpert, Gs, kf, kr);
    calcTB_Concs(density, x, n, tbconcs);
    calcFallOffTerms(Tpert, density, x, n, tbconcs, kf, kr);
    for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); ++i) {
        unsigned int j = i->first;
        kf[j] *= tbconcs[j];
        kr[j] *= tbconcs[j];
    }

    // Recalculate reaction rates-of-progress and molar
    // production rates of species.
    GetRatesOfProgress(density, x, n, kf, kr, rop1);
    wtot1 = GetMolarProdRates(rop1, wdot1);


    // Calculate temperature term dT/dt.
    Tdot1 = 0.0;
    if (!constT) {
        if (constV) {
            // Use internal energies for constant volume.
            thermo.CalcUs_RT(Tpert, Hs);
            Cp = thermo.CalcBulkCv_R(Tpert, x, n);
        } else {
            // Use enthalpies for constant pressure.
            thermo.CalcHs_RT(Tpert, Hs);
            Cp = thermo.CalcBulkCp_R(Tpert, x, n);
        }
        for (unsigned int i=0; i!=m_mech->SpeciesCount(); ++i) {
            Tdot1 += wdot1[i] * Hs[i];
        }
        Tdot1 *= - Tpert / (density * Cp);
    }

    // CALCULATE JACOBIAN ENTRIES.

    // Calculate temperature Jacobian entry.
    J[m_mech->SpeciesCount()][m_mech->SpeciesCount()] = (Tdot1 - Tdot0) * invdconc;

    // Calculate density Jacobian entry.
    if (constV) {
        J[m_mech->SpeciesCount()][m_mech->SpeciesCount()+1] = (wtot1 - wtot0) * invdconc;
    } else {
        J[m_mech->SpeciesCount()][m_mech->SpeciesCount()+1] = 0.0;
    }
}


// PARENT MECHANISM.

// Returns a pointer to the parent mechanism.
const Sprog::Mechanism *const ReactionSet::Mechanism() const
{
    return m_mech;
}

// Sets the parent mechanism.
void ReactionSet::SetMechanism(Sprog::Mechanism &mech)
{
    m_mech = &mech;

    // Set mechanism on all reactions as well.
    for (RxnPtrVector::iterator i=m_rxns.begin(); i!=m_rxns.end(); i++) {
        (*i)->SetMechanism(mech);
    }
}


// READ/WRITE/COPY FUNCTIONS.

// Writes the reaction set to a binary data stream.
void ReactionSet::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialize version to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));
        
        // Write the number of reactions to the stream.
        unsigned int n = m_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write the Reaction objects to the stream.
        for (RxnPtrVector::const_iterator i=m_rxns.begin(); i!=m_rxns.end(); i++) {
            (*i)->Serialize(out);
        }

        // Write number of reactions with explicit reverse parameters.
        n = m_rev_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of reactions with explicit reverse parameters.
        for (RxnMap::const_iterator i=m_rev_rxns.begin(); i!=m_rev_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

        // Write number of third body reactions.
        n = m_tb_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of third body reactions.
        for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

        // Write number of fall-off reactions.
        n = m_fo_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of fall-off reactions.
        for (RxnMap::const_iterator i=m_fo_rxns.begin(); i!=m_fo_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

        // Write number of reactions with Landau Teller parameters.
        n = m_lt_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of Landau Teller reactions.
        for (RxnMap::const_iterator i=m_lt_rxns.begin(); i!=m_lt_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

        // Write number of reactions with explicit reverse LT parameters.
        n = m_revlt_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of reverse Landau Teller reactions.
        for (RxnMap::const_iterator i=m_revlt_rxns.begin(); i!=m_revlt_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

    } else {
        throw invalid_argument("Output stream not ready (Sprog, ReactionSet::Serialize).");
    }
}

// Reads the reaction set data from a binary data stream.
void ReactionSet::Deserialize(std::istream &in)
{
    // Clear the current reaction et.  We do this before checking
    // the stream condition to avoid confusion in the calling code.
    // Even if the possible exception is handled incorrectly, the
    // set will still be empty.
    releaseMemory();

    if (in.good()) {
        // Read the serialized mechanism version.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;

        switch (version) {
            case 0:
                // Read the number of reactions.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the reactions.
                try {
                    for (unsigned int i=0; i<n; i++) {
                        // Create a reaction object using the stream-reading
                        // constructor.
                        Kinetics::Reaction *rxn = NULL;
                        rxn = new Kinetics::Reaction(in);
                        rxn->SetMechanism(*m_mech);

                        // Add the reaction to the vector.
                        m_rxns.push_back(rxn);
                    }
                } catch (exception &e) {
                    // Clear reaction set memory before throwing error to
                    // higher level.
                    releaseMemory();
                    throw;
                }

                // Read number of reactions with explicit reverse parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read map of reactions with explicit reverse parameters.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_rev_rxns[ix] = m_rxns[ix];
                }

                // Write number of third body reactions.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Write map of third body reactions.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_tb_rxns[ix] = m_rxns[ix];
                }

                // Write number of fall-off reactions.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Write map of fall-off reactions.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_fo_rxns[ix] = m_rxns[ix];
                }

                // Write number of reactions with Landau Teller parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Write map of Landau Teller reactions.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_lt_rxns[ix] = m_rxns[ix];
                }

                // Write number of reactions with explicit reverse LT parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Write map of reverse Landau Teller reactions.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_revlt_rxns[ix] = m_rxns[ix];
                }

                break;
            default:
                throw runtime_error("Serialized version number is unsupported (Sprog, ReactionSet::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, ReactionSet::Deserialize).");
    }
}


// MEMORY MANAGEMENT.

// Clears all memory used by the reaction set.
void ReactionSet::releaseMemory()
{
    // Wipe the cross-referencing maps.
    m_rev_rxns.clear();
    m_lt_rxns.clear();
    m_revlt_rxns.clear();
    m_tb_rxns.clear();
    m_fo_rxns.clear();

    // Delete the reactions.
    RxnPtrVector::iterator i;
    for (i=m_rxns.begin(); i!=m_rxns.end(); i++) {
        delete *i; // Remember to delete memory!
    }
    m_rxns.clear();
}