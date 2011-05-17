/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Reaction class declared in the
    gpc_reaction.h header file.

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

#include "gpc_params.h"
#include "gpc_stoich.h"
#include "gpc_rate_params.h"
#include "gpc_reaction.h"
#include "gpc_mech.h"
#include <exception>
#include <stdexcept>
#include <string>
#include <math.h>
#include "string_functions.h"

using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Reaction::Reaction(void)
{
    // Reaction data.
    m_name       = "";
    m_reversible = false;
    m_dstoich    = m_dreac = m_dprod = 0.0;
    m_arrf.A     = m_arrf.n = m_arrf.E = 0.0;
    m_arrr       = NULL;
    m_lt         = NULL;
    m_revlt      = NULL;
    
    // Third-bodies.
    m_usetb = false;
    m_thirdbodies.clear();

    // Fall-off data.
    m_fotype   = None;
    m_foparams = FALLOFF_PARAMS();
    m_fofn     = NULL;

    // Reaction context.
    m_mech = NULL;
}

// Copy constructor.
Reaction::Reaction(const Sprog::Kinetics::Reaction &rxn)
{
    // Nullify pointer members.
    m_arrr    = NULL;
    m_lt      = NULL;
    m_revlt   = NULL;
    m_mech    = NULL;

    // Use assignment operator to copy.
    *this = rxn;
}

// Stream-reading constructor.
Reaction::Reaction(std::istream &in)
{
    // Reaction data.
    m_name       = "";
    m_reversible = false;
    m_dstoich    = m_dreac = m_dprod = 0.0;
    m_arrf.A     = m_arrf.n = m_arrf.E = 0.0;
    m_arrr       = NULL;
    m_lt         = NULL;
    m_revlt      = NULL;
    
    // Third-bodies.
    m_usetb = false;
    m_thirdbodies.clear();

    // Fall-off data.
    m_fotype   = None;
    m_foparams = FALLOFF_PARAMS();
    m_fofn     = NULL;

    // Reaction context.
    m_mech = NULL;

    Deserialize(in);
}

// Destructor.
Reaction::~Reaction(void)
{
    // Clear memory associated with the reaction object.
    releaseMemory();
}


// OPERATOR OVERLOADING.

// Assignment operator.
Reaction &Reaction::operator=(const Sprog::Kinetics::Reaction &rxn)
{
    // Check for self assignment!
    if (this != &rxn) {
        // Clear current memory used.
        releaseMemory();

        // Copy name and reversibility.
        m_name = rxn.m_name;
        m_reversible = rxn.m_reversible;
        
        // Copy products and reactants.
        m_reac.assign(rxn.m_reac.begin(), rxn.m_reac.end());
        m_prod.assign(rxn.m_prod.begin(), rxn.m_prod.end());
        m_freac.assign(rxn.m_freac.begin(), rxn.m_freac.end());
        m_fprod.assign(rxn.m_fprod.begin(), rxn.m_fprod.end());
        
        // Copy total stoichiometries.
        m_dstoich = rxn.m_dstoich;
        m_dreac = rxn.m_dreac;
        m_dprod = rxn.m_dprod;

        // Copy Arrhenius coefficients (forward and reverse).
        m_arrf = rxn.m_arrf;
        if (rxn.m_arrr != NULL) m_arrr = new ARRHENIUS(*rxn.m_arrr);

        // Copy Landau Teller parameters (forward and reverse).
        if (rxn.m_lt != NULL) m_lt = new LTCOEFFS(*rxn.m_lt);
        if (rxn.m_revlt != NULL) m_revlt = new LTCOEFFS(*rxn.m_revlt);

        // Copy third bodies.
        m_usetb = rxn.m_usetb;
        m_thirdbodies.assign(rxn.m_thirdbodies.begin(), rxn.m_thirdbodies.end());

        // Copy fall-off data.
        m_fotype = rxn.m_fotype;
        m_foparams = rxn.m_foparams;
        m_fofn = rxn.m_fofn;

        // Copy pointer to mechanism.
        m_mech = rxn.m_mech;
    }

    return *this;
}


// REACTANTS.

// Adds an integer stoichiometry reactant to the reaction.
void Reaction::AddReactant(const Sprog::Stoich &reac)
{
    // We must check if the reactant is already defined in this reaction.  Need to
    // check both integer and real reactants.

    // Add the new contribution to the total stoichiometry sums.
    m_dstoich -= (real)reac.Mu();
    m_dreac += (real)reac.Mu();

    // Integer reactants.
    vector<Stoich>::iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        if ((*i).Index() == reac.Index()) {
            // We have found this species in the reactants already.  Append
            // the new contribution.
            (*i).IncMu(reac.Mu());
            return;
        }
    }

    // Real reactants.
    vector<Stoichf>::iterator j;
    for (j=m_freac.begin(); j!=m_freac.end(); j++) {
        if ((*j).Index()== reac.Index()) {
            // We have found this species in the reactants already.  Append
            // the new contribution.
            (*j).IncMu((real)reac.Mu());
            return;
        }
    }

    // If we have got here then the reactant is not defined for this reaction, so
    // we must add it.
    m_reac.push_back(reac);
}

// Adds a real stoichiometry reactant to the reaction.
void Reaction::AddReactant(const Sprog::Stoichf &reac)
{
    // We must check if the reactant is already defined in this reaction.  Need to
    // check both integer and real reactants.

    // Add the new contribution to the total stoichiometry sums.
    m_dstoich -= (real)reac.Mu();
    m_dreac += (real)reac.Mu();

    // Integer reactants.
    vector<Stoich>::iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        if ((*i).Index() == reac.Index()) {
            // We have found this species in the reactants already.  It is currently
            // in the integer stoichiometry, and we now want it in the real stoichiometry.
            Stoichf freac(reac);
            freac.IncMu((real)(*i).Mu());
            m_freac.push_back(freac);
            return;
        }
    }

    // Real reactants.
    vector<Stoichf>::iterator j;
    for (j=m_freac.begin(); j!=m_freac.end(); j++) {
        if ((*j).Index() == reac.Index()) {
            // We have found this species in the reactants already.  Append
            // the new contribution.
            (*j).IncMu(reac.Mu());
            return;
        }
    }

    // If we have got here then the reactant is not defined for this reaction, so
    // we must add it.
    m_freac.push_back(reac);
}

// Adds an integer reactant to the reaction given the species name.
void Reaction::AddReactant(const std::string &name, unsigned int stoich)
{
    if (m_mech != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_mech->SpeciesCount(); i++) {
            if (name.compare(m_mech->Species(i)->Name()) == 0) {
                // Found species in the list, now add it as a reactant.
                AddReactant(Stoich(i, stoich));
                return;
            }
        }
        
        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in species "
                           "list (Sprog, Reaction::AddReactant)."));
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Trying to add reactant before assigning "
                          "parent mechanism (Sprog, Reaction::AddReactant).");
    }
}

// Adds a real reactant to the reaction given the species name.
void Reaction::AddReactant(const std::string &name, real stoich)
{
    if (m_mech != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_mech->SpeciesCount(); i++) {
            if (name.compare(m_mech->Species(i)->Name()) == 0) {
                // Found species in the list, now add it as a reactant.
                AddReactant(Stoichf(i, stoich));
                return;
            }
        }
        
        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in species "
                           "list (Sprog, Reaction::AddReactant)."));
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Trying to add reactant before assigning parent "
                          "mechanism (Sprog, Reaction::AddReactant).");
    }
}

// Removes the reactant species with the given name from the reaction.  If the
// species is not a reactant then does nothing.
void Reaction::RemoveReactant(const std::string &name)
{
    if (m_mech != NULL) {
        // Search through the list of species to find that with
        // the given name.
        int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species:  Loop though integer stoichiometry
            // reactant to find that with this index.
            vector<Stoich>::iterator j;
            for (j=m_reac.begin(); j!=m_reac.end(); j++) {
                if ((*j).Index() == i) {
                    // Remove contribution from the total stoichiometry sums.
                    m_dstoich += (real)((*j).Mu());
                    m_dreac -= (real)((*j).Mu());

                    // We have found the species in the list.
                    m_reac.erase(j);
                }
            }        

            // Now search the real stoichiometry reactants.
            vector<Stoichf>::iterator k;
            for (k=m_freac.begin(); k!=m_freac.end(); k++) {
                if ((*k).Index() == i) {
                    // Remove contribution from the total stoichiometry sums.
                    m_dstoich += (*j).Mu();
                    m_dreac -= (*j).Mu();

                    // We have found the species in the real stoich vector.
                    m_freac.erase(k);
                    return;
                }
            }
        }
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Can't search for species before assigning parent "
                          "mechanism (Sprog, Reaction::RemoveReactant).");
    }
}

// Returns the stoichiometry of the kth integer reactant.
const Stoich Reaction::Reactant(unsigned int k) const
{
    if (k < m_reac.size()) {
        return m_reac[k];
    } else {
        return Stoich(0,0);
    }
}

// Returns the stoichiometry of the kth real reactant.
const Stoichf Reaction::FReactant(unsigned int k) const
{
    if (k < m_freac.size()) {
        return m_freac[k];
    } else {
        return Stoichf(0,0.0);
    }
}

// Returns the number of integer reactants.
int Reaction::ReactantCount() const
{
    return m_reac.size();
}

// Returns the number of real reactants.
int Reaction::FReactantCount() const
{
    return m_freac.size();
}


// PRODUCTS.

// Adds an integer stoichiometry product to the reaction.
void Reaction::AddProduct(const Sprog::Stoich &prod)
{
    // We must check if the product is already defined in this reaction.  Need to
    // check both integer and real products.

    // Add the new contribution to the total stoichiometry sums.
    m_dstoich += (real)prod.Mu();
    m_dprod += (real)prod.Mu();

    // Integer products.
    vector<Stoich>::iterator i;
    for (i=m_prod.begin(); i!=m_prod.end(); i++) {
        if ((*i).Index() == prod.Index()) {
            // We have found this species in the products already.  Append
            // the new contribution.
            (*i).IncMu(prod.Mu());
            return;
        }
    }

    // Real products.
    vector<Stoichf>::iterator j;
    for (j=m_fprod.begin(); j!=m_fprod.end(); j++) {
        if ((*j).Index() == prod.Index()) {
            // We have found this species in the products already.  Append
            // the new contribution.
            (*j).IncMu((real)prod.Mu());
            return;
        }
    }

    // If we have got here then the product is not defined for this reaction, so
    // we must add it.
    m_prod.push_back(prod);
}


// Adds a real stoichiometry product to the reaction.
void Reaction::AddProduct(const Sprog::Stoichf &prod)
{
    // We must check if the product is already defined in this reaction.  Need to
    // check both integer and real products.

    // Add the new contribution to the total stoichiometry sums.
    m_dstoich += prod.Mu();
    m_dprod += prod.Mu();

    // Integer products.
    vector<Stoich>::iterator i;
    for (i=m_prod.begin(); i!=m_prod.end(); i++) {
        if ((*i).Index() == prod.Index()) {
            // We have found this species in the products already.  It is currently
            // in the integer stoichiometry, and we now want it in the real stoichiometry.
            Stoichf fprod(prod);
            fprod.IncMu((real)(*i).Mu());
            m_fprod.push_back(fprod);
            return;
        }
    }

    // Real reactants.
    vector<Stoichf>::iterator j;
    for (j=m_fprod.begin(); j!=m_fprod.end(); j++) {
        if ((*j).Index() == prod.Index()) {
            // We have found this species in the products already.  Append
            // the new contribution.
            (*j).IncMu(prod.Mu());
            return;
        }
    }

    // If we have got here then the product is not defined for this reaction, so
    // we must add it.
    m_fprod.push_back(prod);
}

// Adds an integer product to the reaction given the species name.
void Reaction::AddProduct(const std::string &name, unsigned int stoich)
{
    if (m_mech != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_mech->SpeciesCount(); i++) {
            if (name.compare(m_mech->Species(i)->Name()) == 0) {
                // Found species in the list, now add it as a product.
                AddProduct(Stoich(i, stoich));
                return;
            }
        }
        
        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in "
                           "species list (Sprog, Reaction::AddProduct)."));
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Trying to add reactant before assigning "
                          "parent mechanism (Sprog, Reaction::AddProduct).");
    }
}

// Adds a real product to the reaction given the species name.
void Reaction::AddProduct(const std::string &name, real stoich)
{
    if (m_mech != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_mech->SpeciesCount(); i++) {
            if (name.compare(m_mech->Species(i)->Name()) == 0) {
                // Found species in the list, now add it as a product.
                AddProduct(Stoichf(i, stoich));
                return;
            }
        }
        
        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in "
                           "species list (Sprog, Reaction::AddProduct)."));
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Trying to add reactant before assigning "
                          "parent mechanism (Sprog, Reaction::AddProduct).");
    }
}

// Removes the product species with the given name from the reaction.  If the
// species is not a product then does nothing.
void Reaction::RemoveProduct(const std::string &name)
{
    if (m_mech != NULL) {
        // Search through the list of species to find that with
        // the given name.
        int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species:  Loop though integer stoichiometry
            // product to find that with this index.
            vector<Stoich>::iterator j;
            for (j=m_prod.begin(); j!=m_prod.end(); j++) {
                if ((*j).Index() == i) {
                    // Remove contribution from the total stoichiometry sums.
                    m_dstoich -= (real)((*j).Mu());
                    m_dprod -= (real)((*j).Mu());

                    // We have found the species in the list.
                    m_prod.erase(j);
                }
            }        

            // Now search the real stoichiometry products.
            vector<Stoichf>::iterator k;
            for (k=m_fprod.begin(); k!=m_fprod.end(); k++) {
                if ((*k).Index() == i) {
                    // Remove contribution from the total stoichiometry sums.
                    m_dstoich -= (*j).Mu();
                    m_dprod -= (*j).Mu();

                    // We have found the species in the real stoich vector.
                    m_fprod.erase(k);
                    return;
                }
            }
        }
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Can't search for species before assigning "
                          "parent mechanism (Sprog, Reaction::RemoveProduct).");
   }
}

// Returns the stoichiometry of the kth integer product.
const Stoich Reaction::Product(unsigned int k) const
{
    if (k < m_prod.size()) {
        return m_prod[k];
    } else {
        return Stoich(0,0);
    }
}

// Returns the stoichiometry of the kth real product.
const Stoichf Reaction::FProduct(unsigned int k) const
{
    if (k < m_fprod.size()) {
        return m_freac[k];
    } else {
        return Stoichf(0,0.0);
    }
}

// Returns the number of integer products.
int Reaction::ProductCount() const
{
    return m_prod.size();
}

// Returns the number of real products.
int Reaction::FProductCount() const
{
    return m_fprod.size();
}


// ARRHENIUS COEFFICIENTS.

// Sets the forward Arrhenius parameters.
void Reaction::SetArrhenius(const ARRHENIUS &arr)
{
    m_arrf = arr;
}


// Sets the explicit reverse Arrhenius parameters.
void Reaction::SetRevArrhenius(const Sprog::Kinetics::ARRHENIUS &arr)
{
    if (m_arrr == NULL) {
        m_arrr = new ARRHENIUS(arr);
    } else {
        *m_arrr = arr;
    }
}


// LANDAU TELLER COEFFICIENTS.

// Sets the forward Landau Teller rate parameters.
void Reaction::SetLTCoeffs(const Sprog::Kinetics::LTCOEFFS &lt)
{
    if (m_lt == NULL) {
        m_lt = new LTCOEFFS(lt);
    } else {
        *m_lt = lt;
    }
}

// Sets the reverse Landau Teller rate parameters.
void Reaction::SetRevLTCoeffs(const Sprog::Kinetics::LTCOEFFS &lt)
{
    if (m_revlt == NULL) {
        m_revlt = new LTCOEFFS(lt);
    } else {
        *m_revlt = lt;
    }
}


// THIRD-BODIES.

// Returns a constant reference to the vector of third bodies with their
// coefficients.
const vector<Stoichf> &Reaction::ThirdBodies() const
{
    return m_thirdbodies;
}

// Returns the coefficient for the ith third body.
Stoichf Reaction::ThirdBody(unsigned int i) const
{
    if (i < m_thirdbodies.size()) {
        return m_thirdbodies[i];
    } else {
        // Invalid index.
        return Stoichf(0, 1.0);
    }
}

// Returns the number of third body coefficients defined for this reaction.
int Reaction::ThirdBodyCount() const
{
    return m_thirdbodies.size();
}

// Adds a third body to the reaction given the stoichiometric structure.
void Reaction::AddThirdBody(const Sprog::Stoichf &tb)
{
    m_usetb = true;
    m_thirdbodies.push_back(tb);
}

// Adds a third body to the reaction given the species and the
// coefficient.
void Reaction::AddThirdBody(const unsigned int sp, Sprog::real coeff)
{
    // Add a new Stoichf to the array of third bodies.
    m_usetb = true;
    m_thirdbodies.push_back(Stoichf(sp, coeff));
}

// Adds a third body to the reaction given the species name.
void Reaction::AddThirdBody(const std::string &name, Sprog::real coeff)
{
    if (m_mech != NULL) {
        // Find the species in the mechanism with the given name.
        unsigned int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species!
            m_usetb = true;
            m_thirdbodies.push_back(Stoichf(i, coeff));
        }
    } else {
        throw logic_error("Can't add third body before assigning "
                          "parent mechanism (Sprog, Reaction::AddThirdBody).");
    }
}

// Removes a third body, given by name, from the reaction.  If the third
// body is not defined for this reaction then does nothing.
void Reaction::RemoveThirdBody(const std::string &name)
{
    if (m_mech != NULL) {
        // Search through the list of species to find that with
        // the given name.
        int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species:  Loop though third bodies and find
            // that with this index.
            vector<Stoichf>::iterator j;
            for (j=m_thirdbodies.begin(); j!=m_thirdbodies.end(); j++) {
                if ((*j).Index() == i) {
                    // We have found the species in the list.
                    m_thirdbodies.erase(j);
                }
            }        
        }
    } else {
        throw logic_error("Can't find third body name before assigning "
                          "parent mechanism (Sprog, Reaction::RemoveThirdBody).");
    }
}


// LOW PRESSURE LIMIT.

// Returns the Arrhenius parameters for the low pressure limit.
const ARRHENIUS &Reaction::LowPressureLimit() const
{
    return m_foparams.LowP_Limit;
}

// Sets the low pressure limit Arrhenius parameters.
void Reaction::SetLowPressureLimit(const Sprog::Kinetics::ARRHENIUS &lowp)
{
    if (m_fotype == None) m_fotype = Lindemann;
    m_foparams.LowP_Limit = lowp;
}


// FALL-OFF THIRD BODY.

// Returns a pointer to the species used as a third body in the fall-off reaction.
const Species *const Reaction::FallOffThirdBody() const
{
    if ((m_mech != NULL) && (m_foparams.ThirdBody >= 0)) {
        return m_mech->Species(m_foparams.ThirdBody);
    }
    return NULL;
}

// Sets the species to use as a third body for fall-off calculations.
void Reaction::SetFallOffThirdBody(int sp)
{
    m_foparams.ThirdBody = sp;
}

// Sets the species to use as a third body for fall-off calculations
// given the species name.
void Reaction::SetFallOffThirdBody(const std::string &name)
{
    if (m_mech != NULL) {
        // Locate the species by name in the vector.
        unsigned int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species!
            m_foparams.ThirdBody = i;
        } else {
            // If the name is not a valid species then
            // use all species as third body.
            m_foparams.ThirdBody = -1;
        }
    } else {
        throw logic_error("Can't find set fall-off third body before assigning "
                          "parent mechanism (Sprog, Reaction::SetFallOffThirdBody).");
    }
}


// FALL-OFF PARAMETERS.

// Returns the fall-off type.
Kinetics::FALLOFF_FORM Reaction::FallOffType() const
{
    return m_fotype;
}

// Returns the fall-off parameters.
const Kinetics::FALLOFF_PARAMS &Reaction::FallOffParams() const
{
    return m_foparams;
}

// Sets the fall-off type and parameters.
void Reaction::SetFallOffParams(const FALLOFF_FORM form, 
                                const real params[FALLOFF_PARAMS::MAX_FALLOFF_PARAMS])
{
    m_fotype = form;
    
    for (int i=0; i<FALLOFF_PARAMS::MAX_FALLOFF_PARAMS; i++) {
        m_foparams.Params[i] = params[i];
    }
}


// FALL-OFF FUNCTIONAL FORMS.

// 3-parameter Troe fall-off form.
real Reaction::FTROE3(real T, real logpr) const
{
    real fcent, c, n, F;
    const real d = 0.14;

    fcent = log10(((1.0 - m_foparams.Params[0]) * exp(-T / m_foparams.Params[1])) + 
                  (m_foparams.Params[0] * exp(-T / m_foparams.Params[2])));
    c = logpr - 0.4 - (0.67 * fcent);
    n = 0.75 - (1.27 * fcent) - (d * c);
    c = c / n;
    F = pow(10.0, fcent / (1 + (c*c)));
    return F;
}

// 4-parameter Troe fall-off form.
real Reaction::FTROE4(real T, real logpr) const
{
    real fcent, c, n, F;
    const real d = 0.14;

    fcent = log10(((1.0 - m_foparams.Params[0]) * exp(-T / m_foparams.Params[1])) + 
                  (m_foparams.Params[0] * exp(-T / m_foparams.Params[2])) + 
                  exp(-m_foparams.Params[3] / T));
    c = logpr - 0.4 - (0.67 * fcent);
    n = 0.75 - (1.27 * fcent) - (d * c);
    c = c / n;
    F = pow(10.0, fcent / (1 + (c*c)));
    return F;
}

// SRI fall-off form.
real Reaction::FSRI(real T, real logpr) const
{
    real x = 1.0 / (1.0 + (logpr * logpr));
    real F = m_foparams.Params[3] * pow(T, m_foparams.Params[4]) * 
             pow((m_foparams.Params[0]*exp(-m_foparams.Params[1]/T)) + 
                 exp(-T/m_foparams.Params[2]), x);
    return F;
}

// Custom functional form for fall-off.
FallOffFnPtr Reaction::FallOffFn() const
{
    return m_fofn;
}


/*
// SPECIES VECTOR.

// Returns a pointer to the vector of species used to define the reaction.
const SpeciesPtrVector *const Reaction::Species(void) const
{
    return m_species;
}

// Sets the vector of species used to define the reaction.
void Reaction::SetSpecies(const SpeciesPtrVector *const sp)
{
    m_species = sp;
}
*/


// PARENT MECHANISM.

// Returns a pointer to the parent mechanism.
const Mechanism *const Reaction::Mechanism() const
{
    return m_mech;
}

// Sets the parent mechanism.
void Reaction::SetMechanism(Sprog::Mechanism &mech)
{
    m_mech = &mech;
}


// RATE CALCULATION.

// Calculates the rate of progress of this reaction.
real Reaction::RateOfProgress(real density, const real *const x, 
                              unsigned int n, real kforward, 
                              real kreverse) const
{
    int j=0;
    unsigned int k=0;
    real rop=0.0, rev=0.0;

    if (n >= m_mech->Species().size()) {
        // Use rop to store forward rates of production, 
        // and rev to store reverse rates.
        rop = kforward;
        rev = kreverse;

        // Integer reactants.
        for (k=0; k!=m_reac.size(); ++k) {
            // As the stoichiometry is integer, it is more computationally efficient
            // to multiply the values together than to use the pow() function.
            for (j=0; j!=m_reac[k].Mu(); ++j) {
                rop *= density * x[m_reac[k].Index()];
            }
        }

        // Integer products.
        for (k=0; k!=m_prod.size(); ++k) {
            // As the stoichiometry is integer, it is more computationally efficient
            // to multiply the values together than to use the pow() function.
            for (j=0; j!=m_prod[k].Mu(); ++j) {
                rev *= density * x[m_prod[k].Index()];
            }
        }

        // Real reactants.
        for (k=0; k!=m_freac.size(); ++k) {
            // Now the stoichiometry is non-integer, we must use the pow() function.
            rop *= pow(density * x[m_freac[k].Index()], m_freac[k].Mu()); 

        }

        // Real products.
        for (k=0; k!=m_fprod.size(); ++k) {
            // Now the stoichiometry is non-integer, we must use the pow() function.
            rev *= pow(density * x[m_fprod[k].Index()], m_fprod[k].Mu()); 

        }

        // Calculate the net rate of production.
        rop -= rev;
    }

    return rop;
}


// READ/WRITE/COPY FUNCTIONS.

// Creates a copy of the reaction object.
Reaction *Reaction::Clone(void) const
{
    return new Reaction(*this);
}

// Writes the reaction to a binary data stream.
void Reaction::Serialize(std::ostream &out) const
{
    if (out.good()) {
        const unsigned int trueval  = 1;
        const unsigned int falseval = 0;

        // Write the serialisation version number to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write the length of the reaction name to the stream.
        unsigned int n = m_name.length();
        out.write((char*)&n, sizeof(n));

        // Write the reaction name to the stream.
        if (n > 0) {
            out.write(m_name.c_str(), n);
        }

        // Write reaction reversibility.
        if (m_reversible) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write integer reactant stoichiometry count.
        n = m_reac.size();
        out.write((char*)&n, sizeof(n));

        // Write integer reactant stoichiometry.
        for (unsigned int i=0; i<n; i++) {
            // Write species index.
            int ix = m_reac[i].Index();
            out.write((char*)&ix, sizeof(ix));
            
            // Write mu.
            int mu = m_reac[i].Mu();
            out.write((char*)&mu, sizeof(mu));
        }

        // Write integer product stoichiometry count.
        n = m_prod.size();
        out.write((char*)&n, sizeof(n));

        // Write integer product stoichiometry.
        for (unsigned int i=0; i<n; i++) {
            // Write species index.
            int ix = m_prod[i].Index();
            out.write((char*)&ix, sizeof(ix));
            
            // Write mu.
            int mu = m_prod[i].Mu();
            out.write((char*)&mu, sizeof(mu));
        }

        // Write real reactant stoichiometry count.
        n = m_freac.size();
        out.write((char*)&n, sizeof(n));

        // Write real reactant stoichiometry.
        for (unsigned int i=0; i<n; i++) {
            // Write species index.
            int ix = m_freac[i].Index();
            out.write((char*)&ix, sizeof(ix));
            
            // Write mu.
            double mu = m_freac[i].Mu();
            out.write((char*)&mu, sizeof(mu));
        }

        // Write real product stoichiometry count.
        n = m_fprod.size();
        out.write((char*)&n, sizeof(n));

        // Write real product stoichiometry.
        for (unsigned int i=0; i<n; i++) {
            // Write species index.
            int ix = m_fprod[i].Index();
            out.write((char*)&ix, sizeof(ix));
            
            // Write mu.
            double mu = m_fprod[i].Mu();
            out.write((char*)&mu, sizeof(mu));
        }

        // Write stoichiometry changes.
        double mu = (double)m_dstoich;
        out.write((char*)&mu, sizeof(mu));
        mu = (double)m_dreac;
        out.write((char*)&mu, sizeof(mu));
        mu = (double)m_dprod;
        out.write((char*)&mu, sizeof(mu));

        // Write forward Arrhenius coefficients.
        double A  = (double)m_arrf.A;
        double nn = (double)m_arrf.n;
        double E  = (double)m_arrf.E;
        out.write((char*)&A, sizeof(A));
        out.write((char*)&nn, sizeof(nn));
        out.write((char*)&E, sizeof(E));

        // Write reverse Arrhenius coefficients.
        if (m_arrr != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            A  = (double)m_arrr->A;
            nn = (double)m_arrr->n;
            E  = (double)m_arrr->E;
            out.write((char*)&A, sizeof(A));
            out.write((char*)&nn, sizeof(nn));
            out.write((char*)&E, sizeof(E));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write forward LT parameters.
        if (m_lt != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            A = (double)m_lt->B;
            E = (double)m_lt->C;
            out.write((char*)&A, sizeof(A));
            out.write((char*)&E, sizeof(E));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write reverse LT parameters.
        if (m_revlt != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            A = (double)m_revlt->B;
            E = (double)m_revlt->C;
            out.write((char*)&A, sizeof(A));
            out.write((char*)&E, sizeof(E));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write third body flag.
        if (m_usetb) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write third body count.
        n = m_thirdbodies.size();
        out.write((char*)&n, sizeof(n));

        // Write third bodies.
        for (unsigned int i=0; i<n; i++) {
            // Write species index.
            int ix = m_thirdbodies[i].Index();
            out.write((char*)&ix, sizeof(ix));
            
            // Write mu.
            double mu = m_thirdbodies[i].Mu();
            out.write((char*)&mu, sizeof(mu));
        }

        // Write fall-off type.
        n = (unsigned int)m_fotype;
        out.write((char*)&n, sizeof(n));

        // Write fall-off low pressure limit.
        A  = (double)m_foparams.LowP_Limit.A;
        nn = (double)m_foparams.LowP_Limit.n;
        E  = (double)m_foparams.LowP_Limit.E;
        out.write((char*)&A, sizeof(A));
        out.write((char*)&nn, sizeof(nn));
        out.write((char*)&E, sizeof(E));

        // Write fall-off third body.
        int ix = m_foparams.ThirdBody;
        out.write((char*)&ix, sizeof(ix));

        // Write fall-off parameter count.
        n = FALLOFF_PARAMS::MAX_FALLOFF_PARAMS;
        out.write((char*)&n, sizeof(n));

        // Write fall-off parameters.
        for (unsigned int i=0; i<n; i++) {
            A = (double)m_foparams.Params[i];
            out.write((char*)&A, sizeof(A));
        }

        // There is no way of outputting custom fall-off functions at the moment.
        // This needs to corrected in the future.

    } else {
        throw invalid_argument("Output stream not ready (Sprog, Reaction::Serialize).");
    }
}

// Reads the reaction data from a binary data stream.
void Reaction::Deserialize(std::istream &in)
{
    // Clear the reaction of all current data.
    releaseMemory();

    if (in.good()) {
        // Read the serialized species version number.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0; // Need for reading name length.
        char *name = NULL;
        double A = 0.0, nn = 0.0, E = 0.0;
        int itb = 0;

        switch (version) {
            case 0:
                // Read the length of the reaction name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                
                // Read the reaction name.
                if (n > 0) {
                    name = new char[n];
                    in.read(name, n);
                    m_name.assign(name, n);
                    delete [] name;
                } else {
                    m_name = "";
                }

                // Read reaction reversibility.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n==1) {
                    m_reversible = true;
                } else {
                    m_reversible = false;
                }

                // Read integer reactant stoichiometry count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_reac.reserve(n);

                // Read integer reactant stoichiometry.
                for (unsigned int i=0; i<n; i++) {
                    int ix = 0;
                    int mu = 0;

                    // Read species index.
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));
                    
                    // Read mu.
                    in.read(reinterpret_cast<char*>(&mu), sizeof(mu));

                    // Push a new Stoich object into the vector.
                    m_reac.push_back(Stoich(ix, mu));
                }

                // Read integer product stoichiometry count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_prod.reserve(n);

                // Read integer product stoichiometry.
                for (unsigned int i=0; i<n; i++) {
                    int ix = 0;
                    int mu = 0;

                    // Read species index.
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));
                    
                    // Read mu.
                    in.read(reinterpret_cast<char*>(&mu), sizeof(mu));

                    // Push a new Stoich object into the vector.
                    m_prod.push_back(Stoich(ix, mu));
                }

                // Read real reactant stoichiometry count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_freac.reserve(n);

                // Read real reactant stoichiometry.
                for (unsigned int i=0; i<n; i++) {
                    int ix = 0;
                    double mu = 0;

                    // Read species index.
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));
                    
                    // Read mu.
                    in.read(reinterpret_cast<char*>(&mu), sizeof(mu));

                    // Push a new Stoich object into the vector.
                    m_freac.push_back(Stoichf(ix, (real)mu));
                }

                // Read real product stoichiometry count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_fprod.reserve(n);

                // Read real product stoichiometry.
                for (unsigned int i=0; i<n; i++) {
                    int ix = 0;
                    double mu = 0;

                    // Read species index.
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));
                    
                    // Read mu.
                    in.read(reinterpret_cast<char*>(&mu), sizeof(mu));

                    // Push a new Stoich object into the vector.
                    m_fprod.push_back(Stoichf(ix, (real)mu));
                }

                // Read stoichiometry changes.
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                m_dstoich = (real)A;
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                m_dreac = (real)A;
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                m_dprod = (real)A;

                // Read forward Arrhenius coefficients.
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
                in.read(reinterpret_cast<char*>(&E), sizeof(E));
                m_arrf.A = (real)A;
                m_arrf.n = (real)nn;
                m_arrf.E = (real)E;

                // Read reverse Arrhenius coefficients.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    in.read(reinterpret_cast<char*>(&A), sizeof(A));
                    in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
                    in.read(reinterpret_cast<char*>(&E), sizeof(E));
                    m_arrr = new Kinetics::ARRHENIUS();
                    m_arrr->A = (real)A;
                    m_arrr->n = (real)nn;
                    m_arrr->E = (real)E;
                }

                // Read forward LT parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    in.read(reinterpret_cast<char*>(&A), sizeof(A));
                    in.read(reinterpret_cast<char*>(&E), sizeof(E));
                    m_lt = new LTCOEFFS();
                    m_lt->B = (real)A;
                    m_lt->C = (real)E;
                }

                // Read reverse LT parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    in.read(reinterpret_cast<char*>(&A), sizeof(A));
                    in.read(reinterpret_cast<char*>(&E), sizeof(E));
                    m_revlt = new LTCOEFFS();
                    m_revlt->B = (real)A;
                    m_revlt->C = (real)E;
                }

                // Read third body flag.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_usetb = true;
                } else {
                    m_usetb = false;
                }

                // Read third body count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_thirdbodies.reserve(n);

                // Read third bodies.
                for (unsigned int i=0; i<n; i++) {
                    int ix = 0;
                    double mu = 0.0;

                    // Read species index.
                   in.read(reinterpret_cast<char*>(&ix), sizeof(ix));
                    
                    // Read mu.
                   in.read(reinterpret_cast<char*>(&mu), sizeof(mu));

                   // Add third body to vector.
                   m_thirdbodies.push_back(Stoichf(ix,mu));
                }

                // Read fall-off type.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_fotype = (FALLOFF_FORM)n;

                // Read fall-off low pressure limit.
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
                in.read(reinterpret_cast<char*>(&E), sizeof(E));
                m_foparams.LowP_Limit.A = (real)A;
                m_foparams.LowP_Limit.n = (real)nn;
                m_foparams.LowP_Limit.E = (real)E;

                // Read fall-off third body.
                in.read(reinterpret_cast<char*>(&itb), sizeof(itb));
                m_foparams.ThirdBody = itb;

                // Read fall-off parameter count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read fall-off parameters.
                for (unsigned int i=0; i<n; i++) {
                    in.read(reinterpret_cast<char*>(&A), sizeof(A));
                    m_foparams.Params[i] = (real)A;
                }

                break;
            default:
                throw runtime_error("Reaction serialized version number is "
                                    "unsupported (Sprog, Reaction::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, Reaction::Deserialize).");
    }
}

// Prints a diagnostic output file containing all the
// reaction data.  This is used to debug.
void Reaction::WriteDiagnostics(std::ostream &out) const
{
    string data = "";
    real val = 0.0;
    int ival = 0;

    if (out.good()) {
        // Name.
        out.write(string(m_name+" ").c_str(), m_name.length());

        // Reaction reversibility.
        if (m_reversible) {
            data = "R ";
        } else {
            data = "- ";
        }
        out.write(data.c_str(), data.length());

        // Integer reactant stoichiometry.
        for (unsigned int i=0; i!=m_reac.size(); ++i) {
            // Species name.
            data = m_mech->Species(m_reac[i].Index())->Name() + " ";
            out.write(data.c_str(), data.length());
            // Stoichiometry.
            val = m_reac[i].Mu();
            data = cstr(val) + " ";
            out.write(data.c_str(), data.length());
        }

        // Integer product stoichiometry.
        for (unsigned int i=0; i!=m_prod.size(); ++i) {
            // Species name.
            data = m_mech->Species(m_prod[i].Index())->Name() + " ";
            out.write(data.c_str(), data.length());
            // Stoichiometry.
            val = m_prod[i].Mu();
            data = cstr(val) + " ";
            out.write(data.c_str(), data.length());
        }

        // Real reactant stoichiometry.
        for (unsigned int i=0; i!=m_freac.size(); ++i) {
            // Species name.
            data = m_mech->Species(m_freac[i].Index())->Name() + " ";
            out.write(data.c_str(), data.length());
            // Stoichiometry.
            val = m_freac[i].Mu();
            data = cstr(val) + " ";
            out.write(data.c_str(), data.length());
        }

        // Real product stoichiometry.
        for (unsigned int i=0; i!=m_fprod.size(); ++i) {
            // Species name.
            data = m_mech->Species(m_fprod[i].Index())->Name() + " ";
            out.write(data.c_str(), data.length());
            // Stoichiometry.
            val = m_fprod[i].Mu();
            data = cstr(val) + " ";
            out.write(data.c_str(), data.length());
        }

        // Stoichiometry changes.
        data = cstr(m_dstoich) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_dreac) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_dprod) + " ";
        out.write(data.c_str(), data.length());

        // Forward Arrhenius coefficients.
        data = cstr(m_arrf.A) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_arrf.n) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_arrf.E) + " ";
        out.write(data.c_str(), data.length());

        // Reverse Arrhenius coefficients.
        if (m_arrr != NULL) {
            data = cstr(m_arrr->A) + " ";
            out.write(data.c_str(), data.length());
            data = cstr(m_arrr->n) + " ";
            out.write(data.c_str(), data.length());
            data = cstr(m_arrr->E) + " ";
            out.write(data.c_str(), data.length());
        }

        // Forward LT parameters.
        if (m_lt != NULL) {
            data = cstr(m_lt->B) + " ";
            out.write(data.c_str(), data.length());
            data = cstr(m_lt->C) + " ";
            out.write(data.c_str(), data.length());
        }

        // Reverse LT parameters.
        if (m_revlt != NULL) {
            data = cstr(m_revlt->B) + " ";
            out.write(data.c_str(), data.length());
            data = cstr(m_revlt->C) + " ";
            out.write(data.c_str(), data.length());
        }

        // Third body flag.
        if (m_usetb) {
            data = "TB ";
        } else {
            data = "-- ";
        }
        out.write(data.c_str(), data.length());

        // Write third bodies.
        for (unsigned int i=0; i!=m_thirdbodies.size(); ++i) {
            // Write species name.
            data = m_mech->Species(m_thirdbodies[i].Index())->Name() + " ";
            out.write(data.c_str(), data.length());
            // Write mu.
            data = cstr(m_thirdbodies[i].Mu()) + " ";
            out.write(data.c_str(), data.length());
        }

        // Write fall-off type.
        ival = (int)m_fotype;
        data = cstr(ival) + " ";
        out.write(data.c_str(), data.length());

        // Write fall-off low pressure limit.
        data = cstr(m_foparams.LowP_Limit.A) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_foparams.LowP_Limit.n) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_foparams.LowP_Limit.E) + " ";
        out.write(data.c_str(), data.length());

        // Write fall-off third body.
        if (m_foparams.ThirdBody >= 0) {
            data = m_mech->Species(m_foparams.ThirdBody)->Name() + " ";
            out.write(data.c_str(), data.length());
        }

        // Write fall-off parameters.
        for (unsigned int i=0; i!=(unsigned)FALLOFF_PARAMS::MAX_FALLOFF_PARAMS; ++i) {
            data = cstr(m_foparams.Params[i]) + " ";
            out.write(data.c_str(), data.length());
        }

        // New line.
        data = "\n";
        out.write(data.c_str(), data.length());
    }
}

/*!
@param[in]      out             Output stream for file creation.
@param[in]      RejectSpecies   String vector containing the names of LOI rejected species.
*/
void Reaction::WriteReducedMechReacs(std::ostream &out, std::vector<std::string> RejectSpecies) const 
{

    if (out.good()) {

        //Check to see if the reaction contains any of the rejected species.
        //If so, return without writing anything.
        for (unsigned int n = 0; n < RejectSpecies.size(); n++){
            for (unsigned int j = 0; j < m_reac.size(); j++){
                if (RejectSpecies[n] == m_mech->Species(m_reac[j].Index())->Name())
                    return;
            }
            for (unsigned int k = 0; k < m_prod.size(); k++){
                if (RejectSpecies[n] == m_mech->Species(m_prod[k].Index())->Name())
                    return;
            }
            for (unsigned int l = 0; l < m_freac.size(); l++){
                if (RejectSpecies[n] == m_mech->Species(m_freac[l].Index())->Name())
                    return;
            }
            for (unsigned int m = 0; m < m_fprod.size(); m++){
                if (RejectSpecies[n] == m_mech->Species(m_fprod[m].Index())->Name())
                    return;
            }
        }
		
		std::string tbstr;
		if (m_usetb) {
			if (m_foparams.LowP_Limit.A != 0 &&  m_foparams.LowP_Limit.n != 0 && m_foparams.LowP_Limit.E != 0){
				tbstr = "(+M)";
			} else {
                tbstr = "+M";
			};
        }

        // Integer reactant stoichiometry.
        if (m_reac.size() > 0 && m_prod.size() > 0){
            for (unsigned int i = 0; i != m_reac.size(); ++i) {
                // Species name.
                for (int j = 0; j < m_reac[i].Mu(); j++) {
                    out << m_mech->Species(m_reac[i].Index())->Name();
                    if (i < m_reac.size() - 1 || j < m_reac[i].Mu() - 1)
                        out << " + ";
                }
            }


            if (m_usetb){
                out << tbstr;
            }

                // Reaction reversibility.
                if (m_reversible) {
                    out << " =  ";
                } else {
                    out << " => ";
                }
        
            // Integer product stoichiometry.
            for (unsigned int i = 0; i != m_prod.size(); ++i) {
                // Species name.
                for (int j = 0; j < m_prod[i].Mu(); j++) {
                    out << m_mech->Species(m_prod[i].Index())->Name();
                    if (i < m_prod.size() - 1 || j < m_prod[i].Mu() - 1)
                        out << " + ";
                }
            }
            if (m_usetb){
                out << tbstr;
            }

            // Forward Arrhenius coefficients.
            out << " " << m_arrf.A / pow(1.0e-6, ReactantStoich() - 1.0 + (m_usetb?1.0:0.0)) << " " << m_arrf.n << " " << m_arrf.E / 4.184E7 / 1.0e-7 << "\n";

            // Reverse Arrhenius coefficients.
            if (m_arrr != NULL)
                out << "Rev / " << m_arrr->A / pow(1.0e-6, ProductStoich() - 1.0 + (m_usetb?1.0:0.0)) << " " << m_arrr->n << " " << m_arrr->E / 4.184E7 / 1.0e-7 << " /\n";
        }

        // Real reactant stoichiometry.
        if (m_freac.size() > 0 && m_fprod.size() > 0){
            for (unsigned int i = 0; i != m_freac.size(); ++i) {
            // Species name.
                for (unsigned int j = 0; j < m_freac[i].Mu(); j++){
                    out << m_mech->Species(m_freac[i].Index())->Name();
                    if (i < m_freac.size() - 1 || j < m_freac[i].Mu() - 1)
                    out << " + ";
                }
            }

            if (m_usetb){
                out << tbstr;
            }

            // Reaction reversibility.
            if (m_reversible) {
                out << " =  ";
            } else {
                out << " => ";
            }

            // Real product stoichiometry.
            for (unsigned int i = 0; i != m_fprod.size(); ++i) {
            // Species name.
                for (int j = 0; j < m_fprod[i].Mu(); j++){
                   out << m_mech->Species(m_fprod[i].Index())->Name();
                    if (i < m_fprod.size() - 1 || j < m_fprod[i].Mu() - 1)
                        out << " + ";
                }
            }

            if (m_usetb){
                out << tbstr;
            }

            // Forward Arrhenius coefficients.
            out << " " << m_arrf.A / pow(1.0e-6, ReactantStoich() - 1.0 + (m_usetb?1.0:0.0)) << " " << m_arrf.n << " " << m_arrf.E / 4.184E7 / 1.0e-7 << "\n";

            // Reverse Arrhenius coefficients.
            if (m_arrr != NULL)
                out << "Rev / " << m_arrr->A / pow(1.0e-6, ProductStoich() - 1.0 + (m_usetb?1.0:0.0)) << " " << m_arrr->n << " " << m_arrr->E / 4.184E7 / 1.0e-7 << " \n";
        }

        // Forward LT parameters.
        if (m_lt != NULL) {
            out << m_lt->B << " " << m_lt->C << " ";
        }

        // Reverse LT parameters.
        if (m_revlt != NULL) {
            out << m_revlt->B << " " << m_revlt->C << " ";
        }
        
        // Third body flag.
        if (m_usetb) {

            // Write fall-off low pressure limit.
            if (m_foparams.LowP_Limit.A != 0 &&  m_foparams.LowP_Limit.n != 0 && m_foparams.LowP_Limit.E != 0){
                out << "LOW /" << cstr((m_foparams.LowP_Limit.A)/(pow(1.0e-6, ReactantStoich()))) << " " << m_foparams.LowP_Limit.n << " " << (m_foparams.LowP_Limit.E/ 4.184E7 / 1.0e-7) << " ";
                out << " / " << "\n";

                 //Write fall-off third body.
                if (m_foparams.ThirdBody >= 0) {
                    out << m_mech->Species(m_foparams.ThirdBody)->Name() << " ";
                }

                // Write fall-off parameters.
                out << "TROE /";
                for (unsigned int i = 0; i != (unsigned) (FALLOFF_PARAMS::MAX_FALLOFF_PARAMS - 1); ++i) {
                    out << m_foparams.Params[i] << " ";
                }
                out << "/ " << "\n";
            }
        }

        // New line.
        out << "\n";
 
    }
    out.flush();
}


// MEMORY MANAGEMENT.

void Reaction::releaseMemory(void)
{
    m_name.clear();
    m_reac.clear();
    m_prod.clear();
    m_freac.clear();
    m_fprod.clear();
    m_dstoich = m_dreac = m_dprod = 0.0;
    m_arrf.A = m_arrf.n = m_arrf.E = 0.0;
    if (m_arrr != NULL) delete m_arrr;
    m_arrr = NULL;
    if (m_lt != NULL) delete m_lt;
    m_lt = NULL;
    if (m_revlt != NULL) delete m_revlt;
    m_revlt = NULL;
    m_usetb = false;
    m_thirdbodies.clear();
    m_fotype = None;
    m_foparams.LowP_Limit.A = m_foparams.LowP_Limit.n = m_foparams.LowP_Limit.E = 0.0;
    m_foparams.ThirdBody = -1;
    m_fofn = NULL;
    m_mech = NULL;
}
