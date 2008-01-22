#include "gpc_params.h"
#include "gpc_stoich.h"
#include "gpc_rate_params.h"
#include "gpc_reaction.h"
#include "gpc_mech.h"
#include <exception>
#include <stdexcept>
#include <string>
#include <math.h>

using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Reaction::Reaction(void)
{
    // Reaction data.
    m_name = "";
    m_reversible = false;
    m_dstoich = m_dreac = m_dprod = 0.0;
    m_arrf.A = m_arrf.n = m_arrf.E = 0.0;
    m_arrr = NULL;
    m_lt = NULL;
    m_revlt = NULL;
    
    // Third-bodies.
    m_usetb = false;
    m_thirdbodies.clear();

    // Fall-off data.
    m_fotype = None;
    m_foparams = FALLOFF_PARAMS();
    m_fofn = NULL;

    // Reaction context.
    m_species = NULL;
    m_mech = NULL;
}

// Copy constructor.
Reaction::Reaction(const Sprog::Kinetics::Reaction &rxn)
{
    // Nullify pointer members.
    m_arrr = NULL;
    m_lt = NULL;
    m_revlt = NULL;
    m_species = NULL;
    m_mech = NULL;

    // Use assignment operator to copy.
    *this = rxn;
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

        // Copy pointer to species vector.
        m_species = rxn.m_species;

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
    if (m_species != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_species->size(); i++) {
            if (name.compare((*m_species)[i]->Name()) == 0) {
                // Found species in the list, now add it as a reactant.
                AddReactant(Stoich(i, stoich));
                return;
            }
        }
        
        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in species list (Reaction::AddReactant)."));
    } else {
        // Species vector has not been assigned.
        throw exception("Trying to add reactant before assigning species vector (Reaction::AddReactant).");
    }
}

// Adds a real reactant to the reaction given the species name.
void Reaction::AddReactant(const std::string &name, real stoich)
{
    if (m_species != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_species->size(); i++) {
            if (name.compare((*m_species)[i]->Name()) == 0) {
                // Found species in the list, now add it as a reactant.
                AddReactant(Stoichf(i, stoich));
                return;
            }
        }
        
        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in species list (Reaction::AddReactant)."));
    } else {
        // Species vector has not been assigned.
        throw exception("Trying to add reactant before assigning species vector (Reaction::AddReactant).");
    }
}

// Removes the reactant species with the given name from the reaction.  If the
// species is not a reactant then does nothing.
void Reaction::RemoveReactant(const std::string &name)
{
    // Search through the list of species to find that with
    // the given name.
    unsigned int i;
    for (i=0; i<m_species->size(); i++) {
        if (*(*m_species)[i] == name) {
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
    if (m_species != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_species->size(); i++) {
            if (name.compare((*m_species)[i]->Name()) == 0) {
                // Found species in the list, now add it as a product.
                AddProduct(Stoich(i, stoich));
                return;
            }
        }
        
        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in species list (Reaction::AddProduct)."));
    } else {
        // Species vector has not been assigned.
        throw exception("Trying to add reactant before assigning species vector (Reaction::AddProduct).");
    }
}

// Adds a real product to the reaction given the species name.
void Reaction::AddProduct(const std::string &name, real stoich)
{
    if (m_species != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_species->size(); i++) {
            if (name.compare((*m_species)[i]->Name()) == 0) {
                // Found species in the list, now add it as a product.
                AddProduct(Stoichf(i, stoich));
                return;
            }
        }
        
        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in species list (Reaction::AddProduct)."));
    } else {
        // Species vector has not been assigned.
        throw exception("Trying to add reactant before assigning species vector (Reaction::AddProduct).");
    }
}

// Removes the product species with the given name from the reaction.  If the
// species is not a product then does nothing.
void Reaction::RemoveProduct(const std::string &name)
{
    // Search through the list of species to find that with
    // the given name.
    unsigned int i;
    for (i=0; i<m_species->size(); i++) {
        if (*(*m_species)[i] == name) {
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
    // Find the species in the vector with the given name.
    unsigned int i;
    for (i=0; i<m_species->size(); i++) {
        if (*(*m_species)[i] == name) {
            // Found the species!
            m_usetb = true;
            m_thirdbodies.push_back(Stoichf(i, coeff));
            return;
        }
    }
}

// Removes a third body, given by name, from the reaction.  If the third
// body is not defined for this reaction then does nothing.
void Reaction::RemoveThirdBody(const std::string &name)
{
    // Search through the list of species to find that with
    // the given name.
    unsigned int i;
    for (i=0; i<m_species->size(); i++) {
        if (*(*m_species)[i] == name) {
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
    if (m_foparams.ThirdBody >= 0) {
        return (*m_species)[m_foparams.ThirdBody];
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
    // Locate the species by name in the vector.
    unsigned int i;
    for (i=0; i<m_species->size(); i++) {
        if (*(*m_species)[i] == name) {
            // Found the species!
            m_foparams.ThirdBody = i;
            return;
        }
    }

    // If the name is not a valid species then use all species as third body.
    m_foparams.ThirdBody = -1;
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


// PARENT MECHANISM.

// Returns a pointer to the parent mechanism.
Mechanism *const Reaction::Mechanism() const
{
    return m_mech;
}

// Sets the parent mechanism.
void Reaction::SetMechanism(Sprog::Mechanism *const mech)
{
    m_mech = mech;
}


// CLONING.

// Returns a pointer to a copy of the reaction.
Reaction *Reaction::Clone(void) const
{
    return new Reaction(*this);
}


// MEMORY MANAGEMENT.

void Reaction::releaseMemory(void)
{
    m_name.clear();
    if (m_arrr != NULL) delete m_arrr;
    if (m_lt != NULL) delete m_lt;
    if (m_revlt != NULL) delete m_revlt;
    m_thirdbodies.clear();
}