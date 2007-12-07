#include "gpc_params.h"
#include "gpc_stoich.h"
#include "gpc_rate_params.h"
#include "gpc_reaction.h"
#include "gpc_mech.h"
#include <exception>
#include <stdexcept>
#include <string>

using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Reaction::Reaction(void)
{
    m_name = "";
    m_reversible = false;
    m_dstoich = m_dreac = m_dprod = 0.0;
    m_arrf.A = m_arrf.n = m_arrf.E = 0.0;
    m_arrr = NULL;
    m_lt = NULL;
    m_revlt = NULL;
    m_species = NULL;
}

// Copy constructor.
Reaction::Reaction(const Sprog::Kinetics::Reaction &rxn)
{
    // Nullify pointer members.
    m_arrr = NULL;
    m_lt = NULL;
    m_revlt = NULL;
    m_species = NULL;

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

        // Copy pointer to species vector.
        m_species = rxn.m_species;
    }

    return *this;
}


// REACTANTS.

// Adds an integer stoichiometry reactant to the reaction.
void Reaction::AddReactant(const Sprog::Stoich &reac)
{
    // We must check if the reactant is already defined in this reaction.  Need to
    // check both integer and real reactants.

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
                    // We have found the species in the list.
                    m_reac.erase(j);
                }
            }        

            // Now search the real stoichiometry reactants.
            vector<Stoichf>::iterator k;
            for (k=m_freac.begin(); k!=m_freac.end(); k++) {
                if ((*k).Index() == i) {
                    // We have found the species in the real stoich vector.
                    m_freac.erase(k);
                    return;
                }
            }
        }
    }
}


// PRODUCTS.

// Adds an integer stoichiometry product to the reaction.
void Reaction::AddProduct(const Sprog::Stoich &prod)
{
    // We must check if the product is already defined in this reaction.  Need to
    // check both integer and real products.

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
                    // We have found the species in the list.
                    m_prod.erase(j);
                }
            }        

            // Now search the real stoichiometry products.
            vector<Stoichf>::iterator k;
            for (k=m_fprod.begin(); k!=m_fprod.end(); k++) {
                if ((*k).Index() == i) {
                    // We have found the species in the real stoich vector.
                    m_fprod.erase(k);
                    return;
                }
            }
        }
    }
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
}