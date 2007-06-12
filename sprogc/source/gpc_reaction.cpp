#include "gpc_params.h"
#include "gpc_stoich.h"
#include "gpc_rate_params.h"
#include "gpc_reaction.h"

using namespace Sprog;
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
}

// Copy constructor.
Reaction::Reaction(const Sprog::Reaction &rxn)
{
    // Nullify pointer members.
    m_arrr = NULL;
    m_lt = NULL;
    m_revlt = NULL;

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
Reaction &Reaction::operator=(const Sprog::Reaction &rxn)
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
        if ((*i).Species() == reac.Species()) {
            // We have found this species in the reactants already.  Append
            // the new contribution.
            (*i).IncMu(reac.Mu());
            return;
        }
    }

    // Real reactants.
    vector<Stoichf>::iterator j;
    for (j=m_freac.begin(); j!=m_freac.end(); j++) {
        if ((*j).Species()== reac.Species()) {
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
        if ((*i).Species() == reac.Species()) {
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
        if ((*j).Species() == reac.Species()) {
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

// Removes the reactant species with the given name from the reaction.  If the
// species is not a reactant then does nothing.
void Reaction::RemoveReactant(const std::string &name)
{
    // Search the integer stoichiometry reactants.
    vector<Stoich>::iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        if (name.compare((*i).Species()->Name())==0) {
            // We have found the species in the integer stoich vector.
            m_reac.erase(i);
            return;
        }
    }

    // Now search the real stoichiometry reactants.
    vector<Stoichf>::iterator j;
    for (j=m_freac.begin(); j!=m_freac.end(); j++) {
        if (name.compare((*j).Species()->Name())==0) {
            // We have found the species in the real stoich vector.
            m_freac.erase(j);
            return;
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
        if ((*i).Species() == prod.Species()) {
            // We have found this species in the products already.  Append
            // the new contribution.
            (*i).IncMu(prod.Mu());
            return;
        }
    }

    // Real products.
    vector<Stoichf>::iterator j;
    for (j=m_fprod.begin(); j!=m_fprod.end(); j++) {
        if ((*j).Species() == prod.Species()) {
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
        if ((*i).Species() == prod.Species()) {
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
        if ((*j).Species() == prod.Species()) {
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

// Removes the product species with the given name from the reaction.  If the
// species is not a product then does nothing.
void Reaction::RemoveProduct(const std::string &name)
{
    // Search the integer stoichiometry products.
    vector<Stoich>::iterator i;
    for (i=m_prod.begin(); i!=m_prod.end(); i++) {
        if (name.compare((*i).Species()->Name())==0) {
            // We have found the species in the integer stoich vector.
            m_prod.erase(i);
            return;
        }
    }

    // Now search the real stoichiometry products.
    vector<Stoichf>::iterator j;
    for (j=m_fprod.begin(); j!=m_fprod.end(); j++) {
        if (name.compare((*j).Species()->Name())==0) {
            // We have found the species in the real stoich vector.
            m_fprod.erase(j);
            return;
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
void Reaction::SetRevArrhenius(const Sprog::ARRHENIUS &arr)
{
    if (m_arrr == NULL) {
        m_arrr = new ARRHENIUS(arr);
    } else {
        *m_arrr = arr;
    }
}


// LANDAU TELLER COEFFICIENTS.

// Sets the forward Landau Teller rate parameters.
void Reaction::SetLTCoeffs(const Sprog::LTCOEFFS &lt)
{
    if (m_lt == NULL) {
        m_lt = new LTCOEFFS(lt);
    } else {
        *m_lt = lt;
    }
}

// Sets the reverse Landau Teller rate parameters.
void Reaction::SetRevLTCoeffs(const Sprog::LTCOEFFS &lt)
{
    if (m_revlt == NULL) {
        m_revlt = new LTCOEFFS(lt);
    } else {
        *m_revlt = lt;
    }
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