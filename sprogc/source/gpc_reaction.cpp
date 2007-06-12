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


// REACTION NAME/DESRIPTION.

const string &Reaction::Name() const {return m_name;}
void Reaction::SetName(const std::string &name) {m_name = name;}


// REVERSIBLE REACTION?.

bool Reaction::IsReversible() const {return m_reversible;}
void Reaction::SetReversible(const bool isrev) {m_reversible=isrev;}


// REACTANTS.

const std::vector<STOICH> &Reaction::Reactants() const {return m_reac;}
const std::vector<FSTOICH> &Reaction::FReactants() const {return m_freac;}

// Adds an integer stoichiometry reactant to the reaction.
void Reaction::AddReactant(const Sprog::STOICH &reac)
{
    // We must check if the reactant is already defined in this reaction.  Need to
    // check both integer and real reactants.

    // Integer reactants.
    vector<STOICH>::iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        if ((*i).Species.Pointer == reac.Species.Pointer) {
            // We have found this species in the reactants already.  Append
            // the new contribution.
            (*i).Stoich += reac.Stoich;
            return;
        }
    }

    // Real reactants.
    vector<FSTOICH>::iterator j;
    for (j=m_freac.begin(); j!=m_freac.end(); j++) {
        if ((*j).Species.Pointer == reac.Species.Pointer) {
            // We have found this species in the reactants already.  Append
            // the new contribution.
            (*j).Stoich += (real)reac.Stoich;
            return;
        }
    }

    // If we have got here then the reactant is not defined for this reaction, so
    // we must add it.
    m_reac.push_back(reac);
}


// Adds a real stoichiometry reactant to the reaction.
void Reaction::AddReactant(const Sprog::FSTOICH &reac)
{
    // We must check if the reactant is already defined in this reaction.  Need to
    // check both integer and real reactants.

    // Integer reactants.
    vector<STOICH>::iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        if ((*i).Species.Pointer == reac.Species.Pointer) {
            // We have found this species in the reactants already.  It is currently
            // in the integer stoichiometry, and we now want it in the real stoichiometry.
            FSTOICH freac(reac);
            freac.Stoich += (real)(*i).Stoich;
            m_freac.push_back(freac);
            return;
        }
    }

    // Real reactants.
    vector<FSTOICH>::iterator j;
    for (j=m_freac.begin(); j!=m_freac.end(); j++) {
        if ((*j).Species.Pointer == reac.Species.Pointer) {
            // We have found this species in the reactants already.  Append
            // the new contribution.
            (*j).Stoich += reac.Stoich;
            return;
        }
    }

    // If we have got here then the reactant is not defined for this reaction, so
    // we must add it.
    m_freac.push_back(reac);
}


// ARRHENIUS COEFFICIENTS.

// Returns the explicit reverse Arrhenius parameters.
const ARRHENIUS &Reaction::RevArrhenius(void) const
{
    return *m_arrr;
}


// MEMORY MANAGEMENT.

void Reaction::releaseMemory(void)
{
    m_name.clear();
    if (m_arrr != NULL) delete m_arrr;
    if (m_lt != NULL) delete m_lt;
    if (m_revlt != NULL) delete m_revlt;
}