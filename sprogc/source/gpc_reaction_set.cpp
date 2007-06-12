#include "gpc_reaction_set.h"
#include "gpc_third_body_reaction.h"
#include "gpc_fall_off_reaction.h"
#include "gpc_reaction.h"

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ReactionSet::ReactionSet()
{
}

// Copy constructor.

ReactionSet::ReactionSet(const Sprog::ReactionSet &rxn)
{
    *this = rxn;
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
        ThirdBodyRxnMap::const_iterator jtb;
        for (jtb=rxns.m_tb_rxns.begin(); jtb!=rxns.m_tb_rxns.end(); jtb++) {
            j = (*jtb).first;
            // Need to cast the reaction to a ThirdBodyReaction.
            m_tb_rxns.insert(ThirdBodyRxnMap::value_type(j, 
                             dynamic_cast<const ThirdBodyReaction*>(m_rxns[j])));
        }

        // Build fall-off reaction map.
        FallOffRxnMap::const_iterator jfo;
        for (jfo=rxns.m_fo_rxns.begin(); jfo!=rxns.m_fo_rxns.end(); jfo++) {
            j = (*jfo).first;
            // Need to cast the reaction to a FallOffReaction.
            m_fo_rxns.insert(FallOffRxnMap::value_type(j, 
                             dynamic_cast<const FallOffReaction*>(m_rxns[j])));
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
        ThirdBodyRxnMap::const_iterator jtb;
        for (jtb=rxns.m_tb_rxns.begin(); jtb!=rxns.m_tb_rxns.end(); jtb++) {
            j = n + (*jtb).first;
            // Need to cast the reaction to a ThirdBodyReaction.
            m_tb_rxns.insert(ThirdBodyRxnMap::value_type(j, 
                             dynamic_cast<const ThirdBodyReaction*>(m_rxns[j])));
        }

        // Build fall-off reaction map.
        FallOffRxnMap::const_iterator jfo;
        for (jfo=rxns.m_fo_rxns.begin(); jfo!=rxns.m_fo_rxns.end(); jfo++) {
            j = n + (*jfo).first;
            // Need to cast the reaction to a FallOffReaction.
            m_fo_rxns.insert(FallOffRxnMap::value_type(j, 
                             dynamic_cast<const FallOffReaction*>(m_rxns[j])));
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


// SET INFORMATION.

// Returns the number of reactions in the set.
unsigned int ReactionSet::Count(void) const
{
    return m_rxns.size();
}


// REACTIONS.

// Returns the vector of all reactions.
const RxnPtrVector &ReactionSet::Reactions() const
{
    return m_rxns;
}

// Adds a reaction to the set.
void ReactionSet::AddReaction(const Sprog::Reaction &rxn)
{
    // Clone the reaction and add it to the vector.
    Reaction *pr = rxn.Clone();
    m_rxns.push_back(pr);

    // Check for reverse parameters.
    if (rxn.RevArrhenius() != NULL) {
        m_rev_rxns.insert(RxnMap::value_type(m_rxns.size(), pr));
    }

    // Check for forward LT parameters.
    if (rxn.LTCoeffs() != NULL) {
        m_lt_rxns.insert(RxnMap::value_type(m_rxns.size(), pr));
    }

    // Check for reverse LT parameters.
    if (rxn.RevLTCoeffs() != NULL) {
        m_revlt_rxns.insert(RxnMap::value_type(m_rxns.size(), pr));
    }
}

// Adds a third-body reaction to the set.
void ReactionSet::AddReaction(const Sprog::ThirdBodyReaction &rxn)
{
    // Clone the reaction and add it to the vector.
    ThirdBodyReaction *pr = rxn.Clone();
    m_rxns.push_back(pr);

    // Add reaction to third-body reaction map.
    m_tb_rxns.insert(ThirdBodyRxnMap::value_type(m_rxns.size(), pr));
}

// Adds a fall-off reaction to the set.
void ReactionSet::AddReaction(const Sprog::FallOffReaction &rxn)
{
    // Clone the reaction and add it to the vector.
    FallOffReaction *pr = rxn.Clone();
    m_rxns.push_back(pr);

    // Add reaction to third-body reaction map.
    m_fo_rxns.insert(FallOffRxnMap::value_type(m_rxns.size(), pr));
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