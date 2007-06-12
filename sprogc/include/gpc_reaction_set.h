/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a structure for a set of chemical reactions, 
    including code to provide speed enhancements when working with many reactions.
*/

#ifndef GPC_REACTION_SET_H
#define GPC_REACTION_SET_H

#include <vector>
#include <map>
#include "gpc_params.h"
#include "gpc_reaction.h"
#include "gpc_third_body_reaction.h"
#include "gpc_fall_off_reaction.h"

using namespace std;

namespace Sprog
{
/* In the following maps the key is the index in the vector of all reactions. */
typedef map<unsigned int,const Reaction*> RxnMap;
typedef map<unsigned int,const ThirdBodyReaction*> ThirdBodyRxnMap;
typedef map<unsigned int,const FallOffReaction*> FallOffRxnMap;

class ReactionSet
{
protected:
    RxnPtrVector m_rxns;       // Vector of all reactions in the set.
    RxnMap m_rev_rxns;         // Map of reactions which have explicit reverse Arrhenius parameters.
    ThirdBodyRxnMap m_tb_rxns; // Map of third body reactions.
    FallOffRxnMap m_fo_rxns;   // Map of fall-off reactions.
    RxnMap m_lt_rxns;          // Map of reactions with Landau Teller parameters.
    RxnMap m_revlt_rxns;       // Map of reactions with reverse Landau Teller parameters.
public:
    ReactionSet(void);  // Default constructor.
    ~ReactionSet(void); // Default destructor.
    ReactionSet(const ReactionSet &rxn); // Copy constructor.
public:
    ReactionSet &operator=(const ReactionSet &rxns);
    ReactionSet &operator+=(const ReactionSet &rxns);
    const ReactionSet operator+(const ReactionSet &rxns) const;
    Reaction &operator[](unsigned int i);
public:
    /* Returns the number of reactions in the set. */
    unsigned int Count(void) const;
public:
    /* Returns the list of reactions. */
    const vector<Reaction> &Reactions(void) const;
    /* Adds a reaction to the set. */
    void AddReaction(const Reaction &rxn);
    /* Adds a third body reaction to the set. */
    void AddReaction(const ThirdBodyReaction &rxn);
    /* Adds a fall-off reaction to the set. */
    void AddReaction(const FallOffReaction &rxn);
};
};

#endif