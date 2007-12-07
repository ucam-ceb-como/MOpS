/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a structure for a set of chemical reactions, 
    including code to provide speed enhancements when working with many reactions.
*/

#ifndef GPC_REACTION_SET_H
#define GPC_REACTION_SET_H

#include "gpc_params.h"
#include "gpc_reaction.h"
#include "gpc_third_body_reaction.h"
#include "gpc_fall_off_reaction.h"
#include <vector>
#include <map>

namespace Sprog
{
namespace Kinetics
{
class ReactionSet
{
public:
    // In the following maps the key is the index in the vector of all reactions.
    typedef std::map<unsigned int,const Reaction*> RxnMap;
    typedef std::map<unsigned int,const ThirdBodyReaction*> ThirdBodyRxnMap;
    typedef std::map<unsigned int,const FallOffReaction*> FallOffRxnMap;

    // Constructors.
    ReactionSet(void); // Default constructor.
    ReactionSet(const ReactionSet &rxn); // Copy constructor.

    // Destructor.
    virtual ~ReactionSet(void);

    // Operator overloads.
    ReactionSet &operator=(const ReactionSet &rxns);
    ReactionSet &operator+=(const ReactionSet &rxns);
    const ReactionSet operator+(const ReactionSet &rxns) const;
    Reaction *const operator[](unsigned int i);
   
    // Set information.
    unsigned int Count(void) const; // Returns the number of reactions in the set.

    // Reactions.
    const RxnPtrVector &Reactions(void) const;                 // Returns the list of reactions.
    Reaction *const AddReaction(const Reaction &rxn);          // Adds a reaction to the set.
    Reaction *const AddReaction(const ThirdBodyReaction &rxn); // Adds a third body reaction to the set.
    Reaction *const AddReaction(const FallOffReaction &rxn);   // Adds a fall-off reaction to the set.

    // Tidying up.
    void Clear(void); // Clears all reactions from the set.

protected:
    // Reaction set data.
    RxnPtrVector m_rxns;       // Vector of all reactions in the set.
    RxnMap m_rev_rxns;         // Map of reactions which have explicit reverse Arrhenius parameters.
    ThirdBodyRxnMap m_tb_rxns; // Map of third body reactions.
    FallOffRxnMap m_fo_rxns;   // Map of fall-off reactions.
    RxnMap m_lt_rxns;          // Map of reactions with Landau Teller parameters.
    RxnMap m_revlt_rxns;       // Map of reactions with reverse Landau Teller parameters.

    // Memory management.
    virtual void releaseMemory(void); // Clears all memory used by the set.
};
};
};

#endif