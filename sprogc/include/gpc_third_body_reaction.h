/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a reaction with a third body.
*/

#ifndef GPC_THIRD_BODY_REACTION_H
#define GPC_THIRD_BODY_REACTION_H

#include <vector>
#include "gpc_params.h"
#include "gpc_species.h"
#include "gpc_reaction.h"

using namespace std;

namespace Sprog
{
class ThirdBodyReaction : public Reaction
{
protected:
    vector<FSTOICH> m_thirdbodies; // Reaction third bodies and their coefficients.
public:
    ThirdBodyReaction(void);  // Default constructor.
    ~ThirdBodyReaction(void); // Default destructor.
    ThirdBodyReaction(const ThirdBodyReaction &rxn); // Copy constructor.
public:
    ThirdBodyReaction &operator=(const ThirdBodyReaction &rxn);
public:
    /* Returns the vector of third-body coefficients. */
    const vector<FSTOICH> &ThirdBodies(void) const;
    /* Adds a third body to the reaction. */
    void AddThirdBody(const FSTOICH &tb);
    /* Adds a third body to the reaction (real). */
    void AddThirdBody(const Species *sp, real coeff);
    /* Removes a third body, given by name, from the reaction. */
    void RemoveThirdBody(const string &name);

};
};

#endif