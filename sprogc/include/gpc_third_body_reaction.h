/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a reaction with a third body.
*/

#ifndef GPC_THIRD_BODY_REACTION_H
#define GPC_THIRD_BODY_REACTION_H

#include "gpc_params.h"
#include "gpc_species.h"
#include "gpc_reaction.h"
#include "gpc_stoich.h"
#include <vector>

namespace Sprog
{
namespace Kinetics
{
class ThirdBodyReaction : public Reaction
{
public:
    // Constructors.
    ThirdBodyReaction(void);  // Default constructor.
    ThirdBodyReaction(const ThirdBodyReaction &rxn); // Copy constructor.

    // Destructor.
    virtual ~ThirdBodyReaction(void);

    // Operator overloads.
    ThirdBodyReaction &operator=(const ThirdBodyReaction &rxn);

    // Third bodies.
    const std::vector<Stoichf> &ThirdBodies(void) const;    // Returns the vector of third-body coefficients.
    void AddThirdBody(const Stoichf &tb);                   // Adds a third body to the reaction.
    void AddThirdBody(unsigned int sp, real coeff);         // Adds a third body to the reaction.
    void AddThirdBody(const std::string &name, real coeff); // Adds a third body given the species name.
    void RemoveThirdBody(const std::string &name);          // Removes a third body, given by name, from 
                                                            // the reaction.

    // Cloning.
    virtual ThirdBodyReaction* Clone(void) const; // Returns a pointer to a copy of the third-body reaction.

protected:
    std::vector<Stoichf> m_thirdbodies; // Reaction third bodies and their coefficients.
};
};
};
#endif