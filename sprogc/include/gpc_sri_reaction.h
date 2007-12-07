/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of an SRI form pressure fall-off reaction.
*/

#ifndef GPC_SRI_REACTION_H
#define GPC_SRI_REACTION_H

#include "gpc_params.h"
#include "gpc_fall_off_reaction.h"

using namespace std;

namespace Sprog
{
namespace Kinetics
{
class SRIReaction : public FallOffReaction
{
public:
    // A structure to maintain SRI fall-off parameters.
    struct SRI_PARAMS
    {
        real A, B, C, D, E;
    };

    // Constructors.
    SRIReaction(void);  // Default constructor.
    SRIReaction(const SRIReaction &rxn); // Copy constructor.

    // Destructor.
    ~SRIReaction(void);

    // Operator overloads.
    SRIReaction &operator=(const SRIReaction &rxn);

    // SRI parameters.
    const SRI_PARAMS &FallOffParams(void) const;     // Returns the SRI parameters.
    void SetFallOffParams(const SRI_PARAMS &params); // Sets the SRI parameters.

    // Cloning.
    SRIReaction* Clone(void) const; // Returns a pointer to a copy of the SRI reaction.

protected:    
    SRI_PARAMS m_params; // SRI parameters.
};
};
};

#endif