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
struct SRI_PARAMS
{
    real A, B, C, D, E;
};

class SRIReaction : public FallOffReaction
{
protected:
    SRI_PARAMS m_params;
public:
    SRIReaction(void);  // Default constructor.
    ~SRIReaction(void); // Default destructor.
    SRIReaction(const SRIReaction &rxn); // Copy constructor.
public:
    SRIReaction &operator=(const SRIReaction &rxn);
public:
    /* Returns the SRI parameters. */
    const SRI_PARAMS &FallOffParams(void) const;
    /* Sets the SRI parameters. */
    void SetFallOffParams(const SRI_PARAMS &params);
};
};

#endif