/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of Troe form pressure fall-off reaction.
*/

#ifndef GPC_TROE_REACTION_H
#define GPC_TROE_REACTION_H

#include "gpc_params.h"
#include "gpc_fall_off_reaction.h"

using namespace std;

namespace Sprog
{
struct TroeParams
{
    real Alpha;
    real Ts, Tsss;
    real *Tss;
    TroeParams(void);
    ~TroeParams(void);
    TroeParams(const TroeParams &troe);
    TroeParams &operator=(const TroeParams &troe);
};

class TroeReaction : public FallOffReaction
{
protected:
    TroeParams m_params; // TROE form parameters.
public:
    TroeReaction(void);  // Default constructor.
    ~TroeReaction(void); // Default destructor.
    TroeReaction(const TroeReaction &rxn); // Copy constructor.
public:
    TroeReaction &operator=(const TroeReaction &rxn);
public:
    /* Returns the Troe parameters. */
    const TroeParams FallOffParams(void) const;
    /* Sets the Troe parameters. */
    void SetFallOffParams(const TroeParams &params);
};
};

#endif