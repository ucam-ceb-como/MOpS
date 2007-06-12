/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a fall-off pressure dependent reaction
    which defines a low pressure limit.
*/

#ifndef GPC_FALL_OFF_REACTION_H
#define GPC_FALL_OFF_REACTION_H

#include "gpc_params.h"
#include "gpc_species.h"
#include "gpc_reaction.h"
#include "gpc_third_body_reaction.h"

using namespace std;

namespace Sprog
{
class FallOffReaction : public ThirdBodyReaction
{
protected:
    ARRHENIUS m_lowplimit;        // Low pressure limit Arrhenius coefficients.
    IndexedSpecies m_fothirdbody; // Third body for fall-off.
public:
    FallOffReaction(void);  // Default constructor.
    ~FallOffReaction(void); // Default destructor.
    FallOffReaction(const FallOffReaction &rxn); // Copy constructor.
public:
    FallOffReaction &operator=(const FallOffReaction &rxn);
public:
    /* Returns the low pressure limit Arrhenius coefficients. */
    const ARRHENIUS &LowPressureLimit(void) const;
    /* Sets the low pressure limit Arrhenius coefficients. */
    void SetLowPressureLimit(const ARRHENIUS &lowp);
public:
    /* Returns a pointer to the species used as a third body for fall-off calculations. */
    const Species *const FallOffThirdBody(void) const;
    /* Sets the species used as a third body for fall-off calculations. */
    void SetFallOffThirdBody(const IndexedSpecies);
protected:
    /* Calculation of the rate constant of fall-off reactions requires a factor F.  This
       factor is different for different forms of fall-off reactions, for this simple
       form the value is unity (1.0). */
    virtual real F(real T, real logpr);
};

/* Inline function definitions. */

inline real FallOffReaction::F(real T, real logpr) {return 1.0;};

};

#endif