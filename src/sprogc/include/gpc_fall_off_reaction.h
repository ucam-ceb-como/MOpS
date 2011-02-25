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
#include "gpc_rate_params.h"
#include "gpc_reaction.h"
#include "gpc_third_body_reaction.h"

namespace Sprog
{
namespace Kinetics
{
class FallOffReaction : public ThirdBodyReaction
{
public:
    // Constructors.
    FallOffReaction(void); // Default constructor.
    FallOffReaction(const FallOffReaction &rxn); // Copy constructor.

    // Destructor.
    virtual ~FallOffReaction(void);

    // Operator overloads.
    FallOffReaction &operator=(const FallOffReaction &rxn);

    // Low pressure limit.
    const ARRHENIUS &LowPressureLimit(void) const;   // Returns the low pressure limit Arrhenius coefficients.   
    void SetLowPressureLimit(const ARRHENIUS &lowp); // Sets the low pressure limit Arrhenius coefficients.

    // Fall-off third body.   
    const Sprog::Species &FallOffThirdBody(void) const; // Returns a pointer to the species used as a 
                                                        // third body for fall-off calculations.
    void SetFallOffThirdBody(unsigned int sp); // Sets the species used as a third body for 
                                               // fall-off calculations.
    void SetFallOffThirdBody(const std::string &name); // Sets the species used as a third body for 
                                                       // fall-off calculations given the species name.

    // Cloning.
    virtual FallOffReaction* Clone(void) const; // Returns a pointer to a copy of the fall-off reaction.
protected:
    // Fall-off data.
    ARRHENIUS m_lowplimit;   // Low pressure limit Arrhenius coefficients.
    int m_fothirdbody;       // Third body for fall-off.

    // Calculation of the rate constant of fall-off reactions requires a factor F.  This
    // factor is different for different forms of fall-off reactions, for this simple
    // form the value is unity (1.0).
    virtual real F(real T, real logpr);
};

// Inline function definitions.
inline real FallOffReaction::F(real T, real logpr) {return 1.0;};
};
};

#endif