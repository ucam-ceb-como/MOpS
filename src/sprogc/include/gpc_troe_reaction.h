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

namespace Sprog
{
namespace Kinetics
{
class TroeReaction : public FallOffReaction
{
public:
    // A structure to maintain Troe fall-off parameters.
    struct TroeParams
    {
        // Data
        real Alpha;
        real Ts, Tsss;
        real *Tss;
        // Constructors.
        TroeParams(void);
        TroeParams(const TroeParams &troe);
        // Destructor.
        ~TroeParams(void);
        // Operator overloads.
        TroeParams &operator=(const TroeParams &troe);
    };

    // Constructors.
    TroeReaction(void); // Default constructor.
    TroeReaction(const TroeReaction &rxn); // Copy constructor.

    // Destructor.
    ~TroeReaction(void);

    // Operator overloads.
    TroeReaction &operator=(const TroeReaction &rxn);

    // Troe parameters.
    const TroeParams &FallOffParams(void) const;     // Returns the Troe parameters.
    void SetFallOffParams(const TroeParams &params); // Sets the Troe parameters.

    // Cloning.
    TroeReaction* Clone(void) const; // Returns a pointer to a copy of the Troe reaction.

protected:
    TroeParams m_params; // TROE form parameters.
};
};
};
#endif