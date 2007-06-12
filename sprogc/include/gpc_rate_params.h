/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains structures which define reaction rate parameter sets.
*/

#ifndef GPC_RATE_PARAMS_H
#define GPC_RATE_PARAMS_H

#include "gpc_params.h"

namespace Sprog
{
    // Arrhenius parameters.
    struct ARRHENIUS
    {
        real A; // Pre-exponential factor.
        real n; // Temperature exponent.
        real E; // Activation energy.
    };

    // Landau Teller reaction parameters.
    struct LTCOEFFS
    {
        real B, C;
    };
};

#endif