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
namespace Kinetics
{
    // Arrhenius parameters.
    struct ARRHENIUS
    {
        real A; // Pre-exponential factor.
        real n; // Temperature exponent.
        real E; // Activation energy.

        // Constructors.
        ARRHENIUS(void) {A=n=E=0.0;}; // Default constructor.
        ARRHENIUS(real aA, real an, real aE) {A=aA; n=an; E=aE;}; // Initialising constructor.
    };

    // Landau Teller reaction parameters.
    struct LTCOEFFS
    {
        real B, C;

        // Constructors.
        LTCOEFFS(void) {B=C=0.0;}; // Default constructor.
        LTCOEFFS(real aB, real aC) {B=aB; C=aC;}; // Initialising constructor.
    };
};
};
#endif