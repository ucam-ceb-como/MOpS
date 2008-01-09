/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    Defines parameters and typedefs used by sprog.
*/

#ifndef GPC_PARAMS_H
#define GPC_PARAMS_H

namespace Sprog
{
    // Common typedefs.  
    typedef double real;  // typedef for real numbers in Sprog, so they can be easily changed.

    // Mathematical constants.
    const real PI = 3.1415926535897932384626433832795;
    const real ONE_THIRD  = 3.3333334e-01;
    const real TWO_THIRDS = 6.6666667e-01;

    // Physical constants.
    const real NA   = 6.0221367e23;  // Avogadro's number.
    const real R    = 8.31451e0;     // Gas constant SI (J/molK).
    const real R_CGS = 8.31451e7;    // Gas constant CGS (ergs/molK).
    const real RCAL = 1.9872e-3;     // Gas constant Calories (kcal/molK).
};

#endif