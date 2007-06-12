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
    typedef unsigned int index;  // typedef for indexing things.

    // Mathematical constants.
    const real PI = 3.1415926535897932384626433832795;
    const real ONE_THIRD  = 3.3333334e-01;
    const real TWO_THIRDS = 6.6666667e-01;

    // Physical constants.
    const real NA   = 6.022e23;  // Avogadro's number.
    const real R    = 8.3145e0;  // Gas constant SI (J/molK).
    const real RCAL = 1.9872e-3; // Gas constant CGS (kcal/molK).
};

#endif