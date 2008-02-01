/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Defines parameters and typedefs used by sweep.
*/

#ifndef SWEEP_PARAMS_H
#define SWEEP_PARAMS_H

#include <cmath>

namespace Sweep
{
    typedef double real;

    const real PI = 3.1415926535897932384626433832795;
    const real ONE_THIRD  = 3.3333334e-01;
    const real TWO_THIRDS = 6.6666667e-01;

    // Constant term in calculation on Knudsen number.  Actual
    // Kn = K * T [K] / (P [bar] * d [cm]).
    const real KNUDSEN_K = 4.74151636e-8;
 
    // Coagulation kernel parameters.
    const real CFM = 3.23990504e-8;
    const real CSF = 9.2046667e-17;
    const real CFMMAJ = 1.4178;

    // Avogadro's constant.
    const real NA = 6.022e23;

    // Gas constant.
    const real R    = 8.3145E0;  // J/molK.
    const real RCAL = 1.9872E-3; // kcal/molK.

    inline real Viscosity(real T)
    {
        return 14.58e-6 * pow(T,1.5) / (T + 100.4);
    };

};

#endif
