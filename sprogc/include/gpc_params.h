/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    Defines parameters and typedefs used by sprog.
*/

#ifndef GPC_PARAMS_H
#define GPC_PARAMS_H

#include <vector>

namespace Sprog
{
    // COMMON TYPEDEFS.

    // Real numbers in Sprog, so they can be easily changed.
    typedef double real;

    // Real number STL vector.
    typedef std::vector<real> fvector; 


    // MATHEMATICAL CONSTANTS.

    const real PI = 3.1415926535897932384626433832795;
    const real ONE_THIRD  = 3.3333334e-01;
    const real TWO_THIRDS = 6.6666667e-01;


    // PHYSICAL CONSTANTS.

    // Avogadro's number (source = NIST website, physics.nist.gov).
    // Error = 3.0e16 /mol.
    const real NA    = 6.02214179e23; // 1/mol.

    // Gas constant (source = NIST website, physics.nist.gov).
    // Error = 1.5e-5 J/mol/K.
    const real R     = 8.314472e0; // J/mol/K   (SI).
    const real R_CGS = 8.314472e7; // ergs/molK (CGS).
    const real RCAL  = 1.9872e-3;  // kcal/molK (calories).
};

#endif
