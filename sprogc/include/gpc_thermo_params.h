/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a structure for chemical species thermodynamic
    fitting parameters.  Thermo parameters are split into ranges defined by the temperatures
    for which the ranges are valid.
*/

#ifndef GPC_THERMO_PARAMS_H
#define GPC_THERMO_PARAMS_H

#include <map>
#include "gpc_params.h"

namespace Sprog
{
namespace Thermo
{
    // Maximum lengths of thermo parameter arrays (S_PARAM_COUNT).
    const unsigned int CP_PARAM_COUNT    = 5;
    const unsigned int H_PARAM_COUNT     = 6;
    const unsigned int S_PARAM_COUNT     = 7;
    const unsigned int MAX_THERMO_PARAMS = 10;

    struct THERMO_PARAMS
    {
        unsigned int Count;
        real Params[S_PARAM_COUNT];
    };

    typedef std::map<real, THERMO_PARAMS> ThermoMap;
};
};
#endif
