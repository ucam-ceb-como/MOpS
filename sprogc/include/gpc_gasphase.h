/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    The GasPhase class is a child of the Mixture class that also has thermodynamic
    functionality.  In addition to maintaining a description of the mixture composition
    it also allows thermodynamic properties of the gas to be calculated.  This is intended
    as a base class for other classes which will define the equation of state used
    to define the mixture, and hence how the properties should be calculated.
*/

#ifndef GPC_GASPHASE_H
#define GPC_GASPHASE_H

#include "gpc_params.h"
#include "gpc_thermo.h"
#include "gpc_mixture.h"
#include <vector>

namespace Sprog
{
namespace Thermo
{
class GasPhase : public Sprog::Thermo::ThermoInterface, public Mixture
{
public:
    // Constructors.
    GasPhase(void); // Default constructor.

    // Destructors.
    ~GasPhase(void); // Default destructor.
};
};
};

#endif