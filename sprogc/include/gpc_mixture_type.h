/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    Definition of enumeration of different mixture types.  Required
    for mixture serialisation.
*/

#ifndef GPC_MIXTURE_TYPE_H
#define GPC_MIXTURE_TYPE_H

namespace Sprog
{
namespace Thermo
{
    enum Serial_MixtureType {Serial_Mixture, Serial_GasPhase, Serial_IdealGas};
};
};
#endif
