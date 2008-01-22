/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The Mixture class is the base class for gas-phase mixtures used
    by mops.  As mops only solves ideal gas systems, it inherits from
    the Sprog::IdealGas class.
    
    Additionally the Mixture class includes a description of a particle
    population which can be solved with Sweep.
*/

#ifndef MOPS_MIXTURE_H
#define MOPS_MIXTURE_H

#include "mops_params.h"
#include "sprog.h"

namespace Mops
{
class Mixture : public Sprog::Thermo::IdealGas
{
public:
    // Constructors.
    Mixture(void); // Default constructor.
    
    // Destructors.
    ~Mixture(void); // Defaul destructors.

    void Normalise();
};
};

#endif
