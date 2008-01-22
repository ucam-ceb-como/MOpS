/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The Mechanism class combines the gas-phase mechanism of Sprog with
    the particle mechanism of Sweep.
    
    Additionally the Mixture class includes a description of a particle
    population which can be solved with Sweep.
*/

#ifndef MOPS_MECHANISM_H
#define MOPS_MECHANISM_H

#include "sprog.h"
#include "mops_params.h"

namespace Mops
{
class Mechanism : public Sprog::Mechanism
{
public:
    // Constructors.
    Mechanism(void); // Default constructor.

    // Destructors.
    ~Mechanism(void); // Default destructors.
};
};

#endif
