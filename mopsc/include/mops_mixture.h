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
#include "mops_mechanism.h"
#include "sprog.h"
#include "sweep.h"
#include <iostream>

namespace Mops
{
class Mixture : public Sweep::Cell
{
public:
    // Constructors.
    Mixture(const Sprog::SpeciesPtrVector &sp); // Default constructor.
    Mixture(                  // Stream-reading constructor.
        std::istream &in,     //   - Input stream.
        const Mechanism &mech //   - Species which define the mixture.
        );

    // Destructors.
    ~Mixture(void); // Defaul destructors.


    // READ/WRITE/COPY.

    // Creates a clone of the mixture object.
    Mixture *const Clone() const;

protected:
    // As in Sprog, it is meaningless to define a mixture without knowledge
    // of the constitent species, therefore the default constructor is 
    // declared as protected.
    Mixture(void);
};
};

#endif
