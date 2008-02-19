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

#include "mops_params.h"
#include "sprog.h"
#include "sweep.h"

namespace Mops
{
class Mechanism : public Sprog::Mechanism
{
public:
    // Constructors.
    Mechanism(void); // Default constructor.

    // Destructors.
    ~Mechanism(void); // Default destructors.

    // PARTICLE MECHANISM.

    // Returns a reference (non-const) to the particle mechanism.
    Sweep::Mechanism &ParticleMech(void);
    const Sweep::Mechanism &ParticleMech(void) const;


    // READ/WRITE/COPY FUNCTIONS.

    // Writes the mechanism to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the mechanism data from a binary data stream.
    void Deserialize(std::istream &in);

private:
    // The particle mechanism.
    Sweep::Mechanism m_pmech;
};
};

#endif
