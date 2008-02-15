#include "mops_mechanism.h"

using namespace Mops;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Mechanism::Mechanism(void)
{
}

// Default destructor.
Mechanism::~Mechanism(void)
{
}


// PARTICLE MECHANISM.

// Returns a reference (non-const) to the particle mechanism.
Sweep::Mechanism &Mechanism::ParticleMech(void) {return m_pmech;}
const Sweep::Mechanism &Mechanism::ParticleMech(void) const {return m_pmech;}
