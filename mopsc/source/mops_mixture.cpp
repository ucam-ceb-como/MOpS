#include "mops_mixture.h"

using namespace Mops;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Mixture::Mixture(void)
{
}

// Default constructor (public, requires species list).
Mixture::Mixture(const Sprog::SpeciesPtrVector &sp)
: Sweep::Cell(sp)
{
}

// Stream-reading constructor.
Mixture::Mixture(std::istream &in, const Mechanism &mech)
: Sweep::Cell(in, mech.ParticleMech())
{
}

// Default destructor.
Mixture::~Mixture(void)
{
}


// READ/WRITE/COPY.

// Creates a clone of the mixture.
Mixture *const Mixture::Clone() const
{
    return new Mixture(*this);
}