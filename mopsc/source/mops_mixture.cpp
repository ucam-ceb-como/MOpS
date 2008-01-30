#include "mops_mixture.h"

using namespace Mops;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Mixture::Mixture(void)
{
}

// Default constructor (public, requires species list).
Mixture::Mixture(const Sprog::SpeciesPtrVector &sp)
: Sprog::Thermo::IdealGas(sp)
{
}

// Stream-reading constructor.
Mixture::Mixture(std::istream &in, const Sprog::SpeciesPtrVector &sp)
: Sprog::Thermo::IdealGas(in, sp)
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