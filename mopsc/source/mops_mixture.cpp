#include "mops_mixture.h"

using namespace Mops;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Mixture::Mixture(void)
{
}

// Default constructor (public, requires species list).
Mixture::Mixture(const Sweep::ParticleModel &model)
: Sweep::Cell(model)
{
}

// Stream-reading constructor.
Mixture::Mixture(std::istream &in, const Sweep::ParticleModel &model)
: Sweep::Cell(in, model)
{
}

// Default destructor.
Mixture::~Mixture(void)
{
}


// OPERATORS.

// Assignment operator.
Mixture &Mixture::operator =(const Mixture &rhs)
{
    if (this != &rhs) {
        // Invoke operator of base class.
        Sweep::Cell::operator=(rhs);
    }
    return *this;
}


// READ/WRITE/COPY.

// Creates a clone of the mixture.
Mixture *const Mixture::Clone() const
{
    return new Mixture(*this);
}
