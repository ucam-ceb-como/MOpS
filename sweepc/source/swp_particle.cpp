#include "swp_particle.h"
#include "swp_particle_model.h"
#include "swp_submodel.h"
#include <cmath>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Particle::Particle(void)
: m_ensemble(NULL)
{
}

// Initialising constructor.
Particle::Particle(real time, const Sweep::ParticleModel &model)
: SubParticle(time, model), m_ensemble(NULL)
{
}

// Initialising constructor (from Primary particle).
Particle::Particle(Sweep::Primary &pri)
: SubParticle(pri), m_ensemble(NULL)
{
}

// Copy constructor.
Particle::Particle(const Sweep::Particle &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
Particle::Particle(std::istream &in, const Sweep::ParticleModel &model)
: SubParticle(in, model)
{
}

// Default destructor.
Particle::~Particle()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
Particle &Particle::operator=(const Sweep::Particle &rhs)
{
    if (this != &rhs) {
        SubParticle::operator=(rhs);
        m_ensemble = rhs.m_ensemble;
    }
    return *this;
}

// Addition-assignment operator.  This implements coagulation.
Particle &Particle::operator+=(const Sweep::Particle &rhs)
{
    SubParticle::operator+=(rhs);
    return *this;
}

// Addition operator.  This also implements coagulation.
const Particle Particle::operator +(const Sweep::Particle &rhs) const
{
    return Particle(*this) += rhs;
}


// PARENT ENSEMBLE.

// Returns the parent ensemble.
const Sweep::Ensemble *const Particle::Ensemble(void) const {return m_ensemble;}

// Sets the parent ensemble.
void Particle::SetEnsemble(Sweep::Ensemble &ens) {m_ensemble = &ens;}


// READ/WRITE/COPY.

// Creates a clone of the particle.
Particle *const Particle::Clone() const
{
    return new Particle(*this);
}
