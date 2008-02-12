#include "swp_particle.h"
#include "swp_model.h"

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
Particle::Particle(void)
{
}

// Default constructor (public).
Particle::Particle(const CompPtrVector &components, const TrackPtrVector &trackers)
: ParticleData(components, trackers)
{
}

// Copy constructor.
Particle::Particle(const Sweep::Particle &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
Particle::Particle(std::istream &in, const Mechanism &mech)
: ParticleData(in, mech)
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
        ParticleData::operator=(rhs);
        m_ensemble = rhs.m_ensemble;
    }
    return *this;
}

// Addition-assignment operator.  This implements coagulation.
Particle &Particle::operator+=(const Sweep::Particle &rhs)
{
    return Coagulate(rhs);
}

// Addition operator.  This also implements coagulation.
const Particle Particle::operator +(const Sweep::Particle &rhs) const
{
    return Particle(*this) += rhs;
}


// PARTICLE ADJUSTMENT AND PROCESSES.

// Adjusts the particle with the given composition and value
// changes.
void Particle::Adjust(const fvector &dcomp, const fvector &dvalues)
{
	unsigned int i = 0;

	// Add the components.
	for (i=0; i!=m_comp.size(); ++i) {
		m_comp[i] += dcomp[i];
	}

	// Add the tracker values.
	for (i=0; i!=m_values.size(); ++i) {
		m_values[i] += dvalues[i];
	}

	// Now allow the particle models to deal with the additions.
    for (ModelMap::iterator j=m_models.begin(); j!=m_models.end(); ++j) {
        j->second->Model().UpdateParticle(*this, dcomp, dvalues);
	}

	// Allow coagulation model to deal with the addition as well.
	m_coag->Model().UpdateParticle(*this, dcomp, dvalues);

    // Update particle cache.
    UpdateCache();
}

// Adjusts the particle with the given composition and 
// values changes n times.
unsigned int Particle::Adjust(const fvector &dcomp, 
                              const fvector &dvalues, 
                              unsigned int n)
{
	unsigned int i = 0;

	// Add the components.
	for (i=0; i!=m_comp.size(); ++i) {
		m_comp[i] += dcomp[i] * (real)n;
	}

	// Add the tracker values.
	for (i=0; i!=m_values.size(); ++i) {
		m_values[i] += dvalues[i] * (real)n;
	}

	// Now allow the particle models to deal with the additions.
    for (ModelMap::iterator j=m_models.begin(); j!=m_models.end(); ++j) {
        j->second->Model().UpdateParticle(*this, dcomp, dvalues, n);
	}

	// Allow coagulation model to deal with the addition as well.
	m_coag->Model().UpdateParticle(*this, dcomp, dvalues, n);

    // Update particle cache.
    UpdateCache();

    return n;
}

// Combines this particle with another.
Particle &Particle::Coagulate(const Particle &sp)
{
	unsigned int i = 0;

	// Add the components.
	for (i=0; i!=m_comp.size(); ++i) {
		m_comp[i] += sp.m_comp[i];
	}

	// Add the tracker values.
	for (i=0; i!=m_values.size(); ++i) {
		m_values[i] += sp.m_values[i];
	}

	// Now allow the particle models to deal with the coagulation.
    for (ModelMap::iterator j=m_models.begin(); j!=m_models.end(); ++j) {
        j->second->Model().CoagParticles(*this, sp);
	}

	// Allow coagulation model to deal with the coagulation as well.
	m_coag->Model().CoagParticles(*this, sp);

    // Update particle cache.
    UpdateCache();

    return *this;
}


// PARTICLE UPDATE AND CHECKING.

// Recalculates the derived properties from the unique properties.
void Particle::UpdateCache(void)
{
    real m = 0.0;

    // Loop over composition and calculate mass and volume.
    m_mass = m_vol = 0.0;
    for (unsigned int i=0; i!=m_components->size(); ++i) {
        m = (*m_components)[i]->MolWt() * m_comp[i] / NA;
        m_mass += m;
        m_vol  += m / (*m_components)[i]->Density();
    }
    
    // Calculate other properties.
    m_diam = pow(6.0 * m_vol / PI, ONE_THIRD);
    m_dcol = m_diam;
    m_dmob = m_diam;
    m_surf = PI * m_diam * m_diam;

    // Allow models to update themselves and the particle data.
    for (ModelMap::iterator i=m_models.begin(); i!=m_models.end(); ++i) {
        i->second->Model().UpdateCache(*this);
    }

    // Allow coagulation model to update the particle data.
    m_coag->Model().UpdateCache(*this);
}

// Check the that the particle is valid by querying the
// validity conditions of the models and ensuring that it 
// contains any components.
bool Particle::IsValid() const
{
    fvector::const_iterator i;
    for (i=m_comp.begin(); i!=m_comp.end(); ++i) {
        if (*i > 0.0) return true;
    }
    return false;
}


// PARENT ENSEMBLE.

// Returns the parent ensemble.
const Sweep::Ensemble *const Particle::Ensemble(void) const
{
    return m_ensemble;
}

// Sets the parent ensemble.
void Particle::SetEnsemble(Sweep::Ensemble &ens)
{
    m_ensemble = &ens;
}


// READ/WRITE/COPY.

// Creates a clone of the particle.
Particle *const Particle::Clone()
{
    return new Particle(*this);
}
