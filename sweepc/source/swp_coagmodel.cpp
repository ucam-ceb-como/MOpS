#include "swp_coagmodel.h"
#include "swp_particle.h"

using namespace Sweep;

// SINGLETON IMPLEMENTATION.

// Default constructor (private).
CoagModel::CoagModel(void)
{
    // Nothing to do here.
}

// Copy constructor.
CoagModel::CoagModel(const Sweep::CoagModel &copy)
{
    // Nothing to do here.
}

// Default destructor (protected).
CoagModel::~CoagModel(void)
{
    // Nothing special to destruct.
}

// Assignment operator.
CoagModel &CoagModel::operator=(const Sweep::CoagModel &rhs)
{
    // Nothing to do here.
    return *this;
}

// Returns the one and only instance of the CoagModel class.
// If the instance hasn't been created then this is done
// also.
CoagModel &CoagModel::Instance(void)
{
    static CoagModel inst;
    return inst;
}


// PERFORMING COAGULATIONS.

void CoagModel::Perform(Sweep::Particle &p1, Sweep::Particle &p2) const
{
    p1 += p2;
}


// PARTICLE UPDATES.

void CoagModel::UpdateParticle(Sweep::Particle &p, 
                               const Sweep::fvector &dcomp, 
                               const Sweep::fvector &dval) const
{
    // Do nothing.
}

void CoagModel::UpdateParticle(Sweep::Particle &p, 
                               const Sweep::fvector &dcomp, 
                               const Sweep::fvector &dval, 
                               unsigned int n) const
{
    // Do nothing special.
}


// PARTICLE-PARTICLE COAGULATION.

void CoagModel::CoagParticles(Particle &p1, const Particle &p2) const
{
    // Do nothing special.
}


// PROPERTY CALCULATION.

// Recalculates those model properties which are function of
// other particle properties.
void CoagModel::UpdateCache(Sweep::ParticleData &p) const
{
    CoagModelData &cache = *p.CoagModelCache();
    cache.m_dcolsqr      = p.CollDiameter() * p.CollDiameter();
    cache.m_inv_dcol     = 1.0 / p.CollDiameter();
    cache.m_inv_dcolsqr  = cache.m_inv_dcol * cache.m_inv_dcol;
    cache.m_inv_sqrtmass = 1.0 / sqrt(p.Mass());
    cache.m_d2_m_1_2     = cache.m_dcolsqr * cache.m_inv_sqrtmass;
}