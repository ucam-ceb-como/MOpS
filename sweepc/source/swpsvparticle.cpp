#include "swpsvparticle.h"

using namespace Sweep;

SVParticle::SVParticle(void)
{
}

SVParticle::~SVParticle(void)
{
}

void SVParticle::Adjust(const vector<real> &dcomp, const vector<real> &dvalues)
{
    // Save the volume.
    real vold = m_cache[iV];

    // Adjust particle as if it were spherical.
    DefaultParticle::Adjust(dcomp, dvalues);

    if (IsValid()) {
        // Calculate the change in volume and the surface area
        // if this particle were a sphere.
        real dvol = m_cache[iV] - vold;
        real ssph = SphereSurface();

        // Choose radius of gyration based on the direction of component
        // change.
        vector<real>::const_iterator i;
        real dtot=0.0;
        for (i=dcomp.begin(); i!=dcomp.end(); i++) {
            dtot += *i;
        }
        if (dtot > 0.0) {
            m_cache[iS] += 2.0 * dvol / sqrt(m_cache[iS] / (4*PI));
        } else {
            m_cache[iS] += 6.0 * dvol * m_cache[iV] / m_cache[iS];
        }

        // Surface area cannot be less than that of a sphere.
        if (m_cache[iS]<ssph) m_cache[iS] = ssph;
    }
};

void SVParticle::Adjust(const vector<real> &dcomp, const vector<real> &dvalues, const unsigned int n)
{
    // Save the volume.
    real vold = m_cache[iV];

    // Adjust particle as if it were spherical.
    DefaultParticle::Adjust(dcomp, dvalues, n);

    if (IsValid()) {
        // Calculate the change in volume and the surface area
        // if this particle were a sphere.
        real dvol = m_cache[iV] - vold;
        real ssph = SphereSurface();

        // Choose radius of gyration based on the direction of component
        // change.
        vector<real>::const_iterator i;
        real dtot=0.0;
        for (i=dcomp.begin(); i!=dcomp.end(); i++) {
            dtot += *i;
        }
        if (dtot > 0.0) {
            m_cache[iS] += 2.0 * dvol / sqrt(m_cache[iS] / (4*PI));
        } else {
            m_cache[iS] += 6.0 * dvol * m_cache[iV] / m_cache[iS];
        }

        // Surface area cannot be less than that of a sphere.
        if (m_cache[iS]<ssph) m_cache[iS] = ssph;
    }
};

/* Particle actions and interactions. */

DefaultParticle &SVParticle::Coagulate(const DefaultParticle &sp)
{
    // Conserve surface area of particles.
    m_cache[iS] += sp.SurfaceArea();

    // Most of the coagulation process is the same as for
    // spherical particles. So call that routine.
    DefaultParticle::Coagulate(sp);

    // DefaultParticle::Coagulate recalculates the cache, so
    // just return this particle.
    return *this;
}

/* Protected functions. */

void SVParticle::CalcCache(void)
{
    real m;
    vector<Component*>::iterator ic;
    vector<real>::iterator j = m_comp.begin();

    // Recalculate particle mass and volume from components.
    m_cache[iM] = m_cache[iV] = 0.0;
    for (ic=m_components->begin(); ic!=m_components->end(); ic++, j++) {
        m = (*ic)->MolWt() * (*j) / NA;
        m_cache[iM] += m;
        m_cache[iV] += m / (*ic)->Density();
    }

    // Surface area is now tracked, but if it is found to be zero then
    // set it the that of a spherical particle.  In this was the particle
    // is correctly initialised on creation.
    if (m_cache[iS]==0.0) m_cache[iS] = PI * pow(6.0 * m_cache[iV] / PI, TWO_THIRDS);

    // Collision diameter is calculated differently from the
    // spherical particle case.  It is now a mean value calculated using
    // the volume and the surface area.
    m_cache[iD] = (pow(6.0 * m_cache[iV] / PI, ONE_THIRD) + (m_cache[iS]/PI)) * 0.5;

    // All these properties should remain the same as for the spherical
    // particle case.
    m_cache[iD2]      = m_cache[iD] * m_cache[iD];
    m_cache[iD_1]     = 1.0 / m_cache[iD];
    m_cache[iD_2]     = 1.0 / m_cache[iD2];
    m_cache[iM_1_2]   = pow(m_cache[iM], -0.5);
    m_cache[iD2M_1_2] = m_cache[iD2] * m_cache[iM_1_2];
}