#include "swp_surfvolmodel.h"
#include "swp_surfvoldata.h"
#include "swp_particle.h"
#include <string>

using namespace Sweep;
using namespace std;

// SINGLETON IMPLEMENTATION.

// Default constructor.
SurfVolModel::SurfVolModel()
{
    // Do nothing.
}

// Copy constructor.
SurfVolModel::SurfVolModel(const Sweep::SurfVolModel &copy)
{
    // Do nothing.
}

// Assignment operator.
SurfVolModel &SurfVolModel::operator =(const Sweep::SurfVolModel &rhs)
{
    return *this;
}

// Default destructor.
SurfVolModel::~SurfVolModel()
{
    // Do nothing.
}

// Returns the one and only instance.
SurfVolModel &SurfVolModel::Instance()
{
    static SurfVolModel inst;
    return inst;
}


// PERFORMING COAGULATIONS.

// Performs a coagulation event on the two given particles.
// The resultant particle is stored in the location of p1.
void SurfVolModel::Perform(Particle &p1, Particle &p2) const
{
    p1 += p2;
}


// PARTICLE PROPERTY CALCULATION.

void SurfVolModel::UpdateParticle(Particle &p, 
                                  const fvector &dcomp,
                                  const fvector &dval) const
{
    SurfVolData *data = dynamic_cast<SurfVolData*>(p.CoagModelCache());

    // Calculate change in volume.
    real dvol = 0.0;
    for (unsigned int i=0; i!=dcomp.size(); ++i) {
        dvol += dcomp[i] * (*p.Components())[i]->MolWt() / 
                (*p.Components())[i]->Density();
    }
    dvol /= NA;

    // Calculate change in surface area.
    real rad = 0.0;
    if (dvol > 0.0) {
        // Inverse growth radius.
        rad = sqrt(4.0 * PI / data->m_surf);
    } else {
        // Inverse oxidation radius.    
        rad = data->m_surf / (3.0 * p.Volume());
    }

    // Save new surface area.
    data->m_surf += 2.0 * dvol * rad;
}

void SurfVolModel::UpdateParticle(Particle &p, const fvector &dcomp, 
                                  const fvector &dval, unsigned int n) const
{
    SurfVolData *data = dynamic_cast<SurfVolData*>(p.CoagModelCache());

    // Calculate change in volume.
    real dvol = 0.0;
    for (unsigned int i=0; i!=dcomp.size(); ++i) {
        dvol += dcomp[i] * (*p.Components())[i]->MolWt() / 
                (*p.Components())[i]->Density();
    }
    dvol *= (real)n / NA;

    // Calculate change in surface area.
    real rad = 0.0;
    if (dvol > 0.0) {
        // Inverse growth radius.
        rad = sqrt(4.0 * PI / data->m_surf);
    } else {
        // Inverse oxidation radius.    
        rad = data->m_surf / (3.0 * p.Volume());
    }

    // Save new surface area.
    data->m_surf += 2.0 * dvol * rad;
}

// PARTICLE-PARTICLE COAGULATION.

// Updates the first particle according to the model rules for
// coagulation.  The first particle is the recipient of the second.
void SurfVolModel::CoagParticles(Particle &p1, const Particle &p2) const
{
    SurfVolData *data = dynamic_cast<SurfVolData*>(p1.CoagModelCache());
    data->m_surf = p1.SurfaceArea() + p2.SurfaceArea();
}


// PROPERTY CACHE UPDATE.

void SurfVolModel::UpdateCache(Sweep::ParticleData &p) const
{
    // Calculate base class properties.
    CoagModel::UpdateCache(p);

    SurfVolData &cache = dynamic_cast<SurfVolData&>(*p.CoagModelCache());

    // Set correct surface area.
    cache.m_sphsurf = p.SurfaceArea();
    cache.m_surf    = max(cache.m_surf, cache.m_sphsurf);
    p.SetSurfaceArea(cache.m_surf);

    // Set correct collision diameter.
    p.SetCollDiameter((pow(6.0 * p.Volume() / PI, ONE_THIRD) + 
                       sqrt(p.SurfaceArea() / PI)) * 0.5);
}
