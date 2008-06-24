#include "swp_surfvol_cache.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_particle_cache.h"
#include "swp_surfvol_primary.h"
#include <stdexcept>
#include <iostream>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
SurfVolCache::SurfVolCache(void)
: m_sphsurf(0.0), m_surf(0.0), m_ppn(0), m_ppd(0.0)
{
}

// Default constructor (public).
SurfVolCache::SurfVolCache(ParticleCache &parent)
: AggModelCache(parent), m_sphsurf(0.0), m_surf(0.0), m_ppn(0), m_ppd(0.0)
{
}

// Copy constructor.
SurfVolCache::SurfVolCache(const SurfVolCache &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
SurfVolCache::SurfVolCache(std::istream &in, ParticleCache &parent)
{
    Deserialize(in, parent);
}

// Default destructor.
SurfVolCache::~SurfVolCache()
{
    // Nothing special to destruct.
}


// ASSIGNMENT OPERATOR OVERLOADING.

// Assignment operator (SurfVolCache RHS).
SurfVolCache &SurfVolCache::operator=(const SurfVolCache &rhs)
{
    if (this != &rhs) {
        m_sphsurf = rhs.m_sphsurf;
        m_surf    = rhs.m_surf;
        m_ppn     = rhs.m_ppn;
        m_ppd     = rhs.m_ppd;
    }
    return *this;
}

// Assignment operator (SurfVolPrimary RHS).
SurfVolCache &SurfVolCache::operator=(const SurfVolPrimary &rhs)
{
    m_sphsurf = rhs.SphSurfaceArea();
    m_surf    = rhs.SurfaceArea();
    // These are calculated assuming point-contact of identical spheres.
    m_ppn     = (unsigned int)((m_surf * m_surf * m_surf) / 
                               (36.0 * PI * rhs.Volume() * rhs.Volume()));
    m_ppd     = 6.0 * rhs.Volume() / m_surf;
    return *this;
}

// Assignment operator (AggModelCache RHS).
SurfVolCache &SurfVolCache::operator=(const AggModelCache &rhs)
{
    // Attempt to cast the RHS as a SurfVolCache.  This will throw
    // an exception if it isn't possible to cast.
    return operator=(dynamic_cast<const SurfVolCache&>(rhs));
}

// Assignment operator (Primary RHS).
SurfVolCache &SurfVolCache::operator=(const Primary &rhs)
{
    if (rhs.AggID() == Spherical_ID) {
        // If the RHS is actually a spherical primary then copy
        // the required data.
        m_sphsurf = rhs.SurfaceArea();
        m_surf    = rhs.SurfaceArea();
        m_ppn     = 1;
        m_ppd     = rhs.SphDiameter();
    } else {
        // Attempt to cast the RHS as a SurfVolPrimary.  This will throw
        // an exception if it isn't possible to cast.
        operator=(dynamic_cast<const SurfVolPrimary&>(rhs));
    }
    return *this;
}


// COMPOUND ASSIGNMENT OPERATOR OVERLOADING.

// Compound assignment (SurfVolCache RHS).
SurfVolCache &SurfVolCache::operator+=(const SurfVolCache &rhs)
{
    m_sphsurf += rhs.m_sphsurf;
    m_surf    += rhs.m_surf;
    m_ppn     += rhs.m_ppn;
    m_ppd     += rhs.m_ppd;
    return *this;
}

// Compound assignment (SurfVolPrimary RHS).
SurfVolCache &SurfVolCache::operator+=(const SurfVolPrimary &rhs)
{
    m_sphsurf += rhs.SphSurfaceArea();
    m_surf    += rhs.SurfaceArea();
    // These are calculated assuming point-contact of identical spheres.
    m_ppn     += (unsigned int)((m_surf * m_surf * m_surf) / 
                                (36.0 * PI * rhs.Volume() * rhs.Volume()));
    m_ppd     += 6.0 * rhs.Volume() / m_surf;
    return *this;
}

// Compound assignment (AggModelCache RHS).
SurfVolCache &SurfVolCache::operator+=(const AggModelCache &rhs)
{
    // Attempt to cast the RHS as a SurfVolCache.  This will throw
    // an exception if it isn't possible to cast.
    return operator+=(dynamic_cast<const SurfVolCache&>(rhs));
}

// Compound assignment (Primary RHS).
SurfVolCache &SurfVolCache::operator+=(const Primary &rhs)
{
    if (rhs.AggID() == Spherical_ID) {
        m_sphsurf += rhs.SurfaceArea();
        m_surf    += rhs.SurfaceArea();
        m_ppn     += 1;
        m_ppd     += rhs.SphDiameter();
    } else {
        // Attempt to cast the RHS as a SurfVolPrimary.  This will throw
        // an exception if it isn't possible to cast.
        operator+=(dynamic_cast<const SurfVolPrimary&>(rhs));
    }
    return *this;
}


// DATA MANAGEMENT.

// Resets the model data to the default state.
void SurfVolCache::Clear()
{
    m_sphsurf = 0.0;
    m_surf    = 0.0;
    m_ppn     = 0;
    m_ppd     = 0.0;
}


// AGGREGATION MODEL PROPERTIES.

// Returns the equivalent spherical surface area.
real SurfVolCache::SphSurfaceArea(void) const {return m_sphsurf;}

// Sets the equivalent spherical surface area.
void SurfVolCache::SetSphSurfaceArea(real surface) {m_sphsurf = surface;}

// Returns the actual surface area.
real SurfVolCache::SurfaceArea(void) const {return m_surf;}

// Sets the actual surface area.
void SurfVolCache::SetSurfaceArea(real surface) {m_surf = surface;}

// Returns the primary particle count.
unsigned int SurfVolCache::PP_Count(void) const {return m_ppn;}

// Sets the primary particle count.
void SurfVolCache::SetPP_Count(unsigned int n) {m_ppn = n;}

// Returns the average primary particle diameter.
real SurfVolCache::PP_Diameter(void) const {return m_ppd;};

// Sets the average primary particle diameter.
void SurfVolCache::SetPP_Diameter(real d) {m_ppd = d;}


// READ/WRITE/COPY.

// Returns a copy of the data.
SurfVolCache *const SurfVolCache::Clone(void) const {return new SurfVolCache(*this);}

// Returns the model ID.
AggModelType SurfVolCache::ID(void) const {return SurfVol_ID;}

// Writes the object to a binary stream.
void SurfVolCache::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output the sphere surface area.
        double v = (double)m_sphsurf;
        out.write((char*)&v, sizeof(v));

        // Output the true surface area.
        v = (double)m_surf;
        out.write((char*)&v, sizeof(v));

        // Output the primary-particle count.
        unsigned int n = m_ppn;
        out.write((char*)&n, sizeof(n));

        // Output the primary-particle diameter..
        v = (double)m_ppd;
        out.write((char*)&v, sizeof(v));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, SurfVolCache::Serialize).");
    }
}

// Reads the object from a binary stream.
void SurfVolCache::Deserialize(std::istream &in, ParticleCache &parent)
{
    m_parent = &parent;

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double val     = 0.0;

        switch (version) {
            case 0:
                // Read the sphere surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_sphsurf = (real)val;

                // Read the true surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_surf = (real)val;

                // Read the primary-particle count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_ppn = n;

                // Read the primary-particle diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_ppd = (real)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, SurfVolCache::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SurfVolCache::Deserialize).");
    }
}
