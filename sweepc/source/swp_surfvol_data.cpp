#include "swp_surfvoldata.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
SurfVolData::SurfVolData(void)
: m_sphsurf(0.0), m_surf(0.0)
{
}

// Default constructor (public).
SurfVolData::SurfVolData(Sweep::ParticleData &parent)
: CoagModelData(parent), m_sphsurf(0.0), m_surf(0.0)
{
}

// Copy constructor.
SurfVolData::SurfVolData(const Sweep::SurfVolData &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
SurfVolData::SurfVolData(std::istream &in, ParticleData &parent)
{
    Deserialize(in);
    SetParent(parent);
}

// Default destructor.
SurfVolData::~SurfVolData()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
SurfVolData &SurfVolData::operator=(const Sweep::SurfVolData &rhs)
{
    if (this != &rhs) {
        CoagModelData::operator=(rhs);
        m_sphsurf = rhs.m_sphsurf;
        m_surf    = rhs.m_surf;
    }
    return *this;
}

// Compound assignment.
SurfVolData &SurfVolData::operator +=(const Sweep::SurfVolData &rhs)
{
    CoagModelData::operator +=(rhs);
    m_sphsurf += rhs.m_sphsurf;
    m_surf    += rhs.m_surf;
    return *this;
}

// Addition operator.
const SurfVolData SurfVolData::operator +(const Sweep::SurfVolData &rhs) const
{
    return SurfVolData(*this) += rhs;
}

// Resets the model data to the default state.
void SurfVolData::Clear()
{
    m_sphsurf = 0.0;
    m_surf    = 0.0;
}

// COAGULATION MODEL PARTICLE PROPERTIES.

// Returns the equivalent spherical surface area.
real SurfVolData::SphSurfaceArea(void) const {return m_sphsurf;}

// Returns the actual surface area.
real SurfVolData::SurfaceArea(void) const {return m_surf;}


// MODEL WHICH USES THIS DATA.

// Returns the SurfVolModel which operates on this data.
const SurfVolModel &SurfVolData::Model(void) const
{
    return SurfVolModel::Instance();
}


// READ/WRITE/COPY.

// Returns a copy of the data.
SurfVolData *const SurfVolData::Clone(void) const
{
    return new SurfVolData(*this);
}

// Returns the model ID.
ModelType SurfVolData::ID(void) const
{
    return SVModel_ID;
}

// Writes the object to a binary stream.
void SurfVolData::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output base class.
        CoagModelData::Serialize(out);

        // Output the sphere surface area.
        double v = (double)m_sphsurf;
        out.write((char*)&v, sizeof(v));

        // Output the true surface area.
        v = (double)m_surf;
        out.write((char*)&v, sizeof(v));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, SurfVolData::Serialize).");
    }
}

// Reads the object from a binary stream.
void SurfVolData::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;

        switch (version) {
            case 0:
                // Read base class.
                CoagModelData::Deserialize(in);

                // Read the sphere surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_sphsurf = (real)val;

                // Read the true surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_surf = (real)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, SurfVolData::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SurfVolData::Deserialize).");
    }
}
