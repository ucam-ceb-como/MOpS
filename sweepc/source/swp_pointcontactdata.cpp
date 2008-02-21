#include "swp_pointcontactdata.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
PointContactData::PointContactData(void)
: m_sphsurf(0.0), m_surf(0.0)
{
}

// Default constructor (public).
PointContactData::PointContactData(Sweep::ParticleData &parent)
: CoagModelData(parent)
{
}

// Copy constructor.
PointContactData::PointContactData(const Sweep::PointContactData &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
PointContactData::PointContactData(std::istream &in, ParticleData &parent)
{
    Deserialize(in);
    SetParent(parent);
}

// Default destructor.
PointContactData::~PointContactData()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
PointContactData &PointContactData::operator=(const Sweep::PointContactData &rhs)
{
    if (this != &rhs) {
        CoagModelData::operator=(rhs);
        m_sphsurf = rhs.m_sphsurf;
        m_surf    = rhs.m_surf;
    }
    return *this;
}

// Compound assignment.
PointContactData &PointContactData::operator +=(const Sweep::PointContactData &rhs)
{
    CoagModelData::operator +=(rhs);
    m_sphsurf += rhs.m_sphsurf;
    m_surf    += rhs.m_surf;
    return *this;
}

// Addition operator.
const PointContactData PointContactData::operator +(const Sweep::PointContactData &rhs) const
{
    return PointContactData(*this) += rhs;
}

// Resets the model data to the default state.
void PointContactData::Clear()
{
    m_sphsurf = 0.0;
    m_surf    = 0.0;
}

// COAGULATION MODEL PARTICLE PROPERTIES.

// Returns the equivalent spherical surface area.
real PointContactData::SphSurfaceArea(void) const {return m_sphsurf;}

// Returns the actual surface area.
real PointContactData::SurfaceArea(void) const {return m_surf;}


// MODEL WHICH USES THIS DATA.

// Returns the PointContactModel which operates on this data.
const PointContactModel &PointContactData::Model(void) const
{
    return PointContactModel::Instance();
}


// READ/WRITE/COPY.

// Returns a copy of the data.
PointContactData *const PointContactData::Clone(void) const
{
    return new PointContactData(*this);
}

// Returns the model ID.
ModelType PointContactData::ID(void) const
{
    return SVModel_ID;
}

// Writes the object to a binary stream.
void PointContactData::Serialize(std::ostream &out) const
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
                               "(Sweep, PointContactData::Serialize).");
    }
}

// Reads the object from a binary stream.
void PointContactData::Deserialize(std::istream &in)
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
                                    "(Sweep, PointContactData::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PointContactData::Deserialize).");
    }
}
