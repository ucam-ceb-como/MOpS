#include "swp_primary.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Primary::Primary(void)
: m_vol(0.0), m_mass(0.0), m_surf(0.0), m_diam(0.0)
{
}

// Copy constructor.
Primary::Primary(const Primary &copy) 
{
    *this = copy;
}

// Stream-reading constructor.
Primary::Primary(std::istream &in)
{
    Deserialize(in);
}

// Default destructor.
Primary::~Primary()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
Primary &Primary::operator=(const Primary &rhs)
{
    if (this != &rhs) {
        m_vol  = rhs.m_vol;
        m_mass = rhs.m_mass;
        m_diam = rhs.m_diam;
        m_surf = rhs.m_surf;
    }
    return *this;
}

// Compound assignment.
Primary &Primary::operator+=(const Primary &rhs)
{
    m_vol  += rhs.m_vol;
    m_mass += rhs.m_mass;
    m_diam += rhs.m_diam;
    m_surf += rhs.m_surf;
    return *this;
}

// Addition.
const Primary Primary::operator+(const Primary &rhs) const
{
    return Primary(*this) += rhs;
}

// Comparision.
bool Primary::operator==(const Primary &rhs) const
{
    return m_mass == rhs.m_mass;
}

// Inverse comparison.
bool Primary::operator!=(const Primary &rhs) const
{
    return m_mass != rhs.m_mass;
}

// Less than.
bool Primary::operator<(const Primary &rhs) const
{
    return m_mass < rhs.m_mass;
}

// Greater than
bool Primary::operator>(const Primary &rhs) const
{
    return m_mass > rhs.m_mass;
}

// Less than or equal
bool Primary::operator<=(const Primary &rhs) const
{
    return m_mass <= rhs.m_mass;
}

// Greater than or equal
bool Primary::operator>=(const Primary &rhs) const
{
    return m_mass >= rhs.m_mass;
}

// PROPERTIES.

// Returns the primary diameter.
real Primary::Diameter(void) const {return m_diam;}

// Returns the primary volume.
real Primary::Volume(void) const {return m_vol;}

// Sets the primary volume.
void Primary::SetVolume(real vol) {m_vol = vol;}

// Returns the primary mass.
real Primary::Mass(void) const {return m_mass;}

// Sets the primary mass.
void Primary::SetMass(real mass)
{
    m_mass = mass;
    // TODO: Set primary derived properties here.
}

// Adds/removes some mass to/from the primary.
void Primary::ChangeMass(real dm)
{
    m_mass += dm;
    // TODO: Set primary derived properties here.
}

// Returns the primary surface area.
real Primary::SurfaceArea(void) const {return m_surf;}


// READ/WRITE/COPY.

// Returns a copy of the model data.
Primary *const Primary::Clone(void) const
{
    return new Primary(*this);
}

// Writes the object to a binary stream.
void Primary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output the volume
        double v = (double)m_vol;
        out.write((char*)&v, sizeof(v));

        // Output the mass
        v = (double)m_mass;
        out.write((char*)&v, sizeof(v));

        // Output the diameter
        v = (double)m_diam;
        out.write((char*)&v, sizeof(v));

        // Output the surface area
        v = (double)m_surf;
        out.write((char*)&v, sizeof(v));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Primary::Serialize).");
    }
}

// Reads the object from a binary stream.
void Primary::Deserialize(std::istream &in)
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
                // Read the volume.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_vol = (real)val;

                // Read the mass
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_mass = (real)val;

                // Read the diameter
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_diam = (real)val;

                // Read the surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_surf = (real)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Primary::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Primary::Deserialize).");
    }
}
