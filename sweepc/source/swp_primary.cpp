#include "swp_primary.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Primary::Primary(void)
: m_comp(0), m_vol(0.0), m_mass(0.0), m_diam(0.0), m_surf(0.0)
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
        m_comp = rhs.m_comp;
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
    m_comp += rhs.m_comp;
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
    return m_comp == rhs.m_comp;
}

// Inverse comparison.
bool Primary::operator!=(const Primary &rhs) const
{
    return m_comp != rhs.m_comp;
}

// Less than.
bool Primary::operator<(const Primary &rhs) const
{
    return m_comp < rhs.m_comp;
}

// Greater than
bool Primary::operator>(const Primary &rhs) const
{
    return m_comp > rhs.m_comp;
}

// Less than or equal
bool Primary::operator<=(const Primary &rhs) const
{
    return m_comp <= rhs.m_comp;
}

// Greater than or equal
bool Primary::operator>=(const Primary &rhs) const
{
    return m_comp >= rhs.m_comp;
}


// COMPOSITION.

// Returns the composition of the primary (only one component).
unsigned int Primary::Composition(void) const
{
    return m_comp;
}

// Sets the primary composition.
void Primary::SetComposition(unsigned int comp)
{
    m_comp = comp;
}

// Changes the composition by the given amount.
void Primary::ChangeComposition(int dc)
{
    if (dc < 0) {
        if (-dc > m_comp) {
            m_comp = 0;
        } else {
            m_comp += dc;
        }
    } else {
        m_comp += dc;
    }
}


// DERIVED PROPERTIES.

// Returns the primary diameter.
real Primary::Diameter(void) const {return m_diam;}

// Returns the primary volume.
real Primary::Volume(void) const {return m_vol;}

// Sets the primary volume.
void Primary::SetVolume(real vol) {m_vol = vol;}

// Returns the primary mass.
real Primary::Mass(void) const {return m_mass;}

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

        // Output the composition.
        unsigned int n = (unsigned int)m_comp;
        out.write((char*)&n, sizeof(n));

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

        unsigned int n = 0;
        double val = 0.0;

        switch (version) {
            case 0:
                // Read the composition.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_comp = (unsigned int)n;

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
