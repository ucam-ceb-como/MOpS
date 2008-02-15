#include "swp_component.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSRUCTORS AND DESTRUCTORS.

// Default constructor.
Component::Component()
: m_density(0.0), m_molwt(0.0), m_name("")
{
}

// Initialising constructor.
Component::Component(Sweep::real molwt, 
                     Sweep::real dens, 
                     const std::string &name)
{
    // Initialise the component properties.
    m_density = dens;
    m_molwt   = molwt;
    m_name    = name;
}

// Copy constructor.
Component::Component(const Component &copy) 
{
    *this = copy;
}

// Stream-reading constructor.
Component::Component(std::istream &in) 
{
    Deserialize(in);
}

// Default destructor.
Component::~Component()
{
}


// OPERATOR OVERLOADS.
Component &Component::operator=(const Component &rhs)
{
    if (this != &rhs) {
        m_density = rhs.m_density;
        m_molwt   = rhs.m_molwt;
        m_name    = rhs.m_name;
    }
    return *this;
}


// READ/WRITE/COPY.

// Creates a copy of the component.
Component *const Component::Clone(void) const {return new Component(*this);}

// Writes the object to a binary stream.
void Component::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write molecular weight.
        double v = (double)m_molwt;
        out.write((char*)&v, sizeof(v));

        // Write density
        v = (double)m_density;
        out.write((char*)&v, sizeof(v));

        // Write the length of the component name to the stream.
        unsigned int n = m_name.length();
        out.write((char*)&n, sizeof(n));

        // Write the component name to the stream.
        out.write(m_name.c_str(), n);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Component::Serialize).");
    }
}

// Reads the object from a binary stream.
void Component::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;
        unsigned int n = 0;
        char *name = NULL;

        switch (version) {
            case 0:
                // Read molecular weight.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_molwt = (real)val;

                // Read density
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_density = (real)val;

                // Read the length of the species name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                
                // Read the species name.
                name = new char[n];
                in.read(name, n);
                m_name.assign(name, n);
                delete [] name;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Component::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Component::Deserialize).");
    }
}
