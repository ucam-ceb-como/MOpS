#include "mops_mechanism.h"
#include <stdexcept>

using namespace Mops;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Mechanism::Mechanism(void)
{
    m_pmech.SetSpecies(m_species);
}

// Default destructor.
Mechanism::~Mechanism(void)
{
}


// PARTICLE MECHANISM.

// Returns a reference (non-const) to the particle mechanism.
Sweep::Mechanism &Mechanism::ParticleMech(void) {return m_pmech;}
const Sweep::Mechanism &Mechanism::ParticleMech(void) const {return m_pmech;}


// READ/WRITE/COPY FUNCTIONS.

// Writes the mechanism to a binary data stream.
void Mechanism::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialize version to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));
        
        // Write the base class.
        Sprog::Mechanism::Serialize(out);

        // Write the particle mechanism.
        m_pmech.Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Mops, Mechanism::Serialize).");
    }
}

// Reads the mechanism data from a binary data stream.
void Mechanism::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the serialized mechanism version.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                // Read the base class.
                Sprog::Mechanism::Deserialize(in);

                // Read the particle mechanism.
                m_pmech.Deserialize(in);
                m_pmech.SetSpecies(m_species);

                break;
            default:
                throw runtime_error("Mechanism serialized version number "
                                    "is unsupported (Mops, Mechanism::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Mops, Mechanism::Deserialize).");
    }
}
