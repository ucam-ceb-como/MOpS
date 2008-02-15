#include "swp_pripartdata.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
PriPartModelData::PriPartModelData(void)
{
}

// Default constructor (public).
PriPartModelData::PriPartModelData(ParticleData &parent)
: IModelData(parent)
{
}

// Copy constructor.
PriPartModelData::PriPartModelData(const PriPartModelData &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
PriPartModelData::PriPartModelData(std::istream &in, ParticleData &parent)
{
    Deserialize(in);
    SetParent(parent);
}

// Default destructor.
PriPartModelData::~PriPartModelData()
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
PriPartModelData &PriPartModelData::operator=(const Sweep::PriPartModelData &rhs)
{
    if (this != &rhs) {
        m_primaries.assign(rhs.m_primaries.begin(), rhs.m_primaries.end());
    }
    return *this;
}

// Compound assignment.
PriPartModelData &PriPartModelData::operator+=(const Sweep::PriPartModelData &rhs)
{
    return *this;
}

// Addition operator.
const PriPartModelData PriPartModelData::operator+(const Sweep::PriPartModelData &rhs) const
{
    return PriPartModelData(*this) += rhs;
}


// PROPERTIES.

// Returns the vector of primary particles.
std::vector<Primary> &PriPartModelData::Primaries(void) {return m_primaries;}
const std::vector<Primary> &PriPartModelData::Primaries(void) const {return m_primaries;}

// Returns the property with the given ID.
real PriPartModelData::Property(unsigned int id) const
{
    return 0.0;
}


// MODEL WHICH USES THIS DATA.

// Returns the PriPartModel object which operator on this data.
const PriPartModel &PriPartModelData::Model(void) const
{
    return PriPartModel::Instance();
}


// READ/WRITE/COPY.

// Returns a copy of the model data.
PriPartModelData *const PriPartModelData::Clone(void) const
{
    return new PriPartModelData(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
ModelType PriPartModelData::ID(void) const {return PriPartModel_ID;}

// Writes the object to a binary stream.
void PriPartModelData::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output primary count.
        unsigned int n = (unsigned int)m_primaries.size();
        out.write((char*)&n, sizeof(n));

        // Output primaries.
        for (unsigned int i=0; i!=n; ++i) {
            m_primaries[i].Serialize(out);
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PriPartModelData::Serialize).");
    }
}

// Reads the object from a binary stream.
void PriPartModelData::Deserialize(std::istream &in)
{
    m_primaries.clear();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;

        switch (version) {
            case 0:
                // Read the number of primaries.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the primaries.
                for (unsigned int i=0; i!=n; ++i) {
                    m_primaries.push_back(Primary(in));
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, PriPartModelData::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PriPartModelData::Deserialize).");
    }
}
