#include "swp_pripartdata.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
PriPartData::PriPartData(void)
{
}

// Default constructor (public).
PriPartData::PriPartData(ParticleData &parent)
: IModelData(parent)
{
}

// Copy constructor.
PriPartData::PriPartData(const PriPartData &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
PriPartData::PriPartData(std::istream &in, ParticleData &parent)
{
    Deserialize(in);
    SetParent(parent);
}

// Default destructor.
PriPartData::~PriPartData()
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
PriPartData &PriPartData::operator=(const Sweep::PriPartData &rhs)
{
    if (this != &rhs) {
        m_primaries.assign(rhs.m_primaries.begin(), rhs.m_primaries.end());
    }
    return *this;
}

// Compound assignment.
PriPartData &PriPartData::operator+=(const Sweep::PriPartData &rhs)
{
    return *this;
}

// Addition operator.
const PriPartData PriPartData::operator+(const Sweep::PriPartData &rhs) const
{
    return PriPartData(*this) += rhs;
}

// Resets the model data to the default state.
void PriPartData::Clear()
{
    m_primaries.clear();
}


// PROPERTIES.

// Returns the vector of primary particles.
std::vector<Primary> &PriPartData::Primaries(void) {return m_primaries;}
const std::vector<Primary> &PriPartData::Primaries(void) const {return m_primaries;}

// Returns the property with the given ID.
real PriPartData::Property(unsigned int id) const
{
    return 0.0;
}


// MODEL WHICH USES THIS DATA.

// Returns the PriPartModel object which operator on this data.
const PriPartModel &PriPartData::Model(void) const
{
    return PriPartModel::Instance();
}


// READ/WRITE/COPY.

// Returns a copy of the model data.
PriPartData *const PriPartData::Clone(void) const
{
    return new PriPartData(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
ModelType PriPartData::ID(void) const {return PriPartModel_ID;}

// Writes the object to a binary stream.
void PriPartData::Serialize(std::ostream &out) const
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
void PriPartData::Deserialize(std::istream &in)
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
