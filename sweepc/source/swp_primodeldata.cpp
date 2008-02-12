#include "swp_pripartdata.h"

using namespace Sweep;

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
std::vector<Primary> &PriPartModelData::Primaries(void)
{
    return m_primaries;
}

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
