#include "swp_modeldata.h"

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
IModelData::IModelData()
: m_parent(NULL)
{
}

// Default constructor (public).
IModelData::IModelData(Sweep::ParticleData &parent)
{
    m_parent = &parent;
}

// Copy constructor.
IModelData::IModelData(const Sweep::IModelData &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Default destructor.
IModelData::~IModelData()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
IModelData &IModelData::operator =(const Sweep::IModelData &rhs)
{
    m_parent = rhs.m_parent;
    return *this;
}


// PARENT.

// Returns a pointer to the parent particle data.
ParticleData *const IModelData::Parent(void) const
{
    return m_parent;
}

// Sets the parent particle data.
void IModelData::SetParent(ParticleData &parent)
{
    m_parent = &parent;
}
