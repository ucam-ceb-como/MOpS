#include "swp_submodel_cache.h"

using namespace Sweep;
using namespace Sweep::SubModels;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
SubModels::SubModelCache::SubModelCache(void)
: m_parent(NULL)
{
}

// Default constructor (public).
SubModelCache::SubModelCache(Sweep::ParticleCache &parent)
{
    m_parent = &parent;
}

// Copy constructor.
SubModelCache::SubModelCache(const SubModelCache &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Default destructor.
SubModelCache::~SubModelCache()
{
    // Nothing special to destruct.
}

/*
// OPERATOR OVERLOADING.

// Assignment operator.
SubModelCache &SubModelCache::operator=(const SubModelCache &rhs)
{
    return *this;
}

// Compound assignment.
SubModelCache &SubModelCache::operator+=(const SubModelCache &rhs)
{
    return *this;
}
*/

// PARENT.

// Returns a pointer to the parent particle data.
ParticleCache *const SubModelCache::Parent(void) const
{
    return m_parent;
}

// Sets the parent particle data.
void SubModelCache::SetParent(ParticleCache &parent)
{
    m_parent = &parent;
}
