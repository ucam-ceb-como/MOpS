#include "swp_aggmodel_cache.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
AggModelCache::AggModelCache(void)
{
    m_parent = NULL;
}

// Default constructor (public).
AggModelCache::AggModelCache(ParticleCache &parent)
{
    m_parent = &parent;
}

// Copy constructor.
AggModelCache::AggModelCache(const AggModelCache &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
AggModelCache::AggModelCache(std::istream &in, ParticleCache &parent)
{
    m_parent = &parent;
}

// Default destructor.
AggModelCache::~AggModelCache()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator (AggModelCache RHS).
AggModelCache &AggModelCache::operator=(const AggModelCache &rhs)
{
    return *this;
}

/*
// Assignment operator (Primary RHS).
AggModelCache &AggModelCache::operator=(const Primary &rhs)
{
    return (*this=rhs.TypedRef());
}
*/

// Compound assignment operator (AggModelCache RHS).
AggModelCache &AggModelCache::operator+=(const AggModelCache &rhs)
{
    return *this;
}

/*
// Compound assignment operator (Primary RHS).
AggModelCache &AggModelCache::operator+=(const Primary &rhs)
{
    return (*this+=rhs.TypedRef());
}
*/


// PARENT PARTICLE-CACHE.

// Returns a pointer to the parent particle data.
ParticleCache *const AggModelCache::Parent(void) const
{
    return m_parent;
}

// Sets the parent particle data.
void AggModelCache::SetParent(ParticleCache &parent)
{
    m_parent = &parent;
}
