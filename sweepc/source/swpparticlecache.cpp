#include "swpparticlecache.h"
#include <vector>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Constructor.
ParticleCache::ParticleCache(void)
{
	m_fixedcache.assign(NCACHE,0.0);
}

// Copy-constructor.
ParticleCache::ParticleCache(const ParticleCache & ch)
{
	*this = ch;
}

// Destructor.
ParticleCache::~ParticleCache(void)
{
	m_fixedcache.clear();
}


// OPERATOR OVERLOADING.

// Assignment operator.
ParticleCache & ParticleCache::operator=(const ParticleCache &ch)
{
    if (this != &ss) {
		m.fixedcache = ch.m_fixedcache;
    }
    return *this;
}

// Comparison operator:  Returns true if both ParticleCache objects are identical.
bool ParticleCache::operator==(const ParticleCache &ch) const
{
    return (m_fixedcache == ch.m_fixedcache);
}

// Inequality operator:  Returns false if both ParticleCache objects are identical.
bool ParticleCache::operator!=(const ParticleCache &ch) const
{
    return !(*this==ch);
}

// Overloading Subscript operator (for non-constant objects).
real & operator[] (const int i)
{
	return m_fixedcache[i];
}

// Overloading Subscript operator (for constant objects).
const real & operator[] (const int i) const
{
	return m_fixedcache[i];
}

// Overloading += operator.
ParticleCache & ParticleCache::operator+=(const ParticleCache &ch)
{
	m_fixedcache.resize(max(m_fixedcache.size(),ch.m_fixedcache.size()), 0.0);
    vector<real>::iterator i;
    vector<real>::const_iterator j=ch.m_fixedcache.begin();
    for (i=m_fixedcache.begin(); i!=m_fixedcache.end(); i++,j++) *i = *i + *j;
    return *this;
}


