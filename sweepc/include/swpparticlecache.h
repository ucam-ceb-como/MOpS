/*
  Author(s):      Peter Man (plwm2)
  Project:        sweep (population balance solver)

  File purpose:
    Class which holds a vector of cache variables common to all particle types.
	Also holds index of each variable corresponding to the element in the vector.
*/

#ifndef SWEEP_PARTICLECACHE_H
#define SWEEP_PARTICLECACHE_H

#include <vector>

namespace Sweep
{
class ParticleCache
{
public:
	static const unsigned int NCACHE = 11; // Number of properties returned by particle.

	/* Indices in property vectors. */
	static const int iAS      = 0;
	static const int iD       = 1;
	static const int iD2      = 2;
	static const int iD_1     = 3;
	static const int iD_2     = 4;
	static const int iM_1_2   = 5;
	static const int iD2M_1_2 = 6;
	static const int iV       = 7;
	static const int iM       = 8;
	static const int iS       = 9;
	static const int iMOM1    = 10;

	// Constructors.
	ParticleCache(void);                      // Default Constructor.
	ParticleCache(const ParticleCache & ch);  // Copy-constructor.

	// Destructor.
	~ParticleCache(void);


	// OVERLOADING OPERATORS.
	ParticleCache & operator=(const ParticleCache &ch);   // Overloading Assignment operator.
	real & operator[] (const int i);               // Overloading Subscript operator (for non-constant objects).
	const real & operator[] (const int i) const;   // Overloading Subscript operator (for constant objects).
	bool operator==(const ParticleCache &ch) const;       // Overloading Comparison operator.
	bool operator!=(const ParticleCache &ch) const;       // Overloading Inequality operator.
	ParticleCache & operator+=(const ParticleCache &ch);  // Overloading += operator.

	// Vector related functions for m_fixedcache.
	inline real GetElement(const int index) const;   // Get m_fixedcache element corresponding to the given 'index'
	inline std::vector<real> & GetCache(void) const; // Get reference to m_fixedcache vector.
	inline void Clear(void);                         // Clear the m_fixedcache vector.
	inline void Resize(const unsigned int size);     // Resize the m_fixedcache vector.
	inline unsigned int Size(void) const;            // Returns size of m_fixedcache vector.


protected:
	// Vector of Cache variables common to all particle types.
	std::vector<real> m_fixedcache (NCACHE);
};

// Get m_fixedcache element corresponding to the given 'index'
inline real ParticleCache::GetElement(const int index) const { return m_fixedcache[index]; }

// Get reference to m_fixedcache vector.
inline std::vector<real> & GetCache(void) const { return m_fixedcache; }

// Clear the m_fixedcache vector.
inline void Clear(void) { m_fixedcache.assign(m_fixedcache.size(),0.0); }

// Resize the m_fixedcache vector.
inline void Resize(const unsigned int size) { m_fixedcache.resize(size,0.0); }

// Returns size of m_fixedcache vector.
inline unsigned int Size(void) const { return m_fixedcache.size(); }

};

#endif