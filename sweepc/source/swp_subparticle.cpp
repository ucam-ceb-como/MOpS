#include "swp_subparticle.h"
#include "swp_particle_model.h"
#include "swp_submodel.h"
#include "swp_model_factory.h"
#include "swp_cell.h"
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>

/*
 * This SubParticle class is now obsolete.  It was originally used to store the
 * connectivity structure of primary particle aggregates, but this
 * connectivity is now handled in the primary particle classes.
 *
 * This class will eventually removed so that Particle inherits
 * directly from ParticleCache, without having SubParticle in the
 * inheritance hierachy any more.
 *
 * DO NOT add methods or data to this class, since the whole
 * class is due to be removed. (riap2 28 jan 2011)
 *
 */

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
SubParticle::SubParticle(void)
: m_primary(NULL), m_pmodel(NULL),
  m_createt(0.0), m_time(0.0), m_aggcache(NULL)
{
}

// Initialising constructor.
SubParticle::SubParticle(real t, const Sweep::ParticleModel &model)
: m_primary(NULL), m_pmodel(&model),
  m_createt(0.0), m_time(0.0), m_comp(model.ComponentCount(), 0.0),
  m_values(model.TrackerCount(), 0.0)
{
	m_aggcache = ModelFactory::CreateAggCache(model.AggModel(), *this);

}

// Initialising constructor (from Primary particle).
SubParticle::SubParticle(Sweep::Primary &pri)
: m_pmodel(pri.ParticleModel()),
  m_comp(pri.Composition().begin(), pri.Composition().end()),
  m_values(pri.Values().begin(), pri.Values().end())

{
    m_createt = pri.CreateTime();
    m_time = pri.CreateTime();
    m_primary    = &pri;
	m_freesurface=pri.SurfaceArea();
	m_aggcache = pri.CreateAggCache(*this);
}

// Copy constructor.
SubParticle::SubParticle(const SubParticle &copy)
: m_aggcache(NULL)
{
    // Use assignment operator.
    m_primary    = NULL;

	*this        = copy;

}

// Stream-reading constructor.
SubParticle::SubParticle(std::istream &in, const Sweep::ParticleModel &model)
: m_primary(NULL),
  m_pmodel(&model), m_aggcache(NULL)
 {
    init();
    Deserialize(in, model);
}

// Default destructor.
SubParticle::~SubParticle()
{
    releaseMem();
}

// OPERATOR OVERLOADING.

// Assignment operator.
SubParticle &SubParticle::operator=(const SubParticle &rhs)
{
    if (this != &rhs) {
        m_pmodel = rhs.m_pmodel;
        ParticleCache::operator=(rhs);

        // Copy primary.
        if (rhs.m_primary != NULL) {
            if (m_primary == NULL) {
                m_primary = rhs.m_primary->Clone();
            } else {
                *m_primary = *rhs.m_primary;
            }
        } else {
            delete m_primary; m_primary = NULL;
        }

        // Copy composition and tracker data
        m_comp = rhs.m_comp;
        m_values = rhs.m_values;
        m_createt = rhs.m_createt;
        m_time = rhs.m_time;

        if (rhs.m_aggcache != NULL) {
            if ((m_aggcache==NULL) || (m_aggcache->ID() != rhs.m_aggcache->ID())) {
                delete m_aggcache;
                m_aggcache = rhs.m_aggcache->Clone();
                m_aggcache->SetParent(*this);
            } else {
                *m_aggcache = *rhs.m_aggcache;
            }
        } else {
            delete m_aggcache;
            m_aggcache = NULL;
        }
    }


	return *this;
}

real SubParticle::avgeomdiam(double oneovernumsubpart)
{
	Primary::PropID primid=Primary::iD;
	return std::pow(m_primary->Property(primid),oneovernumsubpart);
}

// PRIMARY PARTICLE CHILD.

// Returns a pointer to the child primary particle, NULL if
// this sub-particle has no primary.
Sweep::Primary *const SubParticle::Primary(void) {return m_primary;}
const Sweep::Primary *const SubParticle::Primary(void) const {return m_primary;}

// Sets the pointer to the child primary particle.
void SubParticle::setPrimaryPtr(Sweep::Primary *const pri) {m_primary = pri;}

// Sets the composition vector.
void SubParticle::SetComposition(const fvector &comp)
{
    // Can only set the composition if this is a leaf sub-particle,
    // otherwise it makes no sense.
    if (m_primary != NULL) {
        m_primary->SetComposition(comp);
        UpdateCache();
    }
}

// Sets the values vector.
void SubParticle::SetValues(const fvector &vals)
{
    // Can only set the values if this is a leaf sub-particle,
    // otherwise it makes no sense.
    if (m_primary != NULL) {
        m_primary->SetValues(vals);
        UpdateCache();
    }
}

// Sets the last update time of the particle.
void SubParticle::SetTime(real t)
{
    m_primary->SetTime(t);
    m_time = t;
}

// Sets the spherical particle diameter
void SubParticle::SetSphDiameter(real diam)
{
    if (m_primary != NULL) {
        m_primary->SetSphDiameter(diam);
        UpdateCache();
    }
}

// Sets the collision diameter of the particle.
void SubParticle::SetCollDiameter(real dcol)
{
    if (m_primary != NULL) {
        m_primary->SetCollDiameter(dcol);
        UpdateCache();
    }
}

// Sets the mobility diameter.
void SubParticle::SetMobDiameter(real dmob)
{
    if (m_primary != NULL) {
        m_primary->SetMobDiameter(dmob);
        UpdateCache();
    }
}

// Sets the surface area, subject to minimum spherical area condition.
void SubParticle::SetSurfaceArea(real surf)
{
    if (m_primary != NULL) {
        m_primary->SetVolume(surf);
        UpdateCache();
    }
}

// Sets the volume.
void SubParticle::SetVolume(real vol)
{
    if (m_primary != NULL) {
        m_primary->SetVolume(vol);
        UpdateCache();
    }
}

// Sets the mass.
void SubParticle::SetMass(real m)
{
    if (m_primary != NULL) {
        m_primary->SetMass(m);
        UpdateCache();
    }
}

// PARTICLE ADJUSTMENT AND PROCESSES.

// Adjusts the particle with the given composition and
// values changes n times.
unsigned int SubParticle::Adjust(const fvector &dcomp,
                                 const fvector &dvalues,
                                 unsigned int n
                                 )
{
    unsigned int m = n;

    // This is a leaf-node sub-particle as it contains a
    // primary particle.  The adjustment is applied to
    // the primary.
    m = m_primary->Adjust(dcomp, dvalues, n);

    // Where-ever the adjustment has been applied this sub-particle must
    // now update its cache.
    UpdateCache();

    return m;
}

/*!
 * Combines this particle with another.
 *
 * \param[in]       rhs         Particle to add to current instance
 * \param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 * \param[in,out]   rand_u01    Pointer to function that generates U[0,1] deviates
 *
 * \return      Reference to the current instance after rhs has been added
 */
SubParticle &SubParticle::Coagulate(const SubParticle &rhs, int (*rand_int)(int, int),
                                  Sweep::real(*rand_u01)())
{
    if (rhs.m_aggcache != NULL) {
        if ((m_aggcache==NULL) || (m_aggcache->ID() != rhs.m_aggcache->ID())) {
            delete m_aggcache;
            m_aggcache = rhs.m_aggcache->Clone();
            m_aggcache->SetParent(*this);
        } else {
            *m_aggcache = *rhs.m_aggcache;
        }
    } else {
        delete m_aggcache;
        m_aggcache = NULL;
    }

    m_primary->Coagulate(*rhs.Primary(), rand_int, rand_u01);
    UpdateCache();

    return *this;
}

// Sinters the sub-particle for the given time using the given
// sintering model.
void SubParticle::Sinter(real dt, const Cell &sys,
                         const Processes::SinteringModel &model)
{
    m_primary->Sinter(dt, sys, model);
}
// PARTICLE UPDATE AND CHECKING.

/*!
 *  Recalculates the derived properties from the unique properties.
 *  This function moves down the tree of subparticles from the top to the bottom.
 */
void SubParticle::UpdateCache(void)
{
    // Get cache from primary particle.
    m_primary->UpdateCache();
    ParticleCache::operator=(*m_primary);

    // Also copy composition and tracking data from primary
    m_comp = m_primary->Composition();
    m_values = m_primary->Values();
    m_createt = m_primary->CreateTime();
    m_time    = m_primary->LastUpdateTime();

	m_numsubpart=1;

	// Update the aggregate details from the primary
	if(m_aggcache != NULL)
	    *m_aggcache = *m_primary;

}

/*!
 * Check that this is a valid particle and has not been changed so that it no
 * longer belongs to the particle phase.
 * 
 * It is not quite clear to what extent this method should check the integrity
 * of the data structure and to what extent it should verify the physical
 * meaning.  I think the two purposes may have become a little mixed.
 * (riap2 24 Jan 2011)
 *
 *@return   Particle is valid
 */
bool SubParticle::IsValid() const
{
    return (m_primary != NULL) && (m_primary->IsValid());
}


// READ/WRITE/COPY.

// Creates a clone of the particle.
SubParticle *const SubParticle::Clone() const
{
    return new SubParticle(*this);

}

// Writes the object to a binary stream.
void SubParticle::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output the base class.
        ParticleCache::Serialize(out);

        // Write number of components.
        unsigned int n = (unsigned int)m_comp.size();
        out.write((char*)&n, sizeof(n));

        // Write components.
        real val = 0.0;
        for (unsigned int i=0; i!=n; ++i) {
            val = m_comp[i];
            out.write((char*)&val, sizeof(val));
        }

        // Write number of tracker values.
        n = (unsigned int)m_values.size();
        out.write((char*)&n, sizeof(n));

        // Write values.
        for (unsigned int i=0; i!=n; ++i) {
            val = m_values[i];
            out.write((char*)&val, sizeof(val));
        }

        // Write create time.
        val = (double)m_createt;
        out.write((char*)&val, sizeof(val));

        // Write last update time.
        val = (double)m_time;
        out.write((char*)&val, sizeof(val));

        if (m_primary != NULL) {
            // Write primary particle to stream.
            out.write((char*)&trueval, sizeof(trueval));
            ModelFactory::WritePrimary(*m_primary, out);
        } else {
            throw std::logic_error("Subtrees no longer supported");
        }

        if (m_aggcache != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            ModelFactory::WriteCache(*m_aggcache, out);
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

    } else {
        throw std::invalid_argument("Output stream not ready "
                               "(Sweep, SubParticle::Serialize).");
    }
}

// Reads the object from a binary stream.
void SubParticle::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    releaseMem();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        const unsigned int trueval  = 1;
        real val;

        switch (version) {
            case 0:
                m_pmodel = &model;

                // Read base class.
                ParticleCache::Deserialize(in, model);

                // Read number of components.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read components.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_comp.push_back(val);
                }

                // Read number of tracker values.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read values.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_values.push_back(val);
                }

                // Read create time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_createt = val;

                // Read last update time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_time = val;


                // Read if this sub-particle had a primary particle
                // or two sub-particle children.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                if (n == trueval) {
                    // Read primary particle.
                    m_primary    = ModelFactory::ReadPrimary(in, model);
                } else {
                    throw std::logic_error("Subtrees no longer supported");
                }

                // Read aggregation model.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n==1) {
                    m_aggcache = ModelFactory::ReadAggCache(in, *this);
                } else {
                    m_aggcache = NULL;
                }

                break;
            default:
                throw std::runtime_error("Serialized version number is invalid "
                                    "(Sweep, SubParticle::Deserialize).");
        }
    } else {
        throw std::invalid_argument("Input stream not ready "
                               "(Sweep, SubParticle::Deserialize).");
    }
}


// MEMORY MANAGEMENT.

// Release all memory associated with the SubParticle object.
void SubParticle::releaseMem(void)
{
    delete m_primary;    	m_primary    = NULL;

    m_comp.clear();
    m_values.clear();

    // Clear aggregation model cache.
    delete m_aggcache;
}

// Initialisation routine.
void SubParticle::init(void)
{
    releaseMem();
}

// PARTICLE CREATE TIME.

// Returns the particle create time.
real SubParticle::CreateTime() const {return m_createt;}


// LAST UPDATE TIME.

// Returns the particle last update time.
real SubParticle::LastUpdateTime() const {return m_time;}
// PARTICLE COMPOSITION.

// Returns the composition vector.
const fvector &SubParticle::Composition() const
{
    return m_comp;
}

// Returns the ith component value.  Returns 0.0 if i is invalid.
real SubParticle::Composition(unsigned int i) const
{
    if (i < m_comp.size()) {
        return m_comp[i];
    } else {
        return 0.0;
    }
}

// TRACKER VARIABLE VALUES.

// Returns the tracker value vector.
const fvector &SubParticle::Values() const
{
    return m_values;
}

// Returns the ith tracker variable value.  Returns 0.0 if i is invalid.
real SubParticle::Values(unsigned int i) const
{
    if (i < m_values.size()) {
        return m_values[i];
    } else {
        return 0.0;
    }
}


/*!
 * Pass through to ParticleCache
 */
Sweep::real SubParticle::SphDiameter() const
{
    return ParticleCache::SphDiameter();
}

/*!
 * Pass through to ParticleCache
 */
Sweep::real SubParticle::CollDiameter() const
{
    return ParticleCache::CollDiameter();
}

/*!
 * Pass through to ParticleCache
 */
Sweep::real SubParticle::MobDiameter() const
{
    return ParticleCache::MobDiameter();
}

/*!
 * Pass through to ParticleCache
 */
Sweep::real SubParticle::SurfaceArea() const
{
    return ParticleCache::SurfaceArea();
}

/*!
 * Pass through to ParticleCache
 */
Sweep::real SubParticle::SphSurfaceArea() const
{
    return ParticleCache::SphSurfaceArea();
}

/*!
 * Pass through to ParticleCache
 */
Sweep::real SubParticle::Volume() const
{
    return ParticleCache::Volume();
}

/*!
 * Pass through to ParticleCache
 */
Sweep::real SubParticle::Mass(void) const
{
    return ParticleCache::Mass();
}


