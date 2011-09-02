#include "swp_subparticle.h"
#include "swp_particle_model.h"
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
 * This class will eventually removed so that Particle become self-contained.
 *
 * DO NOT add methods or data to this class, since the whole
 * class is due to be removed. (riap2 28 jan 2011).  New functionality
 * should be placed in Particle or Primary.
 *
 */

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
SubParticle::SubParticle(void)
: m_primary(NULL),
  m_createt(0.0), m_time(0.0), m_aggcache(NULL)
{
}

// Initialising constructor.
SubParticle::SubParticle(real t, const Sweep::ParticleModel &model)
: m_primary(NULL),
  m_createt(0.0), m_time(0.0)
{
	m_aggcache = ModelFactory::CreateAggCache(model.AggModel());

}

// Initialising constructor (from Primary particle).
SubParticle::SubParticle(Sweep::Primary &pri)
{
    m_createt = pri.CreateTime();
    m_time = pri.CreateTime();
    m_primary    = &pri;
	m_aggcache = pri.CreateAggCache();
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
  m_aggcache(NULL)
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

        m_createt = rhs.m_createt;
        m_time = rhs.m_time;

        if (rhs.m_aggcache != NULL) {
            if ((m_aggcache==NULL) || (m_aggcache->ID() != rhs.m_aggcache->ID())) {
                delete m_aggcache;
                m_aggcache = rhs.m_aggcache->Clone();
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

// Sets the last update time of the particle.
void SubParticle::SetTime(real t)
{
    m_primary->SetTime(t);
    m_time = t;
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
void SubParticle::Sinter(real dt, Cell &sys,
                         const Processes::SinteringModel &model,
                         real (*rand_u01)())
{
    m_primary->Sinter(dt, sys, model, rand_u01);
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

    m_createt = m_primary->CreateTime();
    m_time    = m_primary->LastUpdateTime();

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

        real val;

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
                    m_aggcache = ModelFactory::ReadAggCache(in);
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
    return m_primary->Composition();
}

// Returns the ith component value.  Returns 0.0 if i is invalid.
real SubParticle::Composition(unsigned int i) const
{
    if (i < Composition().size()) {
        return Composition()[i];
    } else {
        return 0.0;
    }
}

// TRACKER VARIABLE VALUES.

/*!
 * Returns the tracker value vector.
 */
const fvector &SubParticle::Values() const
{
    return m_primary->Values();
}

// Returns the ith tracker variable value.  Returns 0.0 if i is invalid.
real SubParticle::Values(unsigned int i) const
{
    if (i < Values().size()) {
        return Values()[i];
    } else {
        return 0.0;
    }
}


/*!
 * Pass through to primary particle
 */
Sweep::real SubParticle::SphDiameter() const
{
    return m_primary->SphDiameter();
}

/*!
 * Pass through to primary particle
 */
Sweep::real SubParticle::CollDiameter() const
{
    return m_primary->CollDiameter();
}

/*!
 * Pass through to primary particle
 */
Sweep::real SubParticle::MobDiameter() const
{
    return m_primary->MobDiameter();
}

/*!
 * Pass through to primary particle
 */
Sweep::real SubParticle::SurfaceArea() const
{
    return m_primary->SurfaceArea();
}

/*!
 * Pass through to primary particle
 */
Sweep::real SubParticle::SphSurfaceArea() const
{
    return m_primary->SphSurfaceArea();
}

/*!
 * Pass through to primary particle
 */
Sweep::real SubParticle::Volume() const
{
    return m_primary->Volume();
}

/*!
 * Pass through to primary particle
 */
Sweep::real SubParticle::Mass(void) const
{
    return m_primary->Mass();
}

/*!
 * Provide an interface that allows run time specification of particle properties
 * for use in process rate calculations.  It is currently used for some surface
 * reactions.  Where possible, the use of specific accessors should be preferred.
 *
 *@param[in]    id      Symbolic index of property required
 *
 *@return       Requested particle property
 *
 *@pre          A valid primary particle is present (the data is not stored in SubParticle)
 */
real SubParticle::Property(PropID id) const
{
    switch (id) {
        case iDsph:      // Equivalent sphere diameter.
            return SphDiameter();
        case iDcol:   // Collision diameter.
            return CollDiameter();
        case iDmob:   // Mobility diameter.
            return MobDiameter();
        case iS:      // Surface area.
            return SurfaceArea();
        case iV:      // Volume.
            return Volume();
        case iM:      // Mass.
            return Mass();
        // Collision rate properties:
        case iD2:
            return CollDiameter() * CollDiameter();
        case iD_1:
            return 1.0 / CollDiameter();
        case iD_2:
            return 1.0 / CollDiameter() / CollDiameter();
        case iM_1_2:
            return 1.0 / std::sqrt(Mass());
        case iD2_M_1_2:
            return CollDiameter() * CollDiameter() / std::sqrt(Mass());
        case iFS:
            throw std::logic_error("Free surface no longer supported (SubParticle::Property)");
            return 0.0;
        case iNumCarbon:
            throw std::logic_error("Number of Carbon no longer supported (SubParticle::Property)");
            return 0.0;
        case -1:
            // Special case property, used to select particles
            // uniformly.
            return 1.0;
        default:
            throw std::logic_error("Unrecognised property requested (SubParticle::Property)");
            return 0.0;
    }
}

