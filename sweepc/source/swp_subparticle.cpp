#include "swp_subparticle.h"
#include "swp_particle_model.h"
#include "swp_submodel.h"
#include "swp_model_factory.h"
#include "rng.h"
#include <cmath>
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
SubParticle::SubParticle(void)
: ParticleCache(), m_parent(NULL), m_leftchild(NULL), 
  m_rightchild(NULL), m_primary(NULL)
{
}

// Initialising constructor.
SubParticle::SubParticle(real t, const Sweep::ParticleModel &model)
: ParticleCache(t, model), m_parent(NULL), m_leftchild(NULL), 
  m_rightchild(NULL), m_primary(NULL)
{
}

// Initialising constructor (from Primary particle).
SubParticle::SubParticle(Sweep::Primary &pri)
{
    ParticleCache(pri.CreateTime(), *pri.ParticleModel());
    m_parent     = NULL;
    m_leftchild  = NULL;
    m_rightchild = NULL;
    m_primary    = &pri;
}

// Copy constructor.
SubParticle::SubParticle(const SubParticle &copy)
{
    // Use assignment operator.
    m_parent     = NULL;
    m_leftchild  = NULL;
    m_rightchild = NULL;
    m_primary    = NULL;
    *this        = copy;
}

// Stream-reading constructor.
SubParticle::SubParticle(std::istream &in, const Sweep::ParticleModel &model)
: m_parent(NULL), m_leftchild(NULL), 
  m_rightchild(NULL), m_primary(NULL)
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
        ParticleCache::operator=(rhs);

        // Copy left child.
        if (rhs.m_leftchild != NULL) {
            if (m_leftchild == NULL) {
                m_leftchild = rhs.m_leftchild->Clone();
            } else {
                *m_leftchild = *rhs.m_leftchild;
            }
        } else {
            delete m_leftchild; m_leftchild = NULL;
        }

        // Copy right child.
        if (rhs.m_rightchild != NULL) {
            if (m_rightchild == NULL) {
                m_rightchild = rhs.m_rightchild->Clone();
            } else {
                *m_rightchild = *rhs.m_rightchild;
            }
        } else {
            delete m_rightchild; m_rightchild = NULL;
        }

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
    }
    return *this;
}

// Addition-assignment operator.  This implements coagulation.
SubParticle &SubParticle::operator+=(const SubParticle &rhs)
{
    return Coagulate(rhs);
}

// Addition operator.  This also implements coagulation.
const SubParticle SubParticle::operator+(const SubParticle &rhs) const
{
    return SubParticle(*this) += rhs;
}


// PARENT SUB-PARTICLE.

// Returns the sub-particle's parent.
SubParticle *const SubParticle::Parent(void)
{
    return m_parent;
}

// Sets the sub-particle's parent.
void SubParticle::setParent(SubParticle *const p)
{
    m_parent = p;
}


// PRIMARY PARTICLE CHILD.

// Returns a pointer to the child primary particle, NULL if 
// this sub-particle has no primary.
Sweep::Primary *const SubParticle::Primary(void) {return m_primary;}
const Sweep::Primary *const SubParticle::Primary(void) const {return m_primary;}

// Sets the pointer to the child primary particle.
void SubParticle::setPrimaryPtr(Sweep::Primary *const pri) {m_primary = pri;}


// CHILD SUB-PARTICLES IN SUB-PARTICLE TREE.

// Returns a pointer to the left child.
SubParticle *const SubParticle::Left(void) {return m_leftchild;}
const SubParticle *const SubParticle::Left(void) const {return m_leftchild;}

// Returns a pointer to the right child.
SubParticle *const SubParticle::Right(void) {return m_rightchild;}
const SubParticle *const SubParticle::Right(void) const {return m_rightchild;}

// Sets the pointer to the left child.
void SubParticle::setLeftPtr(SubParticle *const sub) {m_leftchild = sub;}

// Sets the pointer to the right child.
void SubParticle::setRightPtr(SubParticle *const sub) {m_rightchild = sub;}


// PROPERTY SETTING OVERRIDES (FROM PARTICLE CACHE).

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
    if (m_primary != NULL) {
        m_primary->SetTime(t);
        m_time = t;
    } else {
        m_leftchild->SetTime(t);
        m_rightchild->SetTime(t);
    }
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
                                 SubModels::SubModelType model_id,
                                 int id,
                                 unsigned int n
                                 )
{
    unsigned int m = n;

    if (m_primary != NULL) {
        // This is a leaf-node sub-particle as it contains a
        // primary particle.  The adjustment is applied to
        // the primary.
        m = m_primary->Adjust(dcomp, dvalues, n);
    } else {
        // This sub-particle is somewhere in the tree.  The adjustment
        // is passed down the tree until it reaches a primary particle.
        // We just the left or right child sub-particle by using the
        // given weighting variable.

        // Get the sub-particle weightings.
        real left=0.0, right=0.0;
        if (model_id == SubModels::BasicModel_ID) {
            // This is a property is the basic particle cache.
            left  = m_leftchild->Property((ParticleCache::PropID)id);
            right = m_rightchild->Property((ParticleCache::PropID)id);
        } else {
            // This is a sub-model property.
            left  = m_leftchild->SubModel(model_id)->Property(id);
            right = m_rightchild->SubModel(model_id)->Property(id);
        }

        // Randomly select which branch to go down.
        real r = rnd() * (left + right);
        if (r <= left) {
            // Apply adjustment to the left sub-particle.
            m = m_leftchild->adjustLoop(dcomp, dvalues, model_id, id, r, n);
        } else {
            // Apply adjustment to the right sub-particle.
            m = m_rightchild->adjustLoop(dcomp, dvalues, model_id, id, r-=left, n);
        }
    }

    // Where-ever the adjustment has been applied this sub-particle must
    // now update it's cache.
    UpdateCache();

    return m;
}

// Adjusts the particle with the given composition and 
// values changes n times.  This particular function is 
// used internally by the SubParticle class to reuse the 
// random variable used to select the correct leaf which
// is generated by the root sub-particle.
unsigned int SubParticle::adjustLoop(const fvector &dcomp, 
                                     const fvector &dvalues, 
                                     SubModels::SubModelType model_id,
                                     int id, real r,
                                     unsigned int n
                                     )
{
    unsigned int m = n;

    if (m_primary != NULL) {
        // This is a leaf-node sub-particle as it contains a
        // primary particle.  The adjustment is applied to
        // the primary.
        m = m_primary->Adjust(dcomp, dvalues, n);
    } else {
        // This sub-particle is somewhere in the tree.  The adjustment
        // is passed down the tree until it reaches a primary particle.
        // We just the left or right child sub-particle by using the
        // given weighting variable.

        // Get the sub-particle weightings.
        real left=0.0, right=0.0;
        if (model_id == SubModels::BasicModel_ID) {
            // This is a property is the basic particle cache.
            left  = m_leftchild->Property((ParticleCache::PropID)id);
            right = m_rightchild->Property((ParticleCache::PropID)id);
        } else {
            // This is a sub-model property.
            left  = m_leftchild->SubModel(model_id)->Property(id);
            right = m_rightchild->SubModel(model_id)->Property(id);
        }

        // Randomly select which branch to go down.
        if (r <= left) {
            // Apply adjustment to the left sub-particle.
            m = m_leftchild->adjustLoop(dcomp, dvalues, model_id, id, r, n);
        } else {
            // Apply adjustment to the right sub-particle.
            m = m_rightchild->adjustLoop(dcomp, dvalues, model_id, id, r-=left, n);
        }
    }

    // Where-ever the adjustment has been applied this sub-particle must
    // now update it's cache.
    UpdateCache();

    return m;
}


// Combines this particle with another.
SubParticle &SubParticle::Coagulate(const SubParticle &rhs)
{
    if (m_pmodel->UseSubPartTree()) {
        // We are using the sub-particle tree.  This means
        // we branch the tree to incorporate the new RHS
        // sub-particle.

        // Create new sub-particles.
        SubParticle *lsp = new SubParticle(m_time, *m_pmodel);
        SubParticle *rsp = new SubParticle(m_time, *m_pmodel);

        // The new sub-particles' parent is this particle.
        lsp->m_parent = this;
        rsp->m_parent = this;

        // Set left child.   The new sub-particle's children are 
        // the children of this particle.
        lsp->m_primary    = m_primary;
        lsp->m_leftchild  = m_leftchild;
        lsp->m_rightchild = m_rightchild;

        // Set right child.  The children are the children
        // of the right-hand side.
        rsp->m_primary    = rhs.m_primary;
        rsp->m_leftchild  = rhs.m_leftchild;
        rsp->m_rightchild = rhs.m_rightchild;

        // Clear pointers in RHS (so that when it is deleted,
        // it doesn't also delete the sub-tree).
        const_cast<SubParticle&>(rhs).m_primary    = NULL;
        const_cast<SubParticle&>(rhs).m_leftchild  = NULL;
        const_cast<SubParticle&>(rhs).m_rightchild = NULL;

        // This particle now has two sub-particle children.
        m_primary    = NULL;
        m_leftchild  = lsp;
        m_rightchild = rsp;

        // Update the caches on the particle and it's children.
        lsp->UpdateCache();
        rsp->UpdateCache();
        UpdateCache();
    } else {
        // If we are not using the sub-particle tree, then we just need
        // to coagulate the primary particles.
        m_primary->Coagulate(*rhs.Primary());
        UpdateCache();
    }

    return *this;
}


// PARTICLE UPDATE AND CHECKING.

// Recalculates the derived properties from the unique properties.
void SubParticle::UpdateCache(void)
{
    if (m_primary != NULL) {
        // Get cache from primary particle.
        m_primary->UpdateCache();
        ParticleCache::operator=(*m_primary);
    } else {
        // The cache is the sum of the left and right child caches.
        m_leftchild->UpdateCache();
        m_rightchild->UpdateCache();
        ParticleCache::operator=(*m_leftchild);
        ParticleCache::operator+=(*m_rightchild);
    }
}

// Check the that the particle is valid by querying the
// validity conditions of the models and ensuring that it 
// contains any components.
bool SubParticle::IsValid() const
{
    fvector::const_iterator i;
    for (i=m_comp.begin(); i!=m_comp.end(); ++i) {
        if (*i > 0.0) return true;
    }
    return false;
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

        if (m_primary != NULL) {
            // Write primary particle to stream.
            out.write((char*)&trueval, sizeof(trueval));
            ModelFactory::WritePrimary(*m_primary, out);
        } else {
            // Write left and right children to stream.
            out.write((char*)&falseval, sizeof(falseval));
            m_leftchild->Serialize(out);
            m_rightchild->Serialize(out);
        }
    } else {
        throw invalid_argument("Output stream not ready "
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

        switch (version) {
            case 0:
                // Read base class.
                ParticleCache::Deserialize(in, model);

                // Read if this sub-particle had a primary particle
                // or two sub-particle children.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                if (n == trueval) {
                    // Read primary particle.
                    m_primary    = ModelFactory::ReadPrimary(in, model);
                    m_leftchild  = NULL;
                    m_rightchild = NULL;
                } else {
                    // Read sub-particle children.
                    m_primary    = NULL;
                    m_leftchild  = new Sweep::SubParticle(in, model);
                    m_rightchild = new Sweep::SubParticle(in, model);
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, SubParticle::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SubParticle::Deserialize).");
    }
}


// MEMORY MANAGEMENT.

// Release all memory associated with the SubParticle object.
void SubParticle::releaseMem(void)
{
    delete m_leftchild;  m_leftchild  = NULL;
    delete m_rightchild; m_rightchild = NULL;
    delete m_primary;    m_primary    = NULL;
}

// Initialisation routine.
void SubParticle::init(void)
{
    m_parent = NULL;
    releaseMem();
}
