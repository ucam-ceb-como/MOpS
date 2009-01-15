#include "swp_subparticle.h"
#include "swp_particle_model.h"
#include "swp_submodel.h"
#include "swp_model_factory.h"
#include "rng.h"
#include "swp_cell.h"
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>
using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
SubParticle::SubParticle(void)
: ParticleCache(), m_parent(NULL), m_leftchild(NULL), 
  m_rightchild(NULL), m_primary(NULL)
{
	m_sintered=0;
	m_connect_time=0;
	m_avg_diam=0;
	m_leftsinter=NULL;
	m_rightsinter=NULL;
	m_real_surface=0;
	m_sumvol_sinterpart=0;
	dV_left=0;
	dV_right=0;
	m_sinter_level=0;


}

// Initialising constructor.
SubParticle::SubParticle(real t, const Sweep::ParticleModel &model)
: ParticleCache(t, model), m_parent(NULL), m_leftchild(NULL), 
  m_rightchild(NULL), m_primary(NULL)
{
	m_sintered=0;
	m_connect_time=0;
	m_avg_diam=0;
	m_leftsinter=NULL;
	m_rightsinter=NULL;
	m_real_surface=0;
	m_sumvol_sinterpart=0;
	dV_left=0;
	dV_right=0;
	m_sinter_level=0;

}

// Initialising constructor (from Primary particle).
SubParticle::SubParticle(Sweep::Primary &pri)
{  
    ParticleCache(pri.CreateTime(), *pri.ParticleModel());
    m_parent     = NULL;
    m_leftchild  = NULL;
    m_rightchild = NULL;
    m_primary    = &pri;
	m_leftsinter=NULL;
	m_rightsinter=NULL;
	m_sintered=0;
	m_connect_time=0;
	m_avg_diam=pri.CollDiameter();
	m_real_surface=0;
	m_sumvol_sinterpart=0;
	m_freesurface=pri.SurfaceArea();
	vol_sinter=m_primary->Volume();
	dV_left=0;
	dV_right=0;
	m_sinter_level=0;
}

// Copy constructor.
SubParticle::SubParticle(const SubParticle &copy)
{
    // Use assignment operator.
    m_parent     = NULL;
    m_leftchild  = NULL;
    m_rightchild = NULL;
    m_primary    = NULL;
	m_leftsinter=NULL;
	m_rightsinter=NULL;
	m_sintered=0;
	m_connect_time=0;
	m_avg_diam=0;
	m_real_surface=0;
	m_sumvol_sinterpart=0;
	dV_left=0;
	dV_right=0;
	m_sinter_level=0;
	*this        = copy;

}

// Stream-reading constructor.
SubParticle::SubParticle(std::istream &in, const Sweep::ParticleModel &model)
: m_parent(NULL), m_leftchild(NULL), 
  m_rightchild(NULL), m_primary(NULL)
 {
    init();
    Deserialize(in, model);
	/*m_sintered=0;
	m_connect_time=0;
	m_avg_diam=0;
	SubModels::SubModelType model = SubModels::BasicModel_ID;
	ParticleCache::PropID id=ParticleCache::iUniform;
	m_leftsinter=m_leftchild->SelectLeaf(model,id);
	m_rightsinter=m_rightchild->SelectLeaf(model,id);*/
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
			m_leftchild->m_parent=this;
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
			m_rightchild->m_parent=this;
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
	//UpdateSinterParticles();
	UpdatethisSinterParticle(this,&rhs);
	m_real_surface=rhs.m_real_surface;
	m_real_surface_init=rhs.m_real_surface_init;
	m_sintered=rhs.m_sintered;
	m_sph_surfacearea=rhs.m_sph_surfacearea;
	m_sumvol_sinterpart=rhs.m_sumvol_sinterpart;
	m_connect_time=rhs.m_connect_time;
	m_avg_diam=rhs.m_avg_diam;
	m_rightsinterdiam=rhs.m_rightsinterdiam;
	m_leftsinterdiam=rhs.m_leftsinterdiam;
	m_freesurface=rhs.m_freesurface;
	m_sinter_level=rhs.m_sinter_level;
	vol_sinter=rhs.vol_sinter;
	dV_left=rhs.dV_left;
	dV_right=rhs.dV_right;
	m_sinter_level=rhs.m_sinter_level;
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


//Sintering tree
void SubParticle::setSintered(real sintered)
{
	m_sintered=sintered;
}

real SubParticle::Sintered()
{
	return m_sintered;
}


real SubParticle::Diamter()
{
	return m_diam;
}

real SubParticle::avgeomdiam(double oneovernumsubpart)
{
	Primary::PropID primid=Primary::iD;
	if (m_primary!=NULL)
	return pow(m_primary->Property(primid),oneovernumsubpart);
	else 
	{
		return m_leftchild->avgeomdiam(oneovernumsubpart)*m_rightchild->avgeomdiam(oneovernumsubpart);
	}
}

real SubParticle::SubPartSurfaceArea(void) const
{
	return m_real_surface;
}

real SubParticle::SubPartSphSurfaceArea(void) const
{
	return m_sph_surfacearea;
}

real SubParticle::SubPartSumVol(void) const
{
	return m_sumvol_sinterpart;
		   
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

// Selects a leaf from the sub-particle tree below this
// sub-particle using the given property ID to weight the
// primaries.
Sweep::SubParticle *const SubParticle::SelectLeaf(SubModels::SubModelType model_id, 
                                                  int id)
{
    if (m_primary != NULL) {
        // This is a leaf-level sub-particle.
        return this;
    } else {
        // This sub-particle is somewhere in the tree.

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
            // Select from the left sub-particle.
            return m_leftchild->selLeafLoop(model_id, id, r);
        } else {
            // Select from the right sub-particle.
            return m_rightchild->selLeafLoop(model_id, id, r-=left);
        }
    }
}

// Selects a leaf from the tree.  This particular function is 
// used internally by the SubParticle class to reuse the 
// random variable used to select the correct leaf which
// is generated by the root sub-particle.
Sweep::SubParticle *const SubParticle::selLeafLoop(SubModels::SubModelType model_id,
                                                   int id, real r)
{ 
    if (m_primary != NULL) {
        // This is a leaf-level sub-particle.
        return this;
    } else {
        // This sub-particle is somewhere in the tree.

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
            // Select from the left sub-particle.
            return m_leftchild->selLeafLoop(model_id, id, r);
        } else {
            // Select from the right sub-particle.
            return m_rightchild->selLeafLoop(model_id, id, r-=left);
        }
    }
}



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
		m_time = t;
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


// Sets the sintering level
void SubParticle::SetSintering(real sint)
{
	m_sintered=sint;
    
}



real SubParticle::GetSintering()
{
	return m_sintered;
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

		SubParticle *temp;
		

		if (rnd()>0.5)
		{
			// Set left child.   The new sub-particle's children are 
			// the children of this particle.
			lsp->m_primary    = m_primary;
			lsp->m_leftchild  = m_leftchild;
			lsp->m_rightchild = m_rightchild;
			lsp->m_connect_time = m_connect_time;
			lsp->m_sintered		= m_sintered;
			lsp->m_real_surface = m_real_surface;
			lsp->m_real_surface_init = m_real_surface_init;
			lsp->m_sumvol_sinterpart = m_sumvol_sinterpart;
			lsp->m_sph_surfacearea = m_sph_surfacearea;
			lsp->m_avg_diam = m_avg_diam;
			if(lsp->m_primary==NULL)
			{
				lsp->m_rightchild->m_parent=lsp;
				lsp->m_leftchild->m_parent=lsp;
			}
			lsp->m_rightsinter = m_rightsinter;
			lsp->m_leftsinter = m_leftsinter;
			lsp->m_leftsinterdiam = m_leftsinterdiam;
			lsp->m_rightsinterdiam = m_rightsinterdiam;
			lsp->vol_sinter = vol_sinter;
			lsp->m_sinter_level = m_sinter_level;
			lsp->dV_left = dV_left;
			lsp->dV_right = dV_right;
					
//			UpdateTree_sinter(this);
//			UpdateCache_sinter(this);

			// Set right child.  The children are the children
			// of the right-hand side.
			rsp->m_primary    = rhs.m_primary;
			rsp->m_leftchild  = rhs.m_leftchild;
			rsp->m_rightchild = rhs.m_rightchild;
			rsp->m_connect_time = rhs.m_connect_time;
			rsp->m_sintered=rhs.m_sintered;
			rsp->m_real_surface = rhs.m_real_surface;
			rsp->m_real_surface_init = rhs.m_real_surface_init;
			rsp->m_sumvol_sinterpart = rhs.m_sumvol_sinterpart;
			rsp->m_sph_surfacearea = rhs.m_sph_surfacearea;
			rsp->m_avg_diam = rhs.m_avg_diam;
			if(rsp->m_primary==NULL)
			{
				rsp->m_rightchild->m_parent=rsp;
				rsp->m_leftchild->m_parent=rsp;
			}
    		rsp->m_rightsinter = rhs.m_rightsinter;
			rsp->m_leftsinter = rhs.m_leftsinter;
			rsp->m_rightsinterdiam = rhs.m_rightsinterdiam;
			rsp->m_leftsinterdiam = rhs.m_leftsinterdiam;
			temp = &const_cast<SubParticle&>(rhs);
			rsp->vol_sinter = rhs.vol_sinter;
			rsp->m_sinter_level = rhs.m_sinter_level;
			rsp->dV_left = rhs.dV_left;
			rsp->dV_right = rhs.dV_right;
//			UpdateTree_sinter(temp);
//			rsp->UpdateCache_sinter(temp);
		}

		else
		{
			// Set right child.   The new sub-particle's children are 
			// the children of this particle.
			rsp->m_primary    = m_primary;
			rsp->m_leftchild  = m_leftchild;
			rsp->m_rightchild = m_rightchild;
			rsp->m_connect_time = m_connect_time;
			rsp->m_sintered	= m_sintered;
			rsp->m_real_surface = m_real_surface;
			rsp->m_real_surface_init = m_real_surface_init;
			rsp->m_sumvol_sinterpart = m_sumvol_sinterpart;
			rsp->m_sph_surfacearea = m_sph_surfacearea;
			rsp->m_avg_diam = m_avg_diam;
			if(rsp->m_primary==NULL)
			{
				rsp->m_rightchild->m_parent=rsp;
				rsp->m_leftchild->m_parent=rsp;
			}
			rsp->m_rightsinter = m_rightsinter;
			rsp->m_leftsinter = m_leftsinter;
			rsp->m_leftsinterdiam = m_leftsinterdiam;
			rsp->m_rightsinterdiam = m_rightsinterdiam;
			rsp->vol_sinter = vol_sinter;
			rsp->m_sinter_level = m_sinter_level;
			rsp->dV_left = dV_left;
			rsp->dV_right = dV_right;
//			UpdateTree_sinter(this);
//			UpdateCache_sinter(this);

			// Set left child.  The children are the children
			// of the right-hand side.
			lsp->m_primary    = rhs.m_primary;
			lsp->m_leftchild  = rhs.m_leftchild;
			lsp->m_rightchild = rhs.m_rightchild;
			lsp->m_connect_time = rhs.m_connect_time;
			lsp->m_sintered	= rhs.m_sintered;
			lsp->m_real_surface = rhs.m_real_surface;
			lsp->m_real_surface_init = rhs.m_real_surface_init;
			lsp->m_sumvol_sinterpart = rhs.m_sumvol_sinterpart;
			lsp->m_sph_surfacearea = rhs.m_sph_surfacearea;
			lsp->m_avg_diam = rhs.m_avg_diam;
			if(lsp->m_primary==NULL)
			{
				lsp->m_rightchild->m_parent=lsp;
				lsp->m_leftchild->m_parent=lsp;
			}
			lsp->m_rightsinter = rhs.m_rightsinter;
			lsp->m_leftsinter = rhs.m_leftsinter;
			lsp->m_leftsinterdiam = rhs.m_leftsinterdiam;
			lsp->m_rightsinterdiam = rhs.m_rightsinterdiam;
			lsp->vol_sinter = rhs.vol_sinter;
			lsp->m_sinter_level = rhs.m_sinter_level;
			lsp->dV_left = rhs.dV_left;
			lsp->dV_right = rhs.dV_right;
			temp = &const_cast<SubParticle&>(rhs);
		}
		m_connect_time=0;
		m_sintered=0;
        // Clear pointers in RHS (so that when it is deleted,
        // it doesn't also delete the sub-tree).
        const_cast<SubParticle&>(rhs).m_primary    = NULL;
        const_cast<SubParticle&>(rhs).m_leftchild  = NULL;
        const_cast<SubParticle&>(rhs).m_rightchild = NULL;
        const_cast<SubParticle&>(rhs).m_leftsinter  = NULL;
        const_cast<SubParticle&>(rhs).m_rightsinter = NULL;
		const_cast<SubParticle&>(rhs).m_parent = NULL;
		delete &rhs;

        // This particle now has two sub-particle children.
        m_primary    = NULL;
        m_leftchild  = lsp;
        m_rightchild = rsp;

        // Update the caches on the particle and it's children.
        lsp->UpdateCache_thispart();
        rsp->UpdateCache_thispart();
//        UpdateCache_thispart();
		  
		  UpdateCache();
		  UpdateTree();

		SubModels::SubModelType model = SubModels::BasicModel_ID;
		ParticleCache::PropID id=ParticleCache::iFS;
		m_leftsinter=m_leftchild->SelectLeaf(model,id);
		m_rightsinter=m_rightchild->SelectLeaf(model,id);
		m_real_surface=m_leftsinter->m_primary->SurfaceArea()+m_rightsinter->m_primary->SurfaceArea();
		m_real_surface_init=m_real_surface;
		m_sumvol_sinterpart = m_leftsinter->m_primary->Volume()+m_rightsinter->m_primary->Volume();
		m_sph_surfacearea = PI * pow(m_sumvol_sinterpart  * 6.0 / PI, TWO_THIRDS);
		Primary::PropID primid=Primary::iD;
		m_leftsinterdiam = m_leftsinter->m_primary->Property(primid);
		m_rightsinterdiam = m_rightsinter->m_primary->Property(primid);
		real Dfreesurface=min((m_leftsinter->m_diam)*(m_leftsinter->m_diam)*2*PI,(m_rightsinter->m_diam)*(m_rightsinter->m_diam)*2*PI);
		m_rightsinter->m_freesurface-=Dfreesurface;
		if (m_rightsinter->m_freesurface<0) m_rightsinter->m_freesurface=0;
		m_leftsinter->m_freesurface-=Dfreesurface;
		if ( m_leftsinter->m_freesurface<0)	m_leftsinter->m_freesurface=0;

		m_sinter_level=0.;
		m_leftsinter->vol_sinter=m_leftsinter->m_primary->Volume();
		m_rightsinter->vol_sinter=m_rightsinter->m_primary->Volume();
		dV_left=0.;
		dV_right=0.;


/*		UpdateTree_sinter(temp);
		UpdateCache_sinter(temp);
		UpdateTree_sinter(this);
		UpdateCache_sinter(this);*/

//		FindRoot()->CheckTree();
//		cout << "Checktree after coag passed\n";
//		m_real_surface=m_leftchild->m_real_surface+m_rightchild->m_real_surface;
//		m_sumvol_sinterpart = m_leftchild->m_sumvol_sinterpart+m_rightchild->m_sumvol_sinterpart;

    } else {
        // If we are not using the sub-particle tree, then we just need
        // to coagulate the primary particles.
        m_primary->Coagulate(*rhs.Primary());
        UpdateCache();
    }

    return *this;
}
/*
bool FileExists( const char* FileName )
{
    FILE* fp = NULL;

    //will not work if you do not have read permissions

    //to the file, but if you don't have read, it

    //may as well not exist to begin with.

    fp = fopen( FileName, "rb" );
    if( fp != NULL )
    {
        fclose( fp );
        return true;
    }

    return false;
}
*/

// Sinters this particle with another particle
void SubParticle::SinterPart()
{	bool lefttree;
	// Chose the subparticle according to the following property
	//if (!FileExists("beforesinterprim.inp"))
//	{ 
/*		ParticleCache::PropID idsubtreeprint=ParticleCache::iFS;
		ofstream out;
		out.open("beforesinterprim.inp");
	    this->FindRoot()->printSubtree(out,idsubtreeprint);
		out.close();*/
//	}
	if (m_pmodel->UseSubPartTree()) 
	{
		//both childrens are primaries
		if(m_leftchild->m_primary!=NULL && m_rightchild->m_primary!=NULL) 
		{//	FindRoot()->CheckTree();
		//	cout << "Checktree before sinter passed\n";
			//   cout << "sinter";

			//cout << m_leftchild->Property(id);
			m_connect_time=0;
			m_sintered=0;	
			setPrimaryPtr(m_leftchild->m_primary);
			m_primary->Coagulate(*m_rightchild->m_primary);
			m_rightsinter=NULL;
			m_leftsinter=NULL;
			dV_left=0;
			dV_right=0;
	//		out.open("aftersinterprim.inp");
//			printSubtree(out,id);
//			out.close();
			
			UpdateTree_sinter(m_leftchild,this);
			UpdateTree_sinter(m_rightchild,this);

			
			(*m_rightchild).m_primary    = NULL;
			(*m_rightchild).m_parent   = NULL;
			(*m_rightchild).m_leftchild    = NULL;
			(*m_rightchild).m_rightchild    = NULL;
			(*m_rightchild).m_leftsinter  = NULL;
			(*m_rightchild).m_rightsinter    = NULL;
			(*m_leftchild).m_parent   = NULL;
			(*m_leftchild).m_leftchild    = NULL;
			(*m_leftchild).m_rightchild    = NULL;
			(*m_leftchild).m_leftsinter  = NULL;
			(*m_leftchild).m_rightsinter    = NULL;
			(*m_leftchild).m_primary    = NULL;
			vol_sinter=m_leftchild->vol_sinter+m_rightchild->vol_sinter;
			delete m_leftchild;
			delete m_rightchild;
			m_rightchild=NULL;
			m_leftchild=NULL;
	//		UpdateCache();
	//		UpdateTree();
			m_real_surface=m_primary->SurfaceArea();
			m_sumvol_sinterpart = m_primary->Volume();
			m_sph_surfacearea = m_primary->SurfaceArea();
			
			m_leftsinterdiam=0;
			m_rightsinterdiam=0;
			m_sinter_level=0;
			
		

	//		cout <<"break";
	//		FindRoot()->CheckTree();
	//		cout << "Checktree after sinter passed\n";

		}
		else
		{ 
			if (m_rightchild->m_primary==NULL)
			{
				lefttree=true;
			}
			if (m_leftchild->m_primary==NULL)
			{
				lefttree=false;
			}
			if (m_rightchild->m_primary==NULL && m_leftchild->m_primary==NULL)
			{
				if (rnd()>0.5) lefttree=true;
				else lefttree=false;
			}
			

			SubParticle *left=m_leftsinter;
			SubParticle *right=m_rightsinter;

			// if the parent node is the root node we can not delete the node 
			// because the pointer from the particle node points to this node
	
				SubParticle *temp;
				if (lefttree==true)
				{
				//	SubParticle *left=m_leftchild->SelectLeaf(model,id);
				//	SubParticle *right=m_rightchild->SelectLeaf(model,id);
					left->m_primary->Coagulate(*right->m_primary);  
					left->vol_sinter+=right->vol_sinter;
					if (right->Parent()->m_rightchild==right)
						{
							m_leftchild->setParent(right->Parent());
							right->Parent()->setRightPtr(m_leftchild);
						}

					if (right->Parent()->m_leftchild==right)
						{
							m_leftchild->setParent(right->Parent());
							right->Parent()->setLeftPtr(m_leftchild);
						}

					m_leftsinter=m_rightchild->m_leftsinter;
					m_leftsinterdiam=m_rightchild->m_leftsinterdiam;
					m_rightsinter=m_rightchild->m_rightsinter;
					m_rightsinterdiam=m_rightchild->m_rightsinterdiam;
					m_sumvol_sinterpart=m_rightchild->m_sumvol_sinterpart;
					m_real_surface=m_rightchild->m_real_surface;
					m_real_surface_init=m_rightchild->m_real_surface_init;
					m_sph_surfacearea=m_rightchild->m_sph_surfacearea;
					m_sintered=m_rightchild->m_sintered;
					vol_sinter=m_rightchild->vol_sinter;
					m_sinter_level=m_rightchild->m_sinter_level;
					dV_left=m_rightchild->dV_left;
					dV_right=m_rightchild->dV_right;
					m_connect_time=m_rightchild->m_connect_time;
					m_leftchild=m_rightchild->m_leftchild;
					m_leftchild->m_parent=this;
					ParticleCache::operator =(*m_rightchild);
					temp=m_rightchild;
					m_rightchild = m_rightchild->m_rightchild;
					m_rightchild->m_parent=this;
					//todelete=right;
					UpdateCache_sinter(right,left);
					UpdateTree_sinter(right,left);
					
				}

				else		
				
				{	
				//	SubParticle *left=m_leftchild->SelectLeaf(model,id);
				//	SubParticle *right=m_rightchild->SelectLeaf(model,id);
					right->m_primary->Coagulate(*left->m_primary);
					right->vol_sinter+=left->vol_sinter;
					if (left->Parent()->m_rightchild==left)
					{
						m_rightchild->setParent(left->Parent());
						left->Parent()->setRightPtr(m_rightchild);
					}

					if (left->Parent()->m_leftchild==left)
					{
						m_rightchild->setParent(left->Parent());
						left->Parent()->setLeftPtr(m_rightchild);
					}

					m_leftsinter=m_leftchild->m_leftsinter;
					m_leftsinterdiam=m_leftchild->m_leftsinterdiam;
					m_rightsinter=m_leftchild->m_rightsinter;
					m_rightsinterdiam=m_leftchild->m_rightsinterdiam;
					m_sumvol_sinterpart=m_leftchild->m_sumvol_sinterpart;
					m_real_surface=m_leftchild->m_real_surface;
					m_real_surface_init=m_leftchild->m_real_surface_init;
					m_sph_surfacearea=m_leftchild->m_sph_surfacearea;
					vol_sinter=m_leftchild->vol_sinter;
					m_sinter_level=m_leftchild->m_sinter_level;
					dV_left=m_leftchild->dV_left;
					dV_right=m_leftchild->dV_right;
					m_sintered=m_leftchild->m_sintered;
					m_connect_time=m_leftchild->m_connect_time;
					m_rightchild=m_leftchild->m_rightchild;
					m_rightchild->m_parent=this;
					ParticleCache::operator =(*m_leftchild);
					temp=m_leftchild;
					m_leftchild = m_leftchild->m_leftchild;
					m_leftchild->m_parent=this;
					UpdateCache_sinter(left,right);
					UpdateTree_sinter(left,right);
				//	todelete=left;
				}
			
				// we have to change the left and right sinter particle 
				// because the subtree on both sides has changed
				// m_leftsinter=m_leftchild->SelectLeaf(model,id);
				// m_rightsinter=m_rightchild->SelectLeaf(model,id);
				
				//UpdateCache_sinter(todelete);
				//UpdateTree_sinter(todelete);

				(*temp).m_primary    = NULL;
				(*temp).m_leftchild  = NULL;
				(*temp).m_parent  = NULL;
				(*temp).m_pmodel = NULL;
				(*temp).m_rightchild = NULL;
				(*temp).m_leftsinter = NULL;
				(*temp).m_rightsinter = NULL;

				delete temp;

	//			UpdateCache_thispart();
				
	//			UpdateCache();
	//			UpdateTree();
				
			

			//out.open("aftersinter.inp");
			//printSubtree(out,id);
			//out.close();
		
		}
		//if (!FileExists("aftersinter.inp"))
	//	{
	/*		ofstream out;
			out.open("aftersinter.inp");
			this->FindRoot()->printSubtree(out,idsubtreeprint);
			out.close();*/
	//	}
	}

}


// Sinters the sub-particle for the given time using the given
// sintering model.
void SubParticle::Sinter(real dt, const Cell &sys, 
                         const Processes::SinteringModel &model)
{ 	
    if (m_primary != NULL) 
	    {
        m_primary->Sinter(dt, sys, model);
	    } 

	else {	 m_connect_time+=dt;
			// Perform a first order integration method to sinter
			// the subparticle for the given time.
		    
			// Declare time step variables.
			real t1=0.0, t2=0.0, tstop=dt;

			// Define the maximum allowed change in surface
			// area in one internal time step (10% spherical surface).
			real dAmax = 0.1 * SphSurfaceArea();

			// The scale parameter discretises the delta-S when using
			// the Poisson distribution.  This allows a smoother change
			// (smaller scale = higher precision).
			real scale = 0.01;


			m_sumsinterdiameter = 2*min(pow((m_leftsinter->vol_sinter-dV_left)*3/(4*PI),ONE_THIRD),
								  pow((m_rightsinter->vol_sinter-dV_right)*3/(4*PI),ONE_THIRD));	

		//	m_sumsinterdiameter = 2*min(pow((m_leftsinter->vol_sinter)*3/(4*PI),ONE_THIRD),
		//						  pow((m_rightsinter->vol_sinter)*3/(4*PI),ONE_THIRD));
			// Perform integration loop.
				while (t1 < tstop) {
					// Calculate sintering rate.
					real r = model.Rate(m_time+t1, sys, *this);
					// Calculate next time-step end point so that the
					// surface area changes by no more than dAmax.
					t2 = min(t1+(dAmax/max(r,1.0e-300)), tstop); // 1.0e-300 catches DIV ZERO.
					// Approximate sintering by a poisson process.  Calculate
					// number of poisson events.
					int n = ignpoi(r * (t2 - t1) / (scale*dAmax));
					// Adjust the surface area.
					if (n<0) break;
					if (n > 0) {
						m_real_surface -= (real)n * scale * dAmax;
						// Check that primary is not completely sintered.
						if (m_real_surface <= m_sph_surfacearea) {
							m_real_surface = m_sph_surfacearea ;
							m_sintered=1;
						}
					}
					// Set t1 for next time step.
					t1 = t2;
					//Update the volume that is transported between the two particles
					//Reduce the volume by the old sintering level

					m_leftsinter->vol_sinter-=dV_left;
					m_rightsinter->vol_sinter-=dV_right;
					//recalculate the sintering level
					m_sinter_level=m_sph_surfacearea*((1./m_real_surface)-1./m_real_surface_init)/(1-(m_sph_surfacearea/m_real_surface_init));
					//Add the Volume again to the particles according to the new sintering level
					dV_left=m_sinter_level*m_rightsinter->m_primary->Volume();
					dV_right=m_sinter_level*m_leftsinter->m_primary->Volume();
					m_leftsinter->vol_sinter+=dV_left;
					m_rightsinter->vol_sinter+=dV_right;
					//m_leftsinter->surf_sinter=4*PI*pow((3/(4*PI)*m_leftsinter->vol_sinter),2/3);
					//m_rightsinter->surf_sinter=4*PI*pow((3/(4*PI)*m_rightsinter->vol_sinter),2/3);
					if (m_sintered==1) break;

				}

			 if ( m_sinter_level>0.95 || Sintered()==1 || m_leftsinter->m_diam<1e-9 ||  m_rightsinter->m_diam<1e-9 )
			   {	//cout <<"sinter dt="<<dt<<endl;
			   	   	SinterPart();
				    UpdateCache();
					UpdateTree();
					if (m_primary==NULL)
					{
						m_leftchild->Sinter(dt,sys,model);
						m_rightchild->Sinter(dt,sys,model);
					}
			   }
			 else
			 {
				 m_leftchild->Sinter(dt, sys, model);
				 m_rightchild->Sinter(dt, sys, model);
			 }




	}


}
// PARTICLE UPDATE AND CHECKING.

// Recalculates the derived properties from the unique properties.
// This function moves down the tree from the top to the bottom.
void SubParticle::UpdateCache(void)
{	//ParticleCache::PropID id=ParticleCache::iS;
    if (m_primary != NULL) {
        // Get cache from primary particle.
        m_primary->UpdateCache();
        ParticleCache::operator=(*m_primary);
		m_avg_diam=m_primary->CollDiameter();
		m_numsubpart=1;
		ParticleCache::SetCollDiameter((6*m_vol/m_surf)*pow(NumSubPart(),1/1.9));
		

		
    } else {
        // The cache is the sum of the left and right child caches.
        m_leftchild->UpdateCache();
        m_rightchild->UpdateCache();
		m_avg_diam=m_leftchild->m_avg_diam+m_rightchild->m_avg_diam;
        ParticleCache::operator=(*m_leftchild);
        ParticleCache::operator+=(*m_rightchild);
		ParticleCache::SetCollDiameter((6*m_vol/m_surf)*pow(NumSubPart(),1/1.9));
	
    }
}

// Recalculates the derived properties from the unique properties.
// This function calculates only the cache of the current particle 
// assuiming the caches of the childrens do not need to be updated.
void SubParticle::UpdateCache_thispart(void)
{	//ParticleCache::PropID id=ParticleCache::iS;
    if (m_primary != NULL) {
        // Get cache from primary particle.
        m_primary->UpdateCache();
        ParticleCache::operator=(*m_primary);
		m_avg_diam=m_primary->CollDiameter();
		ParticleCache::SetCollDiameter((6*m_vol/m_surf)*pow(NumSubPart(),1/1.9));
		
		
    } else {
        // The cache is the sum of the left and right child caches.
        ParticleCache::operator=(*m_leftchild);
		m_avg_diam=m_leftchild->m_avg_diam+m_rightchild->m_avg_diam;
        ParticleCache::operator+=(*m_rightchild);		
		ParticleCache::SetCollDiameter((6*m_vol/m_surf)*pow(NumSubPart(),1/1.9));
    }
}



//Checks if the two particles that sintered together are not involved in other sintering processes.
/*void SubParticle::UpdateCache_sinter(SubParticle *has_sintered)
{	//ParticleCache::PropID id=ParticleCache::iS;
	bool haschanged=false;
    if (m_primary != NULL) {
		m_leftsinter=NULL;
		m_rightsinter=NULL;			
    } else {
		SubModels::SubModelType model = SubModels::BasicModel_ID;
		ParticleCache::PropID id=ParticleCache::iUniform;
		
		if(m_rightsinter==has_sintered ){
			m_rightsinter=m_rightchild->SelectLeaf(model,id);
			m_sintered=0;
			haschanged=true;
		}
		if(m_leftsinter==has_sintered ){
			m_leftsinter=m_leftchild->SelectLeaf(model,id);
			m_sintered=0;
			haschanged=true;
		}
		if(haschanged){
			m_real_surface=m_leftsinter->m_primary->SurfaceArea()+m_rightsinter->m_primary->SurfaceArea();
		    m_sumvol_sinterpart = m_leftsinter->m_primary->Volume()+m_rightsinter->m_primary->Volume();
		}
        m_leftchild->UpdateCache_sinter(has_sintered);
        m_rightchild->UpdateCache_sinter(has_sintered);
		
    }

	if(m_leftsinter!=NULL)
	if(m_leftsinter->m_primary->Volume()==0) cout <<"test";    // ms785 test
	if(m_rightsinter!=NULL)
	if(m_rightsinter->m_primary->Volume()==0) cout <<"test";    // ms785 test

}
*/
void SubParticle::UpdateCache_sinter(SubParticle *has_sintered, SubParticle *newsinter)
{	//ParticleCache::PropID id=ParticleCache::iS;
    if (m_primary != NULL) {
		m_leftsinter=NULL;
		m_rightsinter=NULL;
			
    } else {
		if(m_rightsinter==has_sintered){
			m_leftsinter->m_freesurface=(m_leftsinter->m_freesurface)+(has_sintered->m_diam)*(has_sintered->m_diam)*2*PI-(newsinter->m_diam)*(newsinter->m_diam)*2*PI;
			if (m_rightsinter->m_freesurface<0) m_rightsinter->m_freesurface=0;
			m_rightsinter=newsinter;
		}
		if(m_leftsinter==has_sintered){
			m_rightsinter->m_freesurface=(m_rightsinter->m_freesurface)+(has_sintered->m_diam)*(has_sintered->m_diam)*2*PI-(newsinter->m_diam)*(newsinter->m_diam)*2*PI;
			if (m_leftsinter->m_freesurface<0) m_leftsinter->m_freesurface=0;
			m_leftsinter=newsinter;
		}

        m_leftchild->UpdateCache_sinter(has_sintered,newsinter);
        m_rightchild->UpdateCache_sinter(has_sintered,newsinter);
		
    }



}


void SubParticle::CheckTree()
{
	if(m_primary==NULL)
	{	
		if(m_leftsinter!=NULL)
		if(m_leftsinter->m_primary->Volume()==0) cout <<"test";    
		if(m_rightsinter!=NULL)
		if(m_rightsinter->m_primary->Volume()==0) cout <<"test";   
		m_leftchild->CheckTree();
		m_rightchild->CheckTree();
	}
}


void SubParticle::Getprimarydistribution(double *distribution)
{
	if(m_primary!=NULL)
	{	
		int intdiam=(int)(CollDiameter()*1e9);
		distribution[intdiam]++;
	}
	else
	{
	    m_leftchild->Getprimarydistribution(distribution);
		m_rightchild->Getprimarydistribution(distribution);
	}
}

/*
void SubParticle::UpdateSinterParticles()
{
    if (m_primary != NULL) {
		m_leftsinter=NULL;
		m_rightsinter=NULL;			
    } else {
		SubModels::SubModelType model = SubModels::BasicModel_ID;
		ParticleCache::PropID id=ParticleCache::iUniform;
		
			m_rightsinter=m_rightchild->SelectLeaf(model,id);
			m_sintered=0;
			m_leftsinter=m_leftchild->SelectLeaf(model,id);
			m_real_surface=m_leftsinter->m_primary->SurfaceArea()+m_rightsinter->m_primary->SurfaceArea();
		    m_sumvol_sinterpart = m_leftsinter->m_primary->Volume()+m_rightsinter->m_primary->Volume();
			Primary::PropID primid=Primary::iD;
			m_leftsinterdiam = m_leftsinter->m_primary->Property(primid);
			m_rightsinterdiam = m_rightsinter->m_primary->Property(primid);
		
    }
}
*/
// This function is used to find the two sinter particles when the entire particle is duplicated
void SubParticle::UpdatethisSinterParticle(SubParticle *target, const SubParticle *original)
{	
    if (original->m_primary != NULL) {
		target->m_leftsinter=NULL;
		target->m_rightsinter=NULL;			
    } else {	
		int depth=0;
		int i;
		bool path[1000];
		SubParticle *orgsinter;
		orgsinter=original->m_leftsinter;
		while (orgsinter!=original)
		{
			if (orgsinter->Parent()->m_leftchild==orgsinter)
				path[depth]=true;
			else
				path[depth]=false;
			orgsinter=orgsinter->m_parent;
			depth++;
		}
		target->m_leftsinter=m_leftchild;
		for (i=depth-2;i>-1;i--)
		{
			if (path[i]==true)
				target->m_leftsinter=target->m_leftsinter->m_leftchild;
			else 
				target->m_leftsinter=target->m_leftsinter->m_rightchild;
		}
		

		depth=0;
		orgsinter=original->m_rightsinter;
		while (orgsinter!=original)
		{
			if (orgsinter->Parent()->m_leftchild==orgsinter)
				path[depth]=true;
			else
				path[depth]=false;
			depth++;
			orgsinter=orgsinter->m_parent;
		}
		target->m_rightsinter=m_rightchild;
		for (i=depth-2;i>-1;i--)
		{
			if (path[i]==true)
				target->m_rightsinter=target->m_rightsinter->m_leftchild;
			else 
				target->m_rightsinter=target->m_rightsinter->m_rightchild;
		}
	}
}




SubParticle *SubParticle::FindRoot()
{
	if (m_parent==NULL) return this;
	else return m_parent->FindRoot();
}




void SubParticle::UpdateFreeSurface()
{
		ResetFreeSurface();
		RecalcFreeSurface();
		UpdateCache();

}

void SubParticle::ResetFreeSurface()
{
	Primary::PropID id=Primary::iS;
	if (m_primary!=NULL)
		m_freesurface=m_primary->Property(id);
	else
	{
		m_freesurface=0;
		m_leftchild->ResetFreeSurface();
		m_rightchild->ResetFreeSurface();
	}
}

void SubParticle::RecalcFreeSurface()
{
	if(m_primary==NULL)
	{
		real Dfreesurface=min((m_leftsinter->m_diam)*(m_leftsinter->m_diam)*0.25*PI,(m_rightsinter->m_diam)*(m_rightsinter->m_diam)*0.25*PI);
		m_rightsinter->m_freesurface-=Dfreesurface;
		if (m_rightsinter->m_freesurface<0) m_rightsinter->m_freesurface=0;
		m_leftsinter->m_freesurface-=Dfreesurface;
		if ( m_leftsinter->m_freesurface<0)	m_leftsinter->m_freesurface=0;
		m_leftchild->RecalcFreeSurface();
		m_rightchild->RecalcFreeSurface();
	}
}

int SubParticle::NumSubPart()
{	return m_numsubpart;
	/*
	if (m_primary!=NULL)
		return 1;
	else
		return m_leftchild->NumSubPart()+m_rightchild->NumSubPart();*/
}

// Tells the parent sub-particle to update its cache of
// derived properties.  This operation is passed up the
// sub-particle tree to the root particle.
void SubParticle::UpdateTree(void)
{   //ParticleCache::PropID id=ParticleCache::iS;
    // Update the cache for this sub-particle.
    if (m_primary != NULL) {
        ParticleCache::operator=(*m_primary); 
		ParticleCache::SetCollDiameter((6*m_vol/m_surf)*pow(NumSubPart(),1/1.9));
		m_avg_diam=m_primary->CollDiameter();
		
    } else {
        ParticleCache::operator=(*m_leftchild);
        ParticleCache::operator+=(*m_rightchild);
	    ParticleCache::SetCollDiameter((6*m_vol/m_surf)*pow(NumSubPart(),1/1.9));
		m_avg_diam=m_leftchild->m_avg_diam+m_rightchild->m_avg_diam;
		
    }

    // Update the tree about this sub-particle.
    if (m_parent != NULL) {
        m_parent->UpdateTree();
    }
/*	//Update the entire tree
	else 
	{
		UpdateCache();
	}*/
}


//Checks that the two sintering children are still existing
/*void SubParticle::UpdateTree_sinter(SubParticle *has_sintered)
{  
	bool haschanged=false;
   
    if (m_primary != NULL) {
		m_leftsinter=NULL;
		m_rightsinter=NULL;
		
    } else {
		SubModels::SubModelType model = SubModels::BasicModel_ID;
		ParticleCache::PropID id=ParticleCache::iUniform;
		if(m_rightsinter==has_sintered){
			m_rightsinter=m_rightchild->SelectLeaf(model,id);
			haschanged=true;
			m_sintered=0;

		}
		if(m_leftsinter==has_sintered){
			m_leftsinter=m_leftchild->SelectLeaf(model,id);
			haschanged=true;
		}
		
    }

	if(haschanged)
	{
			m_real_surface=m_leftsinter->m_primary->SurfaceArea()+m_rightsinter->m_primary->SurfaceArea();
		    m_sumvol_sinterpart = m_leftsinter->m_primary->Volume()+m_rightsinter->m_primary->Volume();
	}

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->UpdateTree_sinter(has_sintered);
    }

}
*/
void SubParticle::UpdateTree_sinter(SubParticle *has_sintered,SubParticle *newsinter)
{    
    if (m_primary != NULL) {
		m_leftsinter=NULL;
		m_rightsinter=NULL;
		
    } else {
		if(m_rightsinter==has_sintered){
			m_rightsinter=newsinter;

		}
		if(m_leftsinter==has_sintered){
			m_leftsinter=newsinter;

		}
		

    }

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->UpdateTree_sinter(has_sintered,newsinter);
    }
	//Update the entire tree
/*	else 
	{
		UpdateCache_sinter(has_sintered,newsinter);
	}*/
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


void SubParticle::printSubtree(std::ostream &out, ParticleCache::PropID id) const
{ out << "digraph unix {"<<endl;
  printSubtreeLoop(out,id);
  out << "}"<<endl;
}

void SubParticle::printSubtreeLoop(std::ostream &out, ParticleCache::PropID id) const
{ 
  if (m_primary==NULL)
  { //out<<"leftchild "<<10E8*m_leftchild->SphDiameter()<<endl;   
	//out<<"rightchild "<<10E8*m_rightchild->SphDiameter()<<endl;    
	out<<"\" "<<this<<"\" "<<" [label = \""<<this->m_sinter_level<<"\"];"<<endl;
	out<<"\" "<<this->m_leftchild<<"\" "<<" [label = \""<<this->m_leftchild->m_sinter_level<<"\"];"<<endl;
	out<<"\" "<<this->m_rightchild<<"\" "<<" [label = \""<<this->m_rightchild->m_sinter_level<<"\"];"<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftchild<<"\"; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightchild<<"\"; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftsinter<<"\"[label=\""<<this->dV_left<<"\",color=\"blue\"]; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightsinter<<"\"[label=\""<<this->dV_right<<"\",color=\"blue\"]; "<<endl;
	m_leftchild->printSubtreeLoop(out,id);
    m_rightchild->printSubtreeLoop(out,id);
  }

  if (m_primary!=NULL)
  {   Primary::PropID idprim=Primary::iD;
      //out<<"\" "<<this<<"\" "<<" [label = \""<<2*pow((this->vol_sinter)*3/(4*PI),ONE_THIRD)<<"\"];"<<endl;
      out<<"\" "<<this<<"\" "<<" [label = \""<<this->m_sinter_level<<"\"];"<<endl;
	  out<<"\" "<<this->m_primary<<"\" "<<" [label = \""<<this->m_primary->Property(idprim)<<"\"];"<<endl;
	  out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_primary<<"\"; "<<endl;  
  }
}



void SubParticle::printSubtreepic(std::ostream &out) const
{ 
  printSubtreepicLoop(out,0,0,0);
}


void SubParticle::printSubtreepicLoop(std::ostream &out,real x, real y, real z) const
{ 
  real x1,z1,y1,random;
  x1=x;
  y1=y;
  z1=z;
  if (m_primary!=NULL)
  { random=rnd();
    out<<(m_primary->SphDiameter())*10E8<<endl;
	if(random<=0.33)
	{   
		 out<<1<<endl;
	}
	if(random>0.33 && random<=0.66)
	{	
		 out<<2<<endl;
	}
	if (random>0.66)
	{
		 out<<3<<endl;
	}
	//  out<<x*10E8<<"   "<<y*10E8<<"   "<<z*10E8<<endl;
	 
  }
  else
  { //out<<"leftchild "<<10E8*m_leftchild->SphDiameter()<<endl;   
	//out<<"rightchild "<<10E8*m_rightchild->SphDiameter()<<endl;
    //out<<"this "<<10E8*SphDiameter()<<endl;
	/*if(random<=0.33)
	{   
		x1=x-0.5*(SphDiameter()-m_leftchild->SphDiameter());
	}
	if(random>0.33 && random<=0.66)
	{	
		y1=y-0.5*(SphDiameter()-m_leftchild->SphDiameter());
	}
	if (random>0.66)
	{
		z1=z-0.5*(SphDiameter()-m_leftchild->SphDiameter());
	}*/
	m_leftchild->printSubtreepicLoop(out, x1, y1, z1);

    x1=x;
    y1=y;
    z1=z;

	/*if(random<=0.33)
	{
		x1=x-0.5*(SphDiameter())+m_leftchild->SphDiameter()+0.5*(m_rightchild->SphDiameter());
	}
	if(random>0.33 && random<=0.66)
	{
		y1=y-0.5*(SphDiameter())+m_leftchild->SphDiameter()+0.5*(m_rightchild->SphDiameter());
	}
	if (random>0.66)
	{
		z1=z-0.5*(SphDiameter())+m_leftchild->SphDiameter()+0.5*(m_rightchild->SphDiameter());
	}*/

    m_rightchild->printSubtreepicLoop(out, x1, y1, z1);
  }
  
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
