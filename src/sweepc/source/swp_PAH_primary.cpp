/*!
 * \file   swp_PAH_primary.cpp
 * \author Markus Sander
 *
 * \brief  Particle Model that stores the individual PAHs of a soot particle
 */
/*
  Author(s):      Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Markus Sander.

  File purpose:
    Implementation of the PAHPrimary class declared in the
    swp_PAH_primary.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#define _USE_MATH_DEFINES //!< First define.
#include <math.h>         //!< Then include so that the pi constant (M_PI) can be used.

#include "swp_primary.h"
#include "swp_PAH_primary.h"
#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"
#include "swp_particle_image.h"
#include "swp_cell.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_pah_structure.h"
#include "swp_PAH.h"
#include "swp_ensemble.h"
#include "swp_particle_model.h"

#include <stdexcept>
#include <cassert>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>

#include "string_functions.h"

int uniquePAHCounter = 0;

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace Sweep::KMC_ARS;

using namespace std;
using namespace Strings;

//used for debugging, testing clone function for PAHStructure.
static unsigned int ID=0; 
static bool m_clone=false;
/*
double PAHPrimary::pow(double a, double b) {
    int tmp = (*(1 + (int *)&a));
    int tmp2 = (int)(b * (tmp - 1072632447) + 1072632447);
    double p = 0.0;
    *(1 + (int * )&p) = tmp2;
    return p;
}
*/

// CONSTRUCTORS AND DESTRUCTORS.
PAHPrimary::PAHPrimary() : Primary(),
    m_numcarbon(0),
    m_numH(0),
    m_numOfEdgeC(0),
    m_numOfRings(0),
	m_numOfRings5(0),
	m_PAHmass(0),
	m_PAHCollDiameter(0),
    m_numPAH(0),
    m_numprimary(0),
    m_primarydiam(0.0),
	m_primaryvol(0.0),
    m_children_radius(0),
    m_children_vol(0),
	m_children_surf(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_roundingLevel(0),
	m_children_sintering(0.0),
	m_children_sumCap(0.0), //used to calculate geometric volume of a particle
	m_sum_cap(0.0), //used to calculate geometric volume of a particle
	m_sph_prim_vol(0.0), //used to calculate geometric volume of a particle
    m_distance_centreToCentre(0.0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0), //if primary coordinates are not tracked
	m_avg_sinter(0.0), //if primary coordinates are tracked
    m_sint_time(0.0),
	m_sint_rate(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL),
	m_free_surf(0.0),
	m_sum_necks(0.0),
	m_r(0.0),
	m_r2(0.0),
	m_r3(0.0)
{
	m_cen_bsph[0] = 0.0;
	m_cen_bsph[1] = 0.0;
	m_cen_bsph[2] = 0.0;

	m_cen_mass[0] = 0.0;
	m_cen_mass[1] = 0.0;
	m_cen_mass[2] = 0.0;
}

/*!
 * @param[in]       time        Time at which particle is being created
 * @param[in]       model       Model which defines the meaning of the primary
 *
 */
PAHPrimary::PAHPrimary(const double time, const Sweep::ParticleModel &model)
	: Primary(time, model),
	m_numcarbon(0),
	m_numH(0),
	m_numOfEdgeC(0),
	m_numOfRings(0),
	m_numOfRings5(0),
	m_PAHmass(0),
	m_PAHCollDiameter(0),
	m_numPAH(0),
	m_numprimary(0),
	m_primarydiam(0.0),
	m_primaryvol(0.0),
	m_children_radius(0),
	m_children_vol(0),
	m_children_surf(0),
	m_leftparticle_vol_old(0),
	m_rightparticle_vol_old(0),
	m_rightparticle_numPAH(0),
	m_leftparticle_numPAH(0),
	m_children_roundingLevel(0),
	m_children_sintering(0.0),
	m_children_sumCap(0.0),
	m_sum_cap(0.0),
	m_sph_prim_vol(0.0),
	m_distance_centreToCentre(0.0),
	m_Rg(0),
	m_fdim(0),
	m_sqrtLW(0),
	m_LdivW(0),
	m_avg_coalesc(0),
	m_avg_sinter(0.0),
	m_sint_time(0.0),
	m_sint_rate(0.0),
	m_leftchild(NULL),
	m_rightchild(NULL),
	m_parent(NULL),
	m_leftparticle(NULL),
	m_rightparticle(NULL),
	m_free_surf(0.0),
	m_sum_necks(0.0),
	m_r(0.0),
	m_r2(0.0),
	m_r3(0.0)
{
	m_cen_bsph[0] = 0.0;
	m_cen_bsph[1] = 0.0;
	m_cen_bsph[2] = 0.0;

	m_cen_mass[0] = 0.0;
	m_cen_mass[1] = 0.0;
	m_cen_mass[2] = 0.0;

	// Other parts of the code check for a non-zero composition
	m_comp[0] = 1;

	AddPAH(time, model);

	//Update the other properties
	UpdateCache();
}

/*!
 * @param[in]       time        Time at which particle is being created
 * @param[in]       position    Position at which particle is being created
 * @param[in]       model       Model which defines the meaning of the primary
 *
 */
PAHPrimary::PAHPrimary(const double time, const double position,
	const Sweep::ParticleModel &model)
	: Primary(time, model),
	m_numcarbon(0),
	m_numH(0),
	m_numOfEdgeC(0),
	m_numOfRings(0),
	m_numOfRings5(0),
	m_PAHmass(0),
	m_PAHCollDiameter(0),
	m_numPAH(0),
	m_numprimary(0),
	m_primarydiam(0.0),
	m_primaryvol(0.0),
	m_children_radius(0),
	m_children_vol(0),
	m_children_surf(0),
	m_leftparticle_vol_old(0),
	m_rightparticle_vol_old(0),
	m_rightparticle_numPAH(0),
	m_leftparticle_numPAH(0),
	m_children_roundingLevel(0),
	m_children_sintering(0.0),
	m_children_sumCap(0.0),
	m_sum_cap(0.0), 
	m_sph_prim_vol(0.0), 
	m_distance_centreToCentre(0.0),
	m_Rg(0),
	m_fdim(0),
	m_sqrtLW(0),
	m_LdivW(0),
	m_avg_coalesc(0),
	m_avg_sinter(0.0),
	m_sint_time(0.0),
	m_sint_rate(0.0), 
	m_leftchild(NULL),
	m_rightchild(NULL),
	m_parent(NULL),
	m_leftparticle(NULL),
	m_rightparticle(NULL),
	m_free_surf(0.0),
	m_sum_necks(0.0), 
	m_r(0.0),
	m_r2(0.0),
	m_r3(0.0)
{
	m_cen_bsph[0] = 0.0; 
	m_cen_bsph[1] = 0.0; 
	m_cen_bsph[2] = 0.0; 

	m_cen_mass[0] = 0.0; 
	m_cen_mass[1] = 0.0; 
	m_cen_mass[2] = 0.0; 
	// Other parts of the code check for a non-zero composition
	m_comp[0] = 1;

	AddPAH(time, model);

	//Update the other properties
	UpdateCache();
}


// Initialising constructor.
PAHPrimary::PAHPrimary(double time, const Sweep::ParticleModel &model, bool noPAH)
	: Primary(time, model),
	m_numcarbon(0),
	m_numH(0),
	m_numOfEdgeC(0),
	m_numOfRings(0),
	m_numOfRings5(0),
	m_PAHmass(0),
	m_PAHCollDiameter(0),
	m_numPAH(0),
	m_numprimary(0),
	m_primarydiam(0.0),
	m_primaryvol(0.0),
	m_children_radius(0),
	m_children_vol(0),
	m_children_surf(0),
	m_leftparticle_vol_old(0),
	m_rightparticle_vol_old(0),
	m_rightparticle_numPAH(0),
	m_leftparticle_numPAH(0),
	m_children_roundingLevel(0),
	m_children_sintering(0.0),
	m_children_sumCap(0.0), 
	m_sum_cap(0.0),
	m_sph_prim_vol(0.0),
	m_distance_centreToCentre(0.0),
	m_Rg(0),
	m_fdim(0),
	m_sqrtLW(0),
	m_LdivW(0),
	m_avg_coalesc(0),
	m_avg_sinter(0.0),
	m_sint_time(0.0),
	m_sint_rate(0.0), 
	m_leftchild(NULL),
	m_rightchild(NULL),
	m_parent(NULL),
	m_leftparticle(NULL),
	m_rightparticle(NULL),
	m_free_surf(0.0), 
	m_sum_necks(0.0), 
	m_r(0.0), 
	m_r2(0.0), 
	m_r3(0.0)
{
	m_cen_bsph[0] = 0.0; 
	m_cen_bsph[1] = 0.0; 
	m_cen_bsph[2] = 0.0; 

	m_cen_mass[0] = 0.0; 
	m_cen_mass[1] = 0.0; 
	m_cen_mass[2] = 0.0; 

	m_comp[0] = 1;
}


/*double PAHPrimary::Fdim() const
{
    return m_fdim;
}*/


/*!
 * Add a PAH to the primary particle
 *
 * @param[in]   time        create time of the PAH

 * @param[in]   model       Particle model containing molecule database
*/
void PAHPrimary::AddPAH(double time,const Sweep::ParticleModel &model)
{
    boost::shared_ptr<PAH> new_PAH (new PAH(time, model.InceptedPAH()));
    new_PAH->PAH_ID=ID;
    m_PAH.push_back(new_PAH);
    ID++;
    // Set the particle mass, diameter etc
    UpdatePrimary();
}

// Copy constructor.
PAHPrimary::PAHPrimary(const PAHPrimary &copy)
{
    *this = copy;
    if (copy.m_leftchild!=NULL)
    {
        CopyTree(&copy);
    }
	//m_clone=false;
}


// Default destructor.
PAHPrimary::~PAHPrimary()
{
    delete m_leftchild;
    delete m_rightchild;
    // it is not necessary to delete m_leftparticle because
    // this is also m_leftchild somewhere down the tree
    if (m_PAH.size()!=0) m_PAH.clear();
    // delete the PAH list
    releaseMem();
    m_clone=false;
}


// Copy constructor.
std::vector<boost::shared_ptr<PAH> > PAHPrimary::GetPAHVector() const
{
	return m_PAH;
}

/*!
 * Recursively copy the tree for non-leaf nodes.
 *
 * @param[in] source Pointer to the primary to be copied
*/
void PAHPrimary::CopyTree( const PAHPrimary *source)
{
    //create the new left and right children with nothing in them
    m_leftchild = new PAHPrimary(source->CreateTime(),*m_pmodel,false);
    m_rightchild = new PAHPrimary(source->CreateTime(),*m_pmodel,false);

    // copy the properties such as the volume, surface area and list of
	// constituent PAH molecules
	m_leftchild->CopyParts(source->m_leftchild);
	m_rightchild->CopyParts(source->m_rightchild);

    //set the pointers to the parent
	m_leftchild->m_parent=this;
	m_rightchild->m_parent=this;

    // the left and right particle are set further down in UpdateAllPointers
	// These are the pointers that specify which primary particles touch each
	// other in the aggregat structure.
	m_leftparticle=NULL;
	m_rightparticle=NULL;

    //recurse to copy the subtrees
	if (source->m_leftchild->m_leftchild!=NULL)
        m_leftchild->CopyTree(source->m_leftchild);

	if (source->m_rightchild->m_leftchild!=NULL)
        m_rightchild->CopyTree(source->m_rightchild);

    //set the leftparticle and rightparticle
	UpdateAllPointers(source);
}

PAHPrimary &PAHPrimary::operator=(const Primary &rhs)
{
    if (this != &rhs) {
        const AggModels::PAHPrimary *pahprimary = NULL;
        pahprimary = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);
        UpdateAllPointers(pahprimary);
    }
    return *this;
}

/*
 * @brief Stream-reading constructor
 *
 * @param[in,out]	 in		       Input binary stream  
 * @param[in]        model	       Particle model defining interpretation of particle data
 * @param[in,out]    duplicates    Addresses of PAHs for use when reading primary particles
 */ 
PAHPrimary::PAHPrimary(std::istream &in, const Sweep::ParticleModel &model, PahDeserialisationMap &pah_duplicates)
{
    Deserialize(in, model, pah_duplicates);
}


/*!
 * This function is like a limited assignment operator, except that the
 * children are not copied and the pointers to the particles may need
 * adjusting after this method has finished.
 *
 * @param[in] source Pointer to the primary to be copied
 */
void PAHPrimary::CopyParts(const PAHPrimary *source)
{
    SetCollDiameter(source->CollDiameter());
    SetMobDiameter(source->MobDiameter());
    SetSphDiameter(source->SphDiameter());

	SetSurfaceArea(source->SurfaceArea());
	SetMass(source->Mass());

	m_values = source->m_values;
	m_comp = source->m_comp;
	m_time = source->m_time;

	m_sint_time = source->m_sint_time;
	m_sint_rate = source->m_sint_rate;
	m_numcarbon = source->m_numcarbon;
	m_numH = source->m_numH;
	m_numOfEdgeC = source->m_numOfEdgeC;
	m_numOfRings = source->m_numOfRings;
	m_numOfRings5 = source->m_numOfRings5;
	m_numPAH = source->m_numPAH;
	m_numprimary = source->m_numprimary;
	m_primarydiam = source->m_primarydiam;
	m_primaryvol = source->m_primaryvol;
    m_PAHCollDiameter=source->m_PAHCollDiameter;
    m_PAHmass=source->m_PAHmass;
    
    m_surf=source->m_surf;
    m_vol=source->m_vol;
	m_free_surf = source->m_free_surf;
    m_distance_centreToCentre = source->m_distance_centreToCentre;
	
	m_cen_bsph = source->m_cen_bsph;
	m_cen_mass = source->m_cen_mass;
	m_r = source->m_r;
	m_r2 = source->m_r2;
	m_r3 = source->m_r3;
    m_children_vol=source->m_children_vol;
    m_children_radius=source->m_children_radius;
	m_children_surf = source->m_children_surf;
	m_children_sintering = source->m_children_sintering;
	m_children_roundingLevel = source->m_children_roundingLevel;
	m_children_sumCap = source->m_children_sumCap;
	m_sum_cap = source->m_sum_cap;
	m_sph_prim_vol = source->m_sph_prim_vol;
	m_sum_necks = source->m_sum_necks;

    m_rightparticle_numPAH=source->m_rightparticle_numPAH;
    m_leftparticle_numPAH=source->m_leftparticle_numPAH;
    m_leftparticle_vol_old=source->m_leftparticle_vol_old;
    m_rightparticle_vol_old=source->m_rightparticle_vol_old;

    m_fdim=source->m_fdim;
    m_Rg=source->m_Rg;
	m_sqrtLW = source->m_sqrtLW;
	m_LdivW = source->m_LdivW;
	m_pmodel = source->m_pmodel;
    m_avg_coalesc=source->m_avg_coalesc;
	m_avg_sinter = source->m_avg_sinter;

	//! Set particles.
	m_parent = source->m_parent;
	m_leftchild = source->m_leftchild;
	m_rightchild = source->m_rightchild;
	m_leftparticle = source->m_leftparticle; 
	m_rightparticle = source->m_rightparticle;

    //! Replace the PAHs with those from the source.
    if (m_clone==true) {
	    m_PAH.resize(source->m_PAH.size());
        m_PAH.clear();
        const int sharedPointers = m_pmodel->Components(0)->SharedPointers();
        for (size_t i=0; i!=source->m_PAH.size();++i) {
		    //! Each time a PAH is cloned 100000 is added to its ID so that we can easily calculate how many times this pah has been cloned.
            if(sharedPointers) {
                //! Create a copy of the shared pointers for PAHs in particles.
		        if (source->m_PAH.size() > 1) {
                    boost::shared_ptr<PAH> new_m_PAH = source->m_PAH[i];
                    new_m_PAH->PAH_ID=source->m_PAH[i]->PAH_ID+100000;
                    m_PAH.push_back(new_m_PAH);
                }
                //! Always create new shared pointers for single PAHs.
                else {
                    boost::shared_ptr<PAH> new_m_PAH (source->m_PAH[i]->Clone());
                    new_m_PAH->PAH_ID=source->m_PAH[i]->PAH_ID+100000;
                    m_PAH.push_back(new_m_PAH);
                }
            }
            else {
                //! Create new shared pointers for single PAHs or PAHs in particles.
                boost::shared_ptr<PAH> new_m_PAH (source->m_PAH[i]->Clone());
                //m_PAH[i]=source->m_PAH[i]->Clone();
                new_m_PAH->PAH_ID=source->m_PAH[i]->PAH_ID+100000;
                m_PAH.push_back(new_m_PAH);
            }
        }
    }
	else {
		m_PAH.assign(source->m_PAH.begin(), source->m_PAH.end());
	}
}
	

/*!
 * Select a primary uniformly at random from this particle
 * and descend the aggregate tree to find the primary.
 * Note that most PAHPrimarys are nodes in a tree representing
 * connectivity within an aggregate, so it is necessary to
 * descend the tree to find a primary that really is a
 * primary.
 *
 * @param[in,out]   rng     Random number generator
 *
 * @return      Pointer to an object representing a phsyical primary
 */
PAHPrimary *PAHPrimary::SelectRandomSubparticle(rng_type &rng)
{
    // We want to choose an integer uniformly on the range [0, m_numprimary - 1]
    typedef boost::uniform_smallint<int> uniform_integer;
    uniform_integer uniformDistrib(0, m_numprimary - 1);

    // Now build an object that can generate a sample and use it
    boost::variate_generator<rng_type&, uniform_integer> uniformGenerator(rng, uniformDistrib);
    return SelectRandomSubparticleLoop(uniformGenerator());

}
/*!
 * @param[in] target The primary to be selected
 *
 * @return      Pointer to the next node down the tree on the path to target primary
*/
PAHPrimary *PAHPrimary::SelectRandomSubparticleLoop(int target)
{
	if (m_leftchild==NULL) return this;
	if (target < m_leftchild->m_numprimary)
	{
		return m_leftchild->SelectRandomSubparticleLoop(target);
	}
	else
	{
		return m_rightchild->SelectRandomSubparticleLoop(target-(m_leftchild->m_numprimary));
	}
}


/*!
 * Each node contains two pointers (m_leftparticle and m_rightparticle)
 * to primary particles that are connected by this node
 * This function is used when the entire particle tree is duplicated.
 * It sets the pointers in the copied node (this), that the connectivity
 * of the primary particles in this node is the same as in the original node.
 *
 * A large number of asserts are provided to check that the pointer network
 * describing the aggregate structure is valid.  This can be re-enabled if
 * there is suspicion of a bug,
 *
 *
 *@todo give this method a more accurate name
 *
 * @param[in] original    Pointer to the primary to be copied
*/
void PAHPrimary::UpdateAllPointers( const PAHPrimary *original)
{
    //the primary has no children => there are no left and right particles
    if (original->m_leftchild == NULL) {
        // Check that both children are missing in the original and the copy
        //assert(original->m_rightchild == NULL);
        //assert(m_leftchild == NULL);
        //assert(m_rightchild == NULL);

        // Check the original really does not have left and right particles
        // since it should be representing a physical primary, not a connection
        // inside an aggregate unit.
        //assert(original->m_leftparticle == NULL);
        //assert(original->m_rightparticle == NULL);

        // Since this is not a connecting node it does not have left and
        // right particles.
		m_leftparticle=NULL;
		m_rightparticle=NULL;
    } else {
//        assert(original->m_rightchild != NULL);
//        assert(m_leftchild != NULL);
//        assert(m_rightchild != NULL);

        // Check the original really does have left and right particles
        // since it should be representing a connection
        // inside an aggregate unit.
//        assert(original->m_leftparticle != NULL);
//        assert(original->m_rightparticle != NULL);

        // m_leftparticle should be a leaf node so assert that it does not have children
//        assert(original->m_leftparticle->m_leftchild == NULL);
//        assert(original->m_leftparticle->m_rightchild == NULL);

        // Find the route to m_leftparticle in the original tree
        std::stack<bool> route = recordPath(original->m_leftparticle, original);

        // Now follow the same route down the new tree to find the new left particle
        m_leftparticle = descendPath(this, route);

        // the new m_leftparticle should be a leaf node so assert
        // that it does not have children
//        assert(m_leftparticle->m_leftchild == NULL);
//        assert(m_leftparticle->m_rightchild == NULL);

        // m_rightparticle should be a leaf node so assert that it does not have children
//        assert(original->m_rightparticle->m_leftchild == NULL);
//        assert(original->m_rightparticle->m_rightchild == NULL);

        // Find the route to m_rightparticle in the original tree
        route = recordPath(original->m_rightparticle, original);

        // Now follow the same route down the new tree to find the new right particle
        m_rightparticle = descendPath(this, route);

        // the new m_rightparticle should be a leaf node so assert
        // that it does not have children
//        assert(m_rightparticle->m_leftchild == NULL);
//        assert(m_rightparticle->m_rightchild == NULL);
	}
}
//void PAHPrimary::ReleasePAH(Primary &rhs)
//{
//    PAHPrimary *rhsparticle = NULL;
//    rhsparticle = dynamic_cast<AggModels::PAHPrimary*>(&rhs);
//    rhsparticle->m_PAH.clear();
//}


/*!
 * Combines this primary with another.
 *
 * \param[in]       rhs         Particle to add to current instance
 * \param[in,out]   rng         Random number generator
 *
 * \return      Reference to the current instance after rhs has been added
 */
PAHPrimary &PAHPrimary::Coagulate(const Primary &rhs, rng_type &rng)
{
	vector<PAH>::const_iterator j;
	const PAHPrimary *rhsparticle = NULL;
	rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);
	double vol_old = 0.0;
	////only one PAH in rhs or this particle -> condensation or inception process.
	if ((rhsparticle->m_numPAH == 1) || (m_numPAH == 1))
	{
		if (rhsparticle->Numprimary() > 1)
		{
			//rhsparticle is a soot particle but this paricle repesents a PAH
			//this paricle will condense on this rhsparticle
			//Get a copy of the rhs ready to add to the current particle
			PAHPrimary copy_rhs(*rhsparticle);
			PAHPrimary *target = copy_rhs.SelectRandomSubparticle(rng);
			vol_old = target->m_vol; //volume before condensation
			target->m_PAH.insert(target->m_PAH.end(), m_PAH.begin(), m_PAH.end());
			target->UpdatePrimary();
			if (m_pmodel->getTrackPrimaryCoordinates() || m_pmodel->getTrackPrimarySeparation())
				target -> Adjust(vol_old); //treat condensation similar to surface growth

			CopyParts(&copy_rhs);
			CopyTree(&copy_rhs);
		}
		else
		{
			//this paricle is a soot particle but rhsparticle repesents a PAH
			//rhsparticle will condense on this particle
			//particle has more then one primary select the primary where
			//the PAH condenses to and add it to the list
			if (m_leftchild != NULL)
			{
				PAHPrimary *target = SelectRandomSubparticle(rng);
				vol_old = target->m_vol; //volume before condensation
				target->m_PAH.insert(target->m_PAH.end(), rhsparticle->m_PAH.begin(), rhsparticle->m_PAH.end());
				target->UpdatePrimary();
				if (m_pmodel->getTrackPrimaryCoordinates() || m_pmodel->getTrackPrimarySeparation())
					target->Adjust(vol_old); //treat condensation similar to surface growth

			}
			else
			{
				//! rhsparticle is a gas-phase PAH but the this pointer may be
				//! pointing to a gas-phase PAH in which case this would be an
				//! inception event, or a single primary particle in which case
				//! it would be a condensation event.
				vol_old = m_vol; //volume before condensation
				m_PAH.insert(m_PAH.end(), rhsparticle->m_PAH.begin(), rhsparticle->m_PAH.end());
				UpdatePrimary();
				if (m_pmodel->getTrackPrimaryCoordinates() || m_pmodel->getTrackPrimarySeparation())
					Adjust(vol_old); //treat condensation similar to surface growth
			}
		}
		UpdateCache();

		//Check the coalescence ratio
		if (m_pmodel->getTrackPrimaryCoordinates() || m_pmodel->getTrackPrimarySeparation())
			CheckSintering();
		else
			CheckRounding();

	}

	else
	{ //comment inception and condensation to test PAH_KMC model and bintree model
		//coagulation process
		//PAHPrimary *newleft = new PAHPrimary(m_time, *m_pmodel);
		//PAHPrimary *newright = new PAHPrimary(m_time, *m_pmodel);
		PAHPrimary *newleft = new PAHPrimary;
		PAHPrimary *newright = new PAHPrimary;
		PAHPrimary copy_rhs(*rhsparticle);
		//rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);

		// bool print=false;
		// if (this->m_numprimary>2 && notconstpah->m_numprimary>2)
		// {
		//    PrintTree("before1");
		//     notconstpah->PrintTree("before2");
		//   print=true;
		//  cout << "printing tree"<<endl;
		//}

		//select where to add the second particle
		boost::bernoulli_distribution<> bernoulliDistrib;
		boost::variate_generator<rng_type &, boost::bernoulli_distribution<> > leftRightChooser(rng, bernoulliDistrib);
		if (leftRightChooser())
		{
			newleft->CopyParts(this);
			newright->CopyParts(&copy_rhs);
		}
		else
		{
			newright->CopyParts(this);
			newleft->CopyParts(&copy_rhs);
		}
		//set the pointers
		m_leftchild = newleft;
		m_rightchild = newright;
		newright->m_parent = this;
		newleft->m_parent = this;

		this->m_leftparticle = NULL;//leftparticle has not been determined yet
		this->m_rightparticle = NULL;//rightparticle has not been determined yet

		//set the pointers to the parent node
		if (newleft->m_leftchild != NULL)
		{
			newleft->m_leftchild->m_parent = newleft;
			newleft->m_rightchild->m_parent = newleft;
		}
		if (newright->m_leftchild != NULL)
		{
			newright->m_leftchild->m_parent = newright;
			newright->m_rightchild->m_parent = newright;
		}

		if (m_pmodel->getTrackPrimaryCoordinates() || m_pmodel->getTrackPrimarySeparation())
			m_children_sintering = 0; //track coordinates
		else
			m_children_roundingLevel = 0; //no track

		UpdateCache();

		//! It is assumed that primary pi from particle Pq and primary pj from
		//! particle Pq are in point contact and by default pi and pj are
		//! uniformly selected. If we track the coordinates of the primaries in
		//! a particle, we can do a smarter selection where pi and pj are
		//! determined by ballistic cluster-cluster aggregation (BCCA):
		//! R. Jullien, Transparency effects in cluster-cluster aggregation with
		//! linear trajectories, J. Phys. A 17 (1984) L771-L776.
		if (m_pmodel->getTrackPrimaryCoordinates() || m_pmodel->getTrackPrimarySeparation()) {
			boost::uniform_01<rng_type&, double> uniformGenerator(rng);

			//! Implementation of Arvo's algorithm, Fast Random Rotation Matrices,
			//! Chapter III.4 in Graphic Gems III edited by David Kirk to generate
			//! a transformation matrix for randomly rotating a particle.
			double theta1 = 2 * PI * uniformGenerator(); //!< Pick a rotation about the pole.
			double phi1 = 2 * PI * uniformGenerator();   //!< Pick a direction to deflect the pole.
			double z1 = uniformGenerator();              //!< Pick the amount of pole deflection.

			//! Construct a vector for performing the reflection.
			fvector V1;
			V1.push_back(cos(phi1) * sqrt(z1));
			V1.push_back(sin(phi1) * sqrt(z1));
			V1.push_back(sqrt(1 - z1));

			//! Rotate centre-of-mass.
			m_leftchild->rotateCOM(theta1, V1);

			//! Implementation of Arvo's algorithm, Fast Random Rotation Matrices,
			//! Chapter III.4 in Graphic Gems III edited by David Kirk to generate
			//! a transformation matrix for randomly rotating a particle.
			double theta2 = 2 * PI * uniformGenerator(); //!< Pick a rotation about the pole.
			double phi2 = 2 * PI * uniformGenerator();   //!< Pick a direction to deflect the pole.
			double z2 = uniformGenerator();              //!< Pick the amount of pole deflection.

			//! Construct a vector for performing the reflection.
			fvector V2;
			V2.push_back(cos(phi2) * sqrt(z2));
			V2.push_back(sin(phi2) * sqrt(z2));
			V2.push_back(sqrt(1 - z2));

			//! Rotate centre-of-mass.
			m_rightchild->rotateCOM(theta2, V2);

			m_leftchild->centreBoundSph();
			m_rightchild->centreBoundSph();

			bool Overlap = false;

			//! Incremental translation.
			while (!Overlap) {
				//! Sphere point picking. This is the random direction step of
				//! Jullien's BCCA algorithm. It is incorrect to select spherical
				//! coordinates theta (polar angle) and phi (azimuthal angle) from
				//! uniform distributions theta E [0, 2 * pi) and phi E [0, pi] as
				//! points picked in this way will be 'bunched' near the poles:
				//! http://mathworld.wolfram.com/SpherePointPicking.html
				double theta = 2.0 * PI * uniformGenerator();
				double phi = acos(2.0 * uniformGenerator() - 1.0);

				//! In terms of Cartesian coordinates.
				double x = cos(theta) * sin(phi);
				double y = sin(theta) * sin(phi);
				double z = cos(phi);

				//! We want to find the rotation matrix which rotates an
				//! arbitrarily chosen unit vector to the randomly chosen unit
				//! vector selected above. Not sure where the formula is from but I
				//! have checked that it works:
				//! https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
				double bx = 0.0;
				double by = 0.0;
				double bz = -1.0;

				//! Cross product of the two vectors.
				double vx = by * z - bz * y;
				double vy = bz * x - bx * z;
				double vz = bx * y - by * x;

				//! Skew-symmetric cross-product matrix of vector v.
				double v[3][3] = { 0 };

				v[0][0] = 0.0;
				v[0][1] = -vz;
				v[0][2] = vy;
				v[1][0] = vz;
				v[1][1] = 0.0;
				v[1][2] = -vx;
				v[2][0] = -vy;
				v[2][1] = vx;
				v[2][2] = 0.0;

				//! Identity.
				double I[3][3] = { 0 };

				I[0][0] = 1.0;
				I[0][1] = 0.0;
				I[0][2] = 0.0;
				I[1][0] = 0.0;
				I[1][1] = 1.0;
				I[1][2] = 0.0;
				I[2][0] = 0.0;
				I[2][1] = 0.0;
				I[2][2] = 1.0;

				//! Rotation matrix.
				double R[3][3] = { 0 };
				double Mult = 1.0 / (1.0 + bx * x + by * y + bz * z); //!< Dot product of the two unit vectors.

				//! Square of the matrix v.
				R[0][0] = Mult * (-vy * vy - vz * vz);
				R[0][1] = Mult * (vx * vy);
				R[0][2] = Mult * (vx * vz);
				R[1][0] = Mult * (vx * vy);
				R[1][1] = Mult * (-vx * vx - vz * vz);
				R[2][2] = Mult * (vy * vz);
				R[2][0] = Mult * (vx * vz);
				R[2][1] = Mult * (-vy * vz);
				R[2][2] = Mult * (-vx * vx - vy * vy);

				R[0][0] = 1.0 + R[0][0];
				R[0][1] = -vz + R[0][1];
				R[0][2] = vy + R[0][2];
				R[1][0] = vz + R[1][0];
				R[1][1] = 1.0 + R[1][1];
				R[1][2] = -vx + R[1][2];
				R[2][0] = -vy + R[2][0];
				R[2][1] = vx + R[2][1];
				R[2][2] = 1.0 + R[2][2];

				//! Disk point picking. This is the random impact step of Jullien's
				//! BCCA algorithm. We first randomly select a point in the x-y
				//! plane over a disk of diameter and at a distance equal to the
				//! sum of the primary particle radii as this is the maximum
				//! distance they can be apart. Then we apply the rotation matrix
				//! obtained above to the point:
				//http://mathworld.wolfram.com/DiskPointPicking.html
				double r = uniformGenerator();
				theta = 2.0 * PI * uniformGenerator();

				double sumr = m_leftchild->Radius() + m_rightchild->Radius();

				double x3 = (sumr / 2.0) * sqrt(r) * cos(theta);
				double y3 = (sumr / 2.0) * sqrt(r) * sin(theta);
				double z3 = -sumr;

				double x4 = R[0][0] * x3 + R[0][1] * y3 + R[0][2] * z3;
				double y4 = R[1][0] * x3 + R[1][1] * y3 + R[1][2] * z3;
				double z4 = R[2][0] * x3 + R[2][1] * y3 + R[2][2] * z3;

				this->m_leftchild->Translate(x4, y4, z4);

				//! The two particles are initially at the origin. Keep doubling
				//! the distance between them so that they do not overlap. This is
				//! more efficient than picking an arbitrarily large distance.
				int numberOfOverlaps = 0;
				int factorApart = 1;
				double Separation = 0.0;

				while (this->checkForOverlap(*m_leftchild, *m_rightchild, numberOfOverlaps, Separation)) {
					this->m_leftchild->Translate(-factorApart * R[0][2] * sumr, -factorApart * R[1][2] * sumr, -factorApart * R[2][2] * sumr);
					factorApart *= 2;
				}

				numberOfOverlaps = 0;

				while (!Overlap) {
					double dx = this->m_leftchild->m_cen_bsph[0];
					double dy = this->m_leftchild->m_cen_bsph[1];
					double dz = this->m_leftchild->m_cen_bsph[2];

					double oldDistance = dx * dx + dy * dy + dz * dz;

					//! Translate particle in 1% increments.
					this->m_leftchild->Translate(-0.01 * x * sumr, -0.01 * y * sumr, -0.01 * z * sumr);

					Overlap = this->checkForOverlap(*m_leftchild, *m_rightchild, numberOfOverlaps, Separation);

					dx = this->m_leftchild->m_cen_bsph[0];
					dy = this->m_leftchild->m_cen_bsph[1];
					dz = this->m_leftchild->m_cen_bsph[2];

					double newDistance = dx * dx + dy * dy + dz * dz;

					//! If the left particle has failed to collide with the
					//! particle at the origin (first condition) or if there are
					//! multiple points of overlap (second condition), the trial is
					//! abandoned and another trajectory is chosen.
					if ((newDistance > oldDistance && newDistance > sumr * sumr) || (numberOfOverlaps > 1)) {
						this->m_leftchild->centreBoundSph();
						this->m_leftchild->centreCOM();
						numberOfOverlaps = 0;
						Overlap = false;
						break;
					}
				}

				//! Newton bisection method to speed up translation.
				//! Needs to be tested.
				//double a = 0.0;
				//double b = 0.0;
				//double c = 0.0;

				//b = this->m_leftchild->m_cen_bsph[2];
				//double Translation = - (b - (a + b) / 2.0);

				//numberOfOverlaps = 0;

				//while (!(abs(Separation - 1) < 0.01 && numberOfOverlaps == 1)) {
				//    this->m_leftchild->Translate(0, 0, Translation);

				//    numberOfOverlaps = 0;

				//    if (this->checkForOverlap(*m_leftchild, *m_rightchild, numberOfOverlaps, Separation)) {
				//        a = this->m_leftchild->m_cen_bsph[2];
				//        Translation = (a + b) / 2.0 - a;
				//    } else {
				//        b = this->m_leftchild->m_cen_bsph[2];
				//        Translation = - (b - (a + b) / 2.0);
				//    }
				//}
			}

			double deltax = m_rightparticle->m_cen_bsph[0] - m_leftparticle->m_cen_bsph[0];
			double deltay = m_rightparticle->m_cen_bsph[1] - m_leftparticle->m_cen_bsph[1];
			double deltaz = m_rightparticle->m_cen_bsph[2] - m_leftparticle->m_cen_bsph[2];

			m_distance_centreToCentre = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);

			//! Calculate properties of this particle.
			this->calcBoundSph();
			this->calcCOM();

			//! If particle contains PAH to be traced write POV-Ray translation
			//! of particle
			this->centreBoundSph();

			//! Particle is centred about its bounding sphere for the purpose
			//! generating its structure. Now that it is complete centre it
			//! about its centre-of-mass.
			this->centreCOM();

		}
		else {
			//! Randomly select the primaries that are touching.
			this->m_leftparticle = m_leftchild->SelectRandomSubparticle(rng);
			this->m_rightparticle = m_rightchild->SelectRandomSubparticle(rng);
		}

		//set the sintertime for the new created primary particle
		SetSinteringTime(std::max(this->m_sint_time, rhsparticle->m_sint_time));
		m_createt = max(m_createt, rhsparticle->m_createt);
		//initialise the variables used to calculate the coalesence ratio
		m_children_vol = m_leftparticle->m_vol + m_rightparticle->m_vol;
		m_children_surf = (m_leftparticle->m_surf + m_rightparticle->m_surf);

		m_leftparticle_vol_old = m_leftparticle->m_vol;
		m_rightparticle_vol_old = m_rightparticle->m_vol;
		m_leftparticle_numPAH = m_leftparticle->m_numPAH;
		m_rightparticle_numPAH = m_rightparticle->m_numPAH;

		m_children_radius = pow(3.0 / (4.0*PI)*(m_children_vol), (ONE_THIRD));

		if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())
			m_children_sintering = SinteringLevel();
		else
			m_children_roundingLevel = CoalescenceLevel();

		m_distance_centreToCentre = m_leftparticle->m_primarydiam / 2.0 + m_rightparticle->m_primarydiam / 2.0;

		if (m_pmodel->getTrackPrimaryCoordinates()) {
			double x1 = m_leftparticle->m_cen_bsph[0];
			double y1 = m_leftparticle->m_cen_bsph[1];
			double z1 = m_leftparticle->m_cen_bsph[2];

			double x2 = m_rightparticle->m_cen_bsph[0];
			double y2 = m_rightparticle->m_cen_bsph[1];
			double z2 = m_rightparticle->m_cen_bsph[2];

			m_distance_centreToCentre = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
		}
		else if (m_pmodel->getTrackPrimarySeparation()) {
			m_distance_centreToCentre = m_leftparticle->m_primarydiam / 2.0 + m_rightparticle->m_primarydiam / 2.0;
		}

		if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())
			CheckSintering();
		else
			CheckRounding();

		//must set all the pointer to NULL otherwise the delete function
		//will also delete the children
		copy_rhs.m_leftchild = NULL;
		copy_rhs.m_rightchild = NULL;
		copy_rhs.m_parent = NULL;
		copy_rhs.m_leftparticle = NULL;
		copy_rhs.m_rightparticle = NULL;
		//	rhsparticle->Clear();
		//    if (fabs(surfacebeforerhs+surfacebefore-m_surf)/m_surf > 1e-6)
		//    {
		//        cout << "error" << surfacebeforerhs<<' ' <<surfacebefore << ' '<<m_surf<<endl;
		////         PrintTree("after");
		//    }
		// if (print)
		//     PrintTree("after");

	}
	return *this;
}

/*!
 * Combines this primary with another.
 *
 * \param[in]       rhs         Particle to add to current instance
 * \param[in,out]   rng         Random number generator
 *
 * \return      Reference to the current instance after rhs has been added
 */
PAHPrimary &PAHPrimary::Fragment(const Primary &rhs, rng_type &rng)
{
    vector<PAH>::const_iterator j;
    const PAHPrimary *rhsparticle = NULL;
    rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);

    //only one PAH in rhs or this particle -> condensation or inception process.
    if ( (rhsparticle->m_numPAH==1) || (m_numPAH==1 ))
    {
        if (rhsparticle->Numprimary()>1)
        {
            // rhsparticle is a soot particle but this paricle repesents a PAH
            // this paricle will condense on this rhsparticle
            // Get a copy of the rhs ready to add to the current particle
            PAHPrimary copy_rhs(*rhsparticle);
            PAHPrimary *target = copy_rhs.SelectRandomSubparticle(rng);

            target->m_PAH.insert(target->m_PAH.end(),m_PAH.begin(),m_PAH.end());
			target->UpdatePrimary();
            CopyParts(&copy_rhs);
            CopyTree(&copy_rhs);
        }
        else
        {
            // this paricle is a soot particle but rhsparticle repesents a PAH
            // rhsparticle will condense on this particle
            // particle has more then one primary select the primary where
            // the PAH condenses to and add it to the list
            if (m_leftchild!=NULL)
            {
                PAHPrimary *target = SelectRandomSubparticle(rng);

                target->m_PAH.insert(target->m_PAH.end(),rhsparticle->m_PAH.begin(),rhsparticle->m_PAH.end());
                target->UpdatePrimary();

            }
            else
            {
                //this particle and rhsparticle are both PAHs, this process should be inception
                m_PAH.insert(m_PAH.end(),rhsparticle->m_PAH.begin(),rhsparticle->m_PAH.end());
                UpdatePrimary();
            }
        }
		UpdateCache();
        //Check the coalescence ratio
        CheckSintering();

	}

    else
    {
            //coagulation process
            PAHPrimary *newleft = new PAHPrimary;
		    PAHPrimary *newright = new PAHPrimary;
            PAHPrimary copy_rhs(*rhsparticle);
		    //rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);

           // bool print=false;
           // if (this->m_numprimary>2 && notconstpah->m_numprimary>2)
           // {
            //    PrintTree("before1");
           //     notconstpah->PrintTree("before2");
             //   print=true;
              //  cout << "printing tree"<<endl;
            //}

            //select where to add the second particle
            boost::bernoulli_distribution<> bernoulliDistrib;
            boost::variate_generator<rng_type &, boost::bernoulli_distribution<> > leftRightChooser(rng, bernoulliDistrib);
		    if (leftRightChooser())
		    {
			    newleft->CopyParts(this);
			    newright->CopyParts(&copy_rhs);
		    }
		    else
		    {
			    newright->CopyParts(this);
			    newleft->CopyParts(&copy_rhs);
		    }
            //set the pointers
		    m_leftchild=newleft;
		    m_rightchild=newright;
		    newright->m_parent=this;
		    newleft->m_parent=this;
            //set the pointers to the parent node
		    if (newleft->m_leftchild!=NULL)
		    {
			    newleft->m_leftchild->m_parent=newleft;
			    newleft->m_rightchild->m_parent=newleft;
		    }
		    if (newright->m_leftchild!=NULL)
		    {
			    newright->m_leftchild->m_parent=newright;
			    newright->m_rightchild->m_parent=newright;
		    }
            m_children_roundingLevel=0;
		    UpdateCache();
            //select the primaries that are touching
		    this->m_leftparticle=m_leftchild->SelectRandomSubparticle(rng);
		    this->m_rightparticle=m_rightchild->SelectRandomSubparticle(rng);

            //set the sintertime for the new created primary particle
            SetSinteringTime(std::max(this->m_sint_time, rhsparticle->m_sint_time));
            //initialise the variables used to calculate the coalesence ratio
            m_children_vol=m_leftparticle->m_vol+m_rightparticle->m_vol;
            m_children_surf=(m_leftparticle->m_surf+m_rightparticle->m_surf);
            m_leftparticle_vol_old=m_leftparticle->m_vol;
            m_rightparticle_vol_old=m_rightparticle->m_vol;
            m_leftparticle_numPAH=m_leftparticle->m_numPAH;
            m_rightparticle_numPAH=m_rightparticle->m_numPAH;
            m_children_radius=pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
            m_children_roundingLevel=CoalescenceLevel();
            CheckSintering();
            //must set all the pointer to NULL otherwise the delete function
            //will also delete the children
	        copy_rhs.m_leftchild=NULL;
	        copy_rhs.m_rightchild=NULL;
	        copy_rhs.m_parent=NULL;
            copy_rhs.m_leftparticle=NULL;
            copy_rhs.m_rightparticle=NULL;
		//	rhsparticle->Clear();
        //    if (fabs(surfacebeforerhs+surfacebefore-m_surf)/m_surf > 1e-6)
        //    {
        //        cout << "error" << surfacebeforerhs<<' ' <<surfacebefore << ' '<<m_surf<<endl;
        ////         PrintTree("after");
        //    }
        // if (print)
    //     PrintTree("after");

	}
    return *this;
}

//from bintree model//
//! Check for the overlap of primary particles.
/*!
* @brief       Copy state-space and derived properties from source *
*  @return Whether there is the overlap of primary particles.
*/
bool PAHPrimary::checkForOverlap(PAHPrimary &target, PAHPrimary &bullet, int &numberOfOverlaps, double &Separation)
{
	bool Overlap = false;

	if (target.isLeaf()) {
		//! Target is a leaf.
		if (bullet.isLeaf()) {
			Overlap = particlesOverlap(target.boundSphCentre(), target.Radius(), bullet.boundSphCentre(), bullet.Radius(), Separation);

			//! Keep a running total of the number of overlaps, and the left
			//! and right particles should not be assigned unless there is
			//! overlap.
			if (Overlap) {
				numberOfOverlaps += 1;

				//! Bullet is a leaf (both leaves).
				this->m_leftparticle = &target;
				this->m_rightparticle = &bullet;
			}

			return Overlap;
		}
		else {
			//! Bullet is not a leaf, call sub-nodes.
			Overlap = checkForOverlap(target, *bullet.m_leftchild, numberOfOverlaps, Separation);
			Overlap = checkForOverlap(target, *bullet.m_rightchild, numberOfOverlaps, Separation) || Overlap;

			return Overlap;
		}
	}
	else {
		//! Target is not a leaf.
		if (bullet.isLeaf()) {
			//! Bullet is a leaf, call target sub-nodes.
			Overlap = checkForOverlap(*target.m_leftchild, bullet, numberOfOverlaps, Separation);
			Overlap = checkForOverlap(*target.m_rightchild, bullet, numberOfOverlaps, Separation) || Overlap;

			return Overlap;
		}
		else {
			//! Bullet is not a leaf (neither is a leaf), check all left/right
			//! collision combinations.
			//! Target left and bullet left.
			Overlap = checkForOverlap(*target.m_leftchild, *bullet.m_leftchild, numberOfOverlaps, Separation);

			//! Target left and bullet right.
			Overlap = checkForOverlap(*target.m_leftchild, *bullet.m_rightchild, numberOfOverlaps, Separation) || Overlap;

			//! Target right and bullet left.
			Overlap = checkForOverlap(*target.m_rightchild, *bullet.m_leftchild, numberOfOverlaps, Separation) || Overlap;

			//! Target right and bullet right.
			Overlap = checkForOverlap(*target.m_rightchild, *bullet.m_rightchild, numberOfOverlaps, Separation) || Overlap;

			return Overlap;
		}
	}
}

//from bintree model.//
//! Determine whether the particles overlap.
/*!
*  @param[in]  p1         Coordinates of sphere 1.
*  @param[in]  r1         Radius of sphere 1.
*  @param[in]  p2         Coordinates of sphere 2.
*  @param[in]  r2         Radius of sphere 2.
*  @param[out] Separation Separation between the centres of the primary
*                         particles for use with the Newton bisection
*                         method.
*
*  @return Do the particles overlap?
*/
bool PAHPrimary::particlesOverlap(const Coords::Vector &p1, double r1,
	const Coords::Vector &p2, double r2, double &Separation)
{
	double sumrsqr;

	sumrsqr = r1 + r2;

	//! Calculate the square of the sum of the radii.
	sumrsqr *= sumrsqr;

	//! Calculate dx, dy and dz.
	double xdev = p2[0] - p1[0];
	double ydev = p2[1] - p1[1];
	double zdev = p2[2] - p1[2];

	//! Calculate dx, dy and dz squared.
	double dxsqr = xdev * xdev;
	double dysqr = ydev * ydev;
	double dzsqr = zdev * zdev;

	//! The particles overlap if the centre-to-centre distance is less than the
	//! sum of the primary radii.
	if (sumrsqr < dxsqr + dysqr + dzsqr) {
		Separation = 0.0;
		return false;
	}
	else {
		Separation = sqrt(dxsqr + dysqr + dzsqr);
		return true;
	}
}

//from bintree model//
//! Calculates the radius of gyration of a particle.
double PAHPrimary::GetRadiusOfGyration() const
{
	double sum = 0;
	double mass;
	double totalmass = 0;
	double r2;
	double Rg;
	double rix, riy, riz, rjx, rjy, rjz, drx, dry, drz;
	vector<fvector> coords;

	if (m_numprimary == 1) {
		//if single primary return the primary radius	
		Rg = m_primarydiam / 2.0;

	}else{

		this->GetPriCoords(coords);

		if (m_pmodel->getTrackPrimaryCoordinates()) {
			//! Calculation is based on Eq. (1) in R. Jullien, Transparency effects
			//! in cluster-cluster aggregation with linear trajectories, J. Phys. A
			//! 17 (1984) L771-L776. 
			for (int i = 0; i != coords.size(); ++i) {
				for (int j = 0; j != coords.size(); ++j) {
					rix = coords[i][0];
					riy = coords[i][1];
					riz = coords[i][2];

					rjx = coords[j][0];
					rjy = coords[j][1];
					rjz = coords[j][2];

					//! Expansion of (r_i - r_j)^2 term. Dot product of r vectors.
					sum += rix * rix + riy * riy + riz * riz +
						rjx * rjx + rjy * rjy + rjz * rjz -
						2 * (rix * rjx + riy * rjy + riz * rjz);
				}
			}

			Rg = sqrt(sum / 2 / coords.size() / coords.size());
		}
		else {
			for (unsigned int i = 0; i != coords.size(); ++i) {
				//! Mass is proportional to the cube of the radius.
				mass = coords[i][3] * coords[i][3] * coords[i][3];
				r2 = coords[i][0] * coords[i][0] + coords[i][1] * coords[i][1] + coords[i][2] * coords[i][2];
				sum += mass * r2;
				totalmass += mass;
			}

			Rg = sqrt(sum / totalmass);
		}	
	}

	return Rg;
}

//from bintree model//
//! Returns a vector of primary coordinates and radius (4D).
/*!
*  @param[in] coords The first three returned values are the cartesian x, y, z
*                    coordinates, the final value is the radius.
*/
void PAHPrimary::GetPriCoords(std::vector<fvector> &coords) const
{
	if (isLeaf()) {
		fvector c(4);
		c[0] = m_cen_mass[0];
		c[1] = m_cen_mass[1];
		c[2] = m_cen_mass[2];
		c[3] = m_r;
		coords.push_back(c);
	}
	else {
		m_leftchild->GetPriCoords(coords);
		m_rightchild->GetPriCoords(coords);
	}
}

//returns the CoalescenceLevel of the two primaries that are connected by this node
//only for coordinates are not tracked
double PAHPrimary::CoalescenceLevel()
{
	if (m_leftparticle != NULL)
	{

		double dV = 0;


		//make sure that the volume growth does not arise from a former coalescence event
		//(number of PAHs must remain const, or increase by 1 in case of a condensation
		// event
		// if the volume changes due to a coalescence event the chilren surface is changed in
		// PAHPrimary::ChangePointer(PAHPrimary *source, PAHPrimary *target)

		// m_leftparticle->m_numPAH is the actual number of PAHs in the left one of the two
		// touching primary particles.

		// m_leftparticle_numPAH is the the number of PAHs in the left one of the two touching
		// particles last time this method was called.

		// The following test checks that at most one PAH has been added to the left one
		// of the two touching primaries since this method was last called.
		if (m_leftparticle->m_numPAH - m_leftparticle_numPAH < 2)
		{
			//calculate the volume change of the left particle in the last timestep
			dV += m_leftparticle->m_vol - m_leftparticle_vol_old;
		}

		//adjust the number of PAHs and vol for the next step
		m_leftparticle_numPAH = m_leftparticle->m_numPAH;
		m_leftparticle_vol_old = m_leftparticle->m_vol;

		//make sure that the volume growth does not arise from a coalescence event
		//(number of PAHs must remain const, or increase by1 in case of a condensation
		// event
		if (m_rightparticle->m_numPAH - m_rightparticle_numPAH < 2)
		{
			//calculate the volume change of the right particle in the last timestep
			dV += m_rightparticle->m_vol - m_rightparticle_vol_old;
		}

		//adjust the number of PAHs and vol for the next step
		m_rightparticle_numPAH = m_rightparticle->m_numPAH;
		m_rightparticle_vol_old = m_rightparticle->m_vol;

		//update the children volume, this value is used to calculate dV n the next timestep
		m_children_vol = m_rightparticle->m_vol + m_leftparticle->m_vol;
		//calculate dS, at the moment it is assumed that the particles always grow
		double ct = m_pmodel->Components(0)->CoalescThresh();
		const double dS = dV*ct / m_children_radius;
		//update the radius for the next event
		m_children_radius = pow(3.0 / (4.0*PI)*(m_children_vol), (ONE_THIRD));
		m_children_surf += dS;

		// const double spherical_surface=4*PI*m_children_radius*m_children_radius;
		//// double two_1_3=pow(2,-1*ONE_THIRD);
		// const double two_1_3=0.79370052231642452;
		// double clevel= ((spherical_surface/m_children_surf)-two_1_3)/(1-two_1_3);
		// if (clevel<0)
		//     return 0;
		// else
		//     return clevel;
		return RoundingLevel();
	}
	else
		return 0;
}

//calculates the fractal dimension of the particle and stores it in m_fdim
//void PAHPrimary::CalcFractalDimension()
//{
//	Sweep::Imaging::ParticleImage img;
//    // construct the particle by colliding the primary particles
//	img.constructSubParttree(this);
//	double L,W;
//    // calculate the length and the width of the particle
//    img.LengthWidth(L,W);
//    // calculate the radius of gyration
//    m_Rg=img.RadiusofGyration();
//	m_sqrtLW=sqrt(L*W);
//	m_LdivW=L/W;
//    m_Rg=m_Rg*1e-9;
//    m_fdim=log((double)m_numprimary)/log(2*m_Rg/(m_primarydiam/m_numprimary));
//  /*  if (m_fdim>0 && m_fdim<3)
//    {
//       string filename;
//       filename=cstr(m_fdim)+".3d";
//       ofstream out;
//       out.open(filename.c_str());
//       img.Write3dout(out,0,0,0);
//       out.close();
//    }*/
//
//}

// sets all the childrenproperties to zero, this function is used after the children are coalesced
void PAHPrimary::ResetChildrenProperties()
{
            m_children_roundingLevel=0.0;
            m_children_surf=0.0;
            m_children_vol=0.0;
            m_distance_centreToCentre=0.0;
            m_rightparticle_vol_old=0.0;
            m_leftparticle_vol_old=0.0;
            m_leftparticle_numPAH=0;
            m_rightparticle_numPAH=0;
            m_avg_coalesc=0; //no track

			if (m_pmodel->getTrackPrimaryCoordinates() || m_pmodel->getTrackPrimarySeparation())
			{
				m_children_sintering = 0.0;
				m_children_radius = 0.0; 
				m_avg_sinter = 0;
				m_sint_rate = 0.0;
				m_children_sumCap = 0.0; //used to calculate geometric volume of a particle
			}
}

//merges the left and the right primary particle to one primary particle
PAHPrimary &PAHPrimary::Merge()
{
	//    if(this->m_numprimary>5)
	//       cout<<"tset";
	//      PrintTree("before.inp");
	//! Declare pointers for coordinate/separation tracking model

	PAHPrimary *small_prim; //!< smaller of the merging primaries
	PAHPrimary *big_prim;	//!< larger of the merging primaries
	PAHPrimary *new_prim;	//!< new (merged) primary

	//initialise parameters
	double r_big, r_small, d_ij, x_ij;

	vector<PAH>::const_iterator j;

	//make sure this primary has children to merge
	if (m_leftchild != NULL)
	{
		if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()){
		
			d_ij = m_distance_centreToCentre;

			//! Update primaries
			m_leftparticle->UpdatePrimary();
			m_rightparticle->UpdatePrimary();

			//! If the centre to centre distance is tracked we need to know which is the smaller primary of the merging pair
			if (m_leftparticle->m_primarydiam > m_rightparticle->m_primarydiam){
				small_prim = m_rightparticle;
				big_prim = m_leftparticle;
			}
			else{
				small_prim = m_leftparticle;
				big_prim = m_rightparticle;
			}

			r_big = big_prim->m_primarydiam / 2.0;
			r_small = small_prim->m_primarydiam / 2.0;
			x_ij = min((d_ij*d_ij - r_small*r_small + r_big*r_big) / (2.0*d_ij), r_big); //!<distance from neck to centre of larger primary (x_ij < r_i)

			double V_prim = small_prim->m_primaryvol; //!< Volume of smaller primary
			
			//! If there is a cap (i.e. the smaller primary isn't completely enclosed)
			if (d_ij + r_small > r_big){
				//! calculate cap volume of larger primary
				double V_cap = 2.0*M_PI*r_big*r_big*r_big / 3.0 + M_PI*x_ij*x_ij*x_ij / 3.0 - M_PI*r_big*r_big*x_ij;
				//! Subtract cap volume from the merging primary's volume
				double dV = max(V_prim - V_cap, 0.0);

				//! Adjust larger primary to incorporate excess volume
				if (dV > 0.0)
					big_prim->AdjustPrimary(dV, d_ij, small_prim);
			}
		}

		if (m_leftchild == m_leftparticle && m_rightchild == m_rightparticle)
		{
			//this node has only two primaries in its subtree
			//it is possible that this node is not the root node and belongs to a bigger particle
			//copy the PAHs of both children to the parent node
			if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())
			{

				new_prim = this; //!< new primary

				new_prim->m_primarydiam = big_prim->m_primarydiam;
			}

			//m_PAH is empty, therefore no need to append
			m_PAH = m_rightparticle->m_PAH;

			m_PAH.insert(this->m_PAH.end(), m_leftparticle->m_PAH.begin(), m_leftparticle->m_PAH.end());

			if (m_pmodel->getTrackPrimaryCoordinates()){
				new_prim->m_cen_bsph = big_prim->m_cen_bsph;
				new_prim->m_cen_mass = big_prim->m_cen_mass;
			}

			//! Update the pointers that pointed to the two former children
			if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()) {
				ChangePointer(m_leftchild, this);
				ChangePointer(m_rightchild, this);
			}
			else{
				//! If coordinates/separation are tracked the pointer to the larger primary is changed first
				//! so that primary properties are correctly calculated when adding neighbours
				ChangePointer(big_prim, this, small_prim, this);
				ChangePointer(small_prim, this, small_prim, this);
			}

			//delete the children (destructor is recursive for this class)
			delete m_leftchild;
			delete m_rightchild;
			m_leftchild = NULL;
			m_rightchild = NULL;
			m_leftparticle = NULL;
			m_rightparticle = NULL;

			//set the children properties to zero, this node has no more children
			ResetChildrenProperties();
			UpdatePrimary();

			// Only update the cache on m_parent if the sintering level of
			// m_parent if the sintering level won't call a merge on the
			// parent node. Otherwise, the *this* memory address could be
			// removed from the tree and segmentation faults will result!
			//if (m_parent != NULL) {
			//	if (!m_parent->MergeCondition()) {
			//		m_parent->UpdateCache();
			//	}
			//}

		}

		else //******m_leftchild != m_leftparticle or m_rightchild != m_rightparticle******//
		{
			//! If primary coordinates or primary separations are not tracked then 
			//! select subtree to keep the tree balanced
			if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()) {

				if (m_leftchild->m_numprimary < m_rightchild->m_numprimary)
				{
					//append to left subtree because there are fewer primaries
					//this is only to keep the tree balanced
					PAHPrimary *oldleftparticle = m_leftparticle;

					//copy the PAHs

					//for (j=oldleftparticle->m_PAH.begin(); j!=oldleftparticle->m_PAH.end(); ++j) {
					//	m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(),PAH(*j));
					//}

					m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(), oldleftparticle->m_PAH.begin(), oldleftparticle->m_PAH.end());
					m_rightparticle->UpdatePrimary();
					//set the pointers from the leftprimary to the rightprimary
					//this will be the new bigger primary
					oldleftparticle->ChangePointer(oldleftparticle, m_rightparticle);
					m_rightparticle->ChangePointer(m_rightparticle, m_rightparticle);

					//set the pointer to the parent node
					if (oldleftparticle->m_parent->m_leftchild == oldleftparticle)
					{
						oldleftparticle->m_parent->m_leftchild = m_rightchild;
					}
					else
					{
						oldleftparticle->m_parent->m_rightchild = m_rightchild;
					}
					m_rightchild->m_parent = oldleftparticle->m_parent;


					PAHPrimary *oldleftchild = m_leftchild;
					PAHPrimary *oldparent = m_parent;

					//copy the properties of the former leftchild to this node
					// so that it can be removed from the aggregate tree structure
					CopyParts(oldleftchild);

					// Now break the links to the tree structure in oldleftchild and free it
					oldleftchild->m_leftchild = NULL;
					oldleftchild->m_rightchild = NULL;
					delete oldleftchild;

					m_parent = oldparent;
					if (m_leftchild != NULL)
					{
						m_rightchild->m_parent = this;
						m_leftchild->m_parent = this;
					}

					delete oldleftparticle;

				}

				else //******m_leftchild->m_numprimary > m_rightchild->m_numprimary******//
				{
					//append to right subtree
					PAHPrimary *oldrightparticle = m_rightparticle;

					//	for (j=oldrightparticle->m_PAH.begin(); j!=oldrightparticle->m_PAH.end(); ++j) {
					//		m_leftparticle->m_PAH.insert(m_leftparticle->m_PAH.end(),PAH(*j));
					//	}

					m_leftparticle->m_PAH.insert(m_leftparticle->m_PAH.end(), oldrightparticle->m_PAH.begin(), oldrightparticle->m_PAH.end());
					m_leftparticle->UpdatePrimary();

					oldrightparticle->ChangePointer(oldrightparticle, m_leftparticle);
					m_leftparticle->ChangePointer(m_leftparticle, m_leftparticle);

					if (oldrightparticle->m_parent->m_leftchild == oldrightparticle)
					{
						oldrightparticle->m_parent->m_leftchild = m_leftchild;
					}
					else
					{
						oldrightparticle->m_parent->m_rightchild = m_leftchild;
					}
					m_leftchild->m_parent = oldrightparticle->m_parent;

					//ReleaseMem(oldrightparticle);
					PAHPrimary *oldrightchild = m_rightchild;
					PAHPrimary *oldparent = m_parent;

					//copy the properties of the former leftchild to this node
					// so that it can be removed from the aggregate tree structure
					CopyParts(oldrightchild);

					// Now break the links to the tree structure in oldrightchild and free it
					oldrightchild->m_leftchild = NULL;
					oldrightchild->m_rightchild = NULL;
					delete oldrightchild;

					m_parent = oldparent;
					if (m_leftchild != NULL)
					{
						m_rightchild->m_parent = this;
						m_leftchild->m_parent = this;
					}
					delete oldrightparticle;
				}
			}
			else{
				//! If the coordinates/separations are tracked then the larger primary becomes the new primary

				//! the new (merged) primary
				new_prim = big_prim;

				//! left/right flag
				bool newleft;
				if (new_prim == m_leftparticle){
					newleft = true;
				}
				else{
					newleft = false;
				}

				PAHPrimary *oldparticle = small_prim;
				//! update composition
				if (newleft)
				{
					new_prim->m_PAH = m_leftparticle->m_PAH;
					new_prim->m_PAH.insert(m_leftparticle->m_PAH.end(), m_rightparticle->m_PAH.begin(), m_rightparticle->m_PAH.end());
					//new_prim->UpdatePrimary();
				}
				else
				{
					new_prim->m_PAH = m_rightparticle->m_PAH;
					new_prim->m_PAH.insert(m_rightparticle->m_PAH.end(), m_leftparticle->m_PAH.begin(), m_leftparticle->m_PAH.end());
					//new_prim->UpdatePrimary();
				}

				//! update pointers to neighbours
				new_prim->ChangePointer(new_prim, new_prim, small_prim, this);
				oldparticle->ChangePointer(oldparticle, new_prim, small_prim, this);

				// Set the pointer to the parent node
				PAHPrimary *oldchild = NULL;

				if (newleft){
					if (oldparticle->m_parent->m_leftchild == oldparticle) {
						oldparticle->m_parent->m_leftchild = m_leftchild;
					}
					else {
						oldparticle->m_parent->m_rightchild = m_leftchild;
					}
					m_leftchild->m_parent = oldparticle->m_parent;

					oldchild = m_rightchild;
				}
				else{
					if (oldparticle->m_parent->m_leftchild == oldparticle) {
						oldparticle->m_parent->m_leftchild = m_rightchild;
					}
					else {
						oldparticle->m_parent->m_rightchild = m_rightchild;
					}
					m_rightchild->m_parent = oldparticle->m_parent;

					oldchild = m_leftchild;
				}

				PAHPrimary *oldparent = m_parent;

				// Copy the properties of the former leftchild to this node
				// so that it can be removed from the aggregate tree structure
				CopyParts(oldchild);

				// Now break the links to the tree structure in oldleftchild
				// in order to free it
				oldchild->m_leftchild = NULL;
				oldchild->m_rightchild = NULL;
				delete oldchild;

				m_parent = oldparent;

				if (m_leftchild != NULL) {
					m_rightchild->m_parent = this;
					m_leftchild->m_parent = this;
				}

				delete oldparticle;

			}
		}

		//if coordinates are tracked then update the tracked radii
		if (m_pmodel->getTrackPrimaryCoordinates()){
			new_prim->setRadius(new_prim->m_primarydiam / 2.0);
		}

		UpdateCache();
		//      PrintTree("after.inp");
	}

	return *this;
}

////from bintree model//
//void PAHPrimary::SumNeighbours(PAHPrimary *prim, double &sumterm) {
//
//	double d_ij = m_parent->m_distance_centreToCentre;
//	double r_i = prim->m_primarydiam / 2.0;
//	double r_j = 0.0;
//	double x_ij = 0.0;
//
//	//check if a neighbour of prim
//	if (m_parent->m_leftparticle == prim) {
//		//right particle is a neighbour
//		r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
//	}
//	else if (m_parent->m_rightparticle == prim) {
//		//left particle is a neighbour
//		r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
//	}
//	else {
//		//not a neighbour
//		r_j = 0.0;
//	}
//
//	//if the node connects a neighbourt then calculate the summation term
//	if (r_j > 0.0){
//		//the volumes and radii of neighbours remain unchanged
//		//the centre to centre separations increase to allow growth of the primary
//		x_ij = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / (2.0*d_ij);
//		sumterm += x_ij - 2.0*r_i + pow(r_i, 2.0) / x_ij;
//	}
//
//	//continue working up the binary tree
//	if (m_parent->m_parent != NULL){
//		m_parent->SumNeighbours(prim, sumterm);
//	}
//}

//! From bintree model.
/*!
* @brief       Identify neighbours and update centre to centre separation ignoring specified neighbour
*
* Works up the binary tree identifying the neighbours of the adjusted primary
* and updates the centre to centre separation. Also sums the contribution from
* neighbours to the free surface area.
*
* @param[in]   prim			Pointer to the primary being adjusted
* @param[in]   delta_r			Change in radius of prim
* @param[in]   sumterm			Sum of contributions from neighbours to the free surface area of prim
* @param[in]   prim_ignore		Pointer to the primary to be ignored
*/
void PAHPrimary::UpdateConnectivity(PAHPrimary *prim, double delta_r, PAHPrimary *prim_ignore){

	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;
	PAHPrimary *neighbour = NULL;

	//! check if a neighbour of prim
	if (m_parent->m_leftparticle == prim && m_parent->m_rightparticle != prim_ignore) {
		//! right particle is a neighbour
		neighbour = m_parent->m_rightparticle;
		r_j = neighbour->m_primarydiam / 2.0;
	}
	else if (m_parent->m_rightparticle == prim &&  m_parent->m_leftparticle != prim_ignore) {
		//! left particle is a neighbour
		neighbour = m_parent->m_leftparticle;
		r_j = neighbour->m_primarydiam / 2.0;
	}
	else {
		//! not a neighbour
		r_j = 0.0;
	}

	if (r_j > 0.0){

		double d_ij_old = d_ij;
		x_ij = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / (2.0*d_ij);
		//! update centre to centre separation
		//making sure centre to centre separation remains smaller than the sum of the radii
		d_ij = min(d_ij + r_i * delta_r / x_ij, r_i + r_j + delta_r);
		m_parent->m_distance_centreToCentre = d_ij;

		//! if primary coordinates are tracked then we need to update the coordinates of the neighbour 
		if (m_pmodel->getTrackPrimaryCoordinates()) {
			//! get (unit) vector separating prim and neighbour
			Coords::Vector u = UnitVector(prim->boundSphCentre(), neighbour->boundSphCentre());
			double delta_d = d_ij - d_ij_old; //!< change in separation (magnitude)
			//! translate the neighbour 
			neighbour->TranslatePrimary(u, delta_d);
			//! translate all neighbours of the neighbour except prim
			neighbour->TranslateNeighbours(neighbour, u, delta_d, prim);
		}
	}

	//continue working up the binary tree
	if (m_parent->m_parent != NULL){
		m_parent->UpdateConnectivity(prim, delta_r, prim_ignore);
	}
}

//! From bintree model.
//! Returns true if this node is a leaf (has no children).
bool PAHPrimary::isLeaf(void) const
{
	return (m_leftchild == NULL) && (m_rightchild == NULL);
}

//! From bintree model.
//! Returns the bounding-sphere centre.
const Coords::Vector &PAHPrimary::boundSphCentre(void) const
{
	return m_cen_bsph;
}

//! From bintree model.
//! Calculates the bounding sphere position and radius using
//! the left and right child node values.
void PAHPrimary::calcBoundSph(void)
{
	if ((m_leftchild != NULL) && (m_rightchild != NULL)) {
		// Calculate bounding spheres of children.
		m_leftchild->calcBoundSph();
		m_rightchild->calcBoundSph();

		// Calculate translation between left and right spheres.
		double dx = m_rightchild->m_cen_bsph[0] - m_leftchild->m_cen_bsph[0];
		double dy = m_rightchild->m_cen_bsph[1] - m_leftchild->m_cen_bsph[1];
		double dz = m_rightchild->m_cen_bsph[2] - m_leftchild->m_cen_bsph[2];

		// Calculate bounding sphere centre.
		m_cen_bsph[0] = m_leftchild->m_cen_bsph[0] + (0.5 * dx);
		m_cen_bsph[1] = m_leftchild->m_cen_bsph[1] + (0.5 * dy);
		m_cen_bsph[2] = m_leftchild->m_cen_bsph[2] + (0.5 * dz);

		// Calculate bounding sphere radius.
		setRadius(sqrt((dx*dx) + (dy*dy) + (dz*dz)));
	}
}

//! From bintree model.
//! Calculates the centre-of-mass using the left and right child node values.
void PAHPrimary::calcCOM(void)
{
	if ((m_leftchild != NULL) && (m_rightchild != NULL)) {
		//! Calculate centres-of-mass of left and right children.
		m_leftchild->calcCOM();
		m_rightchild->calcCOM();

		//! Calculate inverse total mass of left and right children.
		m_mass = m_leftchild->m_mass + m_rightchild->m_mass;
		double invtotmass = 1.0 / m_mass;

		//! Now calculate centre-of-mass.
		for (unsigned int i = 0; i != 3; ++i) {
			m_cen_mass[i] = m_leftchild->m_cen_mass[i] * m_leftchild->m_mass;
			m_cen_mass[i] += m_rightchild->m_cen_mass[i] * m_rightchild->m_mass;
			m_cen_mass[i] *= invtotmass;
		}
	}
	else {
		//! If there are no children, then the centre-of-mass and bounding-
		//! sphere centre are the same.
		m_cen_mass[0] = m_cen_bsph[0];
		m_cen_mass[1] = m_cen_bsph[1];
		m_cen_mass[2] = m_cen_bsph[2];
	}
}

//! From bintree model.
//! Put the bounding-sphere at the origin.
void PAHPrimary::centreBoundSph(void)
{
	Translate(-m_cen_bsph[0], -m_cen_bsph[1], -m_cen_bsph[2]);
}

//! From bintree model.
//! Put the centre-of-mass at the origin.
void PAHPrimary::centreCOM(void)
{
	Translate(-m_cen_mass[0], -m_cen_mass[1], -m_cen_mass[2]);
}

//! From bintree model.
/*!
*  Randomly rotates the aggregate node and child structure about its centre of
*  mass.
*
*  @param[in]    theta    Rotation about the pole.
*  @param[in]    V        Vector for performing the reflection.
*/
void PAHPrimary::rotateCOM(double theta, fvector V)
{
	//! Move the aggregate so that its centre-of-mass is at the origin. Store
	//! the coordinates, so that they can be restored afterwards.
	Coords::Vector D(m_cen_mass);
	Translate(-D.X(), -D.Y(), -D.Z());

	//! Create transformation matrix.
	Coords::Matrix M;
	//M.SetIdentity();
	M.rotateArvo(theta, V);

	//! Rotate child nodes.
	if (m_leftchild != NULL) m_leftchild->transform(M);
	if (m_rightchild != NULL) m_rightchild->transform(M);

	//! Rotate bounding-sphere coordinates.
	m_cen_bsph = M.Mult(m_cen_bsph);

	//! Restore centre-of-mass coordinates.
	Translate(D.X(), D.Y(), D.Z());
}

//! From bintree model.
/*!
*  @brief Sets the radius of the bounding sphere.
*
*  @param[in]    r    Radius of bounding sphere.
*/
void PAHPrimary::setRadius(double r)
{
	m_r = r;
	m_r2 = r * r;
	m_r3 = m_r2 * m_r;
}


//! From bintree model.
//! Returns the bounding sphere radius.
double PAHPrimary::Radius(void) const
{
	return m_r;
}

//! From bintree model.
/*!
*  Transform the primary particle coordinates using the transformation matrix
*  so as to rotate it.
*
*  @param[in]    mat               Transformation matrix.
*  @param[in]    PAHTracerMatch    Flag used to indicate whether the particle
*                                  contains the PAH to be traced.
*/
void PAHPrimary::transform(const Coords::Matrix &mat)
{
	//! Descend binary tree to the leaf nodes, i.e. single primary particles.
	if (m_leftchild != NULL)
		m_leftchild->transform(mat);
	if (m_rightchild != NULL)
		m_rightchild->transform(mat);

	//! Rotate centre-of-mass and bounding sphere coordinates.
	m_cen_mass = mat.Mult(m_cen_mass);
	m_cen_bsph = mat.Mult(m_cen_bsph);
}

//! From bintree model.
/*!
*  Translates (moves) the aggregate node and child structure by the given
*  amounts along the cartesian axes.
*
*  @param[in]    dx    Distance to translate in the x-axis.
*  @param[in]    dy    Distance to translate in the y-axis.
*  @param[in]    dz    Distance to translate in the z-axis.
*/
void PAHPrimary::Translate(double dx, double dy, double dz)
{
	//! Translate child branches.
	if (m_leftchild != NULL) m_leftchild->Translate(dx, dy, dz);
	if (m_rightchild != NULL) m_rightchild->Translate(dx, dy, dz);

	//! Translate bounding sphere centre.
	m_cen_bsph.Translate(dx, dy, dz);

	//! Translate centre-of-mass.
	m_cen_mass.Translate(dx, dy, dz);
}

//! From bintree model.
//! Write the coordinates of the primaries in the particle pointed to by the
//! this pointer. Units of nm.
void PAHPrimary::writePrimaryCoordinatesRadius(void)
{
	if (m_leftchild != NULL)
		m_leftchild->writePrimaryCoordinatesRadius();

	if (m_rightchild != NULL)
		m_rightchild->writePrimaryCoordinatesRadius();

	std::ofstream outfile;

	double r = Radius() * 1.0e9;
	double a = m_cen_mass[0] * 1.0e9;
	double b = m_cen_mass[1] * 1.0e9;
	double c = m_cen_mass[2] * 1.0e9;

	//! There is a need for these two conditions as upon exit of the function
	//! the aggregate node will pass through this part of the code but are only
	//! we are only interested in the coordinates of the primaries.
	if (isLeaf()) {
		outfile.open("Spheres.m", std::ios_base::app);
		outfile << "surf(x*" << r << "+" << a << ",y*" << r << "+" << b << ",z*" << r << "+" << c << ");\n";
		outfile.close();
	}
}

void PAHPrimary::ReleaseMem()
{ 
    m_PAH.clear();
}

//! Overload function. If coordinates are not tracked.
/*!//add explanation
 * @param[in] source Pointer to the original particle
 * @param[in,out] target Pointer to the new particle
*/
void PAHPrimary::ChangePointer(PAHPrimary *source, PAHPrimary *target)
{
		if(m_rightparticle==source){
			m_rightparticle=target;
            double sphericalsurface=
                4*PI*pow(3*(m_leftparticle->Volume()+m_rightparticle->Volume())/(4*PI),TWO_THIRDS);
			m_children_surf = sphericalsurface / (m_children_roundingLevel*0.2063 + 0.7937);    //sphericalsurface/(m_children_roundingLevel*(1-2^(-1/3))+2^(-1/3))
		}
		if(m_leftparticle==source){
			m_leftparticle=target;
            double sphericalsurface=
                4*PI*pow(3*(m_leftparticle->Volume()+m_rightparticle->Volume())/(4*PI),TWO_THIRDS);
			m_children_surf = sphericalsurface / (m_children_roundingLevel*0.2063 + 0.7937);    //sphericalsurface/(m_children_roundingLevel*(1-2^(-1/3))+2^(-1/3))

		}

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->ChangePointer(source,target);
    }

}

//! Overload function. Used if primary coordinates are not tracked.
/*!
 * The actual interval over which the update is carried out on a PAH is from
 * lastupdated to t.
 *
 * @param[in]        t        Time up to which to update.
 * @param[in]        model    Particle model defining interpretation of particle data.
 * @param[in]        sys      Cell containing particle and providing gas phase.
 * @param[in,out]    rng      Random number generator.
 *
 */
void PAHPrimary::UpdatePAHs(const double t, const double dt, const Sweep::ParticleModel &model, Cell &sys, int statweight, 
	int ind, rng_type &rng, PartPtrVector &overflow)
{
    // Either the primary has two children or it is a leaf of the
    // tree
    if (m_leftchild!=NULL)
    {
        // Recurse down to the leaves
		m_leftchild->UpdatePAHs(t, dt, model, sys, statweight, ind, rng, overflow);
		m_rightchild->UpdatePAHs(t, dt, model, sys, statweight, ind, rng, overflow);
    }
    else
    {
        // There are PAHs in this primary so update them, if needed
        // Flag to show if any PAH has been changed
        // this flag store the condition of this cluster(primary particle)
        bool m_PAHclusterchanged = false;
        // this flag store the condition of a particular PAH within this cluster
        bool m_InvalidPAH = false;
        // check whether this updated primary particle is a inceptedPAH, used later to track the num of InceptedPAH in the ensemble
        const int m_InceptedPAH = InceptedPAH();

        //! Loop over each PAH in this primary.
        const std::vector<boost::shared_ptr<PAH> >::iterator itEnd = m_PAH.end();
        for (std::vector<boost::shared_ptr<PAH> >::iterator it = m_PAH.begin(); it != itEnd; ++it) {
            
			bool m_PAHchanged = false;

            //! Model parameter.
			/*!
			 * If a jump process reduces the total number of 6-member rings
			 * (excludes 5-member rings) in a PAH (in a particle) below this
			 * threshold it is removed.
             */
			const int thresholdOxidation = model.Components(0)->ThresholdOxidation();

            //! PAH removal.
			/*!
			 * Allow for the removal of PAHs in particles (2 or more PAHs) when
			 * the total number of 6-member rings (excludes 5-member rings) in
			 * the PAH has fallen below, for example, the minimum number of
			 * rings in a PAH for inception. Particularly relevant when
			 * oxidation is strong.
             */
            if((*it)->m_pahstruct->numofC() == 5){
                m_PAHclusterchanged = true; 
				m_PAHchanged = true;
                m_InvalidPAH = CheckInvalidPAHs(*it);
                continue;
			}			

            //! Model parameter.
			/*!
             * Defines when primary particles contain too many PAHs for them all 
			 * to have full access to the surrounding gas phase.  Once a primary
			 * contains more than minPAH PAHs the growth rate of each PAH is
			 * reduced according to the growth factor.
             */
			const double minPAH = model.Components(0)->MinPAH();
			double growthfact = 1.0;

			if (m_numPAH>=minPAH)
			{
				growthfact = model.Components(0)->GrowthFact();
				//double density = model.Components(0)->Density();
				////! PP mass (kg).
				//double m_mass = m_numcarbon*1.9945e-26 + m_numH*1.6621e-27;

				////! Spherical diameter (nm).
				//double diameter = pow(6.0 * m_mass / density / PI, ONE_THIRD);
				//diameter = diameter * 1.0e9;
				//growthfact = 1.950 / diameter ;
				//if (growthfact > 1.0)
				//{
				//	growthfact = 1.0;
				//}
			}

			//! Time for one particular PAH to grow.
			double growtime = t - (*it)->lastupdated;
			assert(growtime >= 0.0);
			double statweightold = statweight;
			int numloops = 0;
			bool calcrates = true;
			double ratefactor = 1.0;
			const int oldNumCarbon = (*it)->m_pahstruct->numofC(); 
			const int oldNumH = (*it)->m_pahstruct->numofH();

			double updatetime;

			//if this is a particle with a single PAH, it may be weighted. 
			//If it is weighted and IWDSA is being used, we do not want to update the PAH, but rather update a clone of 
			//that PAH and create a new particle
			if (m_PAH.size() == 1 && statweight > 1.0 && ParticleModel()->Components(0)->WeightedPAHs()){ 
				PartPtrVector overflowtemp;

				while (growtime > 0.0 && statweight > 1.0){

					boost::shared_ptr<PAH> new_m_PAH((*it)->Clone());

					new_m_PAH->PAH_ID = ID;
					ID++;

					if (numloops > 0){
						calcrates = false;
					}
					else{
						calcrates = true;
					}

					updatetime = sys.Particles().Simulator()->updatePAH(new_m_PAH->m_pahstruct, (*it)->lastupdated, growtime, 1, 1,
						 rng, growthfact*statweight, new_m_PAH->PAH_ID, calcrates, ratefactor);

					/*updatetime = sys.Particles().Simulator()->updatePAH(new_m_PAH->m_pahstruct, (*it)->lastupdated, growtime, 1, 1,
						rng, growthfact*statweight, new_m_PAH->PAH_ID, true, 1.0);*/

					new_m_PAH->lastupdated = updatetime;
					(*it)->lastupdated = updatetime;

					growtime = t - (*it)->lastupdated;
					numloops++;

					//! Invalidate PAH.
					/*!
					* A PAH is made invalid by setting the number of carbons to 5. The
					* reason for only applying this to particles is that we would not
					* want to remove single PAHs in the gas-phase which fall below the
					* minimum number of rings for nception but are still growing.
					*/
					if (new_m_PAH->m_pahstruct->numofRings() < thresholdOxidation ){
						statweight--;
						(*sys.Particles().At(ind)).setStatisticalWeight(statweight);
						new_m_PAH.reset();
						ID--;
					}

					//! See if anything changed, as this will required a call to UpdatePrimary() below.
					//This will also dictate if the new PAH is used to create a new particle, or if it is just destroyed
					else if (growtime > 0.0)
					//else
					{
						//Reduce statistical weight of the particle being updated
						statweight--;
						(*sys.Particles().At(ind)).setStatisticalWeight(statweight);

						Particle *sp = NULL;
						sp = model.CreateParticle(updatetime);
						AggModels::PAHPrimary *pri =
							dynamic_cast<AggModels::PAHPrimary*>((*sp).Primary());
						pri->m_PAH.clear();
						pri->m_PAH.push_back(new_m_PAH);
						sp->SetTime(updatetime);
						//Update the primary and the cache
						pri->UpdatePrimary();
						sp->UpdateCache();

						overflowtemp.push_back(sp);

					}
					else{
						assert(growtime == 0);
						new_m_PAH.reset();
						ID--;
					}

					ratefactor = statweight / statweightold;
				}
				//Final update after statistical weight reaches 1
				//Now update the PAH one last time
				if (growtime > 0.0){

					updatetime = sys.Particles().Simulator()->updatePAH((*it)->m_pahstruct, (*it)->lastupdated, growtime, 1, 0,
						rng, growthfact, (*it)->PAH_ID, true, 1.0);

					(*it)->lastupdated = t;

					//! Invalidate PAH.
					/*!
					* A PAH is made invalid by setting the number of carbons to 5. The
					* reason for only applying this to particles is that we would not
					* want to remove single PAHs in the gas-phase which fall below the
					* minimum number of rings for inception but are still growing.
					*/
					if ((*it)->m_pahstruct->numofRings() < thresholdOxidation){
						(*it)->m_pahstruct->setnumofC(5);
					}


					//! See if anything changed, as this will required a call to UpdatePrimary() below.
					if (oldNumCarbon != (*it)->m_pahstruct->numofC() || oldNumH != (*it)->m_pahstruct->numofH())
					{
						m_PAHclusterchanged = true;
						m_PAHchanged = true;
					}

				}
				//Now, update all created particles that were added to overflowtemp
				PartPtrVector::iterator it1;
				for (it1 = overflowtemp.begin(); it1 != overflowtemp.end(); ++it1){
					AggModels::PAHPrimary *pri =
						dynamic_cast<AggModels::PAHPrimary*>((*(*it1)).Primary());
					//This new particle must also be updated to time t
					const int oldNumCarbon1 = pri->m_PAH[0]->m_pahstruct->numofC();
					const int oldNumH1 = pri->m_PAH[0]->m_pahstruct->numofH();
					pri->UpdatePAHs(t, t - (*it1)->LastUpdateTime(), model, sys, 1, -1, rng, overflow);
					(*it1)->SetTime(t);

					//Update the primary and the cache if the PAH structure has changed.
					//Place the particle (which is a single PAH) in the overflow vector if it is still a valid PAH
					if (pri->GetPAHVector()[0]->m_pahstruct->numofC() != oldNumCarbon1 &&
						pri->GetPAHVector()[0]->m_pahstruct->numofH() != oldNumH1){
						
						if (pri->m_PAH[0]->m_pahstruct->numofRings() >= thresholdOxidation){
							pri->UpdatePrimary();
							(*it1)->UpdateCache();
							overflow.push_back(*it1);
						}
						else{
							delete *it1;
						}
					}
					else{ //If nothing changed, but it is still a valid PAH, place it in overflow vector
						overflow.push_back(*it1);
					}
				}

			}
			else{
				updatetime = sys.Particles().Simulator()->updatePAH((*it)->m_pahstruct, (*it)->lastupdated, growtime, 1, 0,
					rng, growthfact, (*it)->PAH_ID, true, 1.0);

				(*it)->lastupdated = t;

				//! Invalidate PAH.
				/*!
				* A PAH is made invalid by setting the number of carbons to 5. The
				* reason for only applying this to particles is that we would not
				* want to remove single PAHs in the gas-phase which fall below the
				* minimum number of rings for inception but are still growing.
				*/
				if ((*it)->m_pahstruct->numofRings() < thresholdOxidation && m_numPAH >= minPAH && ind != -1){
					(*it)->m_pahstruct->setnumofC(5);
				}


				//! See if anything changed, as this will required a call to UpdatePrimary() below.
				if (oldNumCarbon != (*it)->m_pahstruct->numofC() || oldNumH != (*it)->m_pahstruct->numofH() && ind != -1)
				{
					m_PAHclusterchanged = true;
					m_PAHchanged = true;
				}

			}

			/*!
             * The second condition ensures that the m_InvalidPAH is modified in the correct way
			 * consider 2 PAH, the first is invalid, the other is valid. If not using the
             * second condition, the m_InvalidPAH will be false enventually, but it should be true.
			 */ 
			if (m_PAHchanged && !m_InvalidPAH && ind != -1)
				m_InvalidPAH = CheckInvalidPAHs(*it);

        } //end of loop for PAHs in a primary
        if (m_InvalidPAH)
        {
            RemoveInvalidPAHs();
        }

        //sys.Particles().Add();
        //this->ParticleModel()->Mode();
        /*!
         * Calculate derived quantities such as collision diameter and surface
         * area by iterating through all the PAHs.  This call is rather expensive.
         */
        if(m_PAHclusterchanged)
            UpdatePrimary();
            // if (m_InceptedPAH!=0) {
            //        sys.Particles().SetNumOfInceptedPAH(-1);
            //        //sys.Particles().NumOfInceptedPAH();
            // }
            //else if (m_InceptedPAH == 0 && InceptedPAH()==1){
            //        sys.Particles().SetNumOfInceptedPAH(1);
            //        //sys.Particles().NumOfInceptedPAH();
            //}
        //! otherwise there is no need to update.
    }
}

//! Overload function. Used if primary coordinates are tracked, then particle free surface area can be used to 
// calculate a free_surf_factor to describe surface growth.
/*!
* The actual interval over which the update is carried out on a PAH is from
* lastupdated to t.
*
* @param[in]        t        Time up to which to update.
* @param[in]        model    Particle model defining interpretation of particle data.
* @param[in]        sys      Cell containing particle and providing gas phase.
* @param[in]        fs       Free surface area of the particle.
* @param[in,out]    rng      Random number generator.
*
*/
void PAHPrimary::UpdatePAHs(const double t, const double dt, const Sweep::ParticleModel &model, Cell &sys, int statweight,
	int ind, rng_type &rng, PartPtrVector &overflow, const double fs)
{
	// Either the primary has two children or it is a leaf of the
	// tree
	if (m_leftchild != NULL)
	{
		// Recurse down to the leaves
		m_leftchild->UpdatePAHs(t, dt, model, sys, statweight, ind, rng, overflow, fs);
		m_rightchild->UpdatePAHs(t, dt, model, sys, statweight, ind, rng, overflow, fs);
	}
	else
	{
		// There are PAHs in this primary so update them, if needed
		// Flag to show if any PAH has been changed
		// this flag store the condition of this cluster(primary particle)
		bool m_PAHclusterchanged = false;
		// this flag store the condition of a particular PAH within this cluster
		bool m_InvalidPAH = false;
		// check whether this updated primary particle is a inceptedPAH, used later to track the num of InceptedPAH in the ensemble
		const int m_InceptedPAH = InceptedPAH();

		double m_vol_old = m_vol; //volume before PAH growth; used for adjust primary after surface growth

		double free_surf_factor = m_free_surf / fs; //if the free surface area of a primary particle is too small, there are less chance for PAHs in it to grow

		//assert(free_surf_factor <= 1.0); //test by hdy

		//! Loop over each PAH in this primary.
		const std::vector<boost::shared_ptr<PAH> >::iterator itEnd = m_PAH.end();
		for (std::vector<boost::shared_ptr<PAH> >::iterator it = m_PAH.begin(); it != itEnd; ++it) {

			bool m_PAHchanged = false;

			//! Model parameter.
			/*!
			* If a jump process reduces the total number of 6-member rings
			* (excludes 5-member rings) in a PAH (in a particle) below this
			* threshold it is removed.
			*/
			const int thresholdOxidation = model.Components(0)->ThresholdOxidation();

			//! PAH removal.
			/*!
			* Allow for the removal of PAHs in particles (2 or more PAHs) when
			* the total number of 6-member rings (excludes 5-member rings) in
			* the PAH has fallen below, for example, the minimum number of
			* rings in a PAH for inception. Particularly relevant when
			* oxidation is strong.
			*/
			if ((*it)->m_pahstruct->numofC() == 5){
				m_PAHclusterchanged = true;
				m_PAHchanged = true;
				m_InvalidPAH = CheckInvalidPAHs(*it);
				continue;
			}

			//! Model parameter.
			/*!
			* Defines when primary particles contain too many PAHs for them all
			* to have full access to the surrounding gas phase.  Once a primary
			* contains more than minPAH PAHs the growth rate of each PAH is
			* reduced according to the growth factor.
			*/
			const double minPAH = model.Components(0)->MinPAH();
			double growthfact = 1.0;
			growthfact = growthfact * free_surf_factor; 

			if (m_numPAH >= minPAH)
			{
				growthfact = model.Components(0)->GrowthFact();
				growthfact = growthfact * free_surf_factor; 
				//double density = model.Components(0)->Density();
				////! PP mass (kg).
				//double m_mass = m_numcarbon*1.9945e-26 + m_numH*1.6621e-27;

				////! Spherical diameter (nm).
				//double diameter = pow(6.0 * m_mass / density / PI, ONE_THIRD);
				//diameter = diameter * 1.0e9;
				//growthfact = 1.950 / diameter ;
				//if (growthfact > 1.0)
				//{
				//	growthfact = 1.0;
				//}
			}

			//! Time for one particular PAH to grow.
			double growtime = t - (*it)->lastupdated;
			assert(growtime >= 0.0);
			double statweightold = statweight;
			int numloops = 0;
			bool calcrates = true;
			double ratefactor = 1.0;
			const int oldNumCarbon = (*it)->m_pahstruct->numofC();
			const int oldNumH = (*it)->m_pahstruct->numofH();

			double updatetime;

			//if this is a particle with a single PAH, it may be weighted. 
			//If it is weighted and IWDSA is being used, we do not want to update the PAH, but rather update a clone of 
			//that PAH and create a new particle
			if (m_PAH.size() == 1 && statweight > 1.0 && ParticleModel()->Components(0)->WeightedPAHs()){
				PartPtrVector overflowtemp;

				while (growtime > 0.0 && statweight > 1.0){

					boost::shared_ptr<PAH> new_m_PAH((*it)->Clone());

					new_m_PAH->PAH_ID = ID;
					ID++;

					if (numloops > 0){
						calcrates = false;
					}
					else{
						calcrates = true;
					}

					updatetime = sys.Particles().Simulator()->updatePAH(new_m_PAH->m_pahstruct, (*it)->lastupdated, growtime, 1, 1,
						rng, growthfact*statweight, new_m_PAH->PAH_ID, calcrates, ratefactor);

					/*updatetime = sys.Particles().Simulator()->updatePAH(new_m_PAH->m_pahstruct, (*it)->lastupdated, growtime, 1, 1,
					rng, growthfact*statweight, new_m_PAH->PAH_ID, true, 1.0);*/

					new_m_PAH->lastupdated = updatetime;
					(*it)->lastupdated = updatetime;

					growtime = t - (*it)->lastupdated;
					numloops++;

					//! Invalidate PAH.
					/*!
					* A PAH is made invalid by setting the number of carbons to 5. The
					* reason for only applying this to particles is that we would not
					* want to remove single PAHs in the gas-phase which fall below the
					* minimum number of rings for nception but are still growing.
					*/
					if (new_m_PAH->m_pahstruct->numofRings() < thresholdOxidation){
						statweight--;
						(*sys.Particles().At(ind)).setStatisticalWeight(statweight);
						new_m_PAH.reset();
						ID--;
					}

					//! See if anything changed, as this will required a call to UpdatePrimary() below.
					//This will also dictate if the new PAH is used to create a new particle, or if it is just destroyed
					else if (growtime > 0.0)
						//else
					{
						//Reduce statistical weight of the particle being updated
						statweight--;
						(*sys.Particles().At(ind)).setStatisticalWeight(statweight);

						Particle *sp = NULL;
						sp = model.CreateParticle(updatetime);
						AggModels::PAHPrimary *pri =
							dynamic_cast<AggModels::PAHPrimary*>((*sp).Primary());
						pri->m_PAH.clear();
						pri->m_PAH.push_back(new_m_PAH);
						sp->SetTime(updatetime);
						//Update the primary and the cache
						pri->UpdatePrimary();
						sp->UpdateCache();

						overflowtemp.push_back(sp);

					}
					else{
						assert(growtime == 0);
						new_m_PAH.reset();
						ID--;
					}

					ratefactor = statweight / statweightold;
				}
				//Final update after statistical weight reaches 1
				//Now update the PAH one last time
				if (growtime > 0.0){

					updatetime = sys.Particles().Simulator()->updatePAH((*it)->m_pahstruct, (*it)->lastupdated, growtime, 1, 0,
						rng, growthfact, (*it)->PAH_ID, true, 1.0);

					(*it)->lastupdated = t;

					//! Invalidate PAH.
					/*!
					* A PAH is made invalid by setting the number of carbons to 5. The
					* reason for only applying this to particles is that we would not
					* want to remove single PAHs in the gas-phase which fall below the
					* minimum number of rings for inception but are still growing.
					*/
					if ((*it)->m_pahstruct->numofRings() < thresholdOxidation){
						(*it)->m_pahstruct->setnumofC(5);
					}


					//! See if anything changed, as this will required a call to UpdatePrimary() below.
					if (oldNumCarbon != (*it)->m_pahstruct->numofC() || oldNumH != (*it)->m_pahstruct->numofH())
					{
						m_PAHclusterchanged = true;
						m_PAHchanged = true;
					}

				}
				//Now, update all created particles that were added to overflowtemp
				PartPtrVector::iterator it1;
				for (it1 = overflowtemp.begin(); it1 != overflowtemp.end(); ++it1){
					AggModels::PAHPrimary *pri =
						dynamic_cast<AggModels::PAHPrimary*>((*(*it1)).Primary());
					//This new particle must also be updated to time t
					pri->UpdatePAHs(t, t - (*it1)->LastUpdateTime(), model, sys, 1, -1, rng, overflow, fs);
					(*it1)->SetTime(t);

					//Update the primary and the cache
					pri->UpdatePrimary();
					(*it1)->UpdateCache();

					//Check if the PAH is still valid after being updated
					if (pri->m_PAH[0]->m_pahstruct->numofRings() >= thresholdOxidation){
						overflow.push_back(*it1);
					}
					else
					{
						delete *it1;
					}
				}

			}
			else{
				updatetime = sys.Particles().Simulator()->updatePAH((*it)->m_pahstruct, (*it)->lastupdated, growtime, 1, 0,
					rng, growthfact, (*it)->PAH_ID, true, 1.0);

				(*it)->lastupdated = t;

				//! Invalidate PAH.
				/*!
				* A PAH is made invalid by setting the number of carbons to 5. The
				* reason for only applying this to particles is that we would not
				* want to remove single PAHs in the gas-phase which fall below the
				* minimum number of rings for inception but are still growing.
				*/
				if ((*it)->m_pahstruct->numofRings() < thresholdOxidation && m_numPAH >= minPAH && ind != -1){
					(*it)->m_pahstruct->setnumofC(5);
				}


				//! See if anything changed, as this will required a call to UpdatePrimary() below.
				if (oldNumCarbon != (*it)->m_pahstruct->numofC() || oldNumH != (*it)->m_pahstruct->numofH() && ind != -1)
				{
					m_PAHclusterchanged = true;
					m_PAHchanged = true;
				}

			}

			/*!
			* The second condition ensures that the m_InvalidPAH is modified in the correct way
			* consider 2 PAH, the first is invalid, the other is valid. If not using the
			* second condition, the m_InvalidPAH will be false enventually, but it should be true.
			*/
			if (m_PAHchanged && !m_InvalidPAH && ind != -1)
				m_InvalidPAH = CheckInvalidPAHs(*it);

		} //end of loop for PAHs in a primary
		if (m_InvalidPAH)
		{
			RemoveInvalidPAHs();
		}

		//sys.Particles().Add();
		//this->ParticleModel()->Mode();
		/*!
		* Calculate derived quantities such as collision diameter and surface
		* area by iterating through all the PAHs.  This call is rather expensive.
		*/
		if (m_PAHclusterchanged) {
			UpdatePrimary();
			if (m_pmodel->getTrackPrimaryCoordinates() || m_pmodel->getTrackPrimarySeparation())
				Adjust(m_vol_old);
			// if (m_InceptedPAH!=0) {
			//        sys.Particles().SetNumOfInceptedPAH(-1);
			//        //sys.Particles().NumOfInceptedPAH();
			// }
			//else if (m_InceptedPAH == 0 && InceptedPAH()==1){
			//        sys.Particles().SetNumOfInceptedPAH(1);
			//        //sys.Particles().NumOfInceptedPAH();
			//}
		}
		//! otherwise there is no need to update.
	}
}


// currently, only A1,A2 and A4 can be specified as InceptedPAH due to the underlying KMC code
bool PAHPrimary::CheckInvalidPAHs(const boost::shared_ptr<PAH> & it) const
{
    ParticleModel::PostProcessStartingStr str = this->ParticleModel()->InceptedPAH();
    int m_control;
    switch (str){
    case ParticleModel::A1:
        m_control=Sweep::KMC_ARS::BENZENE_C;
        break;
    case ParticleModel::A2:
        m_control=10;
        break;
    case ParticleModel::A4:
        m_control=Sweep::KMC_ARS::PYRENE_C;
        break;
    case ParticleModel::A5:
        m_control=Sweep::KMC_ARS::BENZOPYRENE_C;
        break;
    default:
        throw std::runtime_error("no information about the incepted PAH is available (Sweep::PAHPrimary::CheckInvalidPAHs())");
    }
    // if the PAH in the cluster is as the same size as the incepted PAH, it will be released but the current implementation is directly removed which causes mass loss, for a fully coupled model this part should be redesigned.
	return (it->m_pahstruct->numofC() < m_control || (it->m_pahstruct->numofC() <= m_control && NumPAH() != 1));
}

//struct compare_class
//{
//    bool operator () (const boost::shared_ptr<PAH>  & it) const
//    {
//        if (m_pyreneInception) 
//        {
//            if (it->Structure()->numofC() < Sweep::KMC_ARS::PYRENE_C)
//                return true;
//        }
//        else 
//        {
//            if (it->Structure()->numofC() < Sweep::KMC_ARS::BENZENE_C)
//                return true;
//        }
//        return false;
//    }
//};

void PAHPrimary::RemoveInvalidPAHs()
{
    // if the num of PAH in this primary particle is smaller than 2 after removing the InvalidPAHs,
    // FakeRounding() will return true and it will be merged with other part of the aggregate.
    // TODO: thinking suitable solution for moving the InvalidPAH back to gasphase when coupling with the gasphase chemistry.
    // current implementation is only suitable for post-pocessing due to mass loss.
    //remove_if(m_PAH.begin(),m_PAH.end(),compare_class());
    std::vector<boost::shared_ptr<PAH> >::iterator NewEnd = remove_if(m_PAH.begin(),m_PAH.end(),boost::bind(&PAHPrimary::CheckInvalidPAHs, this,_1));
    m_PAH.resize(NewEnd-m_PAH.begin());
}

// Find Xmer (monomer, dimer or trimer) and record its mass which will be used to create mass spectra
// however, this function currently is only implenmented in PAHPrimary class
// this means it is limited to PAH-PP model. 
void PAHPrimary::FindXmer(std::vector<double> &out, int m_xmer) const
{
    if (m_leftchild!=NULL)
        m_leftchild->FindXmer(out, m_xmer);

	// dimer: a partilce with one primary containing 2 PAHs
	if (this->NumPAH()==m_xmer)
		out.push_back(this->MassforXmer());
}

// Find particular primary particle with target number of PAHs ( this can be a range, details please see the implementation) 
// and record num of C and H for each PAH, also push_back a divider (0,0,ID) at end to distinguish 
void PAHPrimary::FindXmer(std::vector<std::vector<double> > &out, int target_num_PAH) const
{
	if (m_leftchild != NULL)
		m_leftchild->FindXmer(out,target_num_PAH);
	if (m_rightchild != NULL)
		m_rightchild->FindXmer(out,target_num_PAH);
    if (target_num_PAH-10<=0) std::cout<<"Warning: (target_num_PAH-10) return 0 or negative value in PAHPrimary::FindXmer()"<<std::endl;
	if (m_PAH.size() <= size_t(target_num_PAH+10) && m_PAH.size() >= size_t(target_num_PAH-10))
		mass_PAH(out);
}

//calculate the mass of Xmer including C and H
double PAHPrimary::MassforXmer() const
{ 
	int sum = 0;
	for (size_t i = 0; i != m_PAH.size(); ++i)
	{
		sum+=m_PAH[i]->m_pahstruct->numofH();
	    sum+=12 * m_PAH[i]->m_pahstruct->numofC();
	}
	return sum;	
}

int PAHPrimary::InceptedPAH() const
{
    // m_parent == NULL is to check whether this primary particle is part of an aggregate.
    if (m_parent == NULL && Numprimary() == 1 && NumPAH() == 1){
        //currently only Num of C and H is used to identify the Pyrene, Naphthalene and benzene
        ParticleModel::PostProcessStartingStr str = ParticleModel()->InceptedPAH();
        switch (str){
        case ParticleModel::A1:
            if (NumCarbon() == BENZENE_C && NumHydrogen() == BENZENE_H)
                return 1;
            else return 0;
            break;
        case ParticleModel::A2:
            if (NumCarbon() == NAPHTHALENE_C && NumHydrogen() == NAPHTHALENE_H)
                return 1;
            else return 0;
            break;
        case ParticleModel::A4:
            if (NumCarbon() == PYRENE_C && NumHydrogen() == PYRENE_H)
                return 1;
            else return 0;
            break;
        case ParticleModel::A5:
            if (NumCarbon() == BENZOPYRENE_C && NumHydrogen() == BENZOPYRENE_H)
                return 1;
            else return 0;
            break;
        default:
            return 0;
        }
    }
    else return 0;
}

// dump information of this Xmer to a vector<vector<double> > 
 void PAHPrimary::mass_PAH(std::vector<std::vector<double> > &out) const
 {
    std::vector<double> temp;
    std::vector<double> divider(2,0);
    for (size_t i = 0; i != m_PAH.size(); ++i)
    {
        temp.push_back(m_PAH[i]->m_pahstruct->numofC());
        temp.push_back(m_PAH[i]->m_pahstruct->numofH());
        m_PAH[i]->saveDOTperLoop((int)ID,(int)i);
        out.push_back(temp);
        temp.clear();
    }
    divider.push_back(ID);
    out.push_back(divider);
    ID++;
}
 // only for num of Primary == 1, carbon only
double PAHPrimary::ReducedMass() const
{
    double val=0;
    for (size_t i = 0; i != m_PAH.size(); ++i)
    {
        int num_C=0;
        num_C=m_PAH[i]->m_pahstruct->numofC();
        val+=1/num_C;
    }
    if (val==0) return 1;
    else 
    return 1/val;
}
bool IsSticked(double val)
{
    return val<2*(32*12+14);
}
bool NotValid(double val)
{
    return val==0;
}

void PAHPrimary::Fragtest(std::vector<double> &out, const int k, std::string mode, double threshold) const
{
    //test8 start
    //int thres=0;
    //if (mode =="MIN"||mode=="MAX")
    //    thres = threshold/2;
    //else if (mode =="COMBINED")
    //    thres = threshold/4;
    //else if (mode =="REDUCED")
    //    thres = threshold;
    //if (m_leftchild!=NULL)
    //m_leftchild->Fragtest(out, k, mode, threshold);

    //std::vector<double> temp;

    ////if (k==0 && ReducedMass()<=thres) {
    ////    for (size_t i = 0; i != m_PAH.size(); i+=1)
    ////        {
    ////            int num_C=0;
    ////            int num_H=0;
    ////            int val=0;
    ////            int val1=0;
    ////            num_C=m_PAH[i]->m_pahstruct->numofC();
    ////            num_H=m_PAH[i]->m_pahstruct->numofH();
    ////            // PAH mass (u)
    ////            val = 12*num_C + num_H;
    ////            temp.push_back(val);
    ////        }
    ////}
    ////else if (k==1 && ReducedMass()>thres) {
    
    //if (k==1 && NumPAH()<=10) {
    //    for (size_t i = 0; i != m_PAH.size(); i+=2)
    //    {
    //        if (i+1>=m_PAH.size())
    //            break;
    //        int num_C=0;
    //        int num_H=0;
    //        int val=0;
    //        int val1=0;
    //        num_C=m_PAH[i]->m_pahstruct->numofC();
    //        num_H=m_PAH[i]->m_pahstruct->numofH();
    //        // PAH mass (u)
    //        val = 12*num_C + num_H;
    //        num_C=m_PAH[i+1]->m_pahstruct->numofC();
    //        num_H=m_PAH[i+1]->m_pahstruct->numofH();
    //        // PAH mass (u)
    //        val1 = 12*num_C + num_H;
    //        val+=val1;
    //        temp.push_back(val);
    //    }
    //}
    //if (temp.size()!=0) {
    //    std::vector<double>::iterator NewEnd = remove_if(temp.begin(),temp.end(),IsSticked);
    //    temp.resize(NewEnd-temp.begin());
    //    out.insert(out.end(),temp.begin(),temp.end());
    //}
    //test8 end
    //test9 start
    std::vector<double> temp;
    if (k==1 && NumPAH()<=5) {
        mass_PAH(temp);
    }
    if (temp.size()!=0)
    {
        for (size_t i = 0; i != static_cast<size_t>(NumPAH()); ++i){
            //c32h14
            if (temp[i]<=398)
                temp[i]=0;
        }
        std::vector<double>::iterator NewEnd = remove_if(temp.begin(),temp.end(),NotValid);
        temp.resize(NewEnd-temp.begin());
        for (size_t i = 0; i !=temp.size(); i+=2) {
            if (i+1>=temp.size()) {
                temp[i]=0;
                break;
            }
            temp[i]=temp[i]+temp[i+1];
            temp[i+1]=0;
        }
        NewEnd = remove_if(temp.begin(),temp.end(),NotValid);
        temp.resize(NewEnd-temp.begin());
        out.insert(out.end(),temp.begin(),temp.end());
    }
    //test9 end

}

/*!
 * @brief Create contents pertaining to PAH specific information to be written to a csv file. 
 *
 * @param[in,out]    out                   Vector (different PAHs) of vectors (different pieces of information about the PAH).
 * @param[in]        index                 Index assigned to particle.
 * @param[in]        density               Density of soot.
 * @param[in,out]    pahUniqueAddresses    Keep a record of which PAHs have been stored in "out" to avoid storing the same memory location twice.
 * @param[in,out]    Mapping               Map of PAH memory locations to PAH in "out".
 * @param[in]        timeStep              Index assigned to time step.
 */ 
void PAHPrimary::OutputPAHPSL(std::vector<std::vector<double> > &out, const int index, const double density, std::set<void*> &pahUniqueAddresses, std::vector<std::string> &Mapping, const double timeStep) const
{
    if (m_leftchild!=NULL)
        m_leftchild->OutputPAHPSL(out, index, density, pahUniqueAddresses, Mapping, timeStep);
    if (m_rightchild!=NULL) m_rightchild->OutputPAHPSL(out, index, density, pahUniqueAddresses, Mapping, timeStep);

	std::stringstream memoryLocation;

    std::vector<double> temp;

    for (size_t i = 0; i != m_PAH.size(); ++i)
    {
        //! Reduce the number of redundant entries.
        /*!
		 * First retrieve the memory location of this PAH. If the memory
		 * location is not in the list of unique memory locations, store the
		 * information for this PAH in "out". Otherwise, we can simply
		 * increment the frequency of the PAH in "out" which shares the same
		 * memory location as this PAH.
		 */
        void *ptr = m_PAH[i].get();
        if(pahUniqueAddresses.find(m_PAH[i].get()) == pahUniqueAddresses.end()) {
            
            //! Initialization of variables.
			int num_C=0;
			int num_H=0;
			double val=0.0;
			double m_mass=0.0;
			double PAHCollDiameter=0.0;
			double diameter=0.0;

            //! If this particle is a single primary with a single PAH it is assigned an index of -1 to distinguish it from PAHs in particles.
			if (this->NumPAH() == 1 && this->Numprimary() == 1) {
				temp.push_back(-1);
			}
			else {
				temp.push_back(index);
			}
				
            //! Number of carbon atoms.
			num_C=m_PAH[i]->m_pahstruct->numofC();
			temp.push_back(num_C);

            //! Number of hydrogen atoms.
			num_H=m_PAH[i]->m_pahstruct->numofH();
			temp.push_back(num_H);

            //! Number of 6-member rings.
			temp.push_back(m_PAH[i]->m_pahstruct->numofRings());
            
            //! Number of 5-member rings.
			temp.push_back(m_PAH[i]->m_pahstruct->numofRings5());
			
            //! Number of carbon atoms on the edge of the PAH.
            temp.push_back(m_PAH[i]->m_pahstruct->numofEdgeC());

			//! PAH mass (u).
			val = 12*num_C + num_H;
			temp.push_back(val);

			//! PAH mass (kg).
			m_mass = num_C*1.9945e-26 + num_H*1.6621e-27;
			temp.push_back(m_mass);

			//! PAH collision diameter (m).
			PAHCollDiameter = sqrt(num_C*2.0/3.);
			PAHCollDiameter *= 2.4162*1e-10;    //! convert from Angstrom to m.
			temp.push_back(PAHCollDiameter);

			//! PAH density (kg/m3).
			temp.push_back(density);

			//! PAH volume (m3).
			val = m_mass / density;
			temp.push_back(val);

			//! Spherical diameter (m).
			diameter = pow(6.0 * val / PI, ONE_THIRD);
			temp.push_back(diameter);

			//! Larger of the spherical or collison diameter (m).
			val = max(diameter, PAHCollDiameter);
			temp.push_back(val);

			//! Time created (s).
			val = m_PAH[i]->time_created;
			temp.push_back(val);

			//! Index of PAH.
			val = m_PAH[i]->PAH_ID;
			temp.push_back(val);

            //! Uncomment the call to saveDOTperLoop to print out the structure of each PAH.
			//m_PAH[i]->saveDOTperLoop(timeStep,uniquePAHCounter);

			uniquePAHCounter = uniquePAHCounter + 1;

			//! Number of PAHs pointing to the same memory location.
			temp.push_back(1);

			out.push_back(temp);

			temp.clear();

			pahUniqueAddresses.insert(m_PAH[i].get());

			memoryLocation << m_PAH[i].get() << std::endl;

			Mapping.push_back(memoryLocation.str());

			memoryLocation.str("");
        }
        else {
			memoryLocation << m_PAH[i].get() << std::endl;

			int pos = find(Mapping.begin(), Mapping.end(), memoryLocation.str()) - Mapping.begin();

			memoryLocation.str("");

			out[pos].back() = out[pos].back() + 1;
        }
    }
}

/*!
* @brief Create contents pertaining to primary particle specific information to be written to a csv file.
*
* @param[in,out]    out                   Vector (different PAHs) of vectors (different pieces of information about the PAH).
* @param[in]        index                 Index assigned to particle.
* @param[in]        density               Density of soot.
* @param[in,out]    pahUniqueAddresses    Keep a record of which PAHs have been stored in "out" to avoid storing the same memory location twice.
* @param[in,out]    Mapping               Map of PAH memory locations to PAH in "out".
* @param[in]        timeStep              Index assigned to time step.
*/
void PAHPrimary::OutputPPPSL(std::vector<std::vector<double> > &out, const int index, const double density, const double timeStep) const
{
	if (m_leftchild != NULL)
		m_leftchild->OutputPPPSL(out, index, density, timeStep);
	if (m_rightchild != NULL) m_rightchild->OutputPPPSL(out, index, density, timeStep);

	std::vector<double> temp;

	//! Initialization of variables.
	double val;
	double m_mass = 0.0;
	double PPCollDiameter = 0.0;
	double diameter = 0.0;

	//! If this particle is a single primary with a single PAH it is assigned an index of -1 to distinguish it from PAHs in particles.
	if (this->NumPAH() == 1 && this->Numprimary() == 1) {
		temp.push_back(-1);
	}
	else {
		temp.push_back(index);
	}

	//! Number of carbon atoms.
	temp.push_back(m_numcarbon);

	//! Number of hydrogen atoms.
	temp.push_back(m_numH);

	//! Number of 6-member rings.
	temp.push_back(m_numOfRings);

	//! Number of 5-member rings. //NICK - TO DO
	temp.push_back(0);

	//! Number of PAHs
	temp.push_back(m_numPAH);

	//! PP mass (kg).
	m_mass = m_numcarbon*1.9945e-26 + m_numH*1.6621e-27;
	temp.push_back(m_mass);

	//! PP volume (m3).
	val = m_mass / density;
	temp.push_back(m_mass / density);

	//! Spherical diameter (nm).
	diameter = pow(6.0 * val / PI, ONE_THIRD);
	temp.push_back(diameter*1.0e9);

	out.push_back(temp);

	temp.clear();
}

// this function is only used to create a vector containing all the mass of individual PAH within this soot particle
void PAHPrimary::mass_PAH(std::vector<double> &out) const
{
    if (m_leftchild != NULL)
        m_leftchild->mass_PAH(out);
    if (m_rightchild != NULL)
        m_rightchild->mass_PAH(out);

    double temp_mass=0.0;
    for (size_t i = 0; i != m_PAH.size(); ++i)
    {
        temp_mass = 12*m_PAH[i]->m_pahstruct->numofC() + m_PAH[i]->m_pahstruct->numofH();
        out.push_back(temp_mass);
    }
}
void PAHPrimary::UpdateCache(void)
{
    UpdateCache(this);
}

bool PAHPrimary::FakeRounding()
{
    // there are two conditions that it should perform FakeRounding, no PAH or only one PAH in this primary particle, thus it should be merged with other primary particle.
    if (m_leftparticle!=NULL)
    {
        if (m_leftparticle->m_numPAH==1||(m_leftparticle->m_numPAH==0&&m_leftparticle->m_numcarbon==0))
            return true;
        else if (m_rightparticle->m_numPAH==1||(m_rightparticle->m_numPAH==0&&m_rightparticle->m_numcarbon==0))
            return true;
        else return false;
    }
    else return false;
}

bool PAHPrimary::CheckRounding()
{
    bool hascoalesced = false;
    bool Condition;
    
    //! The condition for whether a particle has coalesced depends on whether
    //! the distance between the centres of primary particles is tracked. If
    //! tracked, a particle has coalesced if the distance is 0. If not, the
    //! condition depends on whether the rounding level exceeds an arbitrarily
    //! high threshold.
	if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()) {
        Condition = (m_children_roundingLevel > 0.95);
    } else {
        Condition = (m_distance_centreToCentre == 0.0);
    }

    if ((Condition && m_leftparticle != NULL) || FakeRounding()) {
		//if ((Condition && m_leftparticle != NULL)) {
        // PrintTree("before.inp");
        // cout <<"merging"<<m_children_roundingLevel<<endl;
        Merge();
        // PrintTree("after.inp");

        hascoalesced = true;

        //! Check again because this node has changed.
        CheckRounding();
    }

    if (m_leftchild != NULL) {
        hascoalesced = m_leftchild->CheckRounding();
        hascoalesced = m_rightchild->CheckRounding();
    }

    UpdateCache();

    return hascoalesced;
}


//! Update primary particle.
void PAHPrimary::UpdatePrimary(void)
{
	//! If the vector of boost shared pointers to PAHs (m_PAH) is empty, the
	//! primary particle is invalid and the following member variables should
	//! be set to 0.
	if (m_PAH.empty())
	{
		m_mass = 0.0;
		m_numcarbon = 0;
		m_numH = 0;
		m_numOfEdgeC = 0;
		m_numOfRings = 0;
		m_numOfRings5 = 0;
		m_numPAH = m_PAH.size();
		m_PAHmass = 0.0;
		m_PAHCollDiameter = 0.0;
	}
	else
	{
		m_numcarbon = 0;
		m_numH = 0;
		m_numOfEdgeC = 0;
		m_numOfRings = 0;
		m_numOfRings5 = 0;
		m_numPAH = m_PAH.size();
		m_PAHmass = 0.0;
		m_PAHCollDiameter = 0.0;

		//! Initialisation of variables to adjust the primary diameter if the
		//! distance between the centres of primary particles is tracked.
		double d_ij = 0.0;               //!< Distance between the centres of primary particles i and j.
		double r_i = 0.0;                //!< Radius of primary particle i.
		double r_j = 0.0;                //!< Radius of primary particle j.
		double x_i = 0.0;                //!< The distance from the centre of primary particle i to the neck level.
		double A_n = 0.0;                //!< Cross-sectional neck area.
		double A_i = 0.0;                //!< Free surface area of primary particle i.

		//! Initialisation of variables but this is only relevant to particles
		//! with more than one primary.
		//d_ij = m_parent->m_distance_centreToCentre;
		//r_i = m_primarydiam / 2.0;

		//if (m_parent->m_leftparticle == this) {
		//    r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
		//} else {
		//    r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
		//}

		//x_i = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / 2.0 / d_ij; //!< Eq. (3b) of Langmuir 27:6358 (2011).
		//A_n = M_PI * (pow(r_i, 2.0) - pow(x_i, 2.0));                        //!< Eq. (4).
		//A_i = 2.0 * M_PI * (pow(r_i, 2.0) + r_i * x_i);                      //!< Eq. (6).
		//m_primary_diam_old = m_primarydiam;
		//m_vol_old = m_vol;
		//}

		int maxcarbon = 0;

		for (vector<boost::shared_ptr<PAH> >::iterator i = m_PAH.begin(); i != m_PAH.end(); ++i) {
			m_numcarbon += (*i)->m_pahstruct->numofC();
			m_numH += (*i)->m_pahstruct->numofH();
			m_numOfEdgeC += (*i)->m_pahstruct->numofEdgeC();
			m_numOfRings += (*i)->m_pahstruct->numofRings();
			m_numOfRings5 += (*i)->m_pahstruct->numofRings5();
			maxcarbon = max(maxcarbon, (*i)->m_pahstruct->numofC()); //!< Search for the largest PAH-in terms of the number of carbon atoms-in the primary.
		}
		m_numOf6Rings = m_numOfRings;

		m_PAHmass = m_numcarbon * 1.9945e-23 + m_numH * 1.6621e-24;  //!< Units of g.
		m_PAHmass *= 1.0e-3;                                         //!< Units of kg.

		//! Eq. (10.19) in M. Frenklach, H. Wang, Detailed mechanism and
		//! modeling of soot particle formation, in: H. Bockhorn (Ed.), Soot
		//! Formation in Combustion-Mechanisms and Models, Springer, Berlin,
		//! 1994, pp. 165-190.
		//! Note that m_i in Eq. (10.19) refers to the number of carbon atoms.
		m_PAHCollDiameter = 1.395 * sqrt(3.0) * sqrt(maxcarbon * 2.0 / 3.0); //!< Units of Angstroms.
		m_PAHCollDiameter *= 1.0e-10;                             //!< Units of m.

		//! At the moment we have only one component: soot.
		if (m_pmodel->ComponentCount() != 1) {
			throw std::runtime_error("Model contains more then one component. Only soot is supported. (PAHPrimary::UpdatePrimary)");
		}

		m_vol = m_PAHmass / m_pmodel->Components(0)->Density(); //!< Units of m^3.
		m_mass = m_PAHmass;
		m_diam = pow(6.0 * m_vol / PI, ONE_THIRD);
		m_dmob = m_diam;
		m_dcol = max(m_diam, m_PAHCollDiameter); 
		//m_dcol = m_diam; //to test bintree model and PAH-KMC model
		m_surf = PI * m_diam * m_diam;

		//! If the distance between the centres of primary particles is
		//! tracked, the rate of change in the primary diameter is determined
		//! by its neighbour. Therefore, the particle should be made up of more
		//! than one primary.
		//if ((!m_pmodel->getTrackPrimarySeparation()) || (!m_pmodel->getTrackPrimaryCoordinates()) || (m_pmodel->getTrackPrimarySeparation() && m_numprimary == 1) || m_parent == NULL) {
		//m_primarydiam = m_diam; 
		//} //else { 
		//! Differentiating Eq. (3b) with respect to time, and assuming
		//! that r_j and d_ij do not change, the rate of change in x_i with
		//! respect to time can be obtained.
		//! Substituting Eqs. (4) and (6) and the above result into
		//! Eq. (2), the rate of change in r_i with respect to time can be
		//! obtained.
		//! References to equations are to Langmuir 27:6358 (2011).
		//!
		//! @todo Remove derivation and replace with reference to preprint
		//!       or paper if results do get published.
		//m_primarydiam = m_primary_diam_old + 2 * (m_vol - m_vol_old) / (A_i + A_n * r_j / d_ij);
		//} 

		//learn from bintree model//
		if (!(m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()) || m_parent == NULL){
			m_primarydiam = m_diam;
			m_free_surf = m_surf;
			m_primaryvol = m_vol;
			m_sum_necks = 0.0;	
			m_sum_cap = 0.0; //used to calculate geometric volume of a particle.

		}
		else{
			//! Update overlapping primary model
			UpdateOverlappingPrimary();
		}

		m_sph_prim_vol = (1.0 / 6.0) * M_PI * pow(m_primarydiam, 3.0); // used to calculate geometric volume of a particle.
		//m_numprimary = 1;
		m_avg_coalesc = 0.0;

		if (m_pmodel->getTrackPrimaryCoordinates()) {
			setRadius(m_primarydiam / 2.0);
		}
	}
}

void PAHPrimary::Reset()
{
	m_numcarbon = 0;
	m_numH = 0;
	m_numOfEdgeC = 0;
	m_numOfRings = 0;
	m_primarydiam = 0.0;
	m_surf = 0;
	m_vol = 0;
	m_PAH.clear();
	if (m_pmodel->getTrackPrimaryCoordinates() || m_pmodel->getTrackPrimarySeparation())
		m_avg_sinter = 0;
	else
		m_avg_coalesc = 0;
}

//! From bintree model//
/*!
* @brief       Checks the sintering level of the particle
*
* If the sintering level is above 95%, Merge is called and the cache
* is updated.
*
* @return      Boolean telling if the particle has sintered
*/
bool PAHPrimary::CheckSintering()
{
	bool hassintered = false;

	if (m_leftparticle != NULL) {

		// check whether condition for merger is met
		if (MergeCondition() || FakeRounding()) { //add FakeRounding() by hdy
		//if (MergeCondition()) { //used to test bintree model and PAH-KMC model
			Merge();
			UpdateCache();
			hassintered = true;

			// Check again because this node has changed
			CheckSintering();
		}
	}
	if (m_leftchild != NULL) {
		hassintered = m_leftchild->CheckSintering();
		hassintered = m_rightchild->CheckSintering();
	}

	return hassintered;
}

/*!
 * @param[in] root The root node of this particle
 */
void PAHPrimary::UpdateCache(PAHPrimary *root)
{
	//Update the children
	if (m_leftchild != NULL)
	{
		m_leftchild->UpdateCache(root);
		m_rightchild->UpdateCache(root);
		m_numprimary = m_leftchild->m_numprimary + m_rightchild->m_numprimary;
	}
	//this is a primary and the number of primaries below this node is one (this node and no children)
	else
	{
		if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates())
			m_avg_coalesc = 0;

		else{

			if (m_parent == NULL) m_avg_sinter = 1.0;
			else m_avg_sinter = 0.0;

			m_sum_cap = 0.0; //used to calculate geometric volume of a particle.
			UpdatePrimary(); //only update if primary coordinates are tracked; learned from bintree model.
		}

		m_numprimary = 1;

		// Check that the primary is begin kept up to date
		//            const double oldNumCarbons = m_numcarbon;
		//            UpdatePrimary();
		//            if(m_numcarbon != oldNumCarbons)
		//                std::cerr << "UpdatePrimary has changed num carbons inside UpdateCache\n";

	}

	//this is not a primary, sum up the properties
	if (m_leftchild != NULL)
	{
		// remove the PAHs from this node to free memory
		Reset();
		m_surf = m_leftchild->m_surf + m_rightchild->m_surf;
		m_primarydiam = (m_leftchild->m_primarydiam + m_rightchild->m_primarydiam);

		if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())
		{
			m_free_surf = m_leftchild->m_free_surf + m_rightchild->m_free_surf;
			m_primaryvol = m_leftchild->m_primaryvol + m_rightchild->m_primaryvol;
			m_sph_prim_vol = m_leftchild->m_sph_prim_vol + m_rightchild->m_sph_prim_vol; //used to calculate geometric volume of a particle.
		}
		m_vol = m_leftchild->m_vol + m_rightchild->m_vol;
		m_numPAH = m_leftchild->m_numPAH + m_rightchild->m_numPAH;
		m_mass = (m_leftchild->m_mass + m_rightchild->m_mass);
		m_PAHCollDiameter = max(m_leftchild->m_PAHCollDiameter, m_rightchild->m_PAHCollDiameter);

		m_numcarbon = m_leftchild->m_numcarbon + m_rightchild->m_numcarbon;
		m_numH = m_leftchild->m_numH + m_rightchild->m_numH;

		m_numOfEdgeC = m_leftchild->m_numOfEdgeC + m_rightchild->m_numOfEdgeC;
		m_numOfRings = m_leftchild->m_numOfRings + m_rightchild->m_numOfRings;
		m_numOfRings5 = m_leftchild->m_numOfRings5 + m_rightchild->m_numOfRings5;

		if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())
		{
			m_children_sumCap = CalcChildrenSumCap();//used to calculate geometric volume of a particle.
			if ((m_leftchild != NULL) && (m_rightchild != NULL))
			{
				m_sum_cap = m_children_sumCap +
					m_leftchild->m_sum_cap + m_rightchild->m_sum_cap;
			}
			//else m_sum_cap = m_children_sumCap;
		}

		if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates())
		{
			// calculate the coalescence level of the two primaries connected by this node
			m_children_roundingLevel = CoalescenceLevel();
			//sum up the avg coal level
			m_avg_coalesc = m_children_roundingLevel + m_leftchild->m_avg_coalesc + m_rightchild->m_avg_coalesc;
		}
		else{

			//calculate bounding sphere
			calcBoundSph();
			
			// calculate the coalescence level of the two primaries connected by this node
			m_children_sintering = SinteringLevel();

			if (MergeCondition()) CheckSintering();

			// Sum up the avg sintering level (now that sintering is done)
			if ((m_leftchild != NULL) && (m_rightchild != NULL)) {
				m_avg_sinter = m_children_sintering +
					m_leftchild->m_avg_sinter + m_rightchild->m_avg_sinter;		
			}
			else {
				// This should only occur if CheckSintering has merged
				m_avg_sinter = m_children_sintering;
			}

		}

		// calculate the different diameters only for the root node because this goes into the
		// particle tree and gets used by the coagulation kernel
		if (this == root)
		{
			//spherical eqiv radius
			double spherical_radius = pow(3 * m_vol / (4 * PI), ONE_THIRD);
			m_diam = 2 * spherical_radius;

			// there are m_numprimary-1 connections between the primary particles
			if (m_numprimary > 1 && !(m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()))
				// there are m_numprimary-1 connections between the primary particles
				m_avg_coalesc = m_avg_coalesc / (m_numprimary - 1);

			else
				m_avg_sinter = m_avg_sinter / (m_numprimary - 1);

			//approxmiate the surface of the particle

			if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()) {
				// Approxmiate the surface of the particle
				// (same as in ChangePointer)
				const double numprim_1_3 = pow(m_numprimary, -1.0 * ONE_THIRD);

				m_surf = 4 * PI*spherical_radius*spherical_radius /
					(m_avg_coalesc*(1 - numprim_1_3) + numprim_1_3);
			}
			else{
				// if the centre to centre distance is tracked then this is the free surface area
				m_surf = m_free_surf;
			}

			//calculate the surface equivalent radius
			// const double radius_surf=sqrt(m_surf/(4*PI));
			// the average between the surface and voluem equiv diameter
			//const double meandiam=spherical_radius+radius_surf;            //2*0.5*(radius(vol)+radius(sphere))
			const double aggcolldiam = (6 * m_vol / m_surf)*
				pow(pow(m_surf, 3) / (36 * PI*m_vol*m_vol), (1.0 / 1.8));
			// the maximum of the largest PAH diameter and
			// the average between the surface and voluem equiv diameter
			const double cdiam = max(aggcolldiam, m_PAHCollDiameter);
			//const double cdiam = aggcolldiam; //used to test bintree model and PAH-KMC model
			m_dmob = aggcolldiam;
			SetCollDiameter(cdiam);
		}
		else
		{
			m_diam = 0;
			m_dmob = 0;

		}

		//test the difference between geometric volume (geom_vol) and mass-density based volume (m_vol)
		//geom_vol should be very close to m_primaryvol
		//there is some problem to test the geom_vol and m_vol here, maybe because in the function UpdatePAHs, all the primary particles are updated first,
		//then UpdateCache() was called to update the whole particle
		//anyway, geom_vol and m_vol can pass the mass conservation test by post processing the data printed by PrintPrimary() using matlab script.
		//if (this == root && this->m_leftparticle != NULL && this->m_rightparticle != NULL)
		//{
		//	double geom_vol = 0.0;
		//	geom_vol = m_sph_prim_vol - m_sum_cap;

		//	//if (m_primaryvol > 0.0){
		//	//	if (abs((geom_vol - m_primaryvol) / m_primaryvol) > 0.01)
		//	//	{
		//	//		std::cout << "Something wrong with this particle!" << endl;
		//	//		std::cout << "Geometry volume = " << geom_vol << endl;
		//	//		std::cout << "Primary volume = " << m_primaryvol << endl;
		//	//	}
		//	//}

		//	if (abs((geom_vol - m_vol) / m_vol) > 0.01)
		//	{
		//		std::cout << "Something wrong with this particle!" << endl;
		//		std::cout << "Geometry volume = " << geom_vol << endl;
		//		std::cout << "Mass and density volume = " << m_vol << endl;
		//		std::cout << "Num of primary = " << m_numprimary << endl;
		//	}

		//	//assert(abs((geom_vol - m_vol) / m_vol) < 0.01);
		//	//assert(abs((geom_vol - m_primaryvol) / m_primaryvol) < 0.01);
		//}

	}
	
}

/*!
 * @param[in] filename Output filename
*/
void PAHPrimary::PrintTree(string filename)
{
  ofstream out;
  out.open(filename.c_str());
  out << "digraph unix {"<<endl;
  out <<"graph [rankdir = \"LR\"];"<<endl;
  PrintTreeLoop(out);
  out << "}"<<endl;
  out.close();
}

/*!
 * @param[in] out Output stream
*/
void PAHPrimary::PrintTreeLoop(std::ostream &out)
{
  if (m_leftchild!=NULL)
  { //out<<"leftchild "<<10E8*m_leftchild->SphDiameter()<<endl;
	//out<<"rightchild "<<10E8*m_rightchild->SphDiameter()<<endl;
      out<<"\" "<<this<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_sintering<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this->m_leftchild<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_sintering<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this->m_rightchild<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_sintering<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftchild<<"\"; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightchild<<"\"; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftparticle<<"\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightparticle<<"\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
	m_leftchild->PrintTreeLoop(out);
    m_rightchild->PrintTreeLoop(out);
  }

  else
  {
      //out<<"\" "<<this<<"\" "<<" [label = \""<<2*pow((this->vol_sinter)*3/(4*PI),ONE_THIRD)<<"\"];"<<endl;
      out<<"\" "<<this<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_sintering<<"|"<<this<<"\"];"<<endl;
  }
}


const PAHPrimary *PAHPrimary::RightChild() const
{
    return m_rightchild;
}

 const PAHPrimary *PAHPrimary::LeftChild() const
{
    return m_leftchild;
}

double PAHPrimary::PAHCollDiameter() const
{
    return m_PAHCollDiameter;
}

double PAHPrimary::Rg() const
{
    return m_Rg;
}

double PAHPrimary::Fdim() const
{
    return m_fdim;
}

double PAHPrimary::PrimaryDiam() const
{
    return m_primarydiam;
}

double PAHPrimary::LdivW() const
{
    return m_LdivW;
}

int PAHPrimary::Numprimary() const
{
    return m_numprimary;
}

int PAHPrimary::NumCarbon() const
{
    return m_numcarbon;
}

int PAHPrimary::NumHydrogen() const
{
    return m_numH;
}

int PAHPrimary::NumEdgeC() const
{
    return m_numOfEdgeC;
}

int PAHPrimary::NumRings() const
{
    return m_numOfRings;
}

int PAHPrimary::NumRings5() const
{
	return m_numOfRings5;
}

int PAHPrimary::NumPAH() const
{
    return m_numPAH;
}

double PAHPrimary::sqrtLW() const
{
    return m_sqrtLW;
}

double PAHPrimary::AvgCoalesc() const
{
    return m_avg_coalesc;
}

double PAHPrimary::AvgSinter() const
{
	if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())
		return m_avg_sinter;
	else
		return m_avg_coalesc;
}

//! Return distance between the centres of primary particles.
double PAHPrimary::Distance() const
{
    return m_distance_centreToCentre;
}

//! Return free surface area of primary particles.
double PAHPrimary::GetFreeSurfArea() const
{
	return m_free_surf;
}

// READ/WRITE/COPY.

// Returns a copy of the model data.
PAHPrimary *const PAHPrimary::Clone(void) const
{
    m_clone=true;
	PAHPrimary* newPAHPrimary=new PAHPrimary();
	newPAHPrimary->CopyParts(this);
	if (this->m_leftchild!=NULL)
	newPAHPrimary->CopyTree(this);
	m_clone=false;
	return newPAHPrimary;
}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType PAHPrimary::AggID(void) const {return AggModels::PAH_KMC_ID;}

/*
 * @brief Writes a object to a binary stream
 *
 * @param[in,out]    out                 Output binary stream
 * @param[in,out]    duplicates          Addresses of PAHs that have already been serialised
 *
 * @exception		 invalid_argument    Stream not ready
 */
void PAHPrimary::Serialize(std::ostream &out, void *duplicates) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        if (m_pmodel->WriteBinaryTrees()) {
            // Call the binary tree serialiser...
            BinTreeSerializer <PAHPrimary> tree;
            tree.Serialize(out, this, duplicates);
        } else {
            // Just serialise the root node.
            SerializePrimary(out, duplicates);
        }

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PAHPrimary::Serialize).");
    }
}

/*!
 * @brief Writes an individual primary to a binary stream
 *
 * @param[in,out]    out                 Output binary stream
 * @param[in,out]    duplicates          Addresses of PAHs that have already been serialised
 *
 * @exception		 invalid_argument    Stream not ready
 */
void PAHPrimary::SerializePrimary(std::ostream &out, void *duplicates) const
{
	if (out.good()) {

		int  val_int(0);
		double val(0.0);
		// Serialise state space
		val_int = m_numcarbon;
		out.write((char*)&val_int, sizeof(val_int));

		val_int = m_numH;
		out.write((char*)&val_int, sizeof(val_int));

		val_int = m_numOfEdgeC;
		out.write((char*)&val_int, sizeof(val_int));

		val_int = m_numOfRings;
		out.write((char*)&val_int, sizeof(val_int));

		val_int = m_numPAH;
		out.write((char*)&val_int, sizeof(val_int));

		val_int = m_numprimary;
		out.write((char*)&val_int, sizeof(val_int));

		val = m_PAHmass;
		out.write((char*)&val, sizeof(val));

		val = m_PAHCollDiameter;
		out.write((char*)&val, sizeof(val));

		val = m_primarydiam;
		out.write((char*)&val, sizeof(val));

		val = m_children_radius;
		out.write((char*)&val, sizeof(val));

		val = m_children_vol;
		out.write((char*)&val, sizeof(val));

		val = m_leftparticle_vol_old;
		out.write((char*)&val, sizeof(val));

		val = m_rightparticle_vol_old;
		out.write((char*)&val, sizeof(val));

		val_int = m_rightparticle_numPAH;
		out.write((char*)&val_int, sizeof(val_int));

		val_int = m_leftparticle_numPAH;
		out.write((char*)&val_int, sizeof(val_int));

		val = m_children_surf;
		out.write((char*)&val, sizeof(val));

		val = m_free_surf;
		out.write((char*)&val, sizeof(val));

		val = m_sum_necks;
		out.write((char*)&val, sizeof(val));

		val = m_primaryvol;
		out.write((char*)&val, sizeof(val));

		val = m_children_roundingLevel;
		out.write((char*)&val, sizeof(val));

		val = m_distance_centreToCentre;
		out.write((char*)&val, sizeof(val));

		val = m_cen_bsph[0];
		out.write((char*)&val, sizeof(val));

		val = m_cen_bsph[1];
		out.write((char*)&val, sizeof(val));

		val = m_cen_bsph[2];
		out.write((char*)&val, sizeof(val));

		val = m_cen_mass[0];
		out.write((char*)&val, sizeof(val));

		val = m_cen_mass[1];
		out.write((char*)&val, sizeof(val));

		val = m_cen_mass[2];
		out.write((char*)&val, sizeof(val));

		val = m_r;
		out.write((char*)&val, sizeof(val));

		val = m_r2;
		out.write((char*)&val, sizeof(val));

		val = m_r3;
		out.write((char*)&val, sizeof(val));

		val = m_children_sintering;
		out.write((char*)&val, sizeof(val));

		val = m_avg_sinter;
		out.write((char*)&val, sizeof(val));

		/*Imaging properties
		val = (double)m_Rg;
		out.write((char*)&val, sizeof(val));

		val = (double)m_fdim;
		out.write((char*)&val, sizeof(val));

		val = (double)m_sqrtLW;
		out.write((char*)&val, sizeof(val));

		val = (double)m_LdivW;
		out.write((char*)&val, sizeof(val));*/

		val = m_avg_coalesc;
		out.write((char*)&val, sizeof(val));

		val = m_sint_time;
		out.write((char*)&val, sizeof(val));

		val_int = (int)m_PAH.size();
		out.write((char*)&val_int, sizeof(val_int));

		// write the PAH stack (m_PAH)
		PahSerialisationMap *pahDuplicates = reinterpret_cast<PahSerialisationMap*>(duplicates);
		outputPAHs(out, *pahDuplicates);

		// Output base class.
		Primary::Serialize(out);

	}
	else {
		throw invalid_argument("Output stream not ready "
			"(Sweep, PAHPrimary::SerializePrimary).");
	}
}

/*
 * @brief Writes individual PAHs to a binary stream 
 *
 * @param[in]        out               Output binary stream
 * @param[in,out]    pah_duplicates    Addresses of PAHs that have already been serialised
 */
void PAHPrimary::outputPAHs(std::ostream &out, PahSerialisationMap &pah_duplicates) const
{
    if (m_PAH.size() != 0) {
        //count the number of PAH should be serialized
        unsigned count = 0;
		int PAHSize= m_PAH.size();
        
        // Keep a record of which PAHs have been serialised to avoid serialising the same
        // memory location twice.
        //PahSerialisationMap pahUniqueAdresses;

        while (count != m_PAH.size())
        {
            // Serialise the raw memory locations of the PAHs, so that we
            // can recognise PAHs that are referenced multiple times on
            // deserialisation.
            void *ptr = m_PAH[count].get();
            out.write(reinterpret_cast<char*>(&ptr), sizeof(ptr));
            //if(pahUniqueAdresses.find(m_PAH[count].get()) == pahUniqueAdresses.end()) {
            //    m_PAH[count]->Serialize(out);
            //    pahUniqueAdresses.insert(m_PAH[count].get());
            //}
            //else {
            //    std::cout << "PAH at " << m_PAH[count].get() << " has already been serialised\n";
            //}
			if(pah_duplicates.find(m_PAH[count].get()) == pah_duplicates.end()){
				m_PAH[count]->Serialize(out);
				pah_duplicates.insert(m_PAH[count].get());
			}
			else {
                //std::cout << "PAH at " << m_PAH[count].get() << " has already been serialised\n";
            }
            ++count;
        }
    }
}

/* 
 * @brief Reads the object from a binary stream.
 *
 * @param[in,out]	 in		             Stream from which to read
 * @param[in]        model	             Particle model defining interpretation of particle data
 * @param[in,out]    duplicates          Addresses of PAHs for use when reading primary particles
 *
 * @exception		 invalid_argument    Stream not ready
 */
void PAHPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model, PahDeserialisationMap &pah_duplicates)
{
	if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        if (model.WriteBinaryTrees()) {
            // Call the binary tree serialiser...
            BinTreeSerializer <PAHPrimary> tree;
            tree.Deserialize(in, this, model, &pah_duplicates);
        } else {
            // Just deserialise the root node.
            DeserializePrimary(in, model, &pah_duplicates);
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PAHPrimary::Deserialize).");
    }
}

/*
 * @brief Reads an individual primary from the binary stream
 *
 * @param[in,out]	 in		             Input binary stream
 * @param[in]        model	             Particle model defining interpretation of particle data
 * @param[in,out]    duplicates          Addresses of PAHs for use when reading primary particles
 *
 * @exception		 invalid_argument    Stream not ready
 */
void PAHPrimary::DeserializePrimary(std::istream &in, const Sweep::ParticleModel &model, void *duplicates)
{
	if (in.good()) {

		int  val_int(0);
		double val(0.0);

		// Serialise state space
		in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
		m_numcarbon = val_int;

		in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
		m_numH = val_int;

		in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
		m_numOfEdgeC = val_int;

		in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
		m_numOfRings = val_int;

		in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
		m_numPAH = val_int;

		in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
		m_numprimary = val_int;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_PAHmass = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_PAHCollDiameter = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_primarydiam = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_children_radius = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_children_vol = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_leftparticle_vol_old = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_rightparticle_vol_old = val;

		in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
		m_rightparticle_numPAH = val_int;

		in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
		m_leftparticle_numPAH = val_int;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_children_surf = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_free_surf = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_sum_necks = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_primaryvol = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_children_roundingLevel = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_distance_centreToCentre = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_cen_bsph[0] = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_cen_bsph[1] = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_cen_bsph[2] = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_cen_mass[0] = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_cen_mass[1] = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_cen_mass[2] = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_r = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_r2 = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_r3 = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_children_sintering = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_avg_sinter = val;

		/* imaging properies
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_Rg = val;
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_fdim = val;
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_sqrtLW = val;
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_LdivW = val;*/

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_avg_coalesc = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_sint_time = val;

		in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
		// the size of m_PAH
		const int PAHcount = val_int;
		// read the PAH stack (m_PAH)
		PahDeserialisationMap *pahDuplicates = reinterpret_cast<PahDeserialisationMap*>(duplicates);
		if (PAHcount != 0)
			inputPAHs(in, model, PAHcount, *pahDuplicates);

		m_leftchild = NULL;
		m_rightchild = NULL;
		m_leftparticle = NULL;
		m_rightparticle = NULL;
		m_parent = NULL;
		// Output base class.
		Primary::Deserialize(in, model);

	}
	else {
		throw invalid_argument("Input stream not ready "
			"(Sweep, PAHPrimary::DeserializePrimary).");
	}
}

/*
 * @brief Reads individual PAHs from the binary stream 
 *
 * @param[in]        in                Input binary stream
 * @param[in]        model	           Particle model defining interpretation of particle data
 * @param[in]        PAHcount          Number of PAHs in the particle
 * @param[in,out]    pah_duplicates    Addresses of PAHs for use when reading primary particles
 */
void PAHPrimary::inputPAHs(std::istream &in, const Sweep::ParticleModel &model, const int PAHcount, PahDeserialisationMap &pah_duplicates)
{
    int m_count=0;
    
    //PahDeserialisationMap pahPointerMappings;
    while (m_count != PAHcount)
    {
        // Raw memory location of PAH prior to serialisation.
        void* preSerializePointer;
        in.read(reinterpret_cast<char*>(&preSerializePointer), sizeof(preSerializePointer));
        //cout << preSerializePointer;

        // See if this PAH has already been deserialised
        //PahDeserialisationMap::iterator searchResult = pahPointerMappings.find(preSerializePointer);
        //if(searchResult != pahPointerMappings.end())
        //{
        //    // Found PAH in list that has already been deserialised
        //    m_PAH.push_back(searchResult->second);
        //}
        //else {
        //    // create a new PAH to load in its details
        //    boost::shared_ptr<PAH> new_PAH (new PAH());
        //    new_PAH->Deserialize(in);
        //    m_PAH.push_back(new_PAH);
        //    pahPointerMappings.insert(std::make_pair(preSerializePointer, new_PAH));
        //}
        PahDeserialisationMap::iterator searchResult = pah_duplicates.find(preSerializePointer);
        if(searchResult != pah_duplicates.end())
        {
            // Found PAH in list that has already been deserialised
            m_PAH.push_back(searchResult->second);
        }
        else {
            // create a new PAH to load in its details
            boost::shared_ptr<PAH> new_PAH (new PAH());
            new_PAH->Deserialize(in);
            m_PAH.push_back(new_PAH);
            pah_duplicates.insert(std::make_pair(preSerializePointer, new_PAH));
        }
        ++m_count;
    }
}

/*!
 * This is a helper function for ???????
 * It climbs up the tree from bottom to top recording a route
 * suitable for use in call to @see descendPath.
 *
 *@param[in]    bottom      Tree node from which to start climbing
 *@param[in]    top         Tree node at which to stop climbing
 *
 *@pre      top must be above bottom in a tree
 *
 *@return   Stack that can be used to descend the same path by moving to the
 *          left child each time the top of the stack is true.
 */
std::stack<bool> PAHPrimary::recordPath(const PAHPrimary* bottom,
                                        const PAHPrimary* const top) {
    std::stack<bool> wasLeftChild;

    while(bottom != top) {
        // check whether bottom was a left child of its parent
        wasLeftChild.push(bottom == bottom->m_parent->m_leftchild);

        // Climb one level up the tree
        bottom = bottom->m_parent;
    }
    return wasLeftChild;
}

/*!
 *@param[in]        here            Point in tree from which to start descent
 *@param[in,out]    takeLeftBranch  Instructions for which child to move to at each level
 *
 *@return   The node at the bottom of the path
 *
 *@pre  here must be a node of tree in which takeLeftBranch is a valid path
 *@post takeLeftBranch.empty() == true
 */
PAHPrimary* PAHPrimary::descendPath(PAHPrimary *here,
                                    std::stack<bool> &takeLeftBranch) {
    while(!takeLeftBranch.empty()) {
        // Move one step down the tree in the instructed direction
        if(takeLeftBranch.top())
            here = here->m_leftchild;
        else
            here = here->m_rightchild;

        // This instuction has now been processed
        takeLeftBranch.pop();
    }

    return here;
}

/*!
 *  @brief Sinters particles for time dt.
 *
 *  This function only operates on non-leaf nodes. It begins at the root node,
 *  which sinters for time dt. It then descends the tree to sinter nodes below
 *  the root. If the sintering level rises above 95%, or if the distance
 *  between the centres of primary particles is 0 (if tracked), Merge is called
 *  and the particles are combined. 
 *
 *  @param[in] dt    Time for which to sinter.
 *  @param[in] sys   Environment for particles.
 *  @param[in] model Sintering model to apply.
 *  @param[in] rng   Random number generator.
 *  @param[in] wt    Statistical weight.
 */
void PAHPrimary::SinterNode(double dt, Cell &sys, const Processes::SinteringModel &model, rng_type &rng, double wt)
{
	//! Declare sintering rate.
	double r = 0.0;
	// Declare time step variables.
	double t1 = 0.0, delt = 0.0, tstop = dt;

	// The scale parameter discretises the delta-S when using
	// the Poisson distribution.  This allows a smoother change
	// (smaller scale = higher precision).
	double scale = 0.01;

	////! Only update the time on the root node.
	//if (m_parent == NULL) {
	//	m_sint_time += dt;
	//	SetSinteringTime(m_sint_time);
	//}

	//! Do only if there is a particle to sinter.
	if (m_leftparticle != NULL)
	{
		//bool Condition;                          //! Declare variable for condition for complete sintering.

		//! The sintering model depends on whether the distance between the
		//! centres of primary particles is tracked. If tracked, sintering
		//! results in a decrease in the distance between the primaries and an
		//! increase in their diameters. If not, sintering results in a
		//! decrease in the common surface between the primaries.
		if (!(m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())) {
			//! Store the old surface area of particles.
			// double surf_old = m_children_surf;

			//! Calculate the spherical surface.
			const double spherical_surface = 4 * PI * m_children_radius * m_children_radius;

			//! Define the maximum allowed change in surface area in one
			//! internal time step (10% of spherical surface).
			double dAmax = 0.1 * spherical_surface;

			//! Perform integration loop.
			while (t1 < tstop) {
				//! Calculate sintering rate.
				r = model.Rate(m_time + t1, sys, *this);

				if (r > 0) {
					//! Calculate next time-step end point so that the
					//! surface area changes by no more than dAmax.
					delt = dAmax / max(r, 1.0e-300);

					//! Approximate sintering by a poisson process.  Calculate
					//! number of poisson events.
					double mean;

					if (tstop > (t1 + delt)) {
						//! A sub-step, we have changed surface by dAmax, on average.
						mean = 1.0 / scale;
					}
					else {
						//! Step until end. Calculate degree of sintering explicitly.
						mean = r * (tstop - t1) / (scale * dAmax);
					}

					boost::random::poisson_distribution<unsigned, double> repeatDistribution(mean);
					const unsigned n = repeatDistribution(rng);

					//! Adjust the surface area.
					if (n > 0) {
						m_children_surf -= (double)n * scale * dAmax;

						//! Check that primary is not completely sintered.
						if (m_children_surf <= spherical_surface) {
							m_children_surf = spherical_surface;
							break;
						}
					}

					//! Set t1 for next time step.
					t1 += delt;
				}
				else{
					t1 = tstop;
				}
			}
		}
		else {
			//! From bintree model
			//! Declare characteristic sintering time.

			m_leftparticle->UpdatePrimary(); 
			m_rightparticle->UpdatePrimary(); 

			//! Define the maximum allowed change (1%) in the distance between the
			//! centres of primary particles in one internal time step. In the case
			//! of pure sintering it was found that if the allowed change is too
			//! large (~10%) a significant error is incurred in the final spherical
			//! volume determined through comparisons with the mass-derived volume.
			//! Note that the smaller the distance is, the smaller the changes are.
			double dd_ij_Max = m_distance_centreToCentre / 100.0;

			while (t1 < tstop) {

				//! Definition of variables
				double r_i = this->m_leftparticle->m_primarydiam / 2.0;
				double r_j = this->m_rightparticle->m_primarydiam / 2.0;
				double d_ij = m_distance_centreToCentre;

				double d_ij2 = pow(d_ij, 2.0);
				double r_i2 = pow(r_i, 2.0);
				double r_j2 = pow(r_j, 2.0);
				double r_i4 = pow(r_i, 4.0);
				double r_j4 = pow(r_j, 4.0);

				//! Continue if primaries have not coalesced
				if (!MergeCondition() && !FakeRounding()) { // add FakeRounding() by hdy
				//if (!MergeCondition()) { //used to test bintrtee model and PAH_KMC model

					//! Due to rounding, x_i and x_j are sometimes calculated to be larger 
					//! than the respective primary radii resulting in a negative neck area.
					//! Therefore we take the smaller of x_i and r_i.
					double x_i = min((d_ij2 - r_j2 + r_i2) / (2.0 * d_ij), r_i); //!< Eq. (3b).
					double x_j = min((d_ij2 - r_i2 + r_j2) / (2.0 * d_ij), r_j); //!< Eq. (3b).
					double A_n = M_PI * (r_i2 - pow(x_i, 2.0));        //!< Eq. (4).

					//declare more variables
					double dd_ij_dt = 0.0;
					double R_n = 0.0;
					double r4_tau = 0.0;
					double tau = 0.0;
					double gamma_eta = 0.0;

					//! In Section 3.1.2 of Langmuir 27:6358 (2011), it is argued that
					//! the smaller particle dominates the sintering process.

					if (r_i <= r_j) {
						tau = model.SintTime(sys, *this->m_leftparticle); //!< The left particle is smaller than the right. 
					}
					else {
						tau = model.SintTime(sys, *this->m_rightparticle); //!< The right particle is smaller than the left.
					}

					//! Gamma is the surface tension and eta is the viscosity, and the
					//! ratio (gamma/eta) can be related to tau.
					//! J. Colloid Interface Sci. 140:419 (1990).
					gamma_eta = min(r_i, r_j) / tau;

					///////////////////////////////////////////////////////
					/// References to equations in Langmuir 27:6358 (2011).
					///////////////////////////////////////////////////////

					//! Eq. (14a).
					dd_ij_dt = 4.0 * r_i * r_j * d_ij2 * (r_i + r_j) * gamma_eta /
						((r_i + r_j + d_ij) * (r_i4 + r_j4 - 2.0 * r_i2 * r_j2 + 4.0 * d_ij * r_i * r_j *(r_i + r_j) - d_ij2 * (r_i2 + r_j2)));

					//! Get surface area and subtract mutual contribution
					double A_i = m_leftparticle->m_free_surf + m_leftparticle->m_sum_necks - M_PI*(r_i*r_i - x_i*x_i)*r_i / x_i;
					double A_j = m_rightparticle->m_free_surf + m_rightparticle->m_sum_necks - M_PI*(r_j*r_j - x_j*x_j)*r_j / x_j;

					//! @todo Remove derivation and replace with reference to preprint
					//!       or paper if results do get published.
					double B_i = (-r_j*A_n*A_n - x_j*A_j*A_n) / (A_i*A_j*d_ij + r_i*A_j*A_n + r_j*A_i*A_n);
					double B_j = (-r_i*A_n*A_n - x_i*A_i*A_n) / (A_j*A_i*d_ij + r_j*A_i*A_n + r_i*A_j*A_n);

					delt = dd_ij_Max / max(dd_ij_dt, 1.0e-300);
					double mean;

					if (tstop > (t1 + delt)) {
						mean = 1.0 / scale;
					}
					else {
						mean = dd_ij_dt * (tstop - t1) / (scale * dd_ij_Max);
					}

					//! Sinter primaries
					boost::random::poisson_distribution<unsigned, double> repeatDistribution(mean);
					const unsigned n = repeatDistribution(rng);

					double delta_dij = -(double)n * scale * dd_ij_Max; //!< Sintering decreases d_ij hence the negative sign.
					m_distance_centreToCentre += delta_dij;

					//! if coordinates are tracked then we will translate one side of the particle by the change in separation
					//! this is faster than translating both sides by half the change
					if (m_pmodel->getTrackPrimaryCoordinates()) {
						//! get direction of translation (left particle to right particle)
						Coords::Vector vector_change = UnitVector(m_leftparticle->boundSphCentre(), m_rightparticle->boundSphCentre());
						//! translate the leftparticle
						m_leftparticle->TranslatePrimary(vector_change, -delta_dij);
						//! translate all neighbours of the left particle except the right particle
						m_leftparticle->TranslateNeighbours(m_leftparticle, vector_change, -delta_dij, m_rightparticle);
					}

					//! Change in primary radii
					double delta_r_i = -(double)n * scale * B_i * dd_ij_Max;  //!< Eq. (8).
					double delta_r_j = -(double)n * scale * B_j * dd_ij_Max;  //!< Eq. (8).

					//! Adjust separation of neighbours that are not currently sintering
					m_leftparticle->UpdateConnectivity(m_leftparticle, delta_r_i, m_rightparticle);
					m_rightparticle->UpdateConnectivity(m_rightparticle, delta_r_j, m_leftparticle);

					//! Adjust primary radii
					this->m_leftparticle->m_primarydiam += 2.0 * delta_r_i;
					this->m_rightparticle->m_primarydiam += 2.0 * delta_r_j;

					//! update primaries
					m_leftparticle->UpdateOverlappingPrimary();
					m_rightparticle->UpdateOverlappingPrimary();

					t1 += delt;

					//! Return some sintering rate
					r = dd_ij_dt;
				}
				else {
					break; //!do not continue to sinter.
				}
			}
			//! If coordinates are tracked update tracking radius 
			if (m_pmodel->getTrackPrimaryCoordinates()) {
				m_leftparticle->setRadius(m_leftparticle->m_primarydiam / 2.0);
				m_rightparticle->setRadius(m_rightparticle->m_primarydiam / 2.0);
			}

		}

		if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates())
			m_children_roundingLevel = RoundingLevel();
		else m_children_sintering = SinteringLevel();
		
		//m_sint_rate = r;

	}
}

void PAHPrimary::SetSinteringTime(double time) 
{
    m_sint_time = time;
    // Update children
    if (m_leftchild != NULL) {
        m_leftchild->SetSinteringTime(time);
        m_rightchild->SetSinteringTime(time);
    }
    // Update particles
    if (m_leftparticle != NULL) {
        m_leftparticle->SetSinteringTime(time);
        m_rightparticle->SetSinteringTime(time);
    }
}

//returns the CoalescenceLevel of the two primaries that are connected by this node
//used for when primary coordinates are not tracked
double PAHPrimary::RoundingLevel()
{
    if (m_leftparticle!=NULL) {
        // Calculate the spherical surface
        const double spherical_surface=4*PI*m_children_radius*m_children_radius;
        const double two_1_3=0.79370052231642452;
        double slevel;

        if (m_children_surf <= spherical_surface) {
            m_children_surf = spherical_surface;
            return 1.0;
        }

        if (m_children_surf == 0.0) {
            slevel = 0.0;
        } else {
            slevel= ((spherical_surface/m_children_surf)-two_1_3)/(1-two_1_3);
        }

        if (slevel < 0.0) {
            return 0.0;
        } else if (slevel > 1.0) {
            cout << "sweep: PAHPrimary::CoalescenceLevel() > 1.0";
            return 1.0;
        } else {
            return slevel;
        }
    } else {
        // Particle is a primary, should have 0 as default properties
        m_children_surf = 0.0;
        m_children_radius = 0.0;
        m_children_vol = 0.0;
        return 0.0;
    }
}

//from bintree model
/*!
* @brief       Identify neighbours and sum their cap areas and volumes
*
* Works up the binary tree identifying the neighbours of the adjusted primary
* and sums the contribution from neighbours to the free surface area.
* All neighbours of the primary being adjusted are left/rightparticles of nodes directly
* above it.
*
* @param[in]   prim		Pointer to the primary being adjusted
* @param[in]   CapAreas	Sum of cap areas
* @param[in]   CapVolumes	Sum of cap volumes
*/
void PAHPrimary::SumCaps(PAHPrimary *prim, double &CapAreas, double &CapVolumes, double &SumNecks){

	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;

	//! Check if a neighbour of prim
	if (m_parent->m_leftparticle == prim) {

		//! Right primary is a neighbour
		r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
		x_ij = min((pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / (2.0*d_ij), r_i); //ensure -r_i <= x_ij <= r_i
		x_ij = max(x_ij, -r_i);

		//! Calculate cap area and add to sum
		CapAreas += 2 * M_PI*(r_i*r_i - r_i*x_ij);

		//! Calculate cap volume and add to sum
		CapVolumes += M_PI * (2 * pow(r_i, 3.0) + pow(x_ij, 3.0) - 3.0*pow(r_i, 2.0)*x_ij) / 3.0;

		//! Neck area * r_i / x_ij
		SumNecks += abs(M_PI*(r_i*r_i - x_ij*x_ij) * r_i / x_ij);

	}
	else if (m_parent->m_rightparticle == prim) {

		//! Left primary is a neighbour
		r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
		x_ij = min((pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / (2.0*d_ij), r_i);	//ensure -r_i <= x_ij <= r_i
		x_ij = max(x_ij, -r_i);

		//! Calculate cap area and add to sum
		CapAreas += 2 * M_PI*(r_i*r_i - r_i*x_ij);

		//! Calculate cap volume and add to sum
		CapVolumes += M_PI * (2 * pow(r_i, 3.0) + pow(x_ij, 3.0) - 3.0*pow(r_i, 2.0)*x_ij) / 3.0;

		//! Neck area * r_i / x_ij
		SumNecks += abs(M_PI*(r_i*r_i - x_ij*x_ij) * r_i / x_ij);
	}

	//! Continue working up the binary tree
	if (m_parent->m_parent != NULL){
		m_parent->SumCaps(prim, CapAreas, CapVolumes, SumNecks);
	}
}

//! From bintree model
/*! Updates primary free surface area and volume
*
* @param[in]   this		Primary to update
*/
void PAHPrimary::UpdateOverlappingPrimary(){

	//! Get sum of cap areas and volumes
	double CapAreas = 0.0;			//!< Contribution from neighbours to free surface area
	double CapVolumes = 0.0;		//!< Contribution from neighbours to volume
	double SumNecks = 0.0;			//!< Sum of necks * r_i / x_ij
	SumCaps(this, CapAreas, CapVolumes, SumNecks);

	//! Update free surface area
	//if the calculated area is negative (too many overlaps) then set m_free_surf = 0.0
	m_free_surf = max(M_PI*m_primarydiam*m_primarydiam - CapAreas, 0.0);

	//! Update primary volume
	//if calculated volume is negative (too many overlaps) then set m_primaryvol = 0.0
	m_primaryvol = max(M_PI*pow(m_primarydiam, 3.0) / 6.0 - CapVolumes, 0.0);

	//! Update sum of necks 
	m_sum_necks = max(SumNecks, 0.0);
}

//! From bintree model
/*!
* @brief       Checks if condition for merger is met
*
* @return      Boolean telling if the merger condition is met
*/
bool PAHPrimary::MergeCondition()
{
	bool condition = false;

	if (m_leftparticle != NULL) {
		//! The condition for whether a particle has coalesced depends on whether
		//! the distance between the centres of primary particles is tracked. 
		//! If not, the condition depends on whether the rounding 
		//! level exceeds an arbitrarily high threshold.
		if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()) {
			condition = (m_children_roundingLevel > 0.95);
		}
		else {
			//! If tracked, a particle has coalesced when the neck reaches the centre 
			//! of one of the primaries. Condition: min(x_ij,x_ji) <= 0
			//  (If the neck were allowed to move beyond this point and reside outside 
			//  the region between the primary centres any (growth) adjustments would 
			//  fail.) 
			double r_i = m_leftparticle->m_primarydiam / 2.0;
			double r_j = m_rightparticle->m_primarydiam / 2.0;
			double d_ij = m_distance_centreToCentre;

			if (d_ij <= 0.0 ){
			//if (d_ij <= 0.0 || d_ij < 0.5*(r_i + r_j)){
				//ensures that particles are merged if the sintering step overshoots
				condition = true;
			}
			else{
				//! Primaries are merged when the neck radius is 95% of the smaller primary radius.
				//! The second condition ensures that primaries are merged even if sintering overshoots
				//! i.e. the neck crosses the centre of smaller primary
				double x_ij = (d_ij*d_ij - r_j*r_j + r_i*r_i) / (2.0*d_ij);
				double x_ji = (d_ij*d_ij - r_i*r_i + r_j*r_j) / (2.0*d_ij);
				double R_ij = sqrt(r_i*r_i - x_ij*x_ij);	//!neck radius
				condition = R_ij / min(r_i, r_j) >= 0.95 || ((pow(d_ij, 2.0) - pow(max(r_i, r_j), 2.0) + pow(min(r_i, r_j), 2.0)) / (2.0*d_ij)) <= 0.0;
					
			}
		}
	}
	return condition;
}

//! From bintree model
/*!
*  Calculates unit vector between two set coordinates
*
*  @param[in]    x_i	Coordinates
*  @param[in]    x_j	Coordinates
*  @param[out]   unit vector
*/
Coords::Vector PAHPrimary::UnitVector(Coords::Vector x_i, Coords::Vector x_j)
{
	Coords::Vector delta_x;
	double len_delta_x;

	//! Calculate difference x_j - x_i
	delta_x[0] = x_j[0] - x_i[0];
	delta_x[1] = x_j[1] - x_i[1];
	delta_x[2] = x_j[2] - x_i[2];

	//! Calculate the length of the vector
	len_delta_x = sqrt(delta_x[0] * delta_x[0] + delta_x[1] * delta_x[1] + delta_x[2] * delta_x[2]);

	//! Create unit vector
	delta_x[0] /= len_delta_x;
	delta_x[1] /= len_delta_x;
	delta_x[2] /= len_delta_x;

	return delta_x;
}

// !From bintree model
/*!
*  Translates a primary particle along a unit vector
*
*  @param[in]    u			Unit vector (direction)
*  @param[in]    delta_d	Distance to translate
*/
void PAHPrimary::TranslatePrimary(Coords::Vector u, double delta_d)
{
	//! Bounding sphere coordinates
	m_cen_bsph[0] += delta_d * u[0];
	m_cen_bsph[1] += delta_d * u[1];
	m_cen_bsph[2] += delta_d * u[2];
	//! Centre of mass coordinates
	m_cen_mass[0] += delta_d * u[0];
	m_cen_mass[1] += delta_d * u[1];
	m_cen_mass[2] += delta_d * u[2];
}

//! From bintree model
/*!
*  Translates neighbours of a primary particle along a unit vector
*
*  @param[in]    prim			Primary
*  @param[in]    u				Unit vector direction of translation
*  @param[in]    delta_d		Magnitude of translation
*  @param[in]    prim_ignore	Primary to ignore
*/
void PAHPrimary::TranslateNeighbours(PAHPrimary *prim, Coords::Vector u, double delta_d, PAHPrimary *prim_ignore)
{
	PAHPrimary *neighbour = NULL;
	//! Check if a neighbour of prim but not prim_ignore
	if (m_parent->m_leftparticle == prim && m_parent->m_rightparticle != prim_ignore) {
		//! right particle is a neighbour
		neighbour = m_parent->m_rightparticle;
		//! adjust its coordinates
		neighbour->TranslatePrimary(u, delta_d);
		//! adjust its neighbours except for prim
		neighbour->TranslateNeighbours(neighbour, u, delta_d, prim);
	}
	else if (m_parent->m_rightparticle == prim &&  m_parent->m_leftparticle != prim_ignore) {
		//! left particle is a neighbour
		neighbour = m_parent->m_leftparticle;
		//! adjust its coordinates
		neighbour->TranslatePrimary(u, delta_d);
		//! adjust its neighbours except for prim
		neighbour->TranslateNeighbours(neighbour, u, delta_d, prim);
	}

	//! continue working up the binary tree
	if (m_parent->m_parent != NULL){
		m_parent->TranslateNeighbours(prim, u, delta_d, prim_ignore);
	}
}

//! From bintree model
/*!
* @brief    Calculates the sintering level for particles connected by node.
*
* Unlike the original SilicaPrimary model; this model assumes that a single
* primary particle (no tree structure) has a sintering level of 1.0. In the
* case where it is part of a tree, 0.0 is returned so as to not cause
* erroneous calculation of m_children_sintering.
*
* @return   Sintering level
*/
double PAHPrimary::SinteringLevel()
{
	if (m_leftchild != NULL && m_rightchild != NULL) {

		double slevel(0.0);

		//! If the centre-centre separation is not tracked the sintering level is calculated 
		//! as per Shekar et al. (2012)
		if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()) {
			// Calculate the spherical surface
			const double spherical_surface =
				4 * PI * m_children_radius * m_children_radius;

			if (m_children_surf == 0.0) {
				slevel = 0.0;
			}
			else {
				slevel = ((spherical_surface / m_children_surf) - TWO_ONE_THIRD)
					/ (1 - TWO_ONE_THIRD);
			}

			//! if centre-centre separation is tracked the sintering level is calculated as
			//! s = R_ij / min(r_i,r_j)
			//! 0 <= s <= 1 is consistent with the merger condition
			//To do: paper reference here
		}
		else{
			if (m_leftparticle != NULL && m_rightparticle != NULL) {

				double r_i = m_leftparticle->m_primarydiam / 2.0;
				double r_j = m_rightparticle->m_primarydiam / 2.0;
				double d_ij = m_distance_centreToCentre;
				double x_ij = (d_ij*d_ij - r_j*r_j + r_i*r_i) / (2.0*d_ij);
				double R_ij = sqrt(r_i*r_i - x_ij*x_ij);

				slevel = R_ij / min(r_i, r_j);
			}
		}

		if (slevel < 0.0) {
			return 0.0;
		}
		else if (slevel > 1.0) {
			return 1.0;
		}
		else return slevel;

	}
	else {
		// Particle is a primary
		m_children_surf = 0.0;
		m_children_radius = 0.0;
		m_children_vol = 0.0;
		if (m_parent == NULL) return 1.0;        // Single primary case
		else return 0.0;                         // Part of a tree
	}
}

//! From bintree model//
/*!
* @brief       Adjust primary radius after volume addition due to merger.
*				Adjustment is performed in a similar manner to sintering,
*				assuming that the neck radii are unchanged; hence,
*				the neighbours are translated outwwards.
*
* @param[in] V1			Volume to be added
* @param[in] d_ij			Separation of merging primaries
* @param[in] prim_ignore	Primary to be ignored in adjustments
*/
void PAHPrimary::AdjustPrimary(double V1, double d_ij, PAHPrimary *prim_ignore)
{
	double V0 = 0.0;								//!< Variable to track added volume
	double r_i = m_primarydiam / 2.0;				//!< Radius of primary particle i
	double r_j = prim_ignore->m_primarydiam / 2.0;	//!< Radius of primary particle j
	double x_i = (d_ij*d_ij - r_j*r_j + r_i*r_i) / (2.0*d_ij);	//!< Distance to merging neck
	double dr_max = 0.01*r_i;						//!< Maximum change in primary radius during internal step (1% of primary radius)
	double dr_i = 0.0;								//!< Change in radius of i

	while (V0 <= V1){

		r_i = m_primarydiam / 2.0;

		//! change in volume (exclude contribution from merging neck)
		double dV = dr_max * (m_free_surf + 2.0*M_PI*(r_i*r_i - r_i*x_i) + max(m_sum_necks - abs(M_PI*(r_i*r_i - x_i*x_i)*r_i / x_i), 0.0));

		assert(dV > 0.0);

		//! change in radius
		if (V0 + dV > V1){
			dr_i = (V1 - V0)*dr_max / dV;
		}
		else{
			dr_i = dr_max;
		}

		assert(dr_i > 0.0);
		assert(dr_i <= 1.1*dr_max);

		//! Update the particle separations ignoring the smaller primary
		if (m_parent != NULL) UpdateConnectivity(this, dr_i, prim_ignore);

		//! Update primary diameter
		m_primarydiam = 2.0* (r_i + dr_i);

		//! Update primary properties
		this->UpdateOverlappingPrimary();

		V0 += dV;
	}
}

//! From bintree model.//
//! Overload function. Used when primary coordinates are tracked.//
/*!
* @brief       Changes pointer from source to target when centre-centre separation is tracked
*
* If a primary neighbours the smaller of the merging pair the centre to centre separation is
* re-estmated as the smaller of the sum of the separation or the sum of primary radii.
*
* @param[in] source		Pointer to the original particle
* @param[in] target		Pointer to the new particle
* @param[in] small_prim	Smaller of merging primaries
* @param[in] node			Pointer to merging neck (non-leaf node)
*/
void PAHPrimary::ChangePointer(PAHPrimary *source, PAHPrimary *target, PAHPrimary *small_prim, PAHPrimary *node)
{

	if (m_rightparticle == source) {
		//! if the neighbour is the smaller of the merging primaries then add new neighbour
		if (this != node){
			if (source == small_prim){

				double r_j = m_rightparticle->m_primarydiam / 2.0;	//!< radius of smaller merging primary
				double r_k = m_leftparticle->m_primarydiam / 2.0;	//!< radius of neighbour of merging primary
				double d_kj = m_distance_centreToCentre;			//!< primary separation
				double x_kj = (d_kj*d_kj - r_j*r_j + r_k*r_k) / 2.0 / d_kj;	//!< distance form centre of neighbour to neck with smaller merging primary
				double A_n_k = M_PI*(r_k*r_k - x_kj*x_kj);				//!< neck radius

				if (A_n_k > 0.0){
					//! add the neighbouring primary of smaller merging primary as a neighbour of the new merged primary
					double x_ik = target->AddNeighbour(A_n_k, small_prim, node);
					m_distance_centreToCentre = min(x_ik + x_kj, m_rightparticle->m_primarydiam / 2.0 + m_leftparticle->m_primarydiam / 2.0);
				}
				else{
					//! primaries in point contact
					m_distance_centreToCentre = m_rightparticle->m_primarydiam / 2.0 + m_leftparticle->m_primarydiam / 2.0;
				}

				//! adjust coordinates of new neighbour and all its neighbour
				//! this translates the branch along old separation vector d_ik to appropriate separation
				if (m_pmodel->getTrackPrimaryCoordinates()){

					Coords::Vector u_ik = UnitVector(m_leftparticle->boundSphCentre(), target->boundSphCentre());	//!< old separation unit vector
					double d_ik = Separation(m_leftparticle->boundSphCentre(), target->boundSphCentre());			//!< old separation distance 
					//! Translate the neighbour 
					m_leftparticle->TranslatePrimary(u_ik, d_ik - m_distance_centreToCentre);
					//! Translate all neighbours of the neighbour except the old small_prim
					m_leftparticle->TranslateNeighbours(m_leftparticle, u_ik, d_ik - m_distance_centreToCentre, small_prim);
				}
			}

			m_rightparticle = target;

			if (m_distance_centreToCentre < 0.0) m_distance_centreToCentre = -m_distance_centreToCentre; //! this is a length 

		}
		else{
			m_rightparticle = NULL;
		}

	}

	if (m_leftparticle == source){
		//! if the neighbour is the smaller of the merging primaries then add new neighbour
		if (this != node){
			if (source == small_prim){

				double r_j = m_leftparticle->m_primarydiam / 2.0;	//!< radius of smaller merging primary
				double r_k = m_rightparticle->m_primarydiam / 2.0;	//!< radius of neighbour of merging primary
				double d_kj = m_distance_centreToCentre;			//!< primary separation
				double x_kj = (d_kj*d_kj - r_j*r_j + r_k*r_k) / 2.0 / d_kj;	//!< distance form centre of neighbour to neck with smaller merging primary
				double A_n_k = M_PI*(r_k*r_k - x_kj*x_kj);				//!< neck radius

				if (A_n_k > 0.0){
					//! add the neighbouring primary of smaller merging primary as a neighbour of the new merged primary
					double x_ik = target->AddNeighbour(A_n_k, small_prim, node);
					m_distance_centreToCentre = min(x_ik + x_kj, m_rightparticle->m_primarydiam / 2.0 + m_leftparticle->m_primarydiam / 2.0);
				}
				else{
					//! primaries in point contact
					m_distance_centreToCentre = m_rightparticle->m_primarydiam / 2.0 + m_leftparticle->m_primarydiam / 2.0;
				}

				//! adjust coordinates of new neighbour and all its neighbour
				//! this translates the branch along old separation vector d_ik to appropriate separation
				if (m_pmodel->getTrackPrimaryCoordinates()){

					Coords::Vector vector_d_ik = UnitVector(m_rightparticle->boundSphCentre(), target->boundSphCentre());	//!< old separation unit vector
					double d_ik = Separation(m_rightparticle->boundSphCentre(), target->boundSphCentre());					//!< old separation distance 
					//! Translate the neighbour 
					m_rightparticle->TranslatePrimary(vector_d_ik, d_ik - m_distance_centreToCentre);
					//! Translate all neighbours of the neighbour except the old small_prim
					m_rightparticle->TranslateNeighbours(m_rightparticle, vector_d_ik, d_ik - m_distance_centreToCentre, small_prim);
				}
			}

			m_leftparticle = target;

			if (m_distance_centreToCentre < 0.0) m_distance_centreToCentre = -m_distance_centreToCentre; //! this is a length 

		}
		else{
			m_leftparticle = NULL;
		}
	}

	// Update the tree above this sub-particle.
	if (m_parent != NULL) {
		m_parent->ChangePointer(source, target, small_prim, node);
	}

}

//! From bintree model//
/*!
* @brief       Create neck for new neighbours added during merger event
*
* Neighbours of merging primary are added to the new merged primary
* (the larger of the merging pair). The new merged primary is 'sintered'
* to preserve the neck size of the neighbour being added
*
* @param[in]	A_n_k		Neck area of neighbour being added
* @param[in]	small_prim	Small merging primary
* @param[in]	this		New merged primary
* @param[out] centre to neck distance of new merged primary
*/
double PAHPrimary::AddNeighbour(double A_n_k, PAHPrimary *small_prim, PAHPrimary *node)
{

	double dr_i = 0.0;					//!< Change in radius
	double r_i = m_primarydiam / 2.0;		//!< Radius of new (merged) primary
	double dx_max = -0.1*r_i;			//!< Maximum step change in x_ik
	double dx_i = 0.0;					//!< Change in primary centre to neck distance
	double x_ik = 0.999*r_i;			//!< Initialise primary centre to neck distance as 0.999*r_i
	double A_n_i = M_PI*(r_i*r_i - x_ik*x_ik);		//!< initial neck radius

	//! "Sinter" new primary until the neck is the same size as the 
	//! neck on the old neighbour of the old particle.
	//! Loop while the merged primary neck area is less than the desired area A_n_k
	while (A_n_i < A_n_k){

		//! variables for updating surface areas
		double NodeCapAreas = 0.0;
		double NodeCapVolumes = 0.0;
		double NodeSumNecks = 0.0;
		double SmallCapAreas = 0.0;
		double SmallCapVolumes = 0.0;
		double SmallSumNecks = 0.0;
		double CapAreas = 0.0;
		double CapVolumes = 0.0;
		double SumNecks = 0.0;

		//!Update surface areas
		//UpdateOverlappingPrimary ascends from the primary and misses some neighbours
		//becasue they are not on this path during the change pointer update.
		//(The tree is only rearranged after the update.)
		//To include the missed necks we sum the contribution ascending from the 
		//small primary and from the large primary and subtract the contribution 
		//ascending from node to avoid double counting.
		this->SumCaps(this, CapAreas, CapVolumes, SumNecks);
		small_prim->SumCaps(this, SmallCapAreas, SmallCapVolumes, SmallSumNecks);
		if (node->m_parent != NULL) node->SumCaps(this, NodeCapAreas, NodeCapVolumes, NodeSumNecks);

		//! Update free surface area
		//if the calculated area is negative (too many overlaps) then set m_free_surf = 0.0
		m_free_surf = max(M_PI*m_primarydiam*m_primarydiam - CapAreas - SmallCapAreas + NodeCapAreas, 0.0);

		//! Update sum of necks 
		m_sum_necks = max(SumNecks + SmallSumNecks - NodeSumNecks, 0.0);

		r_i = m_primarydiam / 2.0;		//!< Radius of new (merged) primary

		//Add reference to preprint/paper here.
		//Surface areas exclude the contribution from the new neck being added 
		//we must add the contribution to the free surface area,
		//this is to be excluded from the sum over necks anyway.
		double B_ik = -A_n_i / (m_free_surf + m_sum_necks - 2 * M_PI*(r_i*r_i - r_i*x_ik));

		//! Change in radius
		dr_i = B_ik * dx_max;
		//! Save old neck size
		double A_n_i_old = A_n_i;
		//! New neck size
		A_n_i = M_PI*((r_i + dr_i)*(r_i + dr_i) - (x_ik + dx_max)*(x_ik + dx_max));

		//! If desired neck size is exceeded then solve the quadratic for dx
		//! A_n_i = M_PI*((r_i+dr_i)*(r_i+dr_i) - (x_ik+dx)*(x_ik+dx));
		if (A_n_i > A_n_k) {
			//coefficients
			double a = B_ik*B_ik - 1;
			double b = 2.0*r_i*B_ik - 2.0*x_ik;
			double c = r_i*r_i - x_ik*x_ik - A_n_k / M_PI;

			//! dx should be negative
			dx_i = min((-b - sqrt(b*b - 4 * a*c)) / 2.0 / a, (-b + sqrt(b*b - 4 * a*c)) / 2.0 / a);

			//! re-calculate change in radis
			dr_i = B_ik * dx_i;

			assert(dx_i <= 0.0);
			assert(dr_i >= 0.0);
		}
		else{
			dx_i = dx_max;
		}

		//! update connectivity ignoring small (merged) primary
		if (m_parent != NULL) UpdateConnectivity(this, dr_i, small_prim);

		//! update radius and x_ik
		r_i += dr_i;
		m_primarydiam += 2.0*dr_i;
		x_ik += dx_i;

		// if(x_ik <=0.0) break; //break if neck reaches maximum //csl37- necesssary?

	}

	return x_ik;
}

//! From bintree model//
/*!
*  Calculates distance between two points
*
*  @param[in]    x_i	Coordinates
*  @param[in]    x_j	Coordinates
*  @param[out]   Separation
*/
double PAHPrimary::Separation(Coords::Vector x_i, Coords::Vector x_j)
{
	Coords::Vector delta_x;
	double len_delta_x;

	//calculate difference x_j - x_i
	delta_x[0] = x_j[0] - x_i[0];
	delta_x[1] = x_j[1] - x_i[1];
	delta_x[2] = x_j[2] - x_i[2];

	//calculate the length of the vector
	len_delta_x = sqrt(delta_x[0] * delta_x[0] + delta_x[1] * delta_x[1] + delta_x[2] * delta_x[2]);

	return len_delta_x;
}

//! From bintree model//
/*!
* @brief       Sinters particles for time dt
*
* This function only operates on non-leaf nodes. It begins at the root
* node, which sinters for time dt. It then descends the tree to sinter
* nodes below the root. If the sintering level rises above 95%, Merge
* is called and the particles are combined.
*
* @param[in]   dt      Time for which to sinter
* @param[in]   sys     Environment for particles
* @param[in]   model   Sintering model to apply
* @param[in]   rng     Random number generator
* @param[in]   wt      Statistical weight
*/
void PAHPrimary::Sinter(double dt, Cell &sys,
	const Processes::SinteringModel &model,
	rng_type &rng,
	double wt)
{
	// Only update the time on the root node
	if (m_parent == NULL) {
		m_sint_time += dt;
		SetSinteringTime(m_sint_time);
	}

	// Do only if there is a particle to sinter
	if (m_leftparticle != NULL && m_rightparticle != NULL) {

		SinterNode(dt, sys, model, rng, wt);

		// Check if the sintering level is above the threshold, and merge
		if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()){
			bool Condition;
			Condition = (m_children_roundingLevel > 0.95);
			if (Condition) {
				CheckRounding();
				UpdateCache();

				if (m_leftchild != NULL && m_rightchild != NULL) {
					m_leftchild->Sinter(dt, sys, model, rng, wt);
					m_rightchild->Sinter(dt, sys, model, rng, wt);
				}
			}
			else {
				m_leftchild->Sinter(dt, sys, model, rng, wt);
				m_rightchild->Sinter(dt, sys, model, rng, wt);
			}
		}
		else{ //track

			if (MergeCondition()) {
				CheckSintering();
			}

			if (m_leftchild != NULL && m_rightchild != NULL) {
				m_leftchild->Sinter(dt, sys, model, rng, wt);
				m_rightchild->Sinter(dt, sys, model, rng, wt);
			}
		}

		UpdateCache();

		if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())
			m_children_sintering = SinteringLevel();
		else
			m_children_roundingLevel = RoundingLevel();

	}

}

//! Learn from bintree model.
/*!
* @brief       Adjusts one primary of a particle after all PAHs in it has been updated by KMC-ARS code.
*
* Analogous to the implementation in Adjust() in bintree model. However, since only one primary is updated
* here, not the whole particle, so there is no need to climb or descend the tree.
*
* @param[in]   old_vol   Volume of the primary before surface growth.
*/
void PAHPrimary::Adjust(const double old_vol)
{
	double dV(0.0);
	double m_vol_old = old_vol;

	//! If the distance between the centres of primary particles or the
	//! primary coordinates is tracked, the rate of change in the
	//! primary diameter is affected by its neighbours.
	if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()) {

		//! Particle with more than one primary.
		if (m_parent != NULL) {

			while (m_vol_old <= m_vol){

				//! Initialisation of variables to adjust the primary diameter if the
				//! distance between the centres of primary particles is tracked.
				double r_i = m_primarydiam / 2.0;	//!< Radius of primary particle i.
				double dr_max = 0.01*r_i;			//!< Maximum change in primary radius during internal step (1% of primary radius)
				double dr = 0.0;					//!< Change in radius of i

				//Calculate change in volume
				dV = dr_max * m_free_surf;

				//Calculate change in radius
				if (m_vol_old + dV > m_vol){
					dr = (m_vol - m_vol_old)*dr_max / dV;
				}
				else{
					dr = dr_max;
				}

				//Update primary diameter
				m_primarydiam = 2.0* (r_i + dr);

				//! Update free surface area
				this->UpdateOverlappingPrimary();

				m_vol_old += dV;
			}

			//if coordinates are tracked then update coordinate tracking properties
			if (m_pmodel->getTrackPrimaryCoordinates()){
				setRadius(m_primarydiam / 2.0);
				this->calcBoundSph();
				this->calcCOM();
			}

		}
		//! Single primary case: the primary diameter equals the
		//! spherical diameter                
		else {
			m_primarydiam = m_diam;

			if (m_pmodel->getTrackPrimaryCoordinates()) {
				setRadius(m_primarydiam / 2.0);
			}
		}
	}
	
}

//! From bintree model.
/*!
*  Print primary particle details and connectivity
*
*  @param[in]    surface			Primary connectivity
*  @param[in]    primary_diameter	Primary details
*  @param[in]    k					Particle counter
*/
void PAHPrimary::PrintPrimary(vector<fvector> &surface, vector<fvector> &primary_diameter, int k) const
{
	fvector node(10);
	fvector primary(10);

	if ((m_leftchild == NULL) && (m_rightchild == NULL)){
		//if leaf then print diameter
		primary[0] = k + 1;
		primary[1] = m_primarydiam;
		primary[2] = m_diam;
		primary[3] = m_primaryvol;
		primary[4] = m_vol;
		primary[5] = m_free_surf;
		vector<fvector> coords;
		this->GetPriCoords(coords);
		primary[6] = coords[0][0];
		primary[7] = coords[0][1];
		primary[8] = coords[0][2];
		primary[9] = coords[0][3];

		primary_diameter.push_back(primary);

		if (m_parent == NULL){	//single particle case
			node[0] = k + 1;
			node[1] = m_numprimary;
			node[2] = 0.0;
			node[3] = 1.0;
			node[4] = 0.0;
			node[5] = 0.0;
			node[6] = m_primarydiam / 2.0;
			node[7] = 0.0;
			node[8] = reinterpret_cast<uintptr_t>(this);	//print pointer
			node[9] = 0.0;

			surface.push_back(node);
		}
	}
	else {

		double r_i = m_leftparticle->m_primarydiam / 2.0;
		double r_j = m_rightparticle->m_primarydiam / 2.0;
		double d_ij = m_distance_centreToCentre;

		double x_ij = (d_ij*d_ij - r_j*r_j + r_i*r_i) / (2.0*d_ij);
		double R_ij = sqrt(r_i*r_i - x_ij*x_ij);	//!neck radius

		//if non-leaf node then print node and continue down the tree
		node[0] = k + 1;
		node[1] = m_numprimary;
		node[2] = m_children_surf;
		node[3] = m_children_sintering;
		node[4] = d_ij;
		node[5] = R_ij;
		node[6] = r_i;
		node[7] = r_j;
		node[8] = reinterpret_cast<uintptr_t>(m_leftparticle);	//print pointer
		node[9] = reinterpret_cast<uintptr_t>(m_rightparticle);	//print pointer

		surface.push_back(node);

		m_leftchild->PrintPrimary(surface, primary_diameter, k);
		m_rightchild->PrintPrimary(surface, primary_diameter, k);
	}
}

//! Used to calculate the geometric volume of a particle if primary coordinates are tracked.
//! Returns the cap volume of two connected primary parti
double PAHPrimary::CalcChildrenSumCap()
{
	if (m_leftchild != NULL && m_rightchild != NULL) {

		if (m_leftparticle != NULL && m_rightparticle != NULL) {

			double r_i = m_leftparticle->m_primarydiam / 2.0;
			double r_j = m_rightparticle->m_primarydiam / 2.0;
			double d_ij = m_distance_centreToCentre;
			double x_ij = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / (2.0*d_ij);
			double x_ji = (pow(d_ij, 2.0) - pow(r_i, 2.0) + pow(r_j, 2.0)) / (2.0*d_ij);
			double V_cap = 0.0;
			V_cap = V_cap + M_PI*(2 * pow(r_i, 3.0) + pow(x_ij, 3.0) - 3 * r_i*r_i*x_ij) / 3;
			V_cap = V_cap + M_PI*(2 * pow(r_j, 3.0) + pow(x_ji, 3.0) - 3 * r_j*r_j*x_ji) / 3;

			return V_cap;
		}
		else{
			return 0.0;
		}

	}
	else{ //this is a single primary
		
		return 0.0;

	}

}