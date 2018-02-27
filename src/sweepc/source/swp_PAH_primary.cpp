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
    m_numPAH(0),
    m_numprimary(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_roundingLevel(0),
    m_distance_centreToCentre(0.0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL),
	m_sint_rate(0.0), //hdy
	m_r(0.0), //hdy
	m_r2(0.0), //hdy
	m_r3(0.0) //hdy
{
	m_cen_bsph[0] = 0.0; //hdy
	m_cen_bsph[1] = 0.0; //hdy
	m_cen_bsph[2] = 0.0; //hdy

	m_cen_mass[0] = 0.0; //hdy
	m_cen_mass[1] = 0.0; //hdy
	m_cen_mass[2] = 0.0; //hdy
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
    m_numPAH(0),
    m_numprimary(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_roundingLevel(0),
    m_distance_centreToCentre(0.0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL),
	m_r(0.0), //hdy
	m_r2(0.0), //hdy
	m_r3(0.0) //hdy
{
	m_cen_bsph[0] = 0.0; //hdy
	m_cen_bsph[1] = 0.0; //hdy
	m_cen_bsph[2] = 0.0; //hdy

	m_cen_mass[0] = 0.0; //hdy
	m_cen_mass[1] = 0.0; //hdy
	m_cen_mass[2] = 0.0; //hdy

    // Other parts of the code check for a non-zero composition
    m_comp[0]=1;

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
    m_numPAH(0),
    m_numprimary(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_roundingLevel(0),
    m_distance_centreToCentre(0.0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL),
	m_r(0.0), //hdy
	m_r2(0.0), //hdy
	m_r3(0.0) //hdy
{
	m_cen_bsph[0] = 0.0; //hdy
	m_cen_bsph[1] = 0.0; //hdy
	m_cen_bsph[2] = 0.0; //hdy

	m_cen_mass[0] = 0.0; //hdy
	m_cen_mass[1] = 0.0; //hdy
	m_cen_mass[2] = 0.0; //hdy
    // Other parts of the code check for a non-zero composition
    m_comp[0]=1;
   
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
    m_numPAH(0),
    m_numprimary(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_roundingLevel(0),
    m_distance_centreToCentre(0.0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL),
	m_r(0.0), //hdy
	m_r2(0.0), //hdy
	m_r3(0.0) //hdy
{
	m_cen_bsph[0] = 0.0; //hdy
	m_cen_bsph[1] = 0.0; //hdy
	m_cen_bsph[2] = 0.0; //hdy

	m_cen_mass[0] = 0.0; //hdy
	m_cen_mass[1] = 0.0; //hdy
	m_cen_mass[2] = 0.0; //hdy

    m_comp[0]=1;
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
    m_PAHCollDiameter=source->m_PAHCollDiameter;
    SetSurfaceArea(source->SurfaceArea());
    m_time=source->m_time;
    m_PAHmass=source->m_PAHmass;
    m_leftchild=source->m_leftchild;
    m_rightchild=source->m_rightchild;
    m_leftparticle=source->m_leftparticle;
    m_rightparticle=source->m_rightparticle;
    m_parent=source->m_parent;
    SetMass(source->Mass());
    m_numPAH=source->m_numPAH;
    m_numprimary=source->m_numprimary;
    m_primarydiam=source->m_primarydiam;
    m_sqrtLW=source->m_sqrtLW;
    m_LdivW=source->m_LdivW;
    m_pmodel=source->m_pmodel;
    m_surf=source->m_surf;
    m_vol=source->m_vol;
    m_children_surf=source->m_children_surf;
    m_distance_centreToCentre = source->m_distance_centreToCentre;
	m_cen_bsph = source->m_cen_bsph; //hdy
	m_cen_mass = source->m_cen_mass; //hdy
	m_r = source->m_r; //hdy
	m_r2 = source->m_r2; //hdy
	m_r3 = source->m_r3; //hdy
    m_children_vol=source->m_children_vol;
    m_children_radius=source->m_children_radius;
    m_children_roundingLevel=source->m_children_roundingLevel;
    m_rightparticle_numPAH=source->m_rightparticle_numPAH;
    m_leftparticle_numPAH=source->m_leftparticle_numPAH;
    m_leftparticle_vol_old=source->m_leftparticle_vol_old;
    m_rightparticle_vol_old=source->m_rightparticle_vol_old;
    m_fdim=source->m_fdim;
    m_Rg=source->m_Rg;
    m_avg_coalesc=source->m_avg_coalesc;
    m_numcarbon = source->m_numcarbon;
    m_numH = source->m_numH;
    m_numOfEdgeC = source->m_numOfEdgeC;
    m_numOfRings = source->m_numOfRings;
	m_numOfRings5 = source->m_numOfRings5;
    m_values=source->m_values;
    m_comp=source->m_comp;
    m_sint_time=source->m_sint_time;

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

    //only one PAH in rhs or this particle -> condensation or inception process.
 //   if ( (rhsparticle->m_numPAH==1) || (m_numPAH==1 ))
 //   {
 //       if (rhsparticle->Numprimary()>1)
 //       {
 //           // rhsparticle is a soot particle but this paricle repesents a PAH
 //           // this paricle will condense on this rhsparticle
 //           // Get a copy of the rhs ready to add to the current particle
 //           PAHPrimary copy_rhs(*rhsparticle);
 //           PAHPrimary *target = copy_rhs.SelectRandomSubparticle(rng);

 //           target->m_PAH.insert(target->m_PAH.end(),m_PAH.begin(),m_PAH.end());
	//		target->UpdatePrimary();
 //           CopyParts(&copy_rhs);
 //           CopyTree(&copy_rhs);
 //       }
 //       else
 //       {
 //           // this paricle is a soot particle but rhsparticle repesents a PAH
 //           // rhsparticle will condense on this particle
 //           // particle has more then one primary select the primary where
 //           // the PAH condenses to and add it to the list
 //           if (m_leftchild!=NULL)
 //           {
 //               PAHPrimary *target = SelectRandomSubparticle(rng);

 //               target->m_PAH.insert(target->m_PAH.end(),rhsparticle->m_PAH.begin(),rhsparticle->m_PAH.end());
 //               target->UpdatePrimary();

 //           }
 //           else
 //           {
 //               //! rhsparticle is a gas-phase PAH but the this pointer may be
 //               //! pointing to a gas-phase PAH in which case this would be an
 //               //! inception event, or a single primary particle in which case
 //               //! it would be a condensation event.
 //               m_PAH.insert(m_PAH.end(),rhsparticle->m_PAH.begin(),rhsparticle->m_PAH.end());
 //               UpdatePrimary();
 //           }
 //       }
	//	UpdateCache();
 //       //Check the coalescence ratio
 //       CheckRounding();

	//}

    //else
    //{
        //coagulation process
		PAHPrimary *newleft = new PAHPrimary(m_time, *m_pmodel);
		PAHPrimary *newright = new PAHPrimary(m_time, *m_pmodel);
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

		//! It is assumed that primary pi from particle Pq and primary pj from
		//! particle Pq are in point contact and by default pi and pj are
		//! uniformly selected. If we track the coordinates of the primaries in
		//! a particle, we can do a smarter selection where pi and pj are
		//! determined by ballistic cluster-cluster aggregation (BCCA):
		//! R. Jullien, Transparency effects in cluster-cluster aggregation with
		//! linear trajectories, J. Phys. A 17 (1984) L771-L776.
		if (m_pmodel->getTrackPrimaryCoordinates()) {
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
			
		//select the primaries that are touching
		//this->m_leftparticle=m_leftchild->SelectRandomSubparticle(rng); //comment by hdy
		//this->m_rightparticle=m_rightchild->SelectRandomSubparticle(rng); //comment by hdy

		//set the sintertime for the new created primary particle
		SetSinteringTime(std::max(this->m_sint_time, rhsparticle->m_sint_time));
		//initialise the variables used to calculate the coalesence ratio
		m_children_vol = m_leftparticle->m_vol + m_rightparticle->m_vol;
		m_children_surf = (m_leftparticle->m_surf + m_rightparticle->m_surf);
		m_leftparticle_vol_old = m_leftparticle->m_vol;
		m_rightparticle_vol_old = m_rightparticle->m_vol;
		m_leftparticle_numPAH = m_leftparticle->m_numPAH;
		m_rightparticle_numPAH = m_rightparticle->m_numPAH;
		m_children_radius = pow(3.0 / (4.0*PI)*(m_children_vol), (ONE_THIRD));
		m_children_roundingLevel = CoalescenceLevel();
		m_distance_centreToCentre = m_leftparticle->m_primarydiam / 2.0 + m_rightparticle->m_primarydiam / 2.0;

		//****************************************************hdy****************************************************//
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
		//****************************************************hdy****************************************************//

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
			
	//}
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
        CheckRounding();

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
            CheckRounding();
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

//*********************************************hdy************************************************************//
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
//******************************************************hdy**********************************************************//

//******************************************************hdy**********************************************************//
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
//******************************************************hdy**********************************************************//

//******************************************************hdy**********************************************************//
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

	return Rg;
}
//******************************************************hdy**********************************************************//

//******************************************************hdy**********************************************************//
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
//******************************************************hdy**********************************************************//

//returns the CoalescenceLevel of the two primaries that are connected by this node
double PAHPrimary::CoalescenceLevel()
{
    if (m_leftparticle!=NULL)
    {

        double dV=0;


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
        if (m_leftparticle->m_numPAH - m_leftparticle_numPAH<2)
        {
            //calculate the volume change of the left particle in the last timestep
            dV+= m_leftparticle->m_vol - m_leftparticle_vol_old;
        }

        //adjust the number of PAHs and vol for the next step
        m_leftparticle_numPAH=m_leftparticle->m_numPAH;
        m_leftparticle_vol_old=m_leftparticle->m_vol;

        //make sure that the volume growth does not arise from a coalescence event
        //(number of PAHs must remain const, or increase by1 in case of a condensation
        // event
        if (m_rightparticle->m_numPAH-m_rightparticle_numPAH<2)
        {
            //calculate the volume change of the right particle in the last timestep
            dV+= m_rightparticle->m_vol - m_rightparticle_vol_old;
        }

        //adjust the number of PAHs and vol for the next step
        m_rightparticle_numPAH=m_rightparticle->m_numPAH;
        m_rightparticle_vol_old=m_rightparticle->m_vol;

        //update the children volume, this value is used to calculate dV n the next timestep
        m_children_vol=m_rightparticle->m_vol+m_leftparticle->m_vol;
        //calculate dS, at the moment it is assumed that the particles always grow
        double ct=m_pmodel->Components(0)->CoalescThresh();
        const double dS=dV*ct/m_children_radius;
        //update the radius for the next event
        m_children_radius=pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
        m_children_surf+=dS;

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
            m_avg_coalesc=0;
}

//merges the left and the right primary particle to one primary particle
void PAHPrimary::Merge()
{
//    if(this->m_numprimary>5)
 //       cout<<"tset";
 //      PrintTree("before.inp");

	////******************************hdy********************************************//
	////pointer to the new (merged) primary //hdy
	//PAHPrimary *new_prim; //hdy
	////initialise pointer to the smaller of the merging primaries //hdy
	//PAHPrimary *small_prim; //hdy
	////******************************hdy********************************************//

    vector<PAH>::const_iterator j;
      //make sure this primary has children to merge
	  if(m_leftchild!=NULL)
           {
			   
			   ////****************************************************hdy************************************************************//
			   ////if the centre to centre distance is tracked we need to know which is the smaller primary of the merging pair
			   //if (m_leftparticle->m_primarydiam > m_rightparticle->m_primarydiam){
				  // small_prim = m_rightparticle;
			   //}
			   //else{
				  // small_prim = m_leftparticle;
			   //}
			   //
			   //double dV = small_prim->m_vol; //volume of the merging (smaller) primary 
			   ////****************************************************hdy************************************************************//

		if ( m_leftchild==m_leftparticle && m_rightchild==m_rightparticle)
		{
            //this node has only two primaries in its subtree
            //it is possible that this node is not the root node and belongs to a bigger particle
            //copy the PAHs of both children to the parent node

            //m_PAH is empty, therefore no need to append
            m_PAH=m_rightparticle->m_PAH; //comment by hdy

            m_PAH.insert(this->m_PAH.end(),m_leftparticle->m_PAH.begin(),m_leftparticle->m_PAH.end()); //comment by hdy
            //update the pointers that pointed to the two former children
			ChangePointer(m_leftchild,this); //comment by hdy
			ChangePointer(m_rightchild,this); //comment by hdy

			////******************************************************hdy*********************************************************//
			//// This node has only two primaries in its subtree, it is possible
			//// that this node is not the root node and belongs to a bigger
			//// particle.		

			//// Sum up the components first
			//for (size_t i = 0; i != m_comp.size(); i++) {
			//	m_comp[i] = m_leftparticle->Composition(i) +
			//		m_rightparticle->Composition(i);
			//}

			////m_primarydiam of the new particle is the larger of the diameters of the merging particles
			////if the centre to centre separation isn't tracked this will be changed to the spherical equivalent in the call to UpdatePrimary
			//m_primarydiam = max(m_leftparticle->m_primarydiam, m_rightparticle->m_primarydiam);

			////new primary
			//new_prim = this;

			//// Update the pointers that pointed to the two former children
			//if (!m_pmodel->getTrackPrimarySeparation()) {
			//	ChangePointer(m_leftchild, this);
			//	ChangePointer(m_rightchild, this);
			//}
			//else{
			//	// If the priamry separation is tracked then re-estimation of new sintering levels also
			//	// requires the centre-centre separation and the smaller primary
			//	ChangePointer(m_leftchild, this, m_distance_centreToCentre, small_prim);
			//	ChangePointer(m_rightchild, this, m_distance_centreToCentre, small_prim);
			//}
			////******************************************************hdy*********************************************************//

			//delete the children (destructor is recursive for this class)
            delete m_leftchild;
            delete m_rightchild;
			m_leftchild=NULL;
			m_rightchild=NULL;
			m_leftparticle=NULL;
			m_rightparticle=NULL;

            //set the children properties to zero, this node has no more children
            ResetChildrenProperties();
            UpdatePrimary();

			////******************************************************hdy*********************************************************//
			//// Only update the cache on m_parent if the sintering level of
			//// m_parent if the sintering level won't call a merge on the
			//// parent node. Otherwise, the *this* memory address could be
			//// removed from the tree and segmentation faults will result!
			//if (m_parent != NULL) {
			//	if (!m_parent->MergeCondition()) {
			//		m_parent->UpdateCache();
			//	}
			//}
			////******************************************************hdy*********************************************************//

		}


		else
		{
			if (m_leftchild->m_numprimary<m_rightchild->m_numprimary)
			{
				//append to left subtree because there are fewer primaries
                //this is only to keep the tree balanced
				PAHPrimary *oldleftparticle=m_leftparticle;

				////******************************************************hdy*********************************************************//
				//for (size_t i = 0; i != m_comp.size(); i++) {
				//	m_rightparticle->m_comp[i] =
				//		m_leftparticle->Composition(i) +
				//		m_rightparticle->Composition(i);
				//}
				////******************************************************hdy*********************************************************//

                //copy the PAHs

				//for (j=oldleftparticle->m_PAH.begin(); j!=oldleftparticle->m_PAH.end(); ++j) {
				//	m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(),PAH(*j));
				//}

				////******************************************************hdy*********************************************************//
				////m_primarydiam of the new particle is the larger of the diameters of the merging particles
				////if the centre to centre separation isn't tracked this will be changed to the spherical equivalent in the call to UpdatePrimary
				//m_rightparticle->m_primarydiam = max(m_leftparticle->m_primarydiam, m_rightparticle->m_primarydiam);

				////the new (merged) primary
				//new_prim = m_rightparticle;

				//m_rightparticle->UpdatePrimary();

				//// Set the pointers from the leftprimary to the rightprimary
				//// (this will be the new bigger primary)
				//if (!m_pmodel->getTrackPrimarySeparation()) {
				//	oldleftparticle->ChangePointer(oldleftparticle, m_rightparticle);
				//	m_rightparticle->ChangePointer(m_rightparticle, m_rightparticle);
				//}
				//else{
				//	// If the priamry separation is tracked then re-estimation of new sintering levels also
				//	// requires the centre-centre separation and the smaller primary
				//	oldleftparticle->ChangePointer(oldleftparticle, m_rightparticle, m_distance_centreToCentre, small_prim);
				//	m_rightparticle->ChangePointer(m_rightparticle, m_rightparticle, m_distance_centreToCentre, small_prim);
				//} 
				////******************************************************hdy*********************************************************//

                m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(),oldleftparticle->m_PAH.begin(),oldleftparticle->m_PAH.end()); //comment by hdy
                m_rightparticle->UpdatePrimary(); //comment by hdy
                //set the pointers from the leftprimary to the rightprimary
                //this will be the new bigger primary
				oldleftparticle->ChangePointer(oldleftparticle,m_rightparticle); //comment by hdy
                m_rightparticle->ChangePointer(m_rightparticle,m_rightparticle); //comment by hdy

                //set the pointer to the parent node
				if (oldleftparticle->m_parent->m_leftchild==oldleftparticle)
				{
					oldleftparticle->m_parent->m_leftchild=m_rightchild;
				}
				else
				{
					oldleftparticle->m_parent->m_rightchild=m_rightchild;
				}
				m_rightchild->m_parent=oldleftparticle->m_parent;


				PAHPrimary *oldleftchild=m_leftchild;
				PAHPrimary *oldparent=m_parent;

				//copy the properties of the former leftchild to this node
				// so that it can be removed from the aggregate tree structure
				CopyParts(oldleftchild);

				// Now break the links to the tree structure in oldleftchild and free it
				oldleftchild->m_leftchild = NULL;
				oldleftchild->m_rightchild = NULL;
				delete oldleftchild;

				m_parent=oldparent;
                if (m_leftchild!=NULL)
                {
				    m_rightchild->m_parent=this;
				    m_leftchild->m_parent=this;
                }

                delete oldleftparticle;

			}

			else
			{
				//append to right subtree
				PAHPrimary *oldrightparticle=m_rightparticle;

			//	for (j=oldrightparticle->m_PAH.begin(); j!=oldrightparticle->m_PAH.end(); ++j) {
			//		m_leftparticle->m_PAH.insert(m_leftparticle->m_PAH.end(),PAH(*j));
			//	}

				////******************************************************hdy*********************************************************//
				//for (size_t i = 0; i != m_comp.size(); i++) { //hdy
				//	m_leftparticle->m_comp[i] =
				//		m_leftparticle->Composition(i) +
				//		m_rightparticle->Composition(i);
				//}

				////m_primarydiam of the new particle is the larger of the diameters of the merging particles
				////if the centre to centre separation isn't tracked this will be changed to the spherical equivalent in the call to UpdatePrimary
				//m_leftparticle->m_primarydiam = max(m_leftparticle->m_primarydiam, m_rightparticle->m_primarydiam); 

				////the new (merged) primary
				//new_prim = m_leftparticle; 
				////******************************************************hdy*********************************************************//

                m_leftparticle->m_PAH.insert(m_leftparticle->m_PAH.end(),oldrightparticle->m_PAH.begin(),oldrightparticle->m_PAH.end()); //comment by hdy
                m_leftparticle->UpdatePrimary();

				////******************************************************hdy*********************************************************//
				//// All pointers to m_leftparticle now point to oldright particle
				//if (!m_pmodel->getTrackPrimarySeparation()) {
				//	oldrightparticle->ChangePointer(oldrightparticle, m_leftparticle);
				//	m_leftparticle->ChangePointer(m_leftparticle, m_leftparticle);
				//}
				//else{
				//	// If the priamry separation is tracked then re-estimation of new sintering levels also
				//	// requires the centre-centre separation and the smaller primary
				//	oldrightparticle->ChangePointer(oldrightparticle, m_leftparticle, m_distance_centreToCentre, small_prim);
				//	m_leftparticle->ChangePointer(m_leftparticle, m_leftparticle, m_distance_centreToCentre, small_prim);
				//}
				////******************************************************hdy*********************************************************//

                oldrightparticle->ChangePointer(oldrightparticle,m_leftparticle); //comment by hdy
                m_leftparticle->ChangePointer(m_leftparticle,m_leftparticle); //comment by hdy

				if (oldrightparticle->m_parent->m_leftchild==oldrightparticle)
				{
					oldrightparticle->m_parent->m_leftchild=m_leftchild;
				}
				else
				{
					oldrightparticle->m_parent->m_rightchild=m_leftchild;
				}
				m_leftchild->m_parent=oldrightparticle->m_parent;

				//ReleaseMem(oldrightparticle);
				PAHPrimary *oldrightchild=m_rightchild;
				PAHPrimary *oldparent=m_parent;

                //copy the properties of the former leftchild to this node
                // so that it can be removed from the aggregate tree structure
				CopyParts(oldrightchild);

                // Now break the links to the tree structure in oldrightchild and free it
                oldrightchild->m_leftchild = NULL;
                oldrightchild->m_rightchild = NULL;
                delete oldrightchild;

				m_parent=oldparent;
                if (m_leftchild!=NULL)
                {
				    m_rightchild->m_parent=this;
				    m_leftchild->m_parent=this;
                }
                delete oldrightparticle;
			}
		}

		////******************************************************hdy*********************************************************//
		//if (new_prim->m_parent != NULL){
		//	double sumterm = 0.0;
		//	double r_i = new_prim->m_primarydiam / 2.0;
		//	//get contribution of neighours
		//	new_prim->SumNeighbours(new_prim, sumterm);
		//	//change in radius
		//	double delta_r_i = dV / (4 * M_PI*r_i*r_i + M_PI*r_i*sumterm);
		//	//update the particle separations and calculate the new free surface of the particle
		//	double free_surface_term = 0.0;
		//	std::set<void*> duplicates;
		//	new_prim->UpdateConnectivity(new_prim, duplicates, delta_r_i, free_surface_term);

		//	//update primary diameter
		//	new_prim->m_primarydiam = 2.0*(r_i + delta_r_i);
		//	//update the free surface area
		//	new_prim->m_free_surf = max(M_PI*new_prim->m_primarydiam*new_prim->m_primarydiam - 2 * M_PI*free_surface_term, 0.0);
		//}
		////******************************************************hdy*********************************************************//

        UpdateCache();
  //      PrintTree("after.inp");
        }
}

//******************************************hdy*****************************************************//
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
		if (!m_pmodel->getTrackPrimarySeparation()) {
			condition = (m_children_sintering > 0.95);
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

			if (d_ij <= 0.0){
				//ensures that particles are merged if the sintering step overshoots
				condition = true;
			}
			else{
				condition = ((pow(d_ij, 2.0) - pow(max(r_i, r_j), 2.0) + pow(min(r_i, r_j), 2.0)) / (2.0*d_ij) <= 0.0);
			}
		}
	}

	return condition;
}
//******************************************hdy*****************************************************//

//******************************************hdy*****************************************************//
/*!
*
* Sinters the node in the binary tree for time dt
*
* @param[in]   dt      Time for which to sinter
* @param[in]   sys     Environment for particles
* @param[in]   model   Sintering model to apply
* @param[in]   rng     Random number generator
* @param[in]   wt      Statistical weight
*/
void PAHPrimary::SinterNode(
	double dt,
	Cell &sys,
	const Processes::SinteringModel &model,
	rng_type &rng,
	double wt
	) {

	// Declare time step variables.
	double t1 = 0.0, delt = 0.0, tstop = dt;
	double r = 0.0;

	// The scale parameter discretises the delta-S when using
	// the Poisson distribution.  This allows a smoother change
	// (smaller scale = higher precision).
	double scale = 0.01;

	//! The sintering model depends on whether the distance between the centres
	//! of primary particles or the coordinates of the primary particles are
	//! tracked. If tracked, sintering results in a decrease in the distance
	//! between the primaries and an increase in their diameters. If not,
	//! sintering results in a decrease in the common surface between the
	//! primaries.
	if (!(m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())) {

		// Calculate the spherical surface
		const double spherical_surface = 4 * PI*m_children_radius*m_children_radius;

		// Define the maximum allowed change in surface
		// area in one internal time step (10% spherical surface).
		double dAmax = 0.1 * spherical_surface;

		// Perform integration loop.
		while (t1 < tstop)
		{
			// Calculate sintering rate.
			r = model.Rate(m_time + t1, sys, *this);

			if (r > 0) {
				// Calculate next time-step end point so that the
				// surface area changes by no more than dAmax.
				delt = dAmax / max(r, 1.0e-300);

				// Approximate sintering by a poisson process.  Calculate
				// number of poisson events.
				double mean;

				if (tstop > (t1 + delt)) {
					// A sub-step, we have changed surface by dAmax, on average
					mean = 1.0 / scale;
				}
				else {
					// Step until end.  Calculate degree of sintering explicitly.
					mean = r * (tstop - t1) / (scale*dAmax);
				}
				boost::random::poisson_distribution<unsigned, double> repeatDistribution(mean);
				const unsigned n = repeatDistribution(rng);

				// Adjust the surface area.
				if (n > 0) {
					m_children_surf -= (double)n * scale * dAmax;

					// Check that primary is not completely sintered.
					if (m_children_surf <= spherical_surface) {
						m_children_surf = spherical_surface;
						break;
					}
				}

				// Set t1 for next time step.
				t1 += delt;
			}

		}

	}
	else {
		//! Define the maximum allowed change (1%) in the distance between the
		//! centres of primary particles in one internal time step. In the case
		//! of pure sintering it was found that if the allowed change is too
		//! large (~10%) a significant error is incurred in the final spherical
		//! volume determined through comparisons with the mass-derived volume.
		//! Note that the smaller the distance is, the smaller the changes are.
		double dd_ij_Max = m_distance_centreToCentre / 100.0;

		while (t1 < tstop) {
			double r_i = this->m_leftparticle->m_primarydiam / 2.0;
			double r_j = this->m_rightparticle->m_primarydiam / 2.0;
			double d_ij = m_distance_centreToCentre;

			//assert(d_ij<= r_i+r_j);	//debug

			//! Definition of variables for conciseness.
			double d_ij2 = pow(d_ij, 2.0);
			double r_i2 = pow(r_i, 2.0);
			double r_j2 = pow(r_j, 2.0);
			double r_i4 = pow(r_i, 4.0);
			double r_j4 = pow(r_j, 4.0);

			//! Continue if primaries have not coalesced
			if (!MergeCondition()) {

				//! Due to rounding, x_i and x_j are sometimes calculated to be larger 
				//! than the respective primary radii resulting in a negative neck area.
				//! Therefore we take the smaller of x_i and r_i.
				double x_i = min((d_ij2 - r_j2 + r_i2) / (2.0 * d_ij), r_i); //!< Eq. (3b).
				double x_j = min((d_ij2 - r_i2 + r_j2) / (2.0 * d_ij), r_j); //!< Eq. (3b).
				double A_n = M_PI * (r_i2 - pow(x_i, 2.0));        //!< Eq. (4).
				double A_i = 2.0 * M_PI * (r_i2 + r_i * x_i);      //!< Eq. (6).
				double A_j = 2.0 * M_PI * (r_j2 + r_j * x_j);      //!< Eq. (6).

				//declare variables
				double dd_ij_dt = 0.0;
				double R_n = 0.0;
				double r4_tau = 0.0;
				double tau = 0.0;
				double gamma_eta = 0.0;

				//! Sintering model dependent part
				//! Viscous flow model
				if (model.Type() == Processes::SinteringModel::ViscousFlow){

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

					//! Grain boundary diffusion model 
				}
				else if (model.Type() == Processes::SinteringModel::GBD){

					//! If the particles are in point contact set an initial neck radius of 1% 
					//! of the smaller primary radius, otherwise dd_ij_dt is undefined
					if (A_n == 0.0){
						R_n = 0.01*min(r_i, r_j);
						A_n = M_PI * R_n * R_n;
						x_i = sqrt(r_i2 - R_n * R_n);
						x_j = sqrt(r_j2 - R_n * R_n);
					}
					else{
						R_n = sqrt(A_n / M_PI);
					}

					// The primary radius in the numerator cancels with the diameter dependence of tau
					// so we can calculate this for only one of the primaries.
					// In the SintTime the diameter is calculated as 6.0 * m_vol / m_surf
					// so r = 3.0 * m_vol / m_surf
					double r4 = pow(3.0 * m_leftparticle->m_vol / m_leftparticle->m_surf, 4.0);
					r4_tau = r4 / model.SintTime(sys, *this->m_leftparticle);

					//! J Aerosol Sci 46:7-19 (2012) Eq. (A6)
					//! dx_i_dt + dx_j_dt
					//! (this is missing a minus sign, which is accounted for below)
					dd_ij_dt = r4_tau * (1 / (r_i - x_i) + 1 / (r_j - x_j) - 2 / R_n) / A_n;

					//Other model are not coded
				}
				else{
					std::cout << "Sintering model not coded" << endl;
					break;
				}

				//! Account for multiple neighbours
				double sumterm_i = 0.0;
				double sumterm_j = 0.0;
				//get contribution from neighbours working up the binary tree
				m_leftparticle->SumNeighbours(m_leftparticle, sumterm_i);
				m_rightparticle->SumNeighbours(m_rightparticle, sumterm_j);
				//subtract mutual contribution from sumterm
				sumterm_i = max(sumterm_i - (x_i - r_i)*(x_i - r_i) / x_i, 0.0);
				sumterm_j = max(sumterm_j - (x_j - r_j)*(x_j - r_j) / x_j, 0.0);

				//! Modified A_i and A_j
				A_i = M_PI * (2 * r_i*r_i + 2 * r_i*x_i + r_i*sumterm_i);
				A_j = M_PI * (2 * r_j*r_j + 2 * r_j*x_j + r_j*sumterm_j);

				//! The expression for B_i in Eq. (8) is wrong. By combining
				//! Eqs. (5) and (7), we can obtain two equations which are
				//! functions of r_i and r_j. Subsequently combined these two
				//! equations and used Wolfram Alpha to rearrange equation in terms
				//! of r_i (and r_j).
				//!
				//! @todo Remove derivation and replace with reference to preprint
				//!       or paper if results do get published.
				double B_i = (-r_j*A_n*A_n - x_j*A_j*A_n) / (A_i*A_j*d_ij + r_i*A_j*A_n + r_j*A_i*A_n);
				double B_j = (-r_i*A_n*A_n - x_i*A_i*A_n) / (A_j*A_i*d_ij + r_j*A_i*A_n + r_i*A_j*A_n);

				double V_i = 2.0 / 3.0 * M_PI * pow(r_i, 3.0) + M_PI * r_i2 * x_i - 1.0 / 3.0 * M_PI * pow(x_i, 3.0); //!< Eq. (3a).
				double V_j = 2.0 / 3.0 * M_PI * pow(r_j, 3.0) + M_PI * r_j2 * x_j - 1.0 / 3.0 * M_PI * pow(x_j, 3.0); //!< Eq. (3a).

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

				if (m_pmodel->getTrackPrimarySeparation()) {
					m_distance_centreToCentre -= (double)n * scale * dd_ij_Max; //!< Sintering decreases d_ij hence the negative sign.
				}

				//! Needs to be tested.
				else if (m_pmodel->getTrackPrimaryCoordinates()) {
					double dx = this->m_leftparticle->m_cen_bsph[0] - this->m_rightparticle->m_cen_bsph[0];
					double dy = this->m_leftparticle->m_cen_bsph[1] - this->m_rightparticle->m_cen_bsph[1];
					double dz = this->m_leftparticle->m_cen_bsph[2] - this->m_rightparticle->m_cen_bsph[2];

					double Fraction = (m_distance_centreToCentre - (double)n * scale * dd_ij_Max) / m_distance_centreToCentre;

					this->m_rightchild->Translate(Fraction * dx, Fraction * dy, Fraction * dz);

					m_distance_centreToCentre = sqrt(dx * dx + dy * dy + dz * dz);
				}

				//! Change in primary radii
				double delta_r_i = -(double)n * scale * B_i * dd_ij_Max;  //!< Eq. (8).
				double delta_r_j = -(double)n * scale * B_j * dd_ij_Max; //!< Eq. (8).


				//! Adjust separation of neighbours (not currently sintering) and free surface area
				//! due to the increase in primary radius
				if (!MergeCondition()) {
					double free_surf_i = 0.0;
					double free_surf_j = 0.0;

					//adjust separation with neighbours (ignoring p_j) and return the new free surface
					m_leftparticle->UpdateConnectivity(m_leftparticle, delta_r_i, free_surf_i, m_rightparticle);
					//add i,j term
					x_i = (pow(m_distance_centreToCentre, 2.0) - pow(r_j + delta_r_j, 2.0) + pow(r_i + delta_r_i, 2.0)) / (2.0*m_distance_centreToCentre);
					free_surf_i += (r_i + delta_r_i)*(r_i + delta_r_i) - (r_i + delta_r_i)*x_i;

					//adjust separation with neighbours (ignoring p_i) and return the new free surface
					m_rightparticle->UpdateConnectivity(m_rightparticle, delta_r_j, free_surf_j, m_leftparticle);
					//add i,j term
					x_j = (pow(m_distance_centreToCentre, 2.0) - pow(r_i + delta_r_i, 2.0) + pow(r_j + delta_r_j, 2.0)) / (2.0*m_distance_centreToCentre);
					free_surf_j += (r_j + delta_r_j)*(r_j + delta_r_j) - (r_j + delta_r_j)*x_j;

					//update the free surface area 
					m_leftparticle->m_free_surf = max(M_PI*m_leftparticle->m_primarydiam*m_leftparticle->m_primarydiam - 2 * M_PI*free_surf_i, 0.0);
					//update the free surface area
					m_rightparticle->m_free_surf = max(M_PI*m_rightparticle->m_primarydiam*m_rightparticle->m_primarydiam - 2 * M_PI*free_surf_j, 0.0);
				}

				//! Adjust primary radii
				this->m_leftparticle->m_primarydiam += 2.0 * delta_r_i;
				this->m_rightparticle->m_primarydiam += 2.0 * delta_r_j;

				t1 += delt;

				//! Return some sintering rate
				r = dd_ij_dt;
			}
			else {
				break; //!do not continue to sinter.
			}
		}

	}

	m_children_sintering = SinteringLevel();

	m_sint_rate = r;

}
//***************************************************hdy********************************************//

//***************************************************hdy********************************************//
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

		//if the centre-centre separation is not tracked the sintering level is calculated as per
		//Shekar et al. (2012)
		if (!m_pmodel->getTrackPrimarySeparation()) {
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

			//if centre-centre separation is tracked the sintering level is calculated as
			//s = ((d_min/d_ij) - (d_min/d_max))/(1-d_min/d_max)
			//where d_min = d_ij at merger, and d_max = r_i + r_j
		}
		else{
			if (m_leftparticle != NULL && m_rightparticle != NULL) {
				double r_i = m_leftparticle->m_primarydiam / 2.0;
				double r_j = m_rightparticle->m_primarydiam / 2.0;
				double d_ij = m_distance_centreToCentre;

				//calculate the merger condition
				double d_min = sqrt(pow(max(r_i, r_j), 2.0) - pow(min(r_i, r_j), 2.0));

				slevel = (d_min / d_ij - d_min / (r_i + r_j)) / (1 - d_min / (r_i + r_j));
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
//***************************************************hdy********************************************//

//***************************************************hdy********************************************//
void PAHPrimary::SumNeighbours(PAHPrimary *prim, double &sumterm) {

	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;

	//check if a neighbour of prim
	if (m_parent->m_leftparticle == prim) {
		//right particle is a neighbour
		r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
	}
	else if (m_parent->m_rightparticle == prim) {
		//left particle is a neighbour
		r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
	}
	else {
		//not a neighbour
		r_j = 0.0;
	}

	//if the node connects a neighbourt then calculate the summation term
	if (r_j > 0.0){
		//the volumes and radii of neighbours remain unchanged
		//the centre to centre separations increase to allow growth of the primary
		x_ij = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / (2.0*d_ij);
		sumterm += x_ij - 2.0*r_i + pow(r_i, 2.0) / x_ij;
	}

	//continue working up the binary tree
	if (m_parent->m_parent != NULL){
		m_parent->SumNeighbours(prim, sumterm);
	}
}
//***************************************************hdy********************************************//

//***************************************************hdy********************************************//
void PAHPrimary::UpdateConnectivity(PAHPrimary *prim, std::set<void*> &primaryUniqueAddresses, double delta_r, double &sumterm){

	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;
	double dx, dy, dz, Fraction;

	//! Check if a neighbour of prim.
	if (m_parent->m_leftparticle == prim) {
		//! Right particle is a neighbour.
		r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
	}
	else if (m_parent->m_rightparticle == prim) {
		//! Left particle is a neighbour.
		r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
	}
	else {
		//! Not a neighbour.
		r_j = 0.0;
	}

	if (r_j > 0.0){

		x_ij = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / (2.0*d_ij);

		//! Update centre to centre separation making sure centre to centre
		//! separation remains smaller than the sum of the radii.
		//! Needs to be tested.
		if (m_pmodel->getTrackPrimaryCoordinates()) {
			Fraction = min(d_ij + r_i * delta_r / x_ij, r_i + r_j + delta_r) / d_ij - 1.0;

			if (m_parent->m_leftparticle == prim) {
				dx = m_parent->m_rightparticle->m_cen_bsph[0] - m_parent->m_leftparticle->m_cen_bsph[0];
				dy = m_parent->m_rightparticle->m_cen_bsph[1] - m_parent->m_leftparticle->m_cen_bsph[1];
				dz = m_parent->m_rightparticle->m_cen_bsph[2] - m_parent->m_leftparticle->m_cen_bsph[2];

				if (primaryUniqueAddresses.find(m_parent->m_rightparticle) == primaryUniqueAddresses.end()) {
					m_parent->m_rightparticle->Translate(Fraction * dx, Fraction * dy, Fraction * dz);
					dx = m_parent->m_rightparticle->m_cen_bsph[0] - m_parent->m_leftparticle->m_cen_bsph[0];
					dy = m_parent->m_rightparticle->m_cen_bsph[1] - m_parent->m_leftparticle->m_cen_bsph[1];
					dz = m_parent->m_rightparticle->m_cen_bsph[2] - m_parent->m_leftparticle->m_cen_bsph[2];
					m_parent->m_distance_centreToCentre = sqrt(dx * dx + dy * dy + dz * dz);
					primaryUniqueAddresses.insert(m_parent->m_rightparticle);
					m_parent->m_rightparticle->UpdateConnectivity(m_parent->m_rightparticle, primaryUniqueAddresses, delta_r, sumterm);
				}
			}
			else if (m_parent->m_rightparticle == prim) {
				dx = m_parent->m_leftparticle->m_cen_bsph[0] - m_parent->m_rightparticle->m_cen_bsph[0];
				dy = m_parent->m_leftparticle->m_cen_bsph[1] - m_parent->m_rightparticle->m_cen_bsph[1];
				dz = m_parent->m_leftparticle->m_cen_bsph[2] - m_parent->m_rightparticle->m_cen_bsph[2];

				if (primaryUniqueAddresses.find(m_parent->m_leftparticle) == primaryUniqueAddresses.end()) {
					m_parent->m_leftparticle->Translate(Fraction * dx, Fraction * dy, Fraction * dz);
					dx = m_parent->m_leftparticle->m_cen_bsph[0] - m_parent->m_rightparticle->m_cen_bsph[0];
					dy = m_parent->m_leftparticle->m_cen_bsph[1] - m_parent->m_rightparticle->m_cen_bsph[1];
					dz = m_parent->m_leftparticle->m_cen_bsph[2] - m_parent->m_rightparticle->m_cen_bsph[2];
					m_parent->m_distance_centreToCentre = sqrt(dx * dx + dy * dy + dz * dz);
					primaryUniqueAddresses.insert(m_parent->m_leftparticle);
					m_parent->m_leftparticle->UpdateConnectivity(m_parent->m_leftparticle, primaryUniqueAddresses, delta_r, sumterm);
				}
			}
		}
		else if (m_pmodel->getTrackPrimarySeparation()) {
			m_parent->m_distance_centreToCentre = min(d_ij + r_i * delta_r / x_ij, r_i + r_j + delta_r);
		}

		//! Calculate term for the free surface area.
		d_ij = m_parent->m_distance_centreToCentre;
		x_ij = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i + delta_r, 2.0)) / (2.0*d_ij);
		sumterm += (r_i + delta_r)*(r_i + delta_r) - (r_i + delta_r)*x_ij;
	}

	//! Continue working up the binary tree.
	if (m_parent->m_parent != NULL){
		m_parent->UpdateConnectivity(prim, primaryUniqueAddresses, delta_r, sumterm);
	}
}
//***********************************************************hdy*******************************************************//

//***********************************************************hdy*******************************************************//
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
void PAHPrimary::UpdateConnectivity(PAHPrimary *prim, double delta_r, double &sumterm, PAHPrimary *prim_ignore){

	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;

	//check if a neighbour of prim
	if (m_parent->m_leftparticle == prim && m_parent->m_rightparticle != prim_ignore) {
		//right particle is a neighbour
		r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
	}
	else if (m_parent->m_rightparticle == prim &&  m_parent->m_leftparticle != prim_ignore) {
		//left particle is a neighbour
		r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
	}
	else {
		//not a neighbour
		r_j = 0.0;
	}

	if (r_j > 0.0){

		x_ij = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / (2.0*d_ij);
		//update centre to centre separation
		//making sure centre to centre separation remains smaller than some of radii
		m_parent->m_distance_centreToCentre = min(d_ij + r_i * delta_r / x_ij, r_i + r_j + delta_r);

		//calculate term for free surface area 
		d_ij = m_parent->m_distance_centreToCentre;
		x_ij = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i + delta_r, 2.0)) / (2.0*d_ij);
		sumterm += (r_i + delta_r)*(r_i + delta_r) - (r_i + delta_r)*x_ij;
	}

	//continue working up the binary tree
	if (m_parent->m_parent != NULL){
		m_parent->UpdateConnectivity(prim, delta_r, sumterm, prim_ignore);
	}
}

//******************************************hdy************************************************************//
bool PAHPrimary::isLeaf(void) const
{
	return (m_leftchild == NULL) && (m_rightchild == NULL);
}
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
//! Returns the bounding-sphere centre.
const Coords::Vector &PAHPrimary::boundSphCentre(void) const
{
	return m_cen_bsph;
}
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
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
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
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
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
//! Put the bounding-sphere at the origin.
void PAHPrimary::centreBoundSph(void)
{
	Translate(-m_cen_bsph[0], -m_cen_bsph[1], -m_cen_bsph[2]);
}
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
//! Put the centre-of-mass at the origin.
void PAHPrimary::centreCOM(void)
{
	Translate(-m_cen_mass[0], -m_cen_mass[1], -m_cen_mass[2]);
}
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
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
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
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
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
//! Returns the bounding sphere radius.
double PAHPrimary::Radius(void) const
{
	return m_r;
}
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
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
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
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
//******************************************hdy************************************************************//

//******************************************hdy************************************************************//
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
//******************************************hdy************************************************************//

void PAHPrimary::ReleaseMem()
{ 
    m_PAH.clear();
}
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
            m_children_surf=sphericalsurface/(m_children_roundingLevel*0.2063+0.7937);    //sphericalsurface/(m_children_roundingLevel*(1-2^(-1/3))+2^(-1/3))
		}
		if(m_leftparticle==source){
			m_leftparticle=target;
            double sphericalsurface=
                4*PI*pow(3*(m_leftparticle->Volume()+m_rightparticle->Volume())/(4*PI),TWO_THIRDS);
            m_children_surf=sphericalsurface/(m_children_roundingLevel*0.2063+0.7937);    //sphericalsurface/(m_children_roundingLevel*(1-2^(-1/3))+2^(-1/3))

		}

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->ChangePointer(source,target);
    }

}

//******************************************hdy************************************************************//
/*!
* @brief       Changes pointer from source to target when centre-centre separation is tracked
*
* If a primary neighbours the smaller of the merging pair the centre to centre separation is
* re-estmated as the smaller of the sum of the separation or the sum of primary radii.
*
* @param[in] source Pointer to the original particle
* @param[in] target Pointer to the new particle
* @param[in] centre to centre separation of the merging primaries
* @param[in] Pointer to the smaller of merging primaries
*/
void PAHPrimary::ChangePointer(PAHPrimary *source, PAHPrimary *target, double d_ij, PAHPrimary *small_prim)
{
	if (m_rightparticle == source) {
		m_rightparticle = target;
		//if the neighbour is the smaller of the merging primaries then update the centre to centre distance
		if (source == small_prim){
			//estimate new separation as the smaller of the sum of the two separartions or the sum of primary radii
			m_distance_centreToCentre = min(m_distance_centreToCentre + d_ij, m_rightparticle->m_primarydiam / 2.0 + m_leftparticle->m_primarydiam / 2.0);
		}
	}
	if (m_leftparticle == source){
		m_leftparticle = target;
		//if the neighbour is the smaller of the merging primaries then update the centre to centre distance
		if (source == small_prim){
			//estimate new separation as the smaller of the sum of the two separartions or the sum of primary radii
			m_distance_centreToCentre = min(m_distance_centreToCentre + d_ij, m_rightparticle->m_primarydiam / 2.0 + m_leftparticle->m_primarydiam / 2.0);
		}
	}

	// Update the tree above this sub-particle.
	if (m_parent != NULL) {
		m_parent->ChangePointer(source, target, d_ij, small_prim);
	}

}
//******************************************hdy************************************************************//

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
					pri->UpdatePAHs(t, t - (*it1)->LastUpdateTime(), model, sys, 1, -1, rng, overflow);
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
        }
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
        if(m_PAHclusterchanged) {
            UpdatePrimary();
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
	if (!(m_pmodel->getTrackPrimarySeparation() && m_pmodel->getTrackPrimaryCoordinates())){ //modified by hdy
        Condition = (m_children_roundingLevel > 0.95);
    } else {
        Condition = (m_distance_centreToCentre == 0.0);
    }

    //if ((Condition && m_leftparticle != NULL) || FakeRounding()) //comment by hdy
	if ((Condition && m_leftparticle != NULL)) {
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
        m_mass            = 0.0;
        m_numcarbon       = 0;
        m_numH            = 0;
        m_numOfEdgeC      = 0;
        m_numOfRings      = 0;
		m_numOfRings5      = 0;
        m_numPAH          = m_PAH.size();
        m_PAHmass         = 0.0;
        m_PAHCollDiameter = 0.0;
    }
    else 
    {
        m_numcarbon       = 0;
        m_numH            = 0;
        m_numOfEdgeC      = 0;
        m_numOfRings      = 0;
		m_numOfRings5     = 0;
        m_numPAH          = m_PAH.size();
        m_PAHmass         = 0.0;
        m_PAHCollDiameter = 0.0;

        //! Initialisation of variables to adjust the primary diameter if the
        //! distance between the centres of primary particles is tracked.
        double d_ij = 0.0;               //!< Distance between the centres of primary particles i and j.
        double r_i = 0.0;                //!< Radius of primary particle i.
        double r_j = 0.0;                //!< Radius of primary particle j.
        double x_i = 0.0;                //!< The distance from the centre of primary particle i to the neck level.
        double A_n = 0.0;                //!< Cross-sectional neck area.
        double A_i = 0.0;                //!< Free surface area of primary particle i.
        double m_primary_diam_old = 0.0; //!< Old primary diameter.
        double m_vol_old = 0.0;          //!< Old volume.

        //! Initialisation of variables but this is only relevant to particles
        //! with more than one primary.
        if (m_pmodel->getTrackPrimarySeparation() && Numprimary() > 1) {
            d_ij = m_parent->m_distance_centreToCentre;
            r_i = m_primarydiam / 2.0;

            if (m_parent->m_leftparticle == this) {
                r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
            } else {
                r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
            }

            x_i = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / 2.0 / d_ij; //!< Eq. (3b) of Langmuir 27:6358 (2011).
            A_n = M_PI * (pow(r_i, 2.0) - pow(x_i, 2.0));                        //!< Eq. (4).
            A_i = 2.0 * M_PI * (pow(r_i, 2.0) + r_i * x_i);                      //!< Eq. (6).
            m_primary_diam_old = m_primarydiam;
            m_vol_old = m_vol;
        }

        int maxcarbon=0;

        for (vector<boost::shared_ptr<PAH> >::iterator i=m_PAH.begin(); i!=m_PAH.end(); ++i) {
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
        if(m_pmodel->ComponentCount() != 1) {
            throw std::runtime_error("Model contains more then one component. Only soot is supported. (PAHPrimary::UpdatePrimary)");
        }

        m_vol  = m_PAHmass / m_pmodel->Components(0)->Density(); //!< Units of m^3.
        m_mass = m_PAHmass;
        m_diam = pow(6.0 * m_vol / PI, ONE_THIRD);
        m_dmob = m_diam;
        //m_dcol = max(m_diam, m_PAHCollDiameter); comment by hdy to test, will uncomment later
		m_dcol = m_diam; //add by hdy to test, will delete later
        m_surf = PI * m_diam * m_diam;
        
        //! If the distance between the centres of primary particles is
        //! tracked, the rate of change in the primary diameter is determined
        //! by its neighbour. Therefore, the particle should be made up of more
        //! than one primary.
		if (!(m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()) || (m_pmodel->getTrackPrimarySeparation() && m_numprimary == 1) || m_parent == NULL) {//hdy
            m_primarydiam = m_diam;
        } //else { //comment by hdy
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
        //} //comment by hdy

		//***********************************hdy************************************************//
		if (m_pmodel->getTrackPrimaryCoordinates()) {
			setRadius(m_primarydiam / 2.0);
		}
		//***********************************hdy************************************************//

        m_avg_coalesc = 0.0;
    }
}

//***********************************hdy************************************************//
/*!
* @brief       Updates the surface area and sintering level of all parents
*
* @param[in]   dS      Surface area increment to adjust area by
*/
void PAHPrimary::UpdateParents(double dS) {
	if (m_parent != NULL) {
		m_parent->m_children_surf += dS;
		m_parent->m_children_sintering = m_parent->SinteringLevel();
		m_parent->UpdateCache();
	}
}
//***************************************hdy**********************************************//

//***************************************hdy**********************************************//
/*!
* @brief       Updates the BinTreePrimary cache from the root node
*
* Works up the tree to the root node and then calls UpdateCache
*
*/
void PAHPrimary::UpdateCacheRoot(void){
	if (m_parent != NULL){
		m_parent->UpdateCacheRoot();
	}
	else{
		UpdateCache(this);
	}
}
//********************************************hdy*****************************************//

//********************************************hdy*****************************************//
unsigned int PAHPrimary::Adjust(const fvector &dcomp,
	const fvector &dvalues, rng_type &rng, unsigned int n)

{
	if (m_leftchild == NULL && m_rightchild == NULL) {

		double dV(0.0);
		double volOld = m_vol;
		double m_diam_old = m_diam;
		double m_primary_diam_old = m_primarydiam;

		// Call to Primary to adjust the state space
		n = Primary::Adjust(dcomp, dvalues, rng, n);

		// Stop doing the adjustment if n is 0.
		if (n > 0) {
			// Update only the primary
			UpdatePrimary();

			//! If the distance between the centres of primary particles or the
			//! primary coordinates is tracked, the rate of change in the
			//! primary diameter is affected by its neighbours.
			if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()) {

				//! Particle with more than one primary.
				if (m_parent != NULL) {

					while (volOld <= m_vol){

						//! Initialisation of variables to adjust the primary diameter if the
						//! distance between the centres of primary particles is tracked.
						double r_i = m_primarydiam / 2.0;	//!< Radius of primary particle i.
						double dr_max = 0.1*r_i;			//!< Maximum change in primary radius during internal step (10% of primary radius)
						double sumterm = 0.0;				//!< Contribution from neighbours to the change in radius 
						double delta_r_i = 0.0;				//!< Change in radius of i
						double free_surface_term = 0.0;		//!< Contribution from neighbours to free surface area

						//! Get contribution from neighbours working up the
						//! binary tree.
						SumNeighbours(this, sumterm);

						//! Calculate change in volume.
						dV = dr_max * (4 * M_PI*r_i*r_i + M_PI*r_i*sumterm);

						//! Calculate change in radius.
						if (volOld + dV > m_vol){
							delta_r_i = (m_vol - volOld)*dr_max / dV;
						}
						else{
							delta_r_i = dr_max;
						}

						//! Store the memory address of this primary so that it is only updated once. 
						std::set<void*> duplicates;
						duplicates.insert(this);

						//! Update the particle separations and calculate the
						//! new free surface of the particle.
						UpdateConnectivity(this, duplicates, delta_r_i, free_surface_term);

						//! Update primary diameter.
						m_primarydiam = 2.0* (r_i + delta_r_i);

						m_primary_diam_old = m_primarydiam;

						if (m_pmodel->getTrackPrimaryCoordinates()) {
							setRadius(m_primarydiam / 2.0);
						}

						//! Update the free surface area if the calculated area
						//! is negative (too many overlaps) then set
						//! m_free_surf = 0.0
						m_free_surf = max(M_PI*m_primarydiam*m_primarydiam - 2 * M_PI*free_surface_term, 0.0);

						volOld += dV;
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

			//! If the distance between the centres of primary particles or the
			//! primary coordinates is not tracked, just update the surface
			//! area and work up the tree.
			else {
				dV = m_vol - volOld;
				double dS(0.0);

				if (dV > 0.0) {
					//! Surface change due to volume addition.
					dS = dV * 2.0 * m_pmodel->GetBinTreeCoalThresh() / m_diam;
				}
				//! TODO: Implement surface area reduction?

				//! Climb back-up the tree and update the surface area and
				//! sintering of a particle.
				UpdateParents(dS);
			}
		}
	}

	// Else this a non-leaf node (not a primary)
	else
	{
		return SelectRandomSubparticle(rng)->Adjust(dcomp, dvalues, rng, n);

		//Note (csl37): this only picks the left or right particle of the root node 
		//and ignores the rest of the tree
		/*
		// Generate random numbers
		boost::bernoulli_distribution<> leftRightChooser;
		// Select particle
		if(leftRightChooser(rng))
		return m_leftparticle->Adjust(dcomp, dvalues, rng, n);
		else
		return m_rightparticle->Adjust(dcomp, dvalues, rng, n);
		*/

		///////////////////////////////////////////////////////////
		// csl37: new primary selection based on free surface area
		// work down tree selecting left/right child based on sum of primary free surface areas under node
		// generate random number, use bernoulli with p=free_surf(leftchild)/free_surf(this)
		/*
		boost::bernoulli_distribution<> leftRightChooser(m_leftchild->m_free_surf/m_free_surf);
		if(leftRightChooser(rng)){
		return m_leftchild->Adjust(dcomp, dvalues, rng, n);
		}else{
		return m_rightchild->Adjust(dcomp, dvalues, rng, n);
		}
		*/
		///////////////////////////////////////////////////////////
	}

	//csl37:
	/*
	// Update property cache.
	UpdateCache(this);
	*/
	//csl37: update the cache from the root node
	UpdateCacheRoot();

	return n;

}
//********************************************hdy*****************************************//

void PAHPrimary::Reset()
{
    m_numcarbon=0;
    m_numH=0;
    m_numOfEdgeC=0;
    m_numOfRings=0;
    m_primarydiam=0.0;
    m_surf=0;
    m_vol=0;
    m_PAH.clear();
    m_avg_coalesc=0;
}

//********************************************hdy*****************************************//
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
		if (MergeCondition()) {
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
//*******************************hdy*****************************//

/*!
 * @param[in] root The root node of this particle
*/
void PAHPrimary::UpdateCache(PAHPrimary *root)
{
    //Update the children
	if (m_leftchild!=NULL)
	{
		m_leftchild->UpdateCache(root);
		m_rightchild->UpdateCache(root);
		//m_numprimary=m_leftchild->m_numprimary+m_rightchild->m_numprimary; comment by hdy, move to the following block
	}
    //this is a primary and the number of primaries below this node is one (this node and no children)
	else
	{
            m_numprimary=1;

            // Check that the primary is begin kept up to date
//            const double oldNumCarbons = m_numcarbon;
//            UpdatePrimary();
//            if(m_numcarbon != oldNumCarbons)
//                std::cerr << "UpdatePrimary has changed num carbons inside UpdateCache\n";

            m_avg_coalesc=0;
			UpdatePrimary(); //add by hdy, learned from swp_bintree_primary.cpp
	}

    //this is not a primary, sum up the properties
	if (m_leftchild!=NULL)
    {
        // remove the PAHs from this node to free memory
        Reset();
		m_numprimary = m_leftchild->m_numprimary + m_rightchild->m_numprimary; //modified by hdy
		m_surf = m_leftchild->m_surf+m_rightchild->m_surf;
		m_primarydiam = (m_leftchild->m_primarydiam+m_rightchild->m_primarydiam);
        m_vol=m_leftchild->m_vol+m_rightchild->m_vol;
        m_numPAH = m_leftchild->m_numPAH+m_rightchild->m_numPAH;
        m_primarydiam = (m_leftchild->m_primarydiam+m_rightchild->m_primarydiam);
        m_mass=(m_leftchild->m_mass+m_rightchild->m_mass);
        m_PAHCollDiameter=max(m_leftchild->m_PAHCollDiameter,m_rightchild->m_PAHCollDiameter);

        m_numcarbon=m_leftchild->m_numcarbon + m_rightchild->m_numcarbon;
        m_numH=m_leftchild->m_numH + m_rightchild->m_numH;

        m_numOfEdgeC=m_leftchild->m_numOfEdgeC + m_rightchild->m_numOfEdgeC;
        m_numOfRings=m_leftchild->m_numOfRings + m_rightchild->m_numOfRings;
		m_numOfRings5 = m_leftchild->m_numOfRings5 + m_rightchild->m_numOfRings5;
        // calculate the coalescence level of the two primaries connected by this node
        m_children_roundingLevel=CoalescenceLevel();
        //sum up the avg coal level
        m_avg_coalesc=m_children_roundingLevel+m_leftchild->m_avg_coalesc+m_rightchild->m_avg_coalesc;

        // calculate the different diameters only for the root node because this goes into the
        // particle tree and gets used by the coagulation kernel
        if (this==root)
        {
             //spherical eqiv radius
            double spherical_radius=pow(3*m_vol/(4*PI),ONE_THIRD);
            m_diam=2*spherical_radius;
            // there are m_numprimary-1 connections between the primary particles
            m_avg_coalesc=m_avg_coalesc/(m_numprimary-1);
            //approxmiate the surface of the particle
            const double numprim_1_3=pow(m_numprimary,-0.333333);
            m_surf=4*PI*spherical_radius*spherical_radius/
                (m_avg_coalesc*(1-numprim_1_3)+numprim_1_3);

            //calculate the surface equivalent radius
           // const double radius_surf=sqrt(m_surf/(4*PI));
            // the average between the surface and voluem equiv diameter
            //const double meandiam=spherical_radius+radius_surf;            //2*0.5*(radius(vol)+radius(sphere))
            const double aggcolldiam=(6*m_vol/m_surf)*
                                     pow(pow(m_surf,3)/(36*PI*m_vol*m_vol),(1.0/1.8));
            // the maximum of the largest PAH diameter and
            // the average between the surface and voluem equiv diameter
            //const double cdiam=max(aggcolldiam,m_PAHCollDiameter); comment by hdy to test; will uncomment later
			const double cdiam = aggcolldiam; //modified by hdy to test; will delete later
            m_dmob = aggcolldiam;
            SetCollDiameter(cdiam);
        }
        else
        {
            m_diam=0;
            m_dmob=0;

        }
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
      out<<"\" "<<this<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_roundingLevel<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this->m_leftchild<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_roundingLevel<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this->m_rightchild<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_roundingLevel<<"|"<<this<<"\"];"<<endl;
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
      out<<"\" "<<this<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_roundingLevel<<"|"<<this<<"\"];"<<endl;
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

//! Return distance between the centres of primary particles.
double PAHPrimary::Distance() const
{
    return m_distance_centreToCentre;
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

        val_int =  m_numPAH;
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

        val_int =  m_leftparticle_numPAH;
        out.write((char*)&val_int, sizeof(val_int));

        val = m_children_surf;
        out.write((char*)&val, sizeof(val));

        val = m_children_roundingLevel;
        out.write((char*)&val, sizeof(val));

        val = m_distance_centreToCentre;
        out.write((char*)&val, sizeof(val));

		//*************hdy***************//
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
		//*************hdy***************//

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

        val_int = (int) m_PAH.size();
        out.write((char*)&val_int, sizeof(val_int));

        // write the PAH stack (m_PAH)
		PahSerialisationMap *pahDuplicates = reinterpret_cast<PahSerialisationMap*>(duplicates);
        outputPAHs(out, *pahDuplicates);

        // Output base class.
        Primary::Serialize(out);

    } else {
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
        m_children_roundingLevel = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_distance_centreToCentre = val;

		//********************hdy*******************//
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
		//********************hdy*******************//

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

    } else {
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
void PAHPrimary::Sinter(double dt, Cell &sys, const Processes::SinteringModel &model, rng_type &rng, double wt)
{
    //! Only update the time on the root node.
    if (m_parent == NULL) {
        m_sint_time += dt;
        SetSinteringTime(m_sint_time);
    }

    //! Do only if there is a particle to sinter.
    if (m_leftparticle != NULL)
    {
        double t1 = 0.0, delt = 0.0, tstop = dt; //! Declare time step variables.
        bool Condition;                          //! Declare variable for condition for complete sintering.

        //! The scale parameter discretises the delta-S when using
        //! the Poisson distribution.  This allows a smoother change
        //! (smaller scale = higher precision).
        double scale = 0.01;

        //! The sintering model depends on whether the distance between the
        //! centres of primary particles is tracked. If tracked, sintering
        //! results in a decrease in the distance between the primaries and an
        //! increase in their diameters. If not, sintering results in a
        //! decrease in the common surface between the primaries.
        if (!m_pmodel->getTrackPrimarySeparation()) {
            //! Store the old surface area of particles.
            // double surf_old = m_children_surf;

            //! Calculate the spherical surface.
            const double spherical_surface = 4 * PI * m_children_radius * m_children_radius;

            //! Declare sintering rate.
            double r = 0.0;

            //! Define the maximum allowed change in surface area in one
            //! internal time step (10% of spherical surface).
            double dAmax = 0.1 * spherical_surface;

            //! Perform integration loop.
            while (t1 < tstop) {
                //! Calculate sintering rate.
                r = model.Rate(m_time+t1, sys, *this);

                if (r > 0) {
                    //! Calculate next time-step end point so that the
                    //! surface area changes by no more than dAmax.
                    delt = dAmax / max(r, 1.0e-300);

                    //! Approximate sintering by a poisson process.  Calculate
                    //! number of poisson events.
                    double mean;

                    if (tstop > (t1+delt)) {
                        //! A sub-step, we have changed surface by dAmax, on average.
                        mean = 1.0 / scale;
                    } else {
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
        } else {
            //! Declare characteristic sintering time.
            double tau = 0.0;

            //! Define the maximum allowed change (1%) in the distance between the
            //! centres of primary particles in one internal time step. In the case
            //! of pure sintering it was found that if the allowed change is too
            //! large (~10%) a significant error is incurred in the final spherical
            //! volume determined through comparisons with the mass-derived volume.
            //! Note that the smaller the distance is, the smaller the changes are.
            double dd_ij_Max = m_distance_centreToCentre / 100.0;

            while (t1 < tstop) {
                double r_i = this->m_leftparticle->m_primarydiam / 2.0;
                double r_j = this->m_rightparticle->m_primarydiam / 2.0;
                double d_ij = m_distance_centreToCentre;

                //! In Section 3.1.2 of Langmuir 27:6358 (2011), it is argued that
                //! the smaller particle dominates the sintering process.
                if (r_i <= r_j) {
                    tau = model.SintTime(sys, *this->m_leftparticle); //!< The left particle is smaller than the right. 
                } else {
                    tau = model.SintTime(sys, *this->m_rightparticle); //!< The right particle is smaller than the left.
                }

                //! Gamma is the surface tension and eta is the viscosity, and the
                //! ratio (gamma/eta) can be related to tau.
                //! J. Colloid Interface Sci. 140:419 (1990).
                double gamma_eta = min(r_i, r_j) / tau;

                //! Definition of variables for conciseness.
                double d_ij2 = pow(d_ij, 2.0); 
                double r_i2 = pow(r_i, 2.0);
                double r_j2 = pow(r_j, 2.0);
                double r_i4 = pow(r_i, 4.0);
                double r_j4 = pow(r_j, 4.0);

                ///////////////////////////////////////////////////////
                /// References to equations in Langmuir 27:6358 (2011).
                ///////////////////////////////////////////////////////

                //! Eq. (14a).
                double dd_ij_dt = 4.0 * r_i * r_j * d_ij2 * (r_i + r_j) * gamma_eta /
                                ((r_i + r_j + d_ij) * (r_i4  + r_j4 - 2.0 * r_i2 * r_j2 + 4.0 * d_ij * r_i * r_j *(r_i + r_j) - d_ij2 * (r_i2 + r_j2)));

                double x_i = (d_ij2 - r_j2 + r_i2) / (2.0 * d_ij); //!< Eq. (3b).
                double x_j = (d_ij2 - r_i2 + r_j2) / (2.0 * d_ij); //!< Eq. (3b).
                double A_n = M_PI * (r_i2 - pow(x_i, 2.0));        //!< Eq. (4).
                double A_i = 2.0 * M_PI * (r_i2 + r_i * x_i);      //!< Eq. (6).
                double A_j = 2.0 * M_PI * (r_j2 + r_j * x_j);      //!< Eq. (6).

                //! The expression for B_i in Eq. (8) is wrong. By combining
                //! Eqs. (5) and (7), we can obtain two equations which are
                //! functions of r_i and r_j. Subsequently combined these two
                //! equations and used Wolfram Alpha to rearrange equation in terms
                //! of r_i (and r_j).
                //!
                //! @todo Remove derivation and replace with reference to preprint
                //!       or paper if results do get published.
                double B_i = (pow(A_n, 2.0) * r_j * (-2.0 * d_ij + x_j + x_i) / (A_i * d_ij + A_n * r_i) / (d_ij * A_j + A_n * r_j) - A_n * d_ij * A_j * (d_ij - x_i) / (A_i * d_ij + A_n * r_i) / (d_ij * A_j + A_n * r_j)) / 
                             (1.0 - pow(A_n, 2.0) * r_i * r_j / (A_i * d_ij + A_n * r_i) / (d_ij * A_j + A_n * r_j));

                double B_j = (pow(A_n, 2.0) * r_i * (-2.0 * d_ij + x_i + x_j) / (A_j * d_ij + A_n * r_j) / (d_ij * A_i + A_n * r_i) - A_n * d_ij * A_i * (d_ij - x_j) / (A_j * d_ij + A_n * r_j) / (d_ij * A_i + A_n * r_i)) / 
                             (1.0 - pow(A_n, 2.0) * r_j * r_i / (A_j * d_ij + A_n * r_j) / (d_ij * A_i + A_n * r_i));

                double V_i = 2.0 / 3.0 * M_PI * pow(r_i, 3.0) + M_PI * r_i2 * x_i - 1.0 / 3.0 * M_PI * pow(x_i, 3.0); //!< Eq. (3a).
                double V_j = 2.0 / 3.0 * M_PI * pow(r_j, 3.0) + M_PI * r_j2 * x_j - 1.0 / 3.0 * M_PI * pow(x_j, 3.0); //!< Eq. (3a).

                delt = dd_ij_Max / max(dd_ij_dt, 1.0e-300);
                double mean;

                if (tstop > (t1 + delt)) {
                    mean = 1.0 / scale;
                } else {
                    mean = dd_ij_dt * (tstop - t1) / (scale * dd_ij_Max);
                }

                //! Only perform sintering if the distance between the centres
                //! of primary particles i and j is positive.
                if (d_ij > 0.0) {
                    boost::random::poisson_distribution<unsigned, double> repeatDistribution(mean);
                    const unsigned n = repeatDistribution(rng);
                    m_distance_centreToCentre -= (double)n * scale * dd_ij_Max; //!< Sintering decreases d_ij hence the negative sign.

                    //! If m_distance_centreToCentre (or d_ij) is some very small
                    //! value, it was found that B_i and B_j will be 1.#INF (or
                    //! infinite) and will result in an error when adjusting the
                    //! particle diameters below.
                    if (m_distance_centreToCentre < 1.0e-12) {
                        m_distance_centreToCentre = 0.0;
                    }

                    //! The factor of 2 is because Eq. (8) is the rate of change
                    //! in radius.
                    this->m_leftparticle->m_primarydiam -= (double)n * scale * 2.0 * B_i * dd_ij_Max;  //!< Eq. (8).
                    this->m_rightparticle->m_primarydiam -= (double)n * scale * 2.0 * B_j * dd_ij_Max; //!< Eq. (8).

                    t1 += delt;
                } else {
                    break; //! d_ij is 0 so do not continue to sinter.
                }
            }
        }

        m_children_roundingLevel = RoundingLevel();

        //! One can specify a member in PAHPrimary class to store the rate of sintering, but now, it is not useful.
        //m_sint_rate = r;

        //! The condition for whether a particle has sintered depends on
        //! whether the distance between the centres of primary particles is
        //! tracked. If tracked, a particle has sintered if the distance is 0.
        //! If not, the condition depends on whether the rounding level exceeds
        //! an arbitrarily high threshold.
        if (!m_pmodel->getTrackPrimarySeparation()) {
            Condition = (m_children_roundingLevel > 0.95);
        } else {
            Condition = (m_distance_centreToCentre == 0.0);
        }

        //! Check whether particle has sintered. If true, then merge in
        //! CheckRounding function.
        if (Condition) {
            CheckRounding();
            UpdateCache();

            if (m_leftchild != NULL && m_rightchild != NULL) {
                m_leftchild->Sinter(dt, sys, model, rng, wt);
                m_rightchild->Sinter(dt, sys, model, rng, wt);
            }
        } else {
            m_leftchild->Sinter(dt, sys, model, rng, wt);
            m_rightchild->Sinter(dt, sys, model, rng, wt);
        }

        UpdateCache();

        ////! Adjust the gas-phase concentration.
        //fvector dc(sys.GasPhase().Species()->size(), 0.0);

        //double n_NAvol_sint = wt * (double)num_H2O / (NA * sys.SampleVolume());
        //dc[Sprog::Species::Find(string("H2O"),*sys.GasPhase().Species())] += n_NAvol_sint;
        //sys.AdjustConcs(dc);

        m_children_roundingLevel = RoundingLevel();
    } //!< endif m_leftparticle != NULL.
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
