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
#include <math.h>         //!< Then include so that arc cosine or inverse
                          //!< operation of cosine (acos) can be used.
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

//! Generate uniformly distributed random number.
//! To randomly rotate particle and determine point of contact between two
//! particles which collide in the Coagulate function.
#include <boost/random/uniform_01.hpp>

#include "string_functions.h"

int uniquePAHCounter = 0;

//! Stores memory address of pointer to a PAH to trace.
//! For the purpose of creating a video showing the evolution of a particle.
void *PAHTracer = NULL;

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace Sweep::KMC_ARS;

using namespace std;
using namespace Strings;

unsigned int PAHPrimary::Counter = 0;

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
    m_rightparticle(NULL)
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
    m_rightparticle(NULL)
{
    // Other parts of the code check for a non-zero composition
    m_comp[0]=1;

    m_cen_bsph[0] = 0.0;
    m_cen_bsph[1] = 0.0;
    m_cen_bsph[2] = 0.0;

    m_cen_mass[0] = 0.0;
    m_cen_mass[1] = 0.0;
    m_cen_mass[2] = 0.0;

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
    m_rightparticle(NULL)
{
    // Other parts of the code check for a non-zero composition
    m_comp[0]=1;
   
    m_cen_bsph[0] = 0.0;
    m_cen_bsph[1] = 0.0;
    m_cen_bsph[2] = 0.0;

    m_cen_mass[0] = 0.0;
    m_cen_mass[1] = 0.0;
    m_cen_mass[2] = 0.0;

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
    m_rightparticle(NULL)
{
    m_comp[0]=1;

    m_cen_bsph[0] = 0.0;
    m_cen_bsph[1] = 0.0;
    m_cen_bsph[2] = 0.0;

    m_cen_mass[0] = 0.0;
    m_cen_mass[1] = 0.0;
    m_cen_mass[2] = 0.0;
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

    //! Flag to determine whether the pointer PAHTracer has been initialised.
    static bool initialised;

    //! Store the memory address of the first PAH which is added.
    if (!initialised) {
        PAHTracer = new_PAH.get();
        initialised = true;
    }
    
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

    m_PAHCollDiameter         = source->m_PAHCollDiameter;
    m_time                    = source->m_time;
    m_PAHmass                 = source->m_PAHmass;
    m_leftchild               = source->m_leftchild;
    m_rightchild              = source->m_rightchild;
    m_leftparticle            = source->m_leftparticle;
    m_rightparticle           = source->m_rightparticle;
    m_parent                  = source->m_parent;
    m_numPAH                  = source->m_numPAH;
    m_numprimary              = source->m_numprimary;
    m_primarydiam             = source->m_primarydiam;
    m_sqrtLW                  = source->m_sqrtLW;
    m_LdivW                   = source->m_LdivW;
    m_pmodel                  = source->m_pmodel;
    m_surf                    = source->m_surf;
    m_vol                     = source->m_vol;
    m_children_surf           = source->m_children_surf;
    m_children_vol            = source->m_children_vol;
    m_children_radius         = source->m_children_radius;
    m_children_roundingLevel  = source->m_children_roundingLevel;
    m_rightparticle_numPAH    = source->m_rightparticle_numPAH;
    m_leftparticle_numPAH     = source->m_leftparticle_numPAH;
    m_leftparticle_vol_old    = source->m_leftparticle_vol_old;
    m_rightparticle_vol_old   = source->m_rightparticle_vol_old;
    m_fdim                    = source->m_fdim;
    m_Rg                      = source->m_Rg;
    m_avg_coalesc             = source->m_avg_coalesc;
    m_numcarbon               = source->m_numcarbon;
    m_numH                    = source->m_numH;
    m_numOfEdgeC              = source->m_numOfEdgeC;
    m_numOfRings              = source->m_numOfRings;
    m_values                  = source->m_values;
    m_comp                    = source->m_comp;
    m_sint_time               = source->m_sint_time;
    m_cen_bsph                = source->m_cen_bsph;
    m_cen_mass                = source->m_cen_mass;
    m_distance_centreToCentre = source->m_distance_centreToCentre;
    m_r                       = source->m_r;
    m_r2                      = source->m_r2;
    m_r3                      = source->m_r3;

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
    else m_PAH.assign(source->m_PAH.begin(),source->m_PAH.end());
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
	if (target<=m_leftchild->m_numprimary)
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
                //! rhsparticle is a gas-phase PAH but the this pointer may be
                //! pointing to a gas-phase PAH in which case this would be an
                //! inception event, or a single primary particle in which case
                //! it would be a condensation event.
                m_PAH.insert(m_PAH.end(),rhsparticle->m_PAH.begin(),rhsparticle->m_PAH.end());
                UpdatePrimary();
            }
        }
		UpdateCache();
        //Check the coalescence ratio
        CheckRounding();
	} else {       
        //! Coagulation process.
        PAHPrimary *newleft = new PAHPrimary;
		PAHPrimary *newright = new PAHPrimary;
        PAHPrimary copy_rhs(*rhsparticle);
		//rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);

        //bool print=false;
        //if (this->m_numprimary > 2 && notconstpah->m_numprimary > 2)
        //{
        //    PrintTree("before1");
        //    notconstpah->PrintTree("before2");
        //    print = true;
        //    cout << "printing tree" << endl;
        //}

        //! Select where to add the second particle.
        boost::bernoulli_distribution<> bernoulliDistrib;
        boost::variate_generator<rng_type &, boost::bernoulli_distribution<> > leftRightChooser(rng, bernoulliDistrib);
		if (leftRightChooser()) {
			newleft->CopyParts(this);
			newright->CopyParts(&copy_rhs);
		} else {
			newright->CopyParts(this);
			newleft->CopyParts(&copy_rhs);
		}

        //! Set the pointers.
		m_leftchild = newleft;
		m_rightchild = newright;
		newright->m_parent = this;
		newleft->m_parent =this;

        //! Set the pointers to the parent node.
		if (newleft->m_leftchild != NULL) {
			newleft->m_leftchild->m_parent = newleft;
			newleft->m_rightchild->m_parent = newleft;
		}

		if (newright->m_leftchild != NULL) {
			newright->m_leftchild->m_parent = newright;
			newright->m_rightchild->m_parent = newright;
		}

        m_children_roundingLevel = 0;
		UpdateCache();

        std::ofstream outfile;

        //! It is assumed that primary pi from particle Pq and primary pj from
        //! particle Pq are in point contact and by default pi and pj are
        //! uniformly selected. If we track the coordinates of the primaries in
        //! a particle, we can do a smarter selection where pi and pj are
        //! determined by ballistic cluster-cluster aggregation.
        if (m_pmodel->getTrackPrimaryCoordinates()) {
            boost::uniform_01<rng_type&, double> uniformGenerator(rng);

            //! Sphere point picking.
            //! It is incorrect to select spherical coordinates theta (polar
            //! angle) and phi (azimuthal angle) from uniform distributions
            //! theta E [0, 2 * pi) and phi E [0, pi] as was done previously as
            //! points picked in this way will be 'bunched' near the poles.
            //! http://mathworld.wolfram.com/SpherePointPicking.html
            //double phi1 = acos(2 * uniformGenerator() - 1);
            //double theta1 = 2.0 * PI * uniformGenerator();

            double phi1   = uniformGenerator() * 2.0 * PI;
            double theta1 = ((2.0*uniformGenerator())-1.0) * PI;

            //! Write POV-Ray commands only if (1) the left child contains the
            //! PAH to be traced and (2) it is made up of two or more primaries
            //! because rotation of a single primary is not interesting.
			if (false && m_leftchild->m_numprimary > 1) {               
                //! Increment counter to signal the next sequence of events.
				incrementCounter();				

                outfile.open("sphere.pov", std::ios_base::app);

                //! range: the range of values between which the global clock
                //!        variable lies.
                //! clock+counter: will always fall within (0, 1].
                //! theta, phi: declared for use in the POV-Ray transformation
                //!             matrix written in rotateCOM.
                outfile << "\t#range(" << getCounter() << "," << getCounter() + 1 << ")\n"
						<< "\t\t#declare clock" << getCounter() << " = clock - " << getCounter() << ";\n";
						//<< "\t\t#declare theta = " << theta1 << ";\n"
						//<< "\t\t#declare phi = " << phi1 << ";\n";
				outfile.close();

                //! Rotate centre-of-mass and write the correponding POV-Ray
                //! transformation matrix.
                m_leftchild->rotateCOM(theta1, phi1, false);

                //! Terminate range directive.
				outfile.open("sphere.pov", std::ios_base::app);
				outfile << "\t\t#break\n";
				outfile.close();
			} else {
                //! Rotate centre-of-mass but do not write POV-Ray
                //! transformation matrix.
                m_leftchild->rotateCOM(theta1, phi1, false);
            }



            //! Sphere point picking.
            //! It is incorrect to select spherical coordinates theta (polar
            //! angle) and phi (azimuthal angle) from uniform distributions
            //! theta E [0, 2 * pi) and phi E [0, pi] as was done previously as
            //! points picked in this way will be 'bunched' near the poles.
            //! http://mathworld.wolfram.com/SpherePointPicking.html
            //double phi2 = acos(2 * uniformGenerator() - 1);
            //double theta2 = 2.0 * PI * uniformGenerator();

            double phi2   = uniformGenerator() * 2.0 * PI;
            double theta2 = ((2.0*uniformGenerator())-1.0) * PI;

            //! Write POV-Ray commands only if (1) the right child contains the
            //! PAH to be traced and (2) it is made up of two or more primaries
            //! because rotation of a single primary is not interesting.
			if (false && m_rightchild->m_numprimary > 1) {
                //! Increment counter to signal the next sequence of events.
				incrementCounter();

				outfile.open("sphere.pov", std::ios_base::app);
				
                //! range: the range of values between which the global clock
                //!        variable lies.
                //! clock+counter: will always fall within (0, 1].
                //! theta, phi: declared for use in the POV-Ray transformation
                //!             matrix written in rotateCOM.
                outfile << "\t#range(" << getCounter() << "," << getCounter() + 1 << ")\n"
						<< "\t\t#declare clock" << getCounter() << " = clock - " << getCounter() << ";\n"
						<< "\t\t#declare theta = " << theta2 << ";\n"
						<< "\t\t#declare phi = " << phi2 << ";\n";
				outfile.close();

                //! Rotate centre-of-mass and write the correponding POV-Ray
                //! transformation matrix.
                m_rightchild->rotateCOM(theta2, phi2, true);

                //! Terminate range directive.
				outfile.open("sphere.pov", std::ios_base::app);
				outfile << "\t\t#break\n";
				outfile.close();
			} else {
                //! Rotate centre-of-mass but do not write POV-Ray
                //! transformation matrix.
                m_rightchild->rotateCOM(theta2, phi2, false);
                
            }

            //! Do not bother to write POV-Ray translation of nodes as the
            //! movement is likely to be small.
            m_leftchild->centreBoundSph(false);
            m_rightchild->centreBoundSph(false);

            Coords::Vector D;
            double sumr = 0.0;
            bool hit = false;
            while (!hit) {
                //! Primary pi from particle Pq and primary pj from particle Pq
                //! are in point contact.
                sumr = m_leftchild->Radius() + m_rightchild->Radius();
                
                D[0] = ((2.0 * uniformGenerator()) - 1.0) * sumr;
                D[1] = ((2.0 * uniformGenerator()) - 1.0) * sumr;

                //! It is here that the primaries pi and pj are determined
                //! through the ballistic cluster-cluster aggregation.
                hit = this->minCollZ(*m_leftchild, *m_rightchild, D[0], D[1], D[2]);
            }

            //! First check whether the primaries belonging to the node pointed
            //! to by the this pointer contains the PAH to be traced.
            if (PAHTracerCheck()) {
				
                
                //! If Counter is 0 then there is some preliminary information
                //! to be written.
				if (getCounter() == 0) {
                    outfile.open("sphere.pov", std::ios_base::app);
                    outfile << "#include \"colors.inc\"\n\n"
                            
                            << "#declare tem_finish=texture{pigment{color Gray90 transmit 0.6} finish{ambient 0.4 diffuse 0.7 phong 0.1}}\n"
					        << "#declare W=500;\n"
                            << "#declare H=500;\n"
                            << "#declare dLight=10*max(W,H);\n"
                            << "#declare Y=1.5*W*tan(45.0/2.0);\n\n"

                            << "light_source {<0,-dLight,0> White}\n"
                            << "background {color Gray50}\n"
                            << "camera{look_at<0,0,0> location<0,Y,0>}\n\n"

                            << "#switch(clock)\n"
							<< "\t#range(0,1)\n"
							<< "\t\t#declare clock0 = clock;\n";
                    outfile.close();
                    this->m_leftchild->writePrimaryCoordinatesRadius();
                    outfile.open("sphere.pov", std::ios_base::app);
				    outfile << "\t\t#break\n";
				    outfile.close();
                            
                }
    //            else {
    //                //! Increment counter to signal the next sequence of events.
				//	incrementCounter();

				//	outfile << "\t#range(" << getCounter() << "," << getCounter() + 1 << ")\n"
				//			<< "\t\t#declare clock" << getCounter() << " = clock - " << getCounter() << ";\n";
			 //   }
				
			}

			if (false) {
                //! Left child contains PAH to be traced. Since it is already
                //! at the origin just write the coordinates of the primaries
                //! in the left child.
				this->m_leftchild->writePrimaryCoordinatesRadius();

                //! Write POV-Ray translation of the right child.
				this->m_rightchild->Translate(D[0], D[1], D[2], true, true, false);
			} else if (false) {
                //! Right child contains PAH to be traced. Again write the
                //! coordinates of the right child. But this time reverse the
                //! signs on elements of the vector D so that the particle with
                //! the PAH to be traced always remain at the centre.
				this->m_rightchild->writePrimaryCoordinatesRadius();
				this->m_leftchild->Translate(-D[0], -D[1], -D[2], true, true, false);
			} else {
                //! Particle does not contain PAH to be traced, so there is no
                //! need to write POV-Ray translation of the right child.
				//this->m_rightchild->Translate(D[0], D[1], D[2], false, true, false);
			}

            if (this->m_leftchild->PAHTracerCheck()) {
                this->m_rightchild->Translate(D[0], D[1], D[2], false, true, false);
            } else {
                this->m_leftchild->Translate(-D[0], -D[1], -D[2], false, true, false);
            }

            //! If particle contains PAH to be traced terminate range directive.
			if (false) {
				outfile.open("sphere.pov", std::ios_base::app);
				outfile << "\t\t#break\n";
				outfile.close();
			}

            if (this->m_leftchild->PAHTracerCheck()) {
                this->inverseRotateCOM(theta1, phi1, false);
            } else if (this->m_rightchild->PAHTracerCheck()) {
                this->inverseRotateCOM(theta2, phi2, false);
            }

            //! Calculate properties of this particle.
            this->calcBoundSph();
            this->calcCOM();

			if (false) {
                //! Increment counter to signal the next sequence of events.
				incrementCounter();

				outfile.open("sphere.pov", std::ios_base::app);

                //! range: the range of values between which the global clock
                //!        variable lies.
                //! clock+counter: will always fall within (0, 1].				
                outfile << "\t#range(" << getCounter() << "," << getCounter() + 1 << ")\n"
						<< "\t\t#declare clock" << getCounter() << " = clock - " << getCounter() << ";\n";
				outfile.close();
			}

            //! If particle contains PAH to be traced write POV-Ray translation
            //! of particle
            this->centreBoundSph(false);

            //! If particle contains PAH to be traced terminate range directive.
			if (false) {
				outfile.open("sphere.pov", std::ios_base::app);
				outfile << "\t\t#break\n";
				outfile.close();
			}

            //! Particle is centred about its bounding sphere for the purpose
            //! generating its structure. Now that it is complete centre it
            //! about its centre-of-mass.
            centreCOM();

            //outfile.open("sphere.pov", std::ios_base::app);
            //outfile << "#end\n";
            //outfile.close();
        } else {
            //! Randomly select the primaries that are touching.
		    this->m_leftparticle = m_leftchild->SelectRandomSubparticle(rng);
		    this->m_rightparticle = m_rightchild->SelectRandomSubparticle(rng);
        }

        if (PAHTracerCheck()) {
            incrementCounter();
            outfile.open("sphere.pov", std::ios_base::app);			
            outfile << "\t#range(" << getCounter() << "," << getCounter() + 1 << ")\n"
				    << "\t\t#declare clock" << getCounter() << " = clock - " << getCounter() << ";\n";
            outfile.close();
            this->writePrimaryCoordinatesRadius();
            outfile.open("sphere.pov", std::ios_base::app);
			outfile << "\t#break\n";
			outfile.close();
        }

        //! Set the sintertime for the new created primary particle.
        SetSinteringTime(std::max(this->m_sint_time, rhsparticle->m_sint_time));

        //! Initialise the variables used to calculate the coalesence ratio.
        m_children_vol=m_leftparticle->m_vol+m_rightparticle->m_vol;
        m_children_surf=(m_leftparticle->m_surf+m_rightparticle->m_surf);
        m_leftparticle_vol_old=m_leftparticle->m_vol;
        m_rightparticle_vol_old=m_rightparticle->m_vol;
        m_leftparticle_numPAH=m_leftparticle->m_numPAH;
        m_rightparticle_numPAH=m_rightparticle->m_numPAH;
        m_children_radius=pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
        m_children_roundingLevel=CoalescenceLevel();
        m_distance_centreToCentre = m_leftparticle->m_primarydiam / 2.0 + m_rightparticle->m_primarydiam / 2.0;
        CheckRounding();

        //! Must set all the pointer to NULL otherwise the delete function.
        //! will also delete the children.
	    copy_rhs.m_leftchild=NULL;
	    copy_rhs.m_rightchild=NULL;
	    copy_rhs.m_parent=NULL;
        copy_rhs.m_leftparticle=NULL;
        copy_rhs.m_rightparticle=NULL;

        //rhsparticle->Clear();
        //if (fabs(surfacebeforerhs + surfacebefore - m_surf) /m_surf > 1e-6) {
        //    cout << "error" << surfacebeforerhs<<' ' <<surfacebefore << ' '<<m_surf<<endl;
        //    //PrintTree("after");
        //}
        //if (print)
        //    PrintTree("after");
	}
    return *this;
}

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

    vector<PAH>::const_iterator j;
      //make sure this primary has children to merge
	  if(m_leftchild!=NULL)
           {
		if ( m_leftchild==m_leftparticle && m_rightchild==m_rightparticle)
		{
            //this node has only two primaries in its subtree
            //it is possible that this node is not the root node and belongs to a bigger particle
            //copy the PAHs of both children to the parent node

            //m_PAH is empty, therefore no need to append
            m_PAH=m_rightparticle->m_PAH;

            m_PAH.insert(this->m_PAH.end(),m_leftparticle->m_PAH.begin(),m_leftparticle->m_PAH.end());
            //update the pointers that pointed to the two former children
			ChangePointer(m_leftchild,this);
			ChangePointer(m_rightchild,this);

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
		}


		else
		{
			if (m_leftchild->m_numprimary<m_rightchild->m_numprimary)
			{
				//append to left subtree because there are fewer primaries
                //this is only to keep the tree balanced
				PAHPrimary *oldleftparticle=m_leftparticle;

                //copy the PAHs

				//for (j=oldleftparticle->m_PAH.begin(); j!=oldleftparticle->m_PAH.end(); ++j) {
				//	m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(),PAH(*j));
				//}
                m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(),oldleftparticle->m_PAH.begin(),oldleftparticle->m_PAH.end());
                m_rightparticle->UpdatePrimary();
                //set the pointers from the leftprimary to the rightprimary
                //this will be the new bigger primary
				oldleftparticle->ChangePointer(oldleftparticle,m_rightparticle);
                m_rightparticle->ChangePointer(m_rightparticle,m_rightparticle);
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
                m_leftparticle->m_PAH.insert(m_leftparticle->m_PAH.end(),oldrightparticle->m_PAH.begin(),oldrightparticle->m_PAH.end());
                m_leftparticle->UpdatePrimary();
                oldrightparticle->ChangePointer(oldrightparticle,m_leftparticle);
                m_leftparticle->ChangePointer(m_leftparticle,m_leftparticle);
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
        UpdateCache();
  //      PrintTree("after.inp");
        }
}

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
void PAHPrimary::UpdatePAHs(const double t, const Sweep::ParticleModel &model,Cell &sys, rng_type &rng, bool &particleChanged)
{
    // Either the primary has two children or it is a leaf of thef
    // tree
    if (m_leftchild!=NULL)
    {
        // Recurse down to the leaves
        m_leftchild->UpdatePAHs(t, model,sys, rng, particleChanged);
        m_rightchild->UpdatePAHs(t, model,sys, rng, particleChanged);
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
			}

			//! Time for one particular PAH to grow.
			const double growtime = t - (*it)->lastupdated;
			assert(growtime >= 0.0);

			const int oldNumCarbon = (*it)->m_pahstruct->numofC(); 
			const int oldNumH = (*it)->m_pahstruct->numofH();

			/*!
             * Here updatePAH function in KMC_ARS model is called.
             * waitingSteps is set to be 1 by dc516, details seeing KMCSimulator::updatePAH().
			 */
			sys.Particles().Simulator()->updatePAH((*it)->m_pahstruct, (*it)->lastupdated, growtime, 1,
													rng, growthfact, (*it)->PAH_ID);

			(*it)->lastupdated=t;

            //! Invalidate PAH.
            /*!
             * A PAH is made invalid by setting the number of carbons to 5. The
			 * reason for only applying this to particles is that we would not
			 * want to remove single PAHs in the gas-phase which fall below the
			 * minimum number of rings for inception but are still growing.
             */
            if((*it)->m_pahstruct->numofRings() < thresholdOxidation && m_numPAH>=minPAH){
                (*it)->m_pahstruct->setnumofC(5);
			}		

			//! See if anything changed, as this will required a call to UpdatePrimary() below.
			if(oldNumCarbon != (*it)->m_pahstruct->numofC() || oldNumH != (*it)->m_pahstruct->numofH())
			{
				m_PAHclusterchanged = true; 
				m_PAHchanged = true;
                particleChanged = true;
			}

			/*!
             * The second condition ensures that the m_InvalidPAH is modified in the correct way
			 * consider 2 PAH, the first is invalid, the other is valid. If not using the
             * second condition, the m_InvalidPAH will be false enventually, but it should be true.
			 */ 
			if (m_PAHchanged && !m_InvalidPAH)
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
    return (it->m_pahstruct->numofC() < m_control||(it->m_pahstruct->numofC() <= m_control && NumPAH()!=1));
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

//! Returns true if this node is a leaf (has no children).
bool PAHPrimary::isLeaf(void) const
{
    return (m_leftchild == NULL) && (m_rightchild == NULL);
}

/*!
 *  @brief Calculates the minimum collision distance.
 *
 *  Calculates the minimum collision distance between a target and a bullet
 *  node by moving down the binary tree. If the nodes collide then returns
 *  true, otherwise returns false.
 *
 *  @param[in]  target Target node.
 *  @param[in]  bullet Bullet node.
 *  @param[in]  dx     Distance to translate in the x-axis.
 *  @param[in]  dy     Distance to translate in the y-axis.
 *  @param[out] dz     Distance to translate in the z-axis.
 *  
 *  @return            Have the nodes collided?
 */
bool PAHPrimary::minCollZ(PAHPrimary &target, PAHPrimary &bullet, double dx, double dy, double &dz)
{
    bool hit = false, hit1 = false;
    double dz2 = 0.0, dz3 = 0.0, dz4 = 0.0;

    if (target.isLeaf()) {
        //! Target is a leaf.
        if (bullet.isLeaf()) {
            //! Bullet is a leaf (both leaves).

            this->m_leftparticle = &target;
            this->m_rightparticle = &bullet;

            return calcCollZ(target.boundSphCentre(), target.Radius(),
                             bullet.boundSphCentre(), bullet.Radius(),
                             dx, dy, dz, 0.0, false);
        } else {
            //! Bullet is not a leaf, call sub-nodes.
            //! Calculate minimum dz for the target and the bullet left subnode.
            hit1 = calcCollZ(target.boundSphCentre(), target.Radius(),
                             bullet.m_leftchild->boundSphCentre(), bullet.m_leftchild->Radius(),
                             dx, dy, dz, 0.0, false);

            if (hit1) hit = minCollZ(target, *bullet.m_leftchild, dx, dy, dz);
            
            //! Calculate minimum dz for the target and the bullet right subnode.
            hit1 = calcCollZ(target.boundSphCentre(), target.Radius(),
                             bullet.m_rightchild->boundSphCentre(), bullet.m_rightchild->Radius(),
                             dx, dy, dz2, 0.0, false);
            
            if (hit1) hit = minCollZ(target, *bullet.m_rightchild, dx, dy, dz2) || hit;
            
            //! Return minimum dz.
            dz = min(dz, dz2);

            return hit;
        }
    } else {
        //! Target is not a leaf.
        if (bullet.isLeaf()) {
            //! Bullet is a leaf, call target sub-nodes.
            //! Calculate minimum dz for the target left subnode and the bullet.
            hit1 = calcCollZ(target.m_leftchild->boundSphCentre(), target.m_leftchild->Radius(),
                             bullet.boundSphCentre(), bullet.Radius(),
                             dx, dy, dz, 0.0, false);

            if (hit1) hit = minCollZ(*target.m_leftchild, bullet, dx, dy, dz);

            //! Calculate minimum dz for the target right subnode and the bullet.
            hit1 = calcCollZ(target.m_rightchild->boundSphCentre(), target.m_rightchild->Radius(),
                             bullet.boundSphCentre(), bullet.Radius(),
                             dx, dy, dz2, 0.0, false);

            if (hit1) hit = minCollZ(*target.m_rightchild, bullet, dx, dy, dz2) || hit;
            
            //! Return minimum dz.
            dz = min(dz, dz2);

            return hit;
        } else {
            //! Bullet is not a leaf (neither is a leaf), check all left/right
            //! collision combinations.
            //! Target left and bullet left.
            hit1 = calcCollZ(target.m_leftchild->boundSphCentre(), target.m_leftchild->Radius(),
                             bullet.m_leftchild->boundSphCentre(), bullet.m_leftchild->Radius(),
                             dx, dy, dz, 0.0, false);

            if (hit1) hit = minCollZ(*target.m_leftchild, *bullet.m_leftchild, dx, dy, dz);

            //! Target left and bullet right.
            hit1 = calcCollZ(target.m_leftchild->boundSphCentre(), target.m_leftchild->Radius(),
                             bullet.m_rightchild->boundSphCentre(), bullet.m_rightchild->Radius(),
                             dx, dy, dz2, 0.0, false);

            if (hit1) hit = minCollZ(*target.m_leftchild, *bullet.m_rightchild, dx, dy, dz2) || hit;
            
            //! Target right and bullet left.
            hit1 = calcCollZ(target.m_rightchild->boundSphCentre(), target.m_rightchild->Radius(),
                             bullet.m_leftchild->boundSphCentre(), bullet.m_leftchild->Radius(),
                             dx, dy, dz3, 0.0, false);
            
            if (hit1) hit = minCollZ(*target.m_rightchild, *bullet.m_leftchild, dx, dy, dz3) || hit;

            //! Target right and bullet right.
            hit1 = calcCollZ(target.m_rightchild->boundSphCentre(), target.m_rightchild->Radius(),
                             bullet.m_rightchild->boundSphCentre(), bullet.m_rightchild->Radius(),
                             dx, dy, dz4, 0.0, false);

            if (hit1) hit = minCollZ(*target.m_rightchild, *bullet.m_rightchild, dx, dy, dz4) || hit;

            //! Returns minimum dz.
            dz = min(min(dz, dz2), min(dz3, dz4));

            return hit;
        }
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
    if (!m_pmodel->getTrackPrimaryCoordinates()) {
        Condition = (m_children_roundingLevel > 0.95);
    } else {
        Condition = (m_distance_centreToCentre == 0.0);
    }

    if ((Condition && m_leftparticle != NULL) || FakeRounding()) {
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
        if (m_pmodel->getTrackPrimaryCoordinates() && Numprimary() > 1 && m_leftchild != NULL) {
            d_ij = m_parent->m_distance_centreToCentre;
            r_i = m_primarydiam / 2.0;

            if (m_parent->m_leftparticle == this) {
                r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
            } else {
                r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
            }

            x_i = (pow(d_ij, 2.0) - pow(r_j, 2.0) + pow(r_i, 2.0)) / 2.0 / d_ij; //!< Eq. (3b) of Langmuir 27:6358 (2011).
            A_n = PI * (pow(r_i, 2.0) - pow(x_i, 2.0));                        //!< Eq. (4).
            A_i = 2.0 * PI * (pow(r_i, 2.0) + r_i * x_i);                      //!< Eq. (6).
            m_primary_diam_old = m_primarydiam;
            m_vol_old = m_vol;
        }

        int maxcarbon=0;

        for (vector<boost::shared_ptr<PAH> >::iterator i=m_PAH.begin(); i!=m_PAH.end(); ++i) {
            m_numcarbon += (*i)->m_pahstruct->numofC();
            m_numH += (*i)->m_pahstruct->numofH();
            m_numOfEdgeC += (*i)->m_pahstruct->numofEdgeC();
            m_numOfRings += (*i)->m_pahstruct->numofRings();
            maxcarbon = max(maxcarbon, (*i)->m_pahstruct->numofC()); //!< Search for the largest PAH-in terms of the number of carbon atoms-in the primary.
        }

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
        m_dcol = max(m_diam, m_PAHCollDiameter);
        m_surf = PI * m_diam * m_diam;
        
        //! If the distance between the centres of primary particles is
        //! tracked, the rate of change in the primary diameter is determined
        //! by its neighbour. Therefore, the particle should be made up of more
        //! than one primary.
        if (!m_pmodel->getTrackPrimaryCoordinates() || (m_pmodel->getTrackPrimaryCoordinates() && m_numprimary < 2)) {
            m_primarydiam = m_diam;
            setRadius(m_diam / 2.0);
        } else {
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
            m_primarydiam = m_primary_diam_old + 2 * (m_vol - m_vol_old) / (A_i + A_n * r_j / d_ij);
        }
        m_avg_coalesc = 0.0;
    }
}

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
		m_numprimary=m_leftchild->m_numprimary+m_rightchild->m_numprimary;
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
	}

    //this is not a primary, sum up the properties
	if (m_leftchild!=NULL)
    {
        // remove the PAHs from this node to free memory
        Reset();
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
            const double cdiam=max(aggcolldiam,m_PAHCollDiameter);
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
void PAHPrimary::Sinter(double dt, Cell &sys, const Processes::SinteringModel &model, rng_type &rng, double wt, bool &PAHTracerMatch, bool &particleChanged)
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
        if (!m_pmodel->getTrackPrimaryCoordinates()) {
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
                double A_n = PI * (r_i2 - pow(x_i, 2.0));        //!< Eq. (4).
                double A_i = 2.0 * PI * (r_i2 + r_i * x_i);      //!< Eq. (6).
                double A_j = 2.0 * PI * (r_j2 + r_j * x_j);      //!< Eq. (6).

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

                double V_i = 2.0 / 3.0 * PI * pow(r_i, 3.0) + PI * r_i2 * x_i - 1.0 / 3.0 * PI * pow(x_i, 3.0); //!< Eq. (3a).
                double V_j = 2.0 / 3.0 * PI * pow(r_j, 3.0) + PI * r_j2 * x_j - 1.0 / 3.0 * PI * pow(x_j, 3.0); //!< Eq. (3a).

                if (dd_ij_dt > 0.0) {
                    delt = dd_ij_Max / max(dd_ij_dt, 1.0e-300);

                    double mean;

                    if (tstop > (t1 + delt)) {
                        mean = 1.0 / scale;
                    } else {
                        mean = dd_ij_dt * (tstop - t1) / (scale * dd_ij_Max);
                    }

                    boost::random::poisson_distribution<unsigned, double> repeatDistribution(mean);
                    const unsigned n = repeatDistribution(rng);

                    if (n > 0) {
                        m_distance_centreToCentre -= (double)n * scale * dd_ij_Max; //!< Sintering decreases d_ij hence the negative sign.

                        //! The factor of 2 is because Eq. (8) is the rate of change
                        //! in radius.
                        this->m_leftparticle->m_primarydiam -= (double)n * scale * 2.0 * B_i * dd_ij_Max;  //!< Eq. (8).
                        this->m_rightparticle->m_primarydiam -= (double)n * scale * 2.0 * B_j * dd_ij_Max; //!< Eq. (8).

                        particleChanged = true;

                        //! If m_distance_centreToCentre (or d_ij) is some very small
                        //! value, it was found that B_i and B_j will be 1.#INF (or
                        //! infinite) and will result in an error when adjusting the
                        //! particle diameters below.
                        if (m_distance_centreToCentre < 1.0e-12) {
                            m_distance_centreToCentre = 0.0;
                        }

                        //! Only perform sintering if the distance between the centres
                        //! of primary particles i and j is positive.
                        if (m_distance_centreToCentre > 0.0) {
             //               boost::uniform_01<rng_type&, double> uniformGenerator(rng);
             //               Coords::Vector D;
             //               bool hit = false;
             //               while (!hit) {
             //                   D[0] = ((2.0 * uniformGenerator()) - 1.0) * m_distance_centreToCentre;
             //                   D[1] = ((2.0 * uniformGenerator()) - 1.0) * m_distance_centreToCentre;
             //                   hit = calcCollZ(this->m_leftparticle->boundSphCentre(), this->m_leftparticle->PrimaryDiam() / 2.0,
             //                                   this->m_rightparticle->boundSphCentre(), this->m_rightparticle->PrimaryDiam() / 2.0,
             //                                   D[0], D[1], D[2], m_distance_centreToCentre, true);
             //               }
					        //bool PAHTracerMatch = false;

                            Coords::Vector p1 = this->m_leftparticle->boundSphCentre();
                            Coords::Vector p2 = this->m_rightparticle->boundSphCentre();

                            double xdev = p2[0] - p1[0];
                            double ydev = p2[1] - p1[1];
                            double zdev = p2[2] - p1[2];

                            double dxsqr = xdev * xdev;
                            double dysqr = ydev * ydev;
                            double dzsqr = zdev * zdev;

                            double d = sqrt(dxsqr + dysqr + dzsqr);

                            double x = p1[0] + m_distance_centreToCentre * xdev / d;
                            double y = p1[1] + m_distance_centreToCentre * ydev / d;
                            double z = p1[2] + m_distance_centreToCentre * zdev / d;

                            this->m_rightchild->Translate(p2[0] - x, p2[1] - y, p2[2] - z, PAHTracerMatch, false, true);
                        }
                    }
                    t1 += delt;
                } else {
                    break;
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
        if (!m_pmodel->getTrackPrimaryCoordinates()) {
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
                m_leftchild->Sinter(dt, sys, model, rng, wt, PAHTracerMatch, particleChanged);
                m_rightchild->Sinter(dt, sys, model, rng, wt, PAHTracerMatch, particleChanged);
            }
        } else {
            m_leftchild->Sinter(dt, sys, model, rng, wt, PAHTracerMatch, particleChanged);
            m_rightchild->Sinter(dt, sys, model, rng, wt, PAHTracerMatch, particleChanged);
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

//! Returns the bounding-sphere centre.
const Coords::Vector &PAHPrimary::boundSphCentre(void) const
{
    return m_cen_bsph;
}

// Calculates the bounding sphere position and radius using
// the left and right child node values.
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
        setRadius(sqrt((dx*dx)+(dy*dy)+(dz*dz)));
    }
}

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
        for (unsigned int i=0; i!=3; ++i) {
            m_cen_mass[i]  = m_leftchild->m_cen_mass[i] * m_leftchild->m_mass;
            m_cen_mass[i] += m_rightchild->m_cen_mass[i] * m_rightchild->m_mass;
            m_cen_mass[i] *= invtotmass;
        }
    } else {
        //! If there are no children, then the centre-of-mass and bounding-
        //! sphere centre are the same.
        m_cen_mass[0] = m_cen_bsph[0];
        m_cen_mass[1] = m_cen_bsph[1];
        m_cen_mass[2] = m_cen_bsph[2];
    }
}

/*!
 *  @brief Calculates the z-displacement of a sphere.
 *
 *  Calculates the z-displacement of a bullet sphere for a +ve
 *  collision with a target sphere. Returns true if the
 *  spheres collide, otherwise false.
 *
 *  @param[in]  p1                      Coordinates of sphere 1.
 *  @param[in]  r1                      Radius of sphere 1.
 *  @param[in]  p2                      Coordinates of sphere 2.
 *  @param[in]  r2                      Radius of sphere 2
 *  @param[in]  dx                      Bullet x displacement.
 *  @param[in]  dy                      Bullet y displacement.
 *  @param[out] dz                      Bullet z displacement.
 *  @param[in]  distanceCentreToCentre  Distance between the centres of neighbouring primary particles.
 *  @param[in]  trackPrimaryCoordinates Flag used to indicate whether the primary coordinates are tracked.
 *
 *  @return Have the nodes collided?
 */
bool PAHPrimary::calcCollZ(const Coords::Vector &p1, double r1,
                           const Coords::Vector &p2, double r2,
                           double dx, double dy, double &dz,
                           double distanceCentreToCentre, const bool trackPrimaryCoordinates)
{
    double sumrsqr;

    if (!trackPrimaryCoordinates) {
        sumrsqr = r1 + r2;
    } else {
        sumrsqr = distanceCentreToCentre;
    }

    //! Calculate the square of the sum of the radii, or the distance between
    //! the centres of neighbouring primary particles if tracking primary
    //! separation.
    sumrsqr *= sumrsqr;

    //! Calculate dx, dy and dz. Remember to include
    //! argument contributions.
    double xdev = p2[0] - p1[0] + dx;
    double ydev = p2[1] - p1[1] + dy;
    double zdev = p2[2] - p1[2];

    //! Calculate dx, dy and dz squared.
    double dxsqr = xdev * xdev;
    double dysqr = ydev * ydev;
    double dzsqr = zdev * zdev;

    // Calculate quadratic terms.
    double b = 2.0 * zdev;
    double c = dxsqr + dysqr + dzsqr - sumrsqr;

    //! Calculate discriminant.
    double dis = (b*b) - (4.0*c);

    if (dis >= 0.0) {
        //! Spheres intersect.
        dz = - 0.5 * (b + sqrt(dis));
        return true;
    } else {
        //! Spheres do not intersect.
        dz = 1.0e10; //!< A large number.
        return false;
    }
}

/*!
 *  Put the bounding-sphere at the origin.
 *
 *  @param[in] PAHTracerMatch Flag used to indicate whether node contains PAH
 *                            to be traced.
 */
void PAHPrimary::centreBoundSph(bool PAHTracerMatch)
{
    Translate(-m_cen_bsph[0], -m_cen_bsph[1], -m_cen_bsph[2], PAHTracerMatch, false, true);
}

//! Put the centre-of-mass at the origin.
void PAHPrimary::centreCOM(void)
{
    Translate(-m_cen_mass[0], -m_cen_mass[1], -m_cen_mass[2], false, false, false);
}

//! Check whether the particle pointed to by the this pointer contains the PAH
//! to be traced.
bool PAHPrimary::PAHTracerCheck(void)
{
	bool PAHTracerMatch = false;

    if (m_leftchild != NULL)
		PAHTracerMatch = m_leftchild->PAHTracerCheck();

	if (PAHTracerMatch)
		return PAHTracerMatch;

    if (m_rightchild != NULL)
        PAHTracerMatch = m_rightchild->PAHTracerCheck();

	if (PAHTracerMatch)
		return PAHTracerMatch;

    for (size_t i = 0; i != m_PAH.size(); ++i) {
        void *ptr = m_PAH[i].get();
        if(ptr == PAHTracer) {
			PAHTracerMatch = true;
			break;
        }
    }
	return PAHTracerMatch;
}

/*!
 *  Rotates the aggregate node and child structure about its centre
 *  of mass by the given angles (spherical coordinates).
 *
 *  @param[in] dtheta         Change in polar angle.
 *  @param[in] dphi           Change in azimuthal angle.
 *  @param[in] PAHTracerMatch Flag used to indicate whether particle contains
 *                            PAH to be traced. 
 */
void PAHPrimary::rotateCOM(double dtheta, double dphi, bool PAHTracerMatch)
{
    //! Move the aggregate so that its centre-of-mass is at the origin. Store
    //! the coordinates, so that they can be restored afterwards.
    //! Do not bother to write POV-Ray translation of nodes.
    Coords::Vector D(m_cen_mass);
    Translate(-D.X(), -D.Y(), -D.Z(), false, false, false);

    //! Create transformation matrix.
    Coords::Matrix M;
    M.SetIdentity();
    M.Rotate(dtheta, dphi);

    //! Rotate child nodes.
    if (m_leftchild != NULL) m_leftchild->transform(M, PAHTracerMatch);
    if (m_rightchild != NULL) m_rightchild->transform(M, PAHTracerMatch);

    //! Rotate bounding-sphere coordinates.
    m_cen_bsph = M.Mult(m_cen_bsph);

    //! Restore centre-of-mass coordinates.
    //! Do not bother to write POV-Ray translation of nodes.
    Translate(D.X(), D.Y(), D.Z(), false, false, false);
}

void PAHPrimary::inverseRotateCOM(double dtheta, double dphi, bool PAHTracerMatch)
{
    //! Move the aggregate so that its centre-of-mass is at the origin. Store
    //! the coordinates, so that they can be restored afterwards.
    //! Do not bother to write POV-Ray translation of nodes.
    Coords::Vector D(m_cen_mass);
    Translate(-D.X(), -D.Y(), -D.Z(), false, false, false);

    //! Create transformation matrix.
    Coords::Matrix M;
    M.SetIdentity();
    M.Rotate(dtheta, dphi);

    M.Minors();
    M.CoFactors();
    M.Adjugate();
    M.Inverse();

    //! Rotate child nodes.
    if (m_leftchild != NULL) m_leftchild->transform(M, PAHTracerMatch);
    if (m_rightchild != NULL) m_rightchild->transform(M, PAHTracerMatch);

    //! Rotate bounding-sphere coordinates.
    m_cen_bsph = M.Mult(m_cen_bsph);

    //! Restore centre-of-mass coordinates.
    //! Do not bother to write POV-Ray translation of nodes.
    Translate(D.X(), D.Y(), D.Z(), false, false, false);
}

/*!
 *  @brief Sets the radius of the bounding sphere.
 *
 *  @param[in] r Radius of bounding sphere.
 */
void PAHPrimary::setRadius(double r)
{
    m_r  = r;
    m_r2 = r * r;
    m_r3 = m_r2 * m_r;
}

//! Returns the bounding sphere radius.
double PAHPrimary::Radius(void) const
{
    return m_r;
}

/*!
 *  Transform the primary particle coordinates using the transformation matrix
 *  so as to rotate it.
 *
 *  @param[in] mat            Transformation matrix.
 *  @param[in] PAHTracerMatch Flag used to indicate whether the particle
 *                            contains the PAH to be traced.
 */
void PAHPrimary::transform(const Coords::Matrix &mat, bool PAHTracerMatch)
{
    //! Descend binary tree to the leaf nodes, i.e. single primary particles.
    if (m_leftchild != NULL) 
        m_leftchild->transform(mat, PAHTracerMatch);
    if (m_rightchild != NULL) 
        m_rightchild->transform(mat, PAHTracerMatch);
    
    //! If the particle contains the PAH tracer (PAHTracerMatch true) then
    //! write to file the following POV-Ray commands. However, only do this for
    //! leaf nodes, hence, the second and third conditions. As the code ascends
    //! the binary tree to exit this function, the root node will pass through
    //! this part of the code.
	if (PAHTracerMatch && isLeaf()) {
		std::ofstream outfile;

        //! The clock variable runs progressively from 0.0 to 1.0. The counter
        //! is used to indicate the stage of the animation.
		string clock = "*clock" + cstr(getCounter());

        //! Trigonometric functions in the transformation matrix.
		string sin_phi = "sin(phi" + clock + ")";
		string cos_phi = "cos(phi" + clock + ")";
		string sin_theta = "sin(theta" + clock + ")";
		string cos_theta = "cos(theta" + clock + ")";
		
		outfile.open("sphere.pov", std::ios_base::app);

        //! Write the coordinates of the centre-of-mass (units: nm) and radius
        //! (units: nm) before rotation.
		outfile << "\t\tsphere{<" << m_cen_mass[0] * 1.0e9 << "," << m_cen_mass[1] * 1.0e9 << "," << m_cen_mass[2] * 1.0e9 << ">," << Radius() * 1.0e9 << " texture{tem_finish}\n"
		
        //! Given the following transformation matrix:
        //! matrix <Val00, Val01, Val02,
        //!         Val10, Val11, Val12,
        //!         Val20, Val21, Val22,
        //!         Val30, Val31, Val32>
        //!
        //! the point P <px, py, pz> is transformed into the point Q
        //! <qx, qy, qz> by:
        //! qx = Val00 * px + Val10 * py + Val20 * pz + Val30,
        //! qy = Val01 * px + Val11 * py + Val21 * pz + Val31,
        //! qz = Val02 * px + Val12 * py + Val22 * pz + Val32.
        //! 
        //! Note that the matrix multiplication that is performed below is like
        //! so:
        //! qx = cos(phi) * px - sin(phi)cos(theta) * py + sin(phi)sin(theta) * pz,
        //! qy = sin(phi) * px + cos(phi)cos(theta) * py + cos(phi)sin(theta) * pz,
        //! qz =                         sin(theta) * py +         cos(theta) * pz.
        //!
        //! Therefore the transformation matrix in POV-Ray has to be flipped
        //! over the diagonal. The last row of the matrix is 0 as we do not
        //! translate the point P.
                << "\t\t\tmatrix<" << cos_phi << "," << sin_phi << ",0,\n"
				<< "\t\t\t-" << sin_phi << "*" << cos_theta << "," << cos_phi << "*" << cos_theta << sin_theta << ",\n"
				<< "\t\t\t" << sin_phi << "*" << sin_theta << ",-" << cos_phi << "*" << sin_theta << "," << cos_theta << ",\n"
				<< "\t\t\t0,0,0>}\n"; 
		outfile.close();
	}

    //! Rotate centre-of-mass and bounding sphere coordinates.
	m_cen_mass = mat.Mult(m_cen_mass);
    m_cen_bsph = mat.Mult(m_cen_bsph);
}

/*!
 *  Translates (moves) the aggregate node and child structure by the given
 *  amounts along the cartesian axes.
 *
 *  @param[in] dx             Distance to translate in the x-axis.
 *  @param[in] dy             Distance to translate in the y-axis.
 *  @param[in] dz             Distance to translate in the z-axis.
 *  @param[in] PAHTracerMatch Flag used to indicate whether particle contains
 *                            PAH to be traced.
 *  @param[in] Coagulate      Flag used to indicate whether this is the
 *                            translation of a particle to coagulate with
 *                            another particle.
 *  @param[in] Aggregate      Flag used to indicate whether this is the
 *                            translation of an aggregate.
 */
void PAHPrimary::Translate(double dx, double dy, double dz, bool PAHTracerMatch, bool Coagulate, bool Aggregate)
{
    //! Translate child branches.
    if (m_leftchild != NULL) m_leftchild->Translate(dx, dy, dz, PAHTracerMatch, Coagulate, Aggregate);
    if (m_rightchild != NULL) m_rightchild->Translate(dx, dy, dz, PAHTracerMatch, Coagulate, Aggregate);

    //! Translate bounding sphere centre.
    m_cen_bsph.Translate(dx, dy, dz);

    std::ofstream outfile;

    //! If the flag 'Aggregate' is true, this is a call from centreBoundSph.
    //! If the particle contains the PAH tracer (PAHTracerMatch true) then
    //! write to file the following POV-Ray commands. However, only do this for
    //! leaf nodes, hence, the second and third conditions. As the code ascends
    //! the binary tree to exit this function, the root node will pass through
    //! this part of the code.
	if (PAHTracerMatch && Aggregate && isLeaf()) {
        //! The clock variable runs progressively from 0.0 to 1.0. The counter
        //! is used to indicate the stage of the animation.
		string clock = "*clock" + cstr(getCounter());
		
        //! First write the original centre-of-mass and then let POV-Ray
        //! translate the particle. Units of nm.
		outfile.open("sphere.pov", std::ios_base::app);
		outfile << "\t\tsphere{<" << m_cen_mass[0] * 1.0e9 << "," << m_cen_mass[1] * 1.0e9 << "," << m_cen_mass[2] * 1.0e9 << ">," << Radius() * 1.0e9
				<< " translate<" << dx * 1.0e9 << clock << "," << dy * 1.0e9 << clock << "," << dz * 1.0e9 << clock << ">"
				<< " texture{tem_finish}}\n";
		outfile.close();		
	}

    //! Translate centre-of-mass.
    m_cen_mass.Translate(dx, dy, dz);

    //! If the flag 'Coagulate' is true, this is a call from Coagulate.
    //! If the particle contains the PAH tracer (PAHTracerMatch true) then
    //! write to file the following POV-Ray commands. However, only do this for
    //! leaf nodes, hence, the second and third conditions. As the code ascends
    //! the binary tree to exit this function, the root node will pass through
    //! this part of the code.
	if (PAHTracerMatch && Coagulate && isLeaf()) {
        //! The clock variable runs progressively from 0.0 to 1.0. The counter
        //! is used to indicate the stage of the animation.
		string clock = "*clock" + cstr(getCounter());
		
		outfile.open("sphere.pov", std::ios_base::app);

        //! We place the particle at an arbitrarily large distance (10 times
        //! its centre-of-mass) away from its position and allow POV-Ray to
        //! translate it to its new centre-of-mass so that it appears that the 
        //! particle moves into the field of vision and collides and sticks
        //! with a particle.
        outfile << "\t\tsphere{<" << (m_cen_mass[0] + dx * 10.0) * 1.0e9 << "," << (m_cen_mass[1] + dy * 10.0) * 1.0e9 << "," << (m_cen_mass[2] + dz * 10.0) * 1.0e9 << ">," << Radius() * 1.0e9
				<< " translate<" << -dx * 1.0e10 << clock << "," << -dy * 1.0e10 << clock << "," << -dz * 1.0e10 << clock << ">"
				<< " texture{tem_finish}}\n";
		outfile.close();

		//outfile << "\t\tsphere{<" << m_cen_mass[0] * 1.0e10 << "," << m_cen_mass[1] * 1.0e10 << "," << m_cen_mass[2] * 1.0e10 << ">," << Radius() * 1.0e9
		//		<< " translate<" << - m_cen_mass[0] * 9.0e9 << clock << "," << - m_cen_mass[1] * 9.0e9 << clock << "," << - m_cen_mass[2] * 9.0e9 << clock << ">"
		//		<< " texture{tem_finish}}\n";
		//outfile.close();	
	}
}

//! Write the coordinates of the primaries in the particle pointed to by the
//! this pointer. Units of nm.
void PAHPrimary::writePrimaryCoordinatesRadius(void)
{
    if (m_leftchild!=NULL)
        m_leftchild->writePrimaryCoordinatesRadius();

    if (m_rightchild!=NULL)
        m_rightchild->writePrimaryCoordinatesRadius();

    std::ofstream outfile;

    //! There is a need for these two conditions as upon exit of the function
    //! the aggregate node will pass through this part of the code but are only
    //! we are only interested in the coordinates of the primaries.
	if (isLeaf()) {
		outfile.open("sphere.pov", std::ios_base::app);
		outfile << "\t\tsphere{<" << m_cen_mass[0] * 1.0e9 << "," << m_cen_mass[1] * 1.0e9 << "," << m_cen_mass[2] * 1.0e9 << ">," << Radius() * 1.0e9
				<< " texture{tem_finish}}\n";
		outfile.close();
	}
}