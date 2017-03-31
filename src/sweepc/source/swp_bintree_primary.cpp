/*!
 * @file    swp_bintree_primary.cpp
 * @author  William J Menz
 * @brief   Implementation of generic binary tree particle model
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      Implementation of the BinTreePrimary class.
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public License
 *      as published by the Free Software Foundation; either version 2
 *      of the License, or (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this program; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *      02111-1307, USA.
 *
 *   Contact:
 *      Prof Markus Kraft
 *      Dept of Chemical Engineering
 *      University of Cambridge
 *      New Museums Site
 *      Pembroke Street
 *      Cambridge
 *      CB2 3RA, UK
 *
 *      Email:       mk306@cam.ac.uk
 *      Website:     http://como.cheng.cam.ac.uk
*/

//******************************************************csl37
#define _USE_MATH_DEFINES //!< First define.
#include <math.h>         //!< Then include so that the pi constant (M_PI) can be used.
//******************************************************csl37

#include "swp_bintree_primary.h"
#include "swp_bintree_serializer.h"
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <stack>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

/*!
 * @brief       Initialising constructor with no arguments
 *
 * Initialises a particle from no arguments.
 *
 */
BinTreePrimary::BinTreePrimary() : Primary(),
    m_numprimary(0),
    m_primarydiam(0.0),
    m_children_radius(0.0),
    m_children_vol(0.0),
    m_children_surf(0.0),
	m_free_surf(0.0),
    m_distance_centreToCentre(0.0),
    m_children_sintering(0.0),
    m_avg_sinter(0.0),
    m_sint_rate(0.0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL)
{}

/*!
 * @brief       Initialising constructor using time and model
 *
 * Initialises a particle with the default chemical properties. Calls
 * UpdateCache() to calculate the derived properties after creation.
 *
 * @param[in]   time        Time at which particle is being created
 * @param[in]   model       Model which defines the meaning of the primary
 *
 */
BinTreePrimary::BinTreePrimary(const double time,
        const Sweep::ParticleModel &model)
: Primary(time, model),
    m_numprimary(0),
    m_primarydiam(0.0),
    m_children_radius(0.0),
    m_children_vol(0.0),
    m_children_surf(0.0),
	m_free_surf(0.0),
    m_distance_centreToCentre(0.0),
    m_children_sintering(0.0),
    m_avg_sinter(0.0),
    m_sint_rate(0.0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL)

{}

//! Copy constructor.
BinTreePrimary::BinTreePrimary(const BinTreePrimary &copy)
{
    *this = copy;
    if (copy.m_leftchild!=NULL)
    {
        CopyTree(&copy);
    }
}

//! Stream-reading constructor.
BinTreePrimary::BinTreePrimary(std::istream &in,
        const Sweep::ParticleModel &model) :
m_numprimary(0),
m_primarydiam(0.0),
m_children_radius(0.0),
m_children_vol(0.0),
m_children_surf(0.0),
m_free_surf(0.0),
m_distance_centreToCentre(0.0),
m_children_sintering(0.0),
m_avg_sinter(0.0),
m_sint_rate(0.0),
m_sint_time(0.0),
m_leftchild(NULL),
m_rightchild(NULL),
m_parent(NULL),
m_leftparticle(NULL),
m_rightparticle(NULL)
{
    Deserialize(in, model);
}

//! Default destructor.
BinTreePrimary::~BinTreePrimary()
{
    delete m_leftchild;
    delete m_rightchild;

    releaseMem();
}


// Primary operator
BinTreePrimary &BinTreePrimary::operator=(const Primary &rhs)
{
    operator=(dynamic_cast<const BinTreePrimary&>(rhs));

    return *this;
}

//! Return a copy of the model data
BinTreePrimary *const BinTreePrimary::Clone(void) const
{
   return new BinTreePrimary(*this);
}


/*!
 * @brief       Recursively copy the tree for non-leaf nodes.
 *
 * @param[in] source Pointer to the primary to be copied
*/
void BinTreePrimary::CopyTree(const BinTreePrimary *source)
{
    // Create the new left and right children with nothing in them
    m_leftchild     = new BinTreePrimary(source->CreateTime(),*m_pmodel);
    m_rightchild    = new BinTreePrimary(source->CreateTime(),*m_pmodel);

    // Copy the properties such as the volume,
    // surface area and list of constituent chemical units
    m_leftchild->CopyParts(source->m_leftchild);
    m_rightchild->CopyParts(source->m_rightchild);

    // Set the pointers to the parent
    m_leftchild->m_parent=this;
    m_rightchild->m_parent=this;

    // The left and right particle are set further down in UpdateAllPointers
    // These are the pointers that specify which primary particles touch each
    // other in the aggregate structure.
    m_leftparticle=NULL;
    m_rightparticle=NULL;

    // Recursively copy the subtrees
    if (source->m_leftchild->m_leftchild!=NULL)
        m_leftchild->CopyTree(source->m_leftchild);

    if (source->m_rightchild->m_leftchild!=NULL)
        m_rightchild->CopyTree(source->m_rightchild);

    // Set the leftparticle and rightparticle
    UpdateAllPointers(source);

}

/*!
 * @brief       Updates all pointers, ensuring the connectivity of the tree
 *                 is preserved during copying.
 *
 * Each node contains two pointers (m_leftparticle and m_rightparticle)
 * to primary particles that are connected by this node
 * This function is used when the entire particle tree is duplicated.
 * It sets the pointers in the copied node (this), such that the connectivity
 * of the primary particles in this node is the same as in the original node.
 *
 * @param[in] original Pointer to the primary to be copied
*/
void BinTreePrimary::UpdateAllPointers(const BinTreePrimary *original)
{
    // The primary has no children => there are no left and right particles
    if (original->m_leftchild == NULL) {
        m_leftparticle=NULL;
        m_rightparticle=NULL;
     }
    else {
        // Find the route to m_leftparticle in the original tree
        std::stack<bool> route =recordPath(original->m_leftparticle, original);

        // Now follow the same route down the new tree to find
        // the new left particle
        m_leftparticle = descendPath(this, route);

         // Find the route to m_rightparticle in the original tree
        route = recordPath(original->m_rightparticle, original);

        // Now follow the same route down the new tree to find
        // the new right particle
        m_rightparticle = descendPath(this, route);
    }
}


/*!
 *
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
std::stack<bool> BinTreePrimary::recordPath(const BinTreePrimary* bottom,
                                        const BinTreePrimary* const top)
{
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
 *@param[in,out]    takeLeftBranch  Instructions for which child to move to
 *                                   at each level
 *
 *@return   The node at the bottom of the path
 *
 *@pre  here must be a node of tree in which takeLeftBranch is a valid path
 *@post takeLeftBranch.empty() == true
 */
BinTreePrimary* BinTreePrimary::descendPath(BinTreePrimary *here,
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
 * @brief       Coagulates *this* and rhs particle
 *
 * Called by coagulate, this function joins sp1 (this) and sp2 (rhs).
 * Basic properties are first calculated, followed by derived.
 *
 * @param[in] rhs   Pointer to the particle to be coagulated with this
 * @param[in] rng   Random number generator
*/
BinTreePrimary &BinTreePrimary::Coagulate(const Primary &rhs, rng_type &rng)
{
    const BinTreePrimary *rhsparticle = NULL;
    rhsparticle = dynamic_cast<const AggModels::BinTreePrimary*>(&rhs);

    // Don't sum-up components now. Wait for UpdateCache so as to not
    // double-count the values.

    // Create the new particles
    BinTreePrimary *newleft = new BinTreePrimary(m_time, *m_pmodel);
    BinTreePrimary *newright = new BinTreePrimary(m_time, *m_pmodel);
    BinTreePrimary copy_rhs(*rhsparticle);


    //Randomly select where to add the second particle
    boost::bernoulli_distribution<> bernoulliDistrib;
    if (bernoulliDistrib(rng)) {
        newleft->CopyParts(this);
        newright->CopyParts(&copy_rhs);
    }
    else {
        newright->CopyParts(this);
        newleft->CopyParts(&copy_rhs);
    }

    // Set the pointers
    m_leftchild     = newleft;
    m_rightchild    = newright;
    newright->m_parent = this;
    newleft->m_parent  = this;

    // Set the pointers to the parent node
    if (newleft->m_leftchild!=NULL) {
    newleft->m_leftchild->m_parent      = newleft;
    newleft->m_rightchild->m_parent     = newleft;
    }
    if (newright->m_leftchild!=NULL) {
    newright->m_leftchild->m_parent     = newright;
    newright->m_rightchild->m_parent    = newright;
    }
    m_children_sintering=0.0;
    UpdateCache();

    // Select the primaries that are touching
    m_leftparticle      = m_leftchild->SelectRandomSubparticle(rng);
    m_rightparticle     = m_rightchild->SelectRandomSubparticle(rng);

    // Set the sintering times
    SetSinteringTime(std::max(m_sint_time, rhsparticle->m_sint_time));
    m_createt = max(m_createt, rhsparticle->m_createt);

    // Initialise the variables used to calculate the sintering level
    m_children_vol      = m_leftparticle->m_vol + m_rightparticle->m_vol;
    m_children_surf     = m_leftparticle->m_surf + m_rightparticle->m_surf;
    m_children_radius   = pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
    m_children_sintering= SinteringLevel();
	m_free_surf = m_leftparticle->m_free_surf + m_rightparticle->m_free_surf;
	m_distance_centreToCentre = m_leftparticle->m_primarydiam/2.0+m_rightparticle->m_primarydiam/2.0;
	CheckSintering();

    // Must set all the pointer to NULL otherwise the delete function
    // will also delete the children
    copy_rhs.m_leftchild    = NULL;
    copy_rhs.m_rightchild   = NULL;
    copy_rhs.m_parent       = NULL;
    copy_rhs.m_leftparticle = NULL;
    copy_rhs.m_rightparticle= NULL;

    return *this;
}

/*!
 *  @brief Copy state-space and derived properties from source.
 *
 *  This function is like a limited assignment operator, except that the
 *  children are not copied and the pointers to the particles may need
 *  adjusting after this method has finished.
 *
 *  @param[in] source Pointer to the primary to be copied.
 */
void BinTreePrimary::CopyParts(const BinTreePrimary *source)
{
    //! Set primary characteristics.
    SetComposition(source->Composition());
    SetValues(source->Values());
    SetTime(source->LastUpdateTime());
    SetCollDiameter(source->CollDiameter());
    SetMobDiameter(source->MobDiameter());
    SetSphDiameter(source->SphDiameter());
    SetSurfaceArea(source->SurfaceArea());
    SetVolume(source->Volume());
    SetMass(source->Mass());

    //! Set BinTreePrimary model characteristics.
    m_numprimary              = source->m_numprimary;
    m_primarydiam             = source->m_primarydiam;
    m_children_radius         = source->m_children_radius;
    m_children_vol            = source->m_children_vol;
    m_children_surf           = source->m_children_surf;
	m_free_surf				  = source->m_free_surf;
    m_distance_centreToCentre = source->m_distance_centreToCentre;
    m_children_sintering      = source->m_children_sintering;
    m_avg_sinter              = source->m_avg_sinter;
    m_sint_rate               = source->m_sint_rate;
    m_sint_time               = source->m_sint_time;

    //! Set particles.
    m_leftchild     = source->m_leftchild;
    m_rightchild    = source->m_rightchild;
    m_parent        = source->m_parent;
    m_leftparticle  = source->m_leftparticle;
    m_rightparticle = source->m_rightparticle;
}

/*!
 * @brief           Recursively set the sintering time
 *
 * Needed to ensure that every node has an up-to-date value of the
 * sintering time; otherwise coagulation and merge events can lose
 * this information.
 *
 * @param time      Total time for which particles have been sintered
 */
void BinTreePrimary::SetSinteringTime(double time) {
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

/*!
 * @brief       Overload of the SetTime function for BinTreePrimary
 *
 * Sets the LDPA update time throughout the binary tree. This is important
 * as Merge can sometimes delete the originial root node, losing m_time
 * stored in m_primary. This potentially leads to longer dt LDPA action
 * times when UpdateParticle is called.
 *
 * @param t     LDPA update time
 */
void BinTreePrimary::SetTime(double t) {
    m_time = t;

    // Set LDPA time of children
    if (m_leftchild != NULL) {
        m_leftchild->SetTime(t);
        m_rightchild->SetTime(t);
    }

    // Set LDPA time of particles
    if (m_leftparticle != NULL) {
        m_leftparticle->SetTime(t);
        m_rightparticle->SetTime(t);
    }
}

/*!
 * @brief       Randomly selects a primary in the binary tree
 *
 * Select a primary uniformly at random from this particle
 * and descend the aggregate tree to find the primary.
 * Note that most BinTreePrimaries are nodes in a tree representing
 * connectivity within an aggregate, so it is necessary to
 * descend the tree to find a primary that really is a
 * primary.
 *
 * @param[in]   rng     Random number generator
 *
 * @return      Pointer to an object representing a physical primary
 */
BinTreePrimary *BinTreePrimary::SelectRandomSubparticle(rng_type &rng)
{
    // We want to choose an integer uniformly on the range [0, m_numprimary - 1]
    boost::random::uniform_smallint<int> uniformDistrib(0, m_numprimary - 1);

    return SelectRandomSubparticleLoop(uniformDistrib(rng));

}
/*!
 * @brief       Helper function for SelectRandomSubparticle
 *
 * @param[in] target The primary to be selected
 *
 * @return      Pointer to the next node down the tree on the path
*/
BinTreePrimary *BinTreePrimary::SelectRandomSubparticleLoop(int target)
{
    if (m_leftchild==NULL) return this;
    if (target<=m_leftchild->m_numprimary)
    {
        return m_leftchild->SelectRandomSubparticleLoop(target);
    }
    else
    {
        return m_rightchild->SelectRandomSubparticleLoop(
                target-(m_leftchild->m_numprimary));
    }
}

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
double BinTreePrimary::SinteringLevel()
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
			} else {
				slevel= ((spherical_surface/m_children_surf) - TWO_ONE_THIRD)
						/(1 - TWO_ONE_THIRD);
			}

		// if centre-centre separation is tracked
		// s = ((d_min/d_ij) - (d_min/d_mx))/(1-d_min/d_max)
		// where: d_min = d_ij at merger, d_mx = r_i + r_j
		}else{
			if (m_leftparticle != NULL && m_rightparticle != NULL) {
				double r_i = m_leftparticle->m_primarydiam/2.0;
				double r_j = m_rightparticle->m_primarydiam/2.0;
				double d_ij =  m_distance_centreToCentre;

				//calculate the merger condition
				double d_ij_merge = sqrt( pow(max(r_i,r_j),2.0) - pow(min(r_i,r_j),2.0) );
				slevel = ( d_ij_merge/d_ij - d_ij_merge/(r_i + r_j) ) / ( 1 - d_ij_merge/(r_i + r_j) );
			}
		}

		if (slevel < 0.0) {
				return 0.0;
			} else if (slevel > 1.0) {
				return 1.0;
			} else return slevel;

    } else {
        // Particle is a primary
        m_children_surf = 0.0;
        m_children_radius = 0.0;
        m_children_vol = 0.0;
        if (m_parent == NULL) return 1.0;        // Single primary case
        else return 0.0;                         // Part of a tree
    }
}

/*!
 * @brief       Checks the sintering level of the particle
 *
 * If the sintering level is above 95%, Merge is called and the cache
 * is updated.
 *
 * @return      Boolean telling if the particle has sintered
 */
bool BinTreePrimary::CheckSintering()
{
    bool hassintered=false;

	if(m_leftparticle != NULL) {

		// check whether condition for merger is met
		if (MergeCondition()) {
			 Merge();
			 UpdateCache();
			 hassintered=true;

			 // Check again because this node has changed
			 CheckSintering();
		}
	}
    if (m_leftchild!=NULL) {
        hassintered = m_leftchild->CheckSintering();
        hassintered = m_rightchild->CheckSintering();
    }

    return hassintered;
}

//***********************csl37
// function to check the condition for merger
// returns true if the condition is met
bool BinTreePrimary::MergeCondition()
{
	bool condition=false;

	if(m_leftparticle != NULL) {
		//! The condition for whether a particle has coalesced depends on whether
		//! the distance between the centres of primary particles is tracked. 
		//! If not, the condition depends on whether the rounding 
		//! level exceeds an arbitrarily high threshold.
		if (!m_pmodel->getTrackPrimarySeparation()) {
			condition = (m_children_sintering > 0.95);
		} else {
			double r_i = m_leftparticle->m_primarydiam/2.0;
			double r_j = m_rightparticle->m_primarydiam/2.0;
			double d_ij =  m_distance_centreToCentre;
			
			//ensures that particles are merged if sintering step overshoots slightly
			if(d_ij <= 0.0){
				condition = true;
			}else{
				// a particle has coalesced if the neck radius reaches the centre of one
				// primary. The surface growth model fails if the neck resides outside 
				// the region between the primary centres. 
				// min(x_ij,x_ji) <= 0;
				condition = ( (pow(d_ij,2.0) - pow(max(r_i,r_j),2.0) + pow(min(r_i,r_j),2.0) )/(2.0*d_ij) <= 0.0);
			}
		}
	}

	return condition;
}

/*!
 * @brief       Merges the left and right particles
 *
 * Considers two cases of particle structure. If the node only has two
 * primaries in its subtree, they are simply deleted and *this* becomes
 * the new primary. Otherwise, the subtrees are moved in the manner
 * explained in Markus Sander's thesis.
 *
 * @return      Pointer to new merged particle
 */
BinTreePrimary &BinTreePrimary::Merge()
{
	//initialise pointer to the smaller of the merging primaries
	BinTreePrimary *small_prim;
	//pointer to the new (merged) primary
	BinTreePrimary *new_prim;

    // Make sure this primary has children to merge
    if( m_leftchild!=NULL) {

		//if the centre to centre distance is tracked we need to know which is the smaller primary of the merging pair
		if(m_leftparticle->m_primarydiam > m_rightparticle->m_primarydiam){
			small_prim = m_rightparticle;
		}else{
			small_prim = m_leftparticle;
		}

		double dV = small_prim->m_vol; //volume of the merging (smaller) primary 

        if (m_leftchild==m_leftparticle && m_rightchild==m_rightparticle)
        {
            // This node has only two primaries in its subtree, it is possible
            // that this node is not the root node and belongs to a bigger
            // particle.		

            // Sum up the components first
            for (size_t i=0; i != m_comp.size(); i++) {
                m_comp[i] = m_leftparticle->Composition(i) +
                        m_rightparticle->Composition(i);
            }

			//m_primarydiam of the new particle is the larger of the diameters of the merging particles
			//if the centre to centre separation isn't tracked this will be changed to the spherical equivalent in the call to UpdatePrimary
			m_primarydiam = max(m_leftparticle->m_primarydiam,m_rightparticle->m_primarydiam);
			
			//new primary
			new_prim = this;

            // Update the pointers that pointed to the two former children
			if (!m_pmodel->getTrackPrimarySeparation()) {
		        ChangePointer(m_leftchild,this);
			    ChangePointer(m_rightchild,this);
			}else{
				// If the priamry separation is tracked then re-estimation of new sintering levels also
				// requires the centre-centre separation and the smaller primary
				ChangePointer(m_leftchild,this,m_distance_centreToCentre,small_prim);
				ChangePointer(m_rightchild,this,m_distance_centreToCentre,small_prim);
			}

            // Delete the children (destructor is recursive for this class)
            delete m_leftchild;
            delete m_rightchild;
            m_leftchild=NULL;
            m_rightchild=NULL;
            m_leftparticle=NULL;
            m_rightparticle=NULL;

            // Set the children properties to zero, this node has no
            // more children
            ResetChildrenProperties();
            UpdatePrimary();

            // Only update the cache on m_parent if the sintering level of
            // m_parent if the sintering level won't call a merge on the
            // parent node. Otherwise, the *this* memory address could be
            // removed from the tree and segmentation faults will result!
            if(m_parent!=NULL) { 
				if (!m_parent->MergeCondition()) {
                    m_parent->UpdateCache();
                }
            }
        }


        else
        {
            if (m_leftchild->m_numprimary<m_rightchild->m_numprimary)
            {
                // Append to left subtree because there are fewer primaries
                // (this is only to keep the tree balanced)
                BinTreePrimary *oldleftparticle = m_leftparticle;
                for (size_t i=0; i != m_comp.size(); i++) {
                    m_rightparticle->m_comp[i] =
                            m_leftparticle->Composition(i) +
                            m_rightparticle->Composition(i);
                }

				//m_primarydiam of the new particle is the larger of the diameters of the merging particles
				//if the centre to centre separation isn't tracked this will be changed to the spherical equivalent in the call to UpdatePrimary
				m_rightparticle->m_primarydiam = max(m_leftparticle->m_primarydiam,m_rightparticle->m_primarydiam);

				//the new (merged) primary
				new_prim = m_rightparticle;

                m_rightparticle->UpdatePrimary();

                // Set the pointers from the leftprimary to the rightprimary
                // (this will be the new bigger primary)
				if (!m_pmodel->getTrackPrimarySeparation()) {
				    oldleftparticle->ChangePointer(oldleftparticle,m_rightparticle);
					m_rightparticle->ChangePointer(m_rightparticle,m_rightparticle);
				}else{
					// If the priamry separation is tracked then re-estimation of new sintering levels also
					// requires the centre-centre separation and the smaller primary
					oldleftparticle->ChangePointer(oldleftparticle,m_rightparticle,m_distance_centreToCentre,small_prim);
					m_rightparticle->ChangePointer(m_rightparticle,m_rightparticle,m_distance_centreToCentre,small_prim);
				}

                // Set the pointer to the parent node
                if (oldleftparticle->m_parent->m_leftchild==oldleftparticle) {
                    oldleftparticle->m_parent->m_leftchild=m_rightchild;
                }
                else {
                    oldleftparticle->m_parent->m_rightchild=m_rightchild;
                }
                m_rightchild->m_parent=oldleftparticle->m_parent;

                BinTreePrimary *oldleftchild    = m_leftchild;
                BinTreePrimary *oldparent       = m_parent;

                // Copy the properties of the former leftchild to this node
                // so that it can be removed from the aggregate tree structure
                CopyParts(oldleftchild);

                // Now break the links to the tree structure in oldleftchild
                // in order to free it
                oldleftchild->m_leftchild   = NULL;
                oldleftchild->m_rightchild  = NULL;
                delete oldleftchild;

                m_parent = oldparent;

                if (m_leftchild!=NULL) {
                    m_rightchild->m_parent=this;
                    m_leftchild->m_parent=this;
                }

                delete oldleftparticle;

            }

            else
            {
                // Append to right subtree
                BinTreePrimary *oldrightparticle = m_rightparticle;
                for (size_t i=0; i != m_comp.size(); i++) {
                    m_leftparticle->m_comp[i] =
                            m_leftparticle->Composition(i) +
                            m_rightparticle->Composition(i);
                }

				//m_primarydiam of the new particle is the larger of the diameters of the merging particles
				//if the centre to centre separation isn't tracked this will be changed to the spherical equivalent in the call to UpdatePrimary
				m_leftparticle->m_primarydiam = max(m_leftparticle->m_primarydiam,m_rightparticle->m_primarydiam);

				//the new (merged) primary
				new_prim = m_leftparticle;

                m_leftparticle->UpdatePrimary();

                // All pointers to m_leftparticle now point to oldright particle
				if (!m_pmodel->getTrackPrimarySeparation()) {
				    oldrightparticle->ChangePointer(oldrightparticle,m_leftparticle);
					m_leftparticle->ChangePointer(m_leftparticle,m_leftparticle);
				}else{
					// If the priamry separation is tracked then re-estimation of new sintering levels also
					// requires the centre-centre separation and the smaller primary
					oldrightparticle->ChangePointer(oldrightparticle,m_leftparticle,m_distance_centreToCentre,small_prim);
					m_leftparticle->ChangePointer(m_leftparticle,m_leftparticle,m_distance_centreToCentre,small_prim);
				}

                // Set the pointer to the parent node
                if (oldrightparticle->m_parent->m_leftchild==oldrightparticle) {
                    oldrightparticle->m_parent->m_leftchild=m_leftchild;
                }
                else {
                    oldrightparticle->m_parent->m_rightchild=m_leftchild;
                }
                m_leftchild->m_parent=oldrightparticle->m_parent;

                BinTreePrimary *oldrightchild=m_rightchild;
                BinTreePrimary *oldparent=m_parent;

                // Copy the properties of the former leftchild to this node
                // so that it can be removed from the aggregate tree structure
                CopyParts(oldrightchild);

                // Now break the links to the tree structure in oldrightchild
                // in order to free it
                oldrightchild->m_leftchild = NULL;
                oldrightchild->m_rightchild = NULL;
                delete oldrightchild;

                m_parent=oldparent;

                if (m_leftchild!=NULL) {
                    m_rightchild->m_parent=this;
                    m_leftchild->m_parent=this;
                }

                delete oldrightparticle;

            }

        }

		//csl37***********
		//increase the diameter of the new primary to account for the volume of the primary being merged into it 
		//only if it is not a single primary
		//in the single primary case this is correctly updated in the call to updateprimary
		if(new_prim->m_parent != NULL){
			double sumterm = 0.0;
			double r_i = new_prim->m_primarydiam/2.0;
			//get contribution of neighours
			new_prim->SumNeighbours(new_prim,sumterm);
			//change in radius
			double delta_r_i = dV/(4*M_PI*r_i*r_i + M_PI*r_i*sumterm);
			//update the particle separations and calculate the new free surface of the particle
			double free_surface_term =0.0;
			new_prim->UpdateConnectivity(new_prim, delta_r_i, free_surface_term);
			//update primary diameter
			new_prim->m_primarydiam = 2.0*(r_i+delta_r_i);
			//update the free surface area
			new_prim->m_free_surf = max(M_PI*new_prim->m_primarydiam*new_prim->m_primarydiam - 2*M_PI*free_surface_term,0.0);
		}
		//csl37***********

    UpdateCache();

    }
    return *this;
}

/*!
 * @brief       Changes pointer from source to target
 *
 * The children surface is re-estimated using the spherical surface
 * and re-arranging the sintering level formula:
 * Cij = Ssph / (s * (1-2^{1/3}) + 2^{1/3})
 *
 * @param[in] source Pointer to the original particle
 * @param[in] target Pointer to the new particle
*/
void BinTreePrimary::ChangePointer(BinTreePrimary *source, BinTreePrimary *target)
{
	if(m_rightparticle == source) {
        m_rightparticle = target;
        double sphericalsurface =
            4 * PI * pow(3 * (m_leftparticle->Volume() +
                    m_rightparticle->Volume()) /(4 * PI), TWO_THIRDS);
        m_children_surf = sphericalsurface / (m_children_sintering *
                (1.0 - TWO_ONE_THIRD) + TWO_ONE_THIRD);
    }
    if(m_leftparticle == source){
        m_leftparticle = target;
        double sphericalsurface =
            4 * PI * pow(3 * (m_leftparticle->Volume() +
                    m_rightparticle->Volume()) /(4 * PI), TWO_THIRDS);
        m_children_surf = sphericalsurface / (m_children_sintering *
                (1.0 - TWO_ONE_THIRD) + TWO_ONE_THIRD);
    }

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->ChangePointer(source, target);
    }

}

//---------------------------------------------------------------csl37
//ChangePointer for when centre to centre separation is tracked
void BinTreePrimary::ChangePointer(BinTreePrimary *source, BinTreePrimary *target, double d_ij, BinTreePrimary *small_prim)
{
	if(m_rightparticle == source) {
        m_rightparticle = target;
		//if the neighbour is the smaller of the merging primaries then update the centre to centre distance
		if(source == small_prim){
			//estimate new separation as the smaller of the sum of the two separartions or the sum of primary radii
			m_distance_centreToCentre = min(m_distance_centreToCentre +d_ij, m_rightparticle->m_primarydiam/2.0 + m_leftparticle->m_primarydiam/2.0); 
		}
    }
    if(m_leftparticle == source){
        m_leftparticle = target;
		//if the neighbour is the smaller of the merging primaries then update the centre to centre distance
		if(source == small_prim){
			//estimate new separation as the smaller of the sum of the two separartions or the sum of primary radii
			m_distance_centreToCentre = min(m_distance_centreToCentre +d_ij, m_rightparticle->m_primarydiam/2.0 + m_leftparticle->m_primarydiam/2.0); 
		}
    }

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->ChangePointer(source, target,d_ij,small_prim);
    }

}

/*!
 * @brief       Sets all of the properties of the children to zero
 *
 * Used after the children are coalesced
 */
void BinTreePrimary::ResetChildrenProperties()
{
    m_children_radius   = 0.0;
    m_children_vol      = 0.0;
    m_children_surf     = 0.0;
    m_children_sintering= 0.0;
    m_avg_sinter        = 0.0;
    m_sint_rate         = 0.0;
	m_distance_centreToCentre = 0.0;
}

/*!
 * @brief       Updates the properties of a primary
 */
void BinTreePrimary::UpdatePrimary(void)
{
    // Call parent class UpdateCache
    Primary::UpdateCache();

    // Set specific features of this node.
	// if the primary centre-centre separation isn't tracked or this is a single primary 
	// m_primarydiam and m_diam are both the spherical equivalent diameter
	// and the free surface area is the sperhical surface area
	if(!m_pmodel->getTrackPrimarySeparation() || m_parent == NULL){
		m_primarydiam = m_diam;
		m_free_surf = m_surf;
	}
    m_numprimary  = 1;
}

//! UpdateCache helper function
void BinTreePrimary::UpdateCache(void)
{
    UpdateCache(this);
}

/*!
 * @brief       Updates the BinTreePrimary cache
 *
 * Updates the whole particle, from root to the lowest leaf-node.
 * Calculates the sintering level and merges particles (through Check
 * Sintering) if necessary. Only accurately calculates collision
 * diameter and other properties used outside of BinTreePrimary for
 * the root node.
 *
 * @param[in] root The root node of this particle
*/
void BinTreePrimary::UpdateCache(BinTreePrimary *root)
{
    // Update the children
    if (m_leftchild!=NULL) {
        m_leftchild->UpdateCache(root);
        m_rightchild->UpdateCache(root);
    }
    // This is a primary
    else {
        // If it's a single particle, give it a sintering level of 1.0.
        if (m_parent == NULL) m_avg_sinter = 1.0;
        else m_avg_sinter = 0.0;
        m_numprimary    = 1;
        UpdatePrimary();
    }

    // This is not a primary, sum up the properties
    if (m_leftchild!=NULL)
    {
        // Sum up the components first
        for (size_t i=0; i != m_comp.size(); i++) {
            m_comp[i] = m_leftchild->m_comp[i] + m_rightchild->m_comp[i];
        }

        // Now recalculate derived properties
        m_numprimary    = m_leftchild->m_numprimary +
                m_rightchild->m_numprimary;
        m_surf          = m_leftchild->m_surf + m_rightchild->m_surf;
        m_primarydiam   = m_leftchild->m_primarydiam
                + m_rightchild->m_primarydiam;
        m_vol           = m_leftchild->m_vol + m_rightchild->m_vol;
        m_mass          = m_leftchild->m_mass + m_rightchild->m_mass;
		m_free_surf		= m_leftchild->m_free_surf + m_rightchild->m_free_surf;

        // Calculate the sintering level of the two primaries connected by this node
        m_children_sintering = SinteringLevel();
        if (MergeCondition()) CheckSintering();

        // Sum up the avg sintering level (now that sintering is done)
        if((m_leftchild != NULL) && (m_rightchild != NULL)) {
            m_avg_sinter = m_children_sintering +
                    m_leftchild->m_avg_sinter + m_rightchild->m_avg_sinter;
        }
        else {
            // This should only occur if CheckSintering has merged
            m_avg_sinter = m_children_sintering;
        }

        // Calculate the different diameters only for the root node because
        // this is the only part of the tree seen by the other code, for
        // example, the coagulation kernel
        if (this == root) {
             // Get spherical equivalent radius and diameter
            double spherical_radius = pow(3 * m_vol / (4*PI), ONE_THIRD);
            m_diam = 2 * spherical_radius;

            // There are m_numprimary-1 connections between the primary
            // particles
            if (m_numprimary > 1)
                m_avg_sinter = m_avg_sinter / (m_numprimary - 1);

			if (!m_pmodel->getTrackPrimarySeparation()) {
				// Approxmiate the surface of the particle
				// (same as in ChangePointer)
				const double numprim_1_3 = pow(m_numprimary,-1.0 * ONE_THIRD);
				m_surf = 4 * PI * spherical_radius * spherical_radius /
						(m_avg_sinter * (1 - numprim_1_3) + numprim_1_3);
			}else{
				// if the centre to centre distance is tracked then this is the free surface area
				m_surf = m_free_surf;
			}

            // Calculate dcol based-on formula given in Lavvas et al. (2011)
            const double aggcolldiam = (6* m_vol / m_surf) *
                    pow(pow(m_surf, 3) / (36 * PI * m_vol * m_vol),
                            (1.0/m_pmodel->GetFractDim()));
            m_dmob = MobDiameter();
            m_dcol = aggcolldiam;

        }
        else {
            m_diam=0;
            m_dmob=0;
        }
    }
}

/*!
 *@brief        Overload of MobDiameter to return correct dmob
 *
 * This calculation is based on the work of Rogak et al., 1993 Aer. Sci.
 * Tech. 18:25-47, who give the calculation of dmob in the FM, SF and
 * transition regime.
 *
 * Currently, only the FM calculation is used, as the majority of experimental
 * systems which measure dmob run at continuum conditions.
 *
 * dmob,SF = 0.9 * dpri * sqrt(Df/(Df+2)) * npri^(1/Df)
 * dmob,FM = dpri * sqrt(0.802 * (npri - 1) + 1)
 *
 *@return   Mobility diameter of particle
 */
double BinTreePrimary::MobDiameter() const
{
    double dmob(1.0);

    // Is this a single particle?
    if (m_leftchild == NULL && m_parent == NULL) {
        dmob = m_diam;
    } else {
        // Presently hard-coded as T, P are abstracted from this
        // function's view
        // @TODO: somehow give this function access to Kn(T, P, d)

        if (false) {
            // SF regime mobility diameter
            dmob *= 0.9 * m_primarydiam / (double)m_numprimary;
            dmob *= sqrt(m_pmodel->GetFractDim() / (m_pmodel->GetFractDim() + 2));
            dmob *= pow(m_numprimary, (1.0/m_pmodel->GetFractDim()));
        } else {
            // FM regime mobility diameter
            dmob *= m_primarydiam / (double)m_numprimary;
            dmob *= sqrt(0.802*(m_numprimary-1) + 1);
        }

        if (dmob < m_diam) dmob = m_diam;

    }
    return dmob;
}


/*!
 * @brief       Prints a graphical output of the binary tree structure
 *
 * Graph outputted in GraphViz format. Use command 'dot' to convert.
 * e.g. dot -Tpng -o name.png name.dat
 * for PNG terminal, plotting output file name.dat
 *
 * @param[in] filename Output filename
*/
void BinTreePrimary::PrintTree(string filename) const
{
    std::ofstream out;
    out.open(filename.c_str());
    out << "digraph unix {"<<endl;
    out <<"graph [rankdir = \"LR\"];"<<endl;
    PrintTreeLoop(out);
    out << "}"<<endl;
    out.close();
}

/*!
 * @brief       Helper function for printing tree.
 *
 * @param[in] out Output stream
*/
void BinTreePrimary::PrintTreeLoop(std::ostream &out)const
{
    // Non leaf-node case
    if (m_leftchild!=NULL)
    {
        out << "\" " << this << "\" " << " [shape = \"record\" label = \"";
        this->PrintTreeNode(out);
        out <<"\"];"<<endl;

        out << "\" " << this->m_leftchild << "\" " <<
                " [shape = \"record\" label = \"";
        this->PrintTreeNode(out);
        out << "|" << this << "\"];"<<endl;

        out << "\" " << this->m_rightchild << "\" " <<
                " [shape = \"record\" label = \"";
        this->PrintTreeNode(out);
        out << "|" << this << "\"];"<<endl;

        out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftchild<<"\"; "<<endl;
        out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightchild<<"\"; "<<endl;
        out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftparticle<<
                "\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
        out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightparticle<<
                "\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
        m_leftchild->PrintTreeLoop(out);
        m_rightchild->PrintTreeLoop(out);
    }

    // Case when the node is a primary
    else
    {
        out << "\" " << this << "\" " <<
                " [shape = \"record\" color=\"blue\" label = \"";
        this->PrintTreeNode(out);
        out <<"\"];"<<endl;
    }
}

/*!
 * @brief       Prints out each node of the tree
 *
 * @param[in] out Output stream
*/
void BinTreePrimary::PrintTreeNode(std::ostream &out) const
{
    out
        << "|m_surf="             << this->m_surf
        << "|m_vol="             << this->m_vol
        << "|m_numprimary="      << this->m_numprimary
        << "|m_child_sint="      << this->m_children_sintering
        << "|m_child_rad="       << this->m_children_radius
        << "|m_child_surf="      << this->m_children_surf;
        for (size_t i=0; i != m_comp.size(); i++) {
            out << "|" + string(m_pmodel->Components(i)->Name()) + "=" <<
                    m_comp[i];
        }
    out
        << "|m_parent="          << this->m_parent
        << "|this=" << this;
}

void BinTreePrimary::PrintComponents() const
{
    for (size_t i=0; i != m_comp.size(); i++) {
        cout << m_pmodel->Components(i)->Name() << " " << m_comp[i] << " ";
    }
    cout << endl;
}

/*!
 * @brief       Adjusts the particle after a surface rxn event
 *
 * Analogous to the implementation in Primary. The function will
 * however descend the tree to find a primary adjust. The surface area
 * and volume added to the particle is also rippled through the binary
 * tree.
 *
 * @param[in]   dcomp   Vector storing changes in particle composition
 * @param[in]   dvalues Vector storing changes in gas-phase comp
 * @param[in]   rng     Random number generator
 * @param[in]   n       Number of times for adjustment
 */
unsigned int BinTreePrimary::Adjust(const fvector &dcomp,
        const fvector &dvalues, rng_type &rng, unsigned int n)

{
    if (m_leftchild == NULL && m_rightchild == NULL) {
        		
		double dV(0.0);
        double volOld = m_vol;
		double m_diam_old = m_diam;

		//! Initialisation of variables to adjust the primary diameter if the
        //! distance between the centres of primary particles is tracked.
		double m_primary_diam_old = m_primarydiam;
        double r_i = 0.0;                //!< Radius of primary particle i.
		double dr_max = 0.0;			//maximum change in primary radius during internal step
		double sumterm = 0.0;
		double delta_r_i = 0.0;
		double free_surface_term = 0.0;		//sumterm for calculating free surface area

        // Call to Primary to adjust the state space
        n = Primary::Adjust(dcomp, dvalues, rng, n);

        // Stop doing the adjustment if n is 0.
        if (n > 0) {
            // Update only the primary
            UpdatePrimary();

			//! If the distance between the centres of primary particles is
			//! tracked, the rate of change in the primary diameter is affected
			//! by its neighbours.
			if (m_pmodel->getTrackPrimarySeparation()) {
				
				// particle with more than primary
				if (m_parent != NULL) {	
				
				//***************************csl37: new surface growth model
					
					while (volOld <= m_vol){

						r_i = m_primarydiam / 2.0;
						dr_max = 0.1*r_i;					

						//get contribution from neighbours working up the binary tree
						SumNeighbours(this, sumterm);
					
						//calculate change in volume
						dV = dr_max * (4*M_PI*r_i*r_i + M_PI*r_i*sumterm) ;
						assert(dV >= 0.0);	//csl37:debug
						

						if (volOld + dV > m_vol){
							delta_r_i = (m_vol - volOld)*dr_max / dV;
						}else{
							delta_r_i = dr_max;
						}
					
						//update the particle separations and calculate the new free surface of the particle
						if(true){
							UpdateConnectivity(this, delta_r_i, free_surface_term);
						}

						//update primary diameter
						m_primarydiam = 2.0* (r_i + delta_r_i);
						//update the free surface area
						//if the calculated area is negative (too many overlaps) then set m_free_surf = 0.0
						m_free_surf = max(M_PI*m_primarydiam*m_primarydiam - 2*M_PI*free_surface_term, 0.0);

						volOld += dV;

					}

				//single primary case: the primary diameter equals the spherical diameter
				}else{		
					m_primarydiam = m_diam;
				}

			} else {

				dV = m_vol - volOld;
				double dS(0.0);

				if (dV > 0.0) {
					// Surface change due to volume addition
					dS = dV * 2.0 * m_pmodel->GetBinTreeCoalThresh() / m_diam;
				}
				// TODO: Implement surface area reduction?

				// Climb back-up the tree and update the surface area and
				// sintering of a particle
				UpdateParents(dS);
			}
        }
    }

    // Else this a non-leaf node (not a primary)
    else
    {
		
        // Generate random numbers
        boost::bernoulli_distribution<> leftRightChooser;
        // Select particle
        if(leftRightChooser(rng))
            return m_leftparticle->Adjust(dcomp, dvalues, rng, n);
        else
            return m_rightparticle->Adjust(dcomp, dvalues, rng, n);

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
		//csl37

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

//csl37
//function to update particle cache from root node
void BinTreePrimary::UpdateCacheRoot(void){
	if(m_parent != NULL){
		m_parent->UpdateCacheRoot();
	}else{
		UpdateCache(this);
	}
}

//csl37
//function to identify neighbours and sum their contribution to surface growth
//identifies neighbours of primary being adjusted and compute the summation term
//prim: neighbour being adjusted
//sumterm: summation term
void BinTreePrimary::SumNeighbours(BinTreePrimary *prim, double &sumterm) {
	
	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;

	//check if a neighbour of prim
	if (m_parent->m_leftparticle == prim) {
		//right particle is a neighbour
		r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
	} else if(m_parent->m_rightparticle == prim) {
		//left particle is a neighbour
		r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
	} else {
		//not a neighbour
		r_j = 0.0;
	}
	
	//if the node connects a neighbourt then calculate the summation term
	if(r_j > 0.0){
		//the volumes and radii of neighbours remain unchanged
		//the centre to centre separations increase to allow growth of the primary
		x_ij = ( pow(d_ij,2.0) - pow(r_j,2.0) + pow(r_i,2.0) ) / ( 2.0*d_ij );
		sumterm += x_ij - 2.0*r_i + pow(r_i,2.0) / x_ij; 
	}

	//continue working up the binary tree
	if(m_parent->m_parent != NULL){
		m_parent->SumNeighbours(prim, sumterm);
	}
}

//csl37
void BinTreePrimary::UpdateConnectivity(BinTreePrimary *prim, double delta_r, double &sumterm){
	
	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;

	//check if a neighbour of prim
	if (m_parent->m_leftparticle == prim) {
		//right particle is a neighbour
		r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
	} else if(m_parent->m_rightparticle == prim) {
		//left particle is a neighbour
		r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
	} else {
		//not a neighbour
		r_j = 0.0;
	}

	if(r_j > 0.0){

		x_ij = ( pow(d_ij,2.0) - pow(r_j,2.0) + pow(r_i,2.0) ) / ( 2.0*d_ij );
		//update centre to centre separation
		//making sure centre to centre separation remains smaller than some of radii
		m_parent->m_distance_centreToCentre = min(d_ij + r_i * delta_r / x_ij, r_i+r_j);

		//calculate term for free surface area 
		d_ij = m_parent->m_distance_centreToCentre;
		x_ij = ( pow(d_ij,2.0) - pow(r_j,2.0) + pow(r_i+delta_r,2.0) ) / ( 2.0*d_ij );
		sumterm += (r_i+delta_r)*(r_i+delta_r) - (r_i+delta_r)*x_ij;
	}

	//work up binary tree
	if(m_parent->m_parent != NULL){
		m_parent->UpdateConnectivity(prim, delta_r, sumterm);
	}
}

//csl37
//overload of function ignore update to neighbour
void BinTreePrimary::UpdateConnectivity(BinTreePrimary *prim, double delta_r, double &sumterm, BinTreePrimary *prim_ignore){
	
	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;

	//check if a neighbour of prim
	if (m_parent->m_leftparticle == prim && m_parent->m_rightparticle != prim_ignore ) {
		//right particle is a neighbour
		r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
	} else if(m_parent->m_rightparticle == prim &&  m_parent->m_leftparticle != prim_ignore ) {
		//left particle is a neighbour
		r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
	} else {
		//not a neighbour
		r_j = 0.0;
	}

	if(r_j > 0.0){

		x_ij = ( pow(d_ij,2.0) - pow(r_j,2.0) + pow(r_i,2.0) ) / ( 2.0*d_ij );
		//update centre to centre separation
		//making sure centre to centre separation remains smaller than some of radii
		m_parent->m_distance_centreToCentre = min(d_ij + r_i * delta_r / x_ij, r_i+r_j);

		//calculate term for free surface area 
		d_ij = m_parent->m_distance_centreToCentre;
		x_ij = ( pow(d_ij,2.0) - pow(r_j,2.0) + pow(r_i+delta_r,2.0) ) / ( 2.0*d_ij );
		sumterm += (r_i+delta_r)*(r_i+delta_r) - (r_i+delta_r)*x_ij;
	}

	//work up binary tree
	if(m_parent->m_parent != NULL){
		m_parent->UpdateConnectivity(prim, delta_r, sumterm, prim_ignore);
	}
}


/*!
 * @brief       Updates the surface area and sintering level of all parents
 *
 * @param[in]   dS      Surface area increment to adjust area by
 */
void BinTreePrimary::UpdateParents(double dS) {
    if (m_parent != NULL) {
        m_parent->m_children_surf += dS;
        m_parent->m_children_sintering = m_parent->SinteringLevel();
        m_parent->UpdateCache();
    }
}

/*!
 * @brief       Adjusts the particle after an IntP event
 *
 * @param[in]   dcomp   Vector storing changes in particle composition
 * @param[in]   dvalues Vector storing changes in gas-phase comp
 * @param[in]   rng     Random number generator
 * @param[in]   n       Number of times for adjustment
 */
unsigned int BinTreePrimary::AdjustIntPar(const fvector &dcomp,
        const fvector &dvalues, rng_type &rng, unsigned int n)

{
    if(m_numprimary == 1) {
        // Call to Primary to adjust the state space
        n = Primary::AdjustIntPar(dcomp, dvalues, rng, n);
    } else {
        // Generate random numbers
        boost::bernoulli_distribution<> leftRightChooser;

        // Select particle
        if(leftRightChooser(rng))
            return m_leftparticle->AdjustIntPar(dcomp, dvalues, rng, n);
        else
            return m_rightparticle->AdjustIntPar(dcomp, dvalues, rng, n);
    }
    // Update property cache.
    UpdateCache();
    return n;

}

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
void BinTreePrimary::Sinter(double dt, Cell &sys,
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
    if (m_leftparticle!=NULL && m_rightparticle!=NULL) {

        SinterNode(dt, sys, model, rng, wt);
		
		// Check if the sintering level is above the threshold, and merge
		if (MergeCondition()) {
            CheckSintering();
		}

		if (m_leftchild != NULL && m_rightchild != NULL) {
            m_leftchild->Sinter(dt, sys, model, rng, wt);
            m_rightchild->Sinter(dt, sys, model, rng, wt);
        }

		UpdateCache();

		m_children_sintering = SinteringLevel();
    }

}

//! Returns the sintering rate
double BinTreePrimary::GetSintRate() const
{
    double sint_rate = m_sint_rate;
    if(m_leftchild!=NULL)
    {
        sint_rate += (m_leftchild->GetSintRate() +
                m_rightchild->GetSintRate());
    }

    return sint_rate;
}

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
void BinTreePrimary::SinterNode(
        double dt,
        Cell &sys,
        const Processes::SinteringModel &model,
        rng_type &rng,
        double wt
        ) {

	// Declare time step variables.
	double t1=0.0, delt=0.0, tstop=dt;
	double r=0.0;

	// The scale parameter discretises the delta-S when using
	// the Poisson distribution.  This allows a smoother change
	// (smaller scale = higher precision).
	double scale = 0.01;

	//! The sintering model depends on whether the distance between the
	//! centres of primary particles is tracked. If tracked, sintering
	//! results in a decrease in the distance between the primaries and an
	//! increase in their diameters. If not, sintering results in a
	//! decrease in the common surface between the primaries.
	if (!m_pmodel->getTrackPrimarySeparation()) {

		// Calculate the spherical surface
		const double spherical_surface=4*PI*m_children_radius*m_children_radius;	

		// Define the maximum allowed change in surface
		// area in one internal time step (10% spherical surface).
		double dAmax = 0.1 * spherical_surface;

		// Perform integration loop.
		while (t1 < tstop)
		{
			// Calculate sintering rate.
			r = model.Rate(m_time+t1, sys, *this);

			if (r > 0) {
				// Calculate next time-step end point so that the
				// surface area changes by no more than dAmax.
				delt = dAmax / max(r, 1.0e-300);

				// Approximate sintering by a poisson process.  Calculate
				// number of poisson events.
				double mean;

				if (tstop > (t1+delt)) {
					// A sub-step, we have changed surface by dAmax, on average
					mean = 1.0 / scale;
				} else {
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

	} else {
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
		
			assert(d_ij<= r_i+r_j);	//csl37:debug

			//! Definition of variables for conciseness.
            double d_ij2 = pow(d_ij, 2.0); 
            double r_i2 = pow(r_i, 2.0);
            double r_j2 = pow(r_j, 2.0);
            double r_i4 = pow(r_i, 4.0);
            double r_j4 = pow(r_j, 4.0);
			
			//! Continue if primaries have not coalesced
            if (!MergeCondition()) {

				//double x_i = (d_ij2 - r_j2 + r_i2) / (2.0 * d_ij);
				//double x_j = (d_ij2 - r_i2 + r_j2) / (2.0 * d_ij);
				//due to rounding, x_i and x_j sometimes calculated to be larger than respective primary radii resulting in negative neck area
				//therefore we take the smaller of x_i and r_i
				double x_i = min((d_ij2 - r_j2 + r_i2) / (2.0 * d_ij),r_i); //!< Eq. (3b).
				double x_j = min((d_ij2 - r_i2 + r_j2) / (2.0 * d_ij),r_j); //!< Eq. (3b).
				double A_n = M_PI * (r_i2 - pow(x_i, 2.0));        //!< Eq. (4).
				double A_i = 2.0 * M_PI * (r_i2 + r_i * x_i);      //!< Eq. (6).
				double A_j = 2.0 * M_PI * (r_j2 + r_j * x_j);      //!< Eq. (6).

				//declare variables
				double dd_ij_dt=0.0;
				double R_n = 0.0;
				double r4_tau = 0.0;
				double tau = 0.0;
				double gamma_eta = 0.0;

				//Sintering model dependent part
				Processes::SinteringModel::SintType sint_model_type = model.Type();
				switch (sint_model_type){
					//viscous flow model
					case Processes::SinteringModel::SintType::ViscousFlow:
						
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
						gamma_eta = min(r_i, r_j) / tau;
	
						///////////////////////////////////////////////////////
						/// References to equations in Langmuir 27:6358 (2011).
						///////////////////////////////////////////////////////

						//! Eq. (14a).
						dd_ij_dt = 4.0 * r_i * r_j * d_ij2 * (r_i + r_j) * gamma_eta /
										((r_i + r_j + d_ij) * (r_i4  + r_j4 - 2.0 * r_i2 * r_j2 + 4.0 * d_ij * r_i * r_j *(r_i + r_j) - d_ij2 * (r_i2 + r_j2)));
						//****************************************************************************************csl37
						break;

					//grain boundary diffusion model 
					case Processes::SinteringModel::SintType::GBD:

						//if the particles are in point contact set an initial neck radius of 1% of the smaller primary radius
						//otherwise the dd_ij_dt is undefined
						//csl37
						if(A_n == 0.0){
							R_n = 0.01*min(r_i,r_j);
							A_n = M_PI * R_n * R_n;
							x_i = sqrt(r_i2 - R_n * R_n);
							x_j = sqrt(r_j2 - R_n * R_n);
						}else{
							R_n = sqrt(A_n / M_PI);
						}

						// Note: the primary radius in the numerator cancel with the radius dependence of tau
						// so we can calculate this for only one primary 
						r4_tau = r_i4 / model.SintTime(sys, *this->m_leftparticle);

						//J Aerosol Sci 46:7-19 (2012) Eq. (A6)
						//dx_i_dt + dx_j_dt
						// Note: this is missing a minus sign, which is accounted for below
						dd_ij_dt = r4_tau * ( 1/(r_i - x_i) + 1/(r_j - x_j) - 2/R_n ) / A_n;

						break;
						//****************************************************************************************csl37
					default:
						std::cout<<"Sintering model not coded"<<endl;
						break;
				}

				//csl37******************
				//primaries with multiple neighbours
				if(true){
					double sumterm_i = 0.0;
					double sumterm_j = 0.0;
					//get contribution from neighbours working up the binary tree
					m_leftparticle->SumNeighbours(this, sumterm_i);
					m_rightparticle->SumNeighbours(this, sumterm_j);
					//subtract mutual contribution from sumterm
					sumterm_i = max( sumterm_i - (x_i - r_i)*(x_i - r_i)/x_i , 0.0 );
					sumterm_j = max( sumterm_j - (x_j - r_j)*(x_j - r_j)/x_j , 0.0 );

					//modified A_i and A_j
					A_i = M_PI * (2*r_i*r_i + 2*r_i*x_i + r_i*sumterm_i);
					A_j = M_PI * (2*r_j*r_j + 2*r_j*x_j + r_j*sumterm_j);
				}
				//csl37******************

				//! The expression for B_i in Eq. (8) is wrong. By combining
				//! Eqs. (5) and (7), we can obtain two equations which are
				//! functions of r_i and r_j. Subsequently combined these two
				//! equations and used Wolfram Alpha to rearrange equation in terms
				//! of r_i (and r_j).
				//!
				//! @todo Remove derivation and replace with reference to preprint
				//!       or paper if results do get published.
				double B_i = (-r_j*A_n*A_n - x_j*A_j*A_n)/(A_i*A_j*d_ij + r_i*A_j*A_n + r_j*A_i*A_n);
				double B_j = (-r_i*A_n*A_n - x_i*A_i*A_n)/(A_j*A_i*d_ij + r_j*A_i*A_n + r_i*A_j*A_n);
				
				double V_i = 2.0 / 3.0 * M_PI * pow(r_i, 3.0) + M_PI * r_i2 * x_i - 1.0 / 3.0 * M_PI * pow(x_i, 3.0); //!< Eq. (3a).
				double V_j = 2.0 / 3.0 * M_PI * pow(r_j, 3.0) + M_PI * r_j2 * x_j - 1.0 / 3.0 * M_PI * pow(x_j, 3.0); //!< Eq. (3a).

				delt = dd_ij_Max / max(dd_ij_dt, 1.0e-300);
				double mean;

				if (tstop > (t1 + delt)) {
					mean = 1.0 / scale;
				} else {
					mean = dd_ij_dt * (tstop - t1) / (scale * dd_ij_Max);
				}

				//! Sinter primaries
                boost::random::poisson_distribution<unsigned, double> repeatDistribution(mean);
                const unsigned n = repeatDistribution(rng);
                m_distance_centreToCentre -= (double)n * scale * dd_ij_Max; //!< Sintering decreases d_ij hence the negative sign.
				
				//! change in primary radii
				//! The factor of 2 is because Eq. (8) is the rate of change
				//! in radius.
				double delta_r_i = - (double)n * scale * B_i * dd_ij_Max;  //!< Eq. (8).
				double delta_r_j = - (double)n * scale * B_j * dd_ij_Max; //!< Eq. (8).
				

				//****************csl37
				//Adjust separation of neighbours due to increase in radius
				if(true){
				if (!MergeCondition()) {
					double free_surf_i=0.0;
					double free_surf_j=0.0;

					//adjust separation with neighbours (ignoring p_j) and return the new free surface
					m_leftparticle->UpdateConnectivity(m_leftparticle, delta_r_i, free_surf_i, m_rightparticle);
					//add i,j term
					x_i = ( pow(m_distance_centreToCentre,2.0) - pow(r_j+delta_r_j,2.0) + pow(r_i+delta_r_i,2.0) ) / ( 2.0*m_distance_centreToCentre );
					free_surf_i += (r_i+delta_r_i)*(r_i+delta_r_i) - (r_i+delta_r_i)*x_i;

					//adjust separation with neighbours (ignoring p_i) and return the new free surface
					m_rightparticle->UpdateConnectivity(m_rightparticle, delta_r_j, free_surf_j, m_leftparticle);					
					//add i,j term
					x_j = ( pow(m_distance_centreToCentre,2.0) - pow(r_i+delta_r_i,2.0) + pow(r_j+delta_r_j,2.0) ) / ( 2.0*m_distance_centreToCentre );
					free_surf_j += (r_j+delta_r_j)*(r_j+delta_r_j) - (r_j+delta_r_j)*x_j;

					//update the free surface area 
					m_leftparticle->m_free_surf = max(M_PI*m_leftparticle->m_primarydiam*m_leftparticle->m_primarydiam - 2*M_PI*free_surf_i, 0.0);
					//update the free surface area
					m_rightparticle->m_free_surf = max(M_PI*m_rightparticle->m_primarydiam*m_rightparticle->m_primarydiam - 2*M_PI*free_surf_j, 0.0);
				}
				}
				//****************csl37

				//!Adjust primary radii
				this->m_leftparticle->m_primarydiam += 2.0 * delta_r_i;				
				this->m_rightparticle->m_primarydiam += 2.0 * delta_r_j;

				t1 += delt;

				//return some sintering rate
				r = dd_ij_dt;
            } else {
                break; //!do not continue to sinter.
            }
        }

	}

    m_children_sintering = SinteringLevel();

    m_sint_rate = r;

}


/*!
 * Get the number of units of a component
 *
 * @param name  Name of component
 * @return      Number of units of a component
 */
double BinTreePrimary::GetComponent(std::string name) const
{
	return m_comp[m_pmodel->ComponentIndex(name)];
}

/*!
 * Set the number of units of a component
 *
 * @param name  Name of component
 * @param val   Value to set
 */
void BinTreePrimary::SetComponent(std::string name, double val)
{
    try {
        m_comp[m_pmodel->ComponentIndex(name)] = val;
    } catch(std::exception& e) {
        throw e.what();
    }
}

/*!
 * @return  Arithmetic standard dev. of primary diameter
 */
double BinTreePrimary::GetPrimaryAStdDev() const
{
    // Get a list of all primary diameters first
    fvector diams;
    GetAllPrimaryDiameters(diams);
    double dpri = GetPrimaryDiam() / ((double) GetNumPrimary());

    // If binary tree writing hasn't been enabled, an erroneous
    // value will be returned. So just give zero.
    if (!m_pmodel->WriteBinaryTrees())
        return 0.0;

    // Loop over diameters to get stdev
    double stdev(0.0), dev(0.0);
    for (size_t i = 0; i != diams.size(); i++) {
        dev = (diams[i] - dpri);
        stdev += dev * dev;
    }

    stdev = sqrt(stdev / diams.size());
    return stdev;
}

/*!
 * @return  Geometric mean primary diameter
 */
double BinTreePrimary::GetPrimaryGMean() const
{
    // Get a list of all primary diameters first
    fvector diams;
    GetAllPrimaryDiameters(diams);

    // If binary tree writing hasn't been enabled, an erroneous
    // value will be returned. So just give zero.
    if (!m_pmodel->WriteBinaryTrees())
        return 0.0;

    // Calculate the geometric mean diameter
    double dpri(1.0);
    double inv_n = 1.0 / (double) diams.size();
    for (size_t i = 0; i != diams.size(); i++) {
        dpri *= pow(diams[i], inv_n);
    }

    return dpri;
}

/*!
 *
 * @return  Geometric stdev of primary diameter
 */
double BinTreePrimary::GetPrimaryGStdDev() const
{
    // Get a list of all primary diameters first
    fvector diams;
    GetAllPrimaryDiameters(diams);
    double dpri = GetPrimaryGMean();

    // If binary tree writing hasn't been enabled, an erroneous
    // value will be returned. So just give zero.
    if (!m_pmodel->WriteBinaryTrees())
        return 0.0;

    // Loop over diameters to get stdev
    double stdev(0.0), dev(0.0);
    for (size_t i = 0; i != diams.size(); i++) {
        dev = log(diams[i] / dpri); // (natural log)
        stdev += dev * dev;
    }

    stdev = exp(sqrt(stdev / diams.size()));
    return stdev;
}

/*!
 * Loop through the particles and collect a list of primary
 * particle diameters.
 *
 * @param diams     Vector of diameters
 */
void BinTreePrimary::GetAllPrimaryDiameters(fvector &diams) const
{
    // Only add diameter to list if it's a primary
    if (m_leftchild == NULL && m_rightchild == NULL)
        //diams.push_back(m_diam);
		diams.push_back(m_primarydiam);		//csl37: use m_primarydiam for consistency

    if (m_leftchild != NULL && m_rightchild != NULL) {
        m_leftchild->GetAllPrimaryDiameters(diams);
        m_rightchild->GetAllPrimaryDiameters(diams);
    }
}

/*!
 * @brief Writes a particle to a binary stream
 *
 * @param[in,out]    out                 Output binary stream
 *
 * @exception        invalid_argument    Stream not ready
 */
void BinTreePrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        if (m_pmodel->WriteBinaryTrees()) {
            // Call the binary tree serialiser...
            BinTreeSerializer <BinTreePrimary> tree;
            tree.Serialize(out, this, NULL);
        } else {
            // Just serialise the root node.
            SerializePrimary(out, NULL);
        }

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, BinTreePrimary::Serialize).");
    }
}

/*!
 *  @brief Writes an individual primary to a binary stream.
 *
 *  @param[in,out] out  Output binary stream.
 *  @param         void It is a pointer, but the type that it points to is not known.
 *
 *  @exception invalid_argument Stream not ready.
 */
void BinTreePrimary::SerializePrimary(std::ostream &out, void*) const
{
    if (out.good()) {

        int  val_int(0);
        double val(0.0);

        val_int = m_numprimary;
        out.write((char*)&val_int, sizeof(val_int));

        val = m_primarydiam;
        out.write((char*)&val, sizeof(val));

        val = m_children_radius;
        out.write((char*)&val, sizeof(val));

        val = m_children_vol;
        out.write((char*)&val, sizeof(val));

        val = m_children_surf;
        out.write((char*)&val, sizeof(val));

		val = m_free_surf;
        out.write((char*)&val, sizeof(val));

        val = m_distance_centreToCentre;
        out.write((char*)&val, sizeof(val));

        val = m_children_sintering;
        out.write((char*)&val, sizeof(val));

        val = m_avg_sinter;
        out.write((char*)&val, sizeof(val));

        val = m_sint_rate;
        out.write((char*)&val, sizeof(val));

        val = m_sint_time;
        out.write((char*)&val, sizeof(val));

        // Output base class.
        Primary::Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, BinTreePrimary::SerializePrimary).");
    }
}



/*!
 * @brief Deserialise the binary tree
 *
 * Will only deserialise the full binary tree if the particle model has
 * writing/reading of trees activated. Otherwise, it will just load the
 * root node data.
 *
 * @param[in,out]    in                  Input binary stream
 * @param[in]        model               Particle model defining interpretation of particle data
 *
 * @exception        invalid_argument    Stream not ready
 */
void BinTreePrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    //UpdateCache();
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        if (model.WriteBinaryTrees()) {
            // Call the binary tree serialiser...
            BinTreeSerializer <BinTreePrimary> tree;
            tree.Deserialize(in, this, model, NULL);
        } else {
            // Just deserialise the root node.
            DeserializePrimary(in, model, NULL);
        }


    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, BinTreePrimary::Deserialize).");
    }
}

/*!
 *  @brief Deserialise attributes of a single particle node.
 *
 *  @param[in,out] in    Input binary stream.
 *  @param[in]     model Particle model defining interpretation of particle data.
 *  @param         void  It is a pointer, but the type that it points to is not known.
 *
 *  @exception invalid_argument Stream not ready.
 */
void BinTreePrimary::DeserializePrimary(std::istream &in,
        const Sweep::ParticleModel &model,
        void*)
{
    if (in.good()) {

        int  val_int(0);
        double val(0.0);

        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        m_numprimary = val_int;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_primarydiam = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_children_radius = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_children_vol = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_children_surf = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_free_surf = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_distance_centreToCentre = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_children_sintering = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_avg_sinter = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sint_rate = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sint_time = val;

        // Input base class.
        Primary::Deserialize(in, model);

    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, BinTreePrimary::DeserializePrimary).");
    }
}
