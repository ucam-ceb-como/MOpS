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
 *      Implementation of the BintreePrimary class.
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
BintreePrimary::BintreePrimary() : Primary(),
    m_numprimary(0),
    m_primarydiam(0.0),
    m_children_radius(0.0),
    m_children_vol(0.0),
    m_children_surf(0.0),
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
BintreePrimary::BintreePrimary(const real time,
        const Sweep::ParticleModel &model)
: Primary(time, model),
    m_numprimary(0),
    m_primarydiam(0.0),
    m_children_radius(0.0),
    m_children_vol(0.0),
    m_children_surf(0.0),
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
BintreePrimary::BintreePrimary(const BintreePrimary &copy)
{
    *this = copy;
    if (copy.m_leftchild!=NULL)
    {
        CopyTree(&copy);
    }
}

//! Stream-reading constructor.
BintreePrimary::BintreePrimary(std::istream &in,
        const Sweep::ParticleModel &model) :
m_numprimary(0),
m_primarydiam(0.0),
m_children_radius(0.0),
m_children_vol(0.0),
m_children_surf(0.0),
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
BintreePrimary::~BintreePrimary()
{
    delete m_leftchild;
    delete m_rightchild;
}


//! Operator definition
BintreePrimary &BintreePrimary::operator=(const Primary &rhs)
{
    return *this;
}

//! Return a copy of the model data
BintreePrimary *const BintreePrimary::Clone(void) const
{
   return new BintreePrimary(*this);
}


/*!
 * @brief       Recursively copy the tree for non-leaf nodes.
 *
 * @param[in] source Pointer to the primary to be copied
*/
void BintreePrimary::CopyTree(const BintreePrimary *source)
{
    // Create the new left and right children with nothing in them
    m_leftchild     = new BintreePrimary(source->CreateTime(),*m_pmodel);
    m_rightchild    = new BintreePrimary(source->CreateTime(),*m_pmodel);

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
void BintreePrimary::UpdateAllPointers(const BintreePrimary *original)
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
std::stack<bool> BintreePrimary::recordPath(const BintreePrimary* bottom,
                                        const BintreePrimary* const top)
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
BintreePrimary* BintreePrimary::descendPath(BintreePrimary *here,
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
BintreePrimary &BintreePrimary::Coagulate(const Primary &rhs, rng_type &rng)
{
    const BintreePrimary *rhsparticle = NULL;
    rhsparticle = dynamic_cast<const AggModels::BintreePrimary*>(&rhs);

    // Don't sum-up components now. Wait for UpdateCache so as to not
    // double-count the values.

    // Create the new particles
    BintreePrimary *newleft = new BintreePrimary(m_time, *m_pmodel);
    BintreePrimary *newright = new BintreePrimary(m_time, *m_pmodel);
    BintreePrimary copy_rhs(*rhsparticle);


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
 * @brief       Copy state-space and derived properties from source
 *
 * This function is like a limited assignment operator, except that the
 * children are not copied and the pointers to the particles may need
 * adjusting after this method has finished.
 *
 * @param[in] source Pointer to the primary to be copied
 */
void BintreePrimary::CopyParts(const BintreePrimary *source)
{
    // Set primary characteristics
    SetComposition(source->Composition());
    SetValues(source->Values());
    SetTime(source->LastUpdateTime());
    SetCollDiameter(source->CollDiameter());
    SetMobDiameter(source->MobDiameter());
    SetSphDiameter(source->SphDiameter());
    SetSurfaceArea(source->SurfaceArea());
    SetVolume(source->Volume());
    SetMass(source->Mass());

    // Set BintreePrimary model characteristics
    m_numprimary        = source->m_numprimary;
    m_primarydiam       = source->m_primarydiam;
    m_children_radius   = source->m_children_radius;
    m_children_vol      = source->m_children_vol;
    m_children_surf     = source->m_children_surf;
    m_children_sintering= source->m_children_sintering;
    m_avg_sinter        = source->m_avg_sinter;
    m_sint_rate         = source->m_sint_rate;
    m_sint_time         = source->m_sint_time;

    // Set particles
    m_leftchild         = source->m_leftchild;
    m_rightchild        = source->m_rightchild;
    m_parent            = source->m_parent;
    m_leftparticle      = source->m_leftparticle;
    m_rightparticle     = source->m_rightparticle;
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
void BintreePrimary::SetSinteringTime(real time) {
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
 * @brief       Randomly selects a primary in the binary tree
 *
 * Select a primary uniformly at random from this particle
 * and descend the aggregate tree to find the primary.
 * Note that most BintreePrimaries are nodes in a tree representing
 * connectivity within an aggregate, so it is necessary to
 * descend the tree to find a primary that really is a
 * primary.
 *
 * @param[in]   rng     Random number generator
 *
 * @return      Pointer to an object representing a physical primary
 */
BintreePrimary *BintreePrimary::SelectRandomSubparticle(rng_type &rng)
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
BintreePrimary *BintreePrimary::SelectRandomSubparticleLoop(int target)
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
real BintreePrimary::SinteringLevel()
{
    if (m_leftchild != NULL && m_rightchild != NULL) {
        // Calculate the spherical surface
        const real spherical_surface =
                4 * PI * m_children_radius * m_children_radius;
        real slevel(0.0);

        if (m_children_surf <= spherical_surface) {
            m_children_surf = spherical_surface;
            return 1.0;
        }

        if (m_children_surf == 0.0) {
            slevel = 0.0;
        } else {
            slevel= ((spherical_surface/m_children_surf) - TWO_ONE_THIRD)
                    /(1 - TWO_ONE_THIRD);
        }

        if (slevel < 0.0) {
            return 0.0;
        } else if (slevel > 1.0) {
            cout << "sweep: BintreePrimary::SinteringLevel() s.l. > 1.0";
            return 1.0;
        } else {
            return slevel;
        }
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
bool BintreePrimary::CheckSintering()
{
    bool hassintered=false;
    if (m_children_sintering > 0.95 && m_leftparticle != NULL) {
         Merge();
         UpdateCache();
         hassintered=true;

         // Check again because this node has changed
         CheckSintering();
    }
    if (m_leftchild!=NULL) {
        hassintered = m_leftchild->CheckSintering();
        hassintered = m_rightchild->CheckSintering();
    }

    return hassintered;
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
BintreePrimary &BintreePrimary::Merge()
{

    // Make sure this primary has children to merge
    if( m_leftchild!=NULL) {
        if (m_leftchild==m_leftparticle && m_rightchild==m_rightparticle)
        {
            // This node has only two primaries in its subtree, it is possible
            // that this node is not the root node and belongs to a bigger
            // particle.

            // Sum up the components first
            for (unsigned int i=0; i != m_comp.size(); i++) {
                m_comp[i] = m_leftparticle->Composition(i) +
                        m_rightparticle->Composition(i);
            }

            // Update the pointers that pointed to the two former children
            ChangePointer(m_leftchild,this);
            ChangePointer(m_rightchild,this);


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
                if (m_parent->SinteringLevel() <= 0.95) {
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
                BintreePrimary *oldleftparticle = m_leftparticle;
                for (unsigned int i=0; i != m_comp.size(); i++) {
                    m_rightparticle->m_comp[i] =
                            m_leftparticle->Composition(i) +
                            m_rightparticle->Composition(i);
                }

                m_rightparticle->UpdatePrimary();

                // Set the pointers from the leftprimary to the rightprimary
                // (this will be the new bigger primary)
                oldleftparticle->ChangePointer(oldleftparticle,m_rightparticle);
                m_rightparticle->ChangePointer(m_rightparticle,m_rightparticle);

                // Set the pointer to the parent node
                if (oldleftparticle->m_parent->m_leftchild==oldleftparticle) {
                    oldleftparticle->m_parent->m_leftchild=m_rightchild;
                }
                else {
                    oldleftparticle->m_parent->m_rightchild=m_rightchild;
                }
                m_rightchild->m_parent=oldleftparticle->m_parent;

                BintreePrimary *oldleftchild    = m_leftchild;
                BintreePrimary *oldparent       = m_parent;

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
                BintreePrimary *oldrightparticle = m_rightparticle;
                for (unsigned int i=0; i != m_comp.size(); i++) {
                    m_leftparticle->m_comp[i] =
                            m_leftparticle->Composition(i) +
                            m_rightparticle->Composition(i);
                }

                m_leftparticle->UpdatePrimary();
                // All pointers to m_leftparticle now point to oldright particle
                oldrightparticle->ChangePointer(oldrightparticle,m_leftparticle);
                m_leftparticle->ChangePointer(m_leftparticle,m_leftparticle);

                // Set the pointer to the parent node
                if (oldrightparticle->m_parent->m_leftchild==oldrightparticle) {
                    oldrightparticle->m_parent->m_leftchild=m_leftchild;
                }
                else {
                    oldrightparticle->m_parent->m_rightchild=m_leftchild;
                }
                m_leftchild->m_parent=oldrightparticle->m_parent;

                BintreePrimary *oldrightchild=m_rightchild;
                BintreePrimary *oldparent=m_parent;

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
void BintreePrimary::ChangePointer(BintreePrimary *source, BintreePrimary *target)
{
    if(m_rightparticle == source) {
        m_rightparticle = target;
        real sphericalsurface =
            4 * PI * pow(3 * (m_leftparticle->Volume() +
                    m_rightparticle->Volume()) /(4 * PI), TWO_THIRDS);
        m_children_surf = sphericalsurface / (m_children_sintering *
                (1.0 - TWO_ONE_THIRD) + TWO_ONE_THIRD);
    }
    if(m_leftparticle == source){
        m_leftparticle = target;
        real sphericalsurface =
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

/*!
 * @brief       Sets all of the properties of the children to zero
 *
 * Used after the children are coalesced
 */
void BintreePrimary::ResetChildrenProperties()
{
    m_children_radius   = 0.0;
    m_children_vol      = 0.0;
    m_children_surf     = 0.0;
    m_children_sintering= 0.0;
    m_avg_sinter        = 0.0;
    m_sint_rate         = 0.0;
}

/*!
 * @brief       Updates the properties of a primary
 */
void BintreePrimary::UpdatePrimary(void)
{
    // Call parent class UpdateCache
    Primary::UpdateCache();

    // Set specific features of this node.
    m_primarydiam = m_diam;
    m_numprimary  = 1;
}

//! UpdateCache helper function
void BintreePrimary::UpdateCache(void)
{
    UpdateCache(this);
}

/*!
 * @brief       Updates the BintreePrimary cache
 *
 * Updates the whole particle, from root to the lowest leaf-node.
 * Calculates the sintering level and merges particles (through Check
 * Sintering) if necessary. Only accurately calculates collision
 * diameter and other properties used outside of BintreePrimary for
 * the root node.
 *
 * @param[in] root The root node of this particle
*/
void BintreePrimary::UpdateCache(BintreePrimary *root)
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
        for (unsigned int i=0; i != m_comp.size(); i++) {
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

        // Calculate the sintering level of the two primaries connected by this node
        m_children_sintering = SinteringLevel();
        if (m_children_sintering > 0.95) CheckSintering();

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
            real spherical_radius = pow(3 * m_vol / (4*PI), ONE_THIRD);
            m_diam = 2 * spherical_radius;

            // There are m_numprimary-1 connections between the primary
            // particles
            m_avg_sinter = m_avg_sinter / (m_numprimary - 1);

            // Approxmiate the surface of the particle
            // (same as in ChangePointer)
            const real numprim_1_3 = pow(m_numprimary,-1.0 * ONE_THIRD);
            m_surf = 4 * PI * spherical_radius * spherical_radius /
                    (m_avg_sinter * (1 - numprim_1_3) + numprim_1_3);

            // Calculate dcol based-on formula given in Lavvas et al. (2011)
            // assume fractal dimension Df = 1.8
            const real aggcolldiam = (6* m_vol / m_surf) *
                    pow(pow(m_surf, 3) / (36 * PI * m_vol * m_vol), (1.0/1.8));
            m_dmob = aggcolldiam;
            m_dcol = aggcolldiam;

        }
        else {
            m_diam=0;
            m_dmob=0;
        }
    }
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
void BintreePrimary::PrintTree(string filename) const
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
void BintreePrimary::PrintTreeLoop(std::ostream &out)const
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
void BintreePrimary::PrintTreeNode(std::ostream &out) const
{
    out
        << "|m_surf="             << this->m_surf
        << "|m_vol="             << this->m_vol
        << "|m_numprimary="      << this->m_numprimary
        << "|m_child_sint="      << this->m_children_sintering
        << "|m_child_rad="       << this->m_children_radius
        << "|m_child_surf="      << this->m_children_surf;
        for (unsigned int i=0; i != m_comp.size(); i++) {
            out << "|" + string(m_pmodel->Components(i)->Name()) + "=" <<
                    m_comp[i];
        }
    out
        << "|m_parent="          << this->m_parent
        << "|this=" << this;
}

void BintreePrimary::PrintComponents() const
{
    for (unsigned int i=0; i != m_comp.size(); i++) {
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
unsigned int BintreePrimary::Adjust(const fvector &dcomp,
        const fvector &dvalues, rng_type &rng, unsigned int n)

{
    if(m_numprimary == 1)
    {
        unsigned int i = 0;

        real dV;
        real m_vol_old = m_vol;

        // Add the components.
        for (i=0; i!=min(m_comp.size(),dcomp.size()); ++i)
        {
            m_comp[i] += dcomp[i] * (real)n;
            // Set component to zero if too many are removed.
            if (m_comp[i] < 1.0) m_comp[i] = 0.0;
        }

        // Add the tracker values.
        for (i=0; i!=min(m_values.size(),dvalues.size()); ++i)
        {
            m_values[i] += dvalues[i] * (real)n;
        }

        // Update only the primary
        UpdatePrimary();

        dV = m_vol - m_vol_old;

        // Surface change due to volume addition
        real dS = dV * 2.0 * m_pmodel->GetBintreeCoalThresh() / m_diam;

        // Climb back-up the tree and update the surface area and
        // sintering of a particle
        UpdateParents(dS);

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
    }

    // Update property cache.
    UpdateCache();
    return n;

}

/*!
 * @brief       Updates the surface area and sintering level of all parents
 *
 * @param[in]   dS      Surface area increment to adjust area by
 */
void BintreePrimary::UpdateParents(real dS) {
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
unsigned int BintreePrimary::AdjustIntPar(const fvector &dcomp,
        const fvector &dvalues, rng_type &rng, unsigned int n)

{
    if(m_numprimary == 1)
    {
        unsigned int i = 0;

        // Add the components.
        for (i=0; i!=min(m_comp.size(),dcomp.size()); ++i)
        {
            m_comp[i] += dcomp[i] * (real)n;
            // Set component to zero if too many are removed.
            if (m_comp[i] < 1.0) m_comp[i] = 0.0;
        }

        // Add the tracker values.
        for (i=0; i!=min(m_values.size(),dvalues.size()); ++i)
        {
            m_values[i] += dvalues[i] * (real)n;
        }
    }
    else
    {
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
void BintreePrimary::Sinter(real dt, Cell &sys,
                            const Processes::SinteringModel &model,
                            rng_type &rng,
                            real wt)
{
    //Only update the time on the root node
    if (m_parent == NULL) {
        m_sint_time += dt;
        SetSinteringTime(m_sint_time);
    }

    //Do only if there is a particle to sinter
    if (m_leftparticle!=NULL)
    {
        // Store the old surface area of particles
        //real surf_old = m_children_surf;

        // Calculate the spherical surface
        const real spherical_surface=4*PI*m_children_radius*m_children_radius;

        // Declare time step variables.
        real t1=0.0, delt=0.0, tstop=dt;
        real r=0.0;

        // Define the maximum allowed change in surface
        // area in one internal time step (10% spherical surface).
        real dAmax = 0.1 * spherical_surface;

        // The scale parameter discretises the delta-S when using
        // the Poisson distribution.  This allows a smoother change
        // (smaller scale = higher precision).
        real scale = 0.01;

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
                real mean;

                if (tstop > (t1+delt)) {
                    // A sub-step, we have changed surface by dAmax, on average
                    mean = 1.0 / scale;
                } else {
                    // Step until end.  Calculate degree of sintering explicitly.
                    mean = r * (tstop - t1) / (scale*dAmax);
                }
                boost::random::poisson_distribution<unsigned, real> repeatDistribution(mean);
                const unsigned n = repeatDistribution(rng);

                // Adjust the surface area.
                if (n > 0) {
                    m_children_surf -= (real)n * scale * dAmax;

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

        m_children_sintering = SinteringLevel();
        m_sint_rate = r;

        // Check if the sintering level is above the threshold, and merge
        if(m_children_sintering > 0.95)
          {
                CheckSintering();
                UpdateCache();
                if (m_leftchild!=NULL && m_rightchild!=NULL)
                {
                    m_leftchild->Sinter(dt,sys,model,rng,wt);
                    m_rightchild->Sinter(dt,sys,model,rng,wt);
                }
           }
         else {
             m_leftchild->Sinter(dt, sys, model,rng,wt);
             m_rightchild->Sinter(dt, sys, model,rng,wt);
         }

        UpdateCache();

        m_children_sintering = SinteringLevel();
    }  // endif m_leftparticle != NULL

}

//! Returns the sintering rate
real BintreePrimary::GetSintRate() const
{
    real sint_rate = m_sint_rate;
    if(m_leftchild!=NULL)
    {
        sint_rate += (m_leftchild->GetSintRate() +
                m_rightchild->GetSintRate());
    }

    return sint_rate;
}




/*!
 * @brief       Writes a particle to the binary stream
 *
 * @param out   Output stream
 */
void BintreePrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        if (m_pmodel->WriteBinaryTrees()) {
            // Call the binary tree serialiser...
            BinTreeSerializer <BintreePrimary> tree;
            tree.Serialize(out, this);
        } else {
            // Just serialise the root node.
            SerializePrimary(out);
        }

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, BintreePrimary::Serialize).");
    }
}

/*!
 * @brief       Writes an individual primary to binary stream
 *
 * @param out   Output binary stream
 */
void BintreePrimary::SerializePrimary(std::ostream &out) const
{
    if (out.good()) {

        int  val_int(0);
        real val(0.0);

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
                               "(Sweep, BintreePrimary::SerializePrimary).");
    }
}



/*!
 * @brief       Deserialise the binary tree
 *
 * Will only deserialise the full binary tree if the particle model has
 * writing/reading of trees activated. Otherwise, it will just load the
 * root node data.
 *
 * @param in        Input binary stream
 * @param model     Particle model being used
 */
void BintreePrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
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
            BinTreeSerializer <BintreePrimary> tree;
            tree.Deserialize(in, this, model);
        } else {
            // Just deserialise the root node.
            DeserializePrimary(in, model);
        }


    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, BintreePrimary::Deserialize).");
    }
}

/*!
 * @brief           Deserialise attributes of a single particle node
 *
 * @param in        Input binary stream
 * @param model     Particle model being used
 */
void BintreePrimary::DeserializePrimary(std::istream &in,
        const Sweep::ParticleModel &model)
{
    if (in.good()) {

        int  val_int(0);
        real val(0.0);

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
                               "(Sweep, BintreePrimary::DeserializePrimary).");
    }
}