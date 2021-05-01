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

#define _USE_MATH_DEFINES //!< First define.
#include <math.h>         //!< Then include so that the pi constant (M_PI) can be used.

#include "swp_bintree_primary.h"
#include "swp_bintree_serializer.h"
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/math/tools/roots.hpp>
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
	m_sum_necks(0.0),
	m_primaryvol(0.0),
    m_distance_centreToCentre(0.0),
    m_children_sintering(0.0),
    m_avg_sinter(0.0),
    m_sint_rate(0.0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL),
    m_r(0.0),
    m_r2(0.0),
    m_r3(0.0),
	m_Rg(0.0),
	m_tracked(false)
{
    m_cen_bsph[0] = 0.0;
    m_cen_bsph[1] = 0.0;
    m_cen_bsph[2] = 0.0;

    m_cen_mass[0] = 0.0;
    m_cen_mass[1] = 0.0;
    m_cen_mass[2] = 0.0;

	m_frame_orient_z[0] = 0.0;
	m_frame_orient_z[1] = 0.0;
	m_frame_orient_z[2] = 1.0e-9;

	m_frame_orient_x[0] = 1.0e-9;
	m_frame_orient_x[1] = 0.0;
	m_frame_orient_x[2] = 0.0;
}

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
	m_sum_necks(0.0),
	m_primaryvol(0.0),
    m_distance_centreToCentre(0.0),
    m_children_sintering(0.0),
    m_avg_sinter(0.0),
    m_sint_rate(0.0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL),
    m_r(0.0),
    m_r2(0.0),
    m_r3(0.0),
	m_Rg(0.0),
	m_tracked(false)
{
    m_cen_bsph[0] = 0.0;
    m_cen_bsph[1] = 0.0;
    m_cen_bsph[2] = 0.0;

    m_cen_mass[0] = 0.0;
    m_cen_mass[1] = 0.0;
    m_cen_mass[2] = 0.0;

	m_frame_orient_z[0] = 0.0;
	m_frame_orient_z[1] = 0.0;
	m_frame_orient_z[2] = 1.0e-9;

	m_frame_orient_x[0] = 1.0e-9;
	m_frame_orient_x[1] = 0.0;
	m_frame_orient_x[2] = 0.0;
}

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
m_sum_necks(0.0),
m_primaryvol(0.0),
m_distance_centreToCentre(0.0),
m_children_sintering(0.0),
m_avg_sinter(0.0),
m_sint_rate(0.0),
m_sint_time(0.0),
m_leftchild(NULL),
m_rightchild(NULL),
m_parent(NULL),
m_leftparticle(NULL),
m_rightparticle(NULL),
m_tracked(false)
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
		//if both particles are being tracked remove tracking from the one that will be deleted
		if (newleft->m_tracked == true && newright->m_tracked == true) newright->removeTracking();
    }
    else {
        newright->CopyParts(this);
        newleft->CopyParts(&copy_rhs);
		//if both particles are being tracked remove tracking from the one that will be deleted
		if (newright->m_tracked == true && newleft->m_tracked == true) newleft->removeTracking();
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

    //! It is assumed that primary pi from particle Pq and primary pj from
    //! particle Pq are in point contact and by default pi and pj are
    //! uniformly selected. If we track the coordinates of the primaries in
    //! a particle, we can do a smarter selection where pi and pj are
    //! determined by ballistic cluster-cluster aggregation (BCCA):
    //! R. Jullien, Transparency effects in cluster-cluster aggregation with
    //! linear trajectories, J. Phys. A 17 (1984) L771-L776.
    if (m_pmodel->getTrackPrimaryCoordinates()) {
        boost::uniform_01<rng_type&, double> uniformGenerator(rng);
        
		//! Calculate centre of mass and bounding sphere
		m_leftchild->calcBoundSph();
		m_leftchild->calcCOM();
		m_rightchild->calcBoundSph();
		m_rightchild->calcCOM();

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

		if (true){	//Select between BCCA and DLCA (currently only using BCCA)
			///////////////////////////////////////////////////////////////////////////////
			//	Ballistic Cluster Cluster Aggregation (BCCA)
			//	Gives D_f ~ 1.91 and k_f ~ 1.333
			//	Lindberg et al., J. Comp. Phys. 397, 108799, (2019)
			///////////////////////////////////////////////////////////////////////////////
        
			//! Incremental translation.
			while (!Overlap) {
				//! Sphere point picking. This is the random direction step of
				//! Jullien's BCCA algorithm. It is incorrect to select spherical
				//! coordinates theta (azimuthal angle) and phi (polar angle) from
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
				double v[3][3] = {0};

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
				double I[3][3] = {0};

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
				double R[3][3] = {0};
				double Mult = 1.0 / (1.0 + bx * x + by * y + bz * z); //!< Dot product of the two unit vectors.

				//! Square of the matrix v.
				R[0][0] = Mult * (-vy * vy - vz * vz);
				R[0][1] = Mult * (vx * vy);
				R[0][2] = Mult * (vx * vz);
				R[1][0] = Mult * (vx * vy);
				R[1][1] = Mult * (-vx * vx - vz * vz);
				R[1][2] = Mult * (vy * vz);
				R[2][0] = Mult * (vx * vz);
				R[2][1] = Mult * (vy * vz);
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
				// http://mathworld.wolfram.com/DiskPointPicking.html
				double r = uniformGenerator();
				theta  = 2.0 * PI * uniformGenerator();

				double sumr = m_leftchild->Radius() + m_rightchild->Radius();

				double x3 = sumr * sqrt(r) * cos(theta);
				double y3 = sumr * sqrt(r) * sin(theta);
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
						//this->m_leftchild->centreCOM();
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
			///////////////////////////////////////////////////////////////////////////////
		}
		else{	//DLCA disabled
			///////////////////////////////////////////////////////////////////////////////
			//	Diffusion Limited Cluster Aggregation (DLCA)
			//	Gives D_f ~ 1.8 and k_f ~ 1.36
			//	Instead of ballistic trajectories as with BCCA, DLCA models Brownian motion. 
			//	The "bullet" particle takes steps in random directions with step size 
			//	equal to the average primary diameter.
			///////////////////////////////////////////////////////////////////////////////
			while (!Overlap) {
				//! Generate random point on sphere
				double theta = 2.0 * PI * uniformGenerator();
				double phi = acos(2.0 * uniformGenerator() - 1.0);

				//! In terms of Cartesian coordinates.
				double x = cos(theta) * sin(phi);
				double y = sin(theta) * sin(phi);
				double z = cos(phi);

				double sumr = m_leftchild->Radius() + m_rightchild->Radius();

				// translate the left child
				this->m_leftchild->Translate(sumr*x, sumr*y, sumr*z);

				//! The two particles are initially at the origin. Keep doubling
				//! the distance between them so that they do not overlap. This is
				//! more efficient than picking an arbitrarily large distance.
				int numberOfOverlaps = 0;
				int factorApart = 1;
				double Separation = 0.0;

				while (this->checkForOverlap(*m_leftchild, *m_rightchild, numberOfOverlaps, Separation)) {
					this->m_leftchild->Translate(factorApart * x * sumr, factorApart * y * sumr, factorApart * z * sumr);
					factorApart *= 2;
				}

				double dx = this->m_leftchild->m_cen_bsph[0];
				double dy = this->m_leftchild->m_cen_bsph[1];
				double dz = this->m_leftchild->m_cen_bsph[2];

				double initialDistance = dx * dx + dy * dy + dz * dz;

				numberOfOverlaps = 0;

				//Loop over Brownian steps 
				while (!Overlap) {

					//Generate a random direction
					double theta2 = 2.0 * PI * uniformGenerator();
					double phi2 = acos(2.0 * uniformGenerator() - 1.0);

					//! In terms of Cartesian coordinates.
					double x2 = cos(theta2) * sin(phi2);
					double y2 = sin(theta2) * sin(phi2);
					double z2 = cos(phi2);

					//Brownian step size: this is the average primary size
					double step_size = m_leftchild->m_primarydiam / m_leftchild->m_numprimary;

					//First take an entire Brownian step
					this->m_leftchild->Translate(step_size * x2, step_size * y2, step_size * z2);
					Overlap = this->checkForOverlap(*m_leftchild, *m_rightchild, numberOfOverlaps, Separation);

					//If particles overlap then reverse the step and retake it in smaller increments 
					//until the particles are approximately in point contact
					if (Overlap == true){

						//Reverse the step
						this->m_leftchild->Translate(-step_size * x2, -step_size * y2, -step_size * z2);
						Overlap = false;
						numberOfOverlaps = 0;

						//Take the Brownian step in small increments
						unsigned int increment = 0;
						unsigned int tot_increments = 10;
					
						while (Overlap == false && increment < tot_increments){

							//! Translate particle in 10% increments of Brownian step
							this->m_leftchild->Translate(step_size * x2 / tot_increments, step_size * y2 / tot_increments, step_size * z2 / tot_increments);

							Overlap = this->checkForOverlap(*m_leftchild, *m_rightchild, numberOfOverlaps, Separation);

							increment++;
						}
					}

					// Get the boundins sphere separation
					dx = this->m_leftchild->m_cen_bsph[0];
					dy = this->m_leftchild->m_cen_bsph[1];
					dz = this->m_leftchild->m_cen_bsph[2];

					double newDistance = dx * dx + dy * dy + dz * dz;

					//! If the left particle has failed to collide with the particle at the
					//! origin (first condition) i.e. travelled too far in the wrong direction
					//! or if there are multiple points of overlap (second condition), 
					//! the trial is abandoned and another trajectory is chosen.
					if ((newDistance > initialDistance*2) || (numberOfOverlaps > 1)) {
						this->m_leftchild->centreBoundSph();
						numberOfOverlaps = 0;
						Overlap = false;
						break;
					}
				}
			}
		}
		///////////////////////////////////////////////////////////////////////////////

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

		//! Rotate back to original orientation of one particle 
		this->rotateCOM(-theta1, V1);

    } else {
        //! Randomly select the primaries that are touching.
        this->m_leftparticle = m_leftchild->SelectRandomSubparticle(rng);
        this->m_rightparticle = m_rightchild->SelectRandomSubparticle(rng);
    }

    // Set the sintering times
    SetSinteringTime(std::max(m_sint_time, rhsparticle->m_sint_time));
    m_createt = max(m_createt, rhsparticle->m_createt);

    // Initialise the variables used to calculate the sintering level
    m_children_vol      = m_leftparticle->m_vol + m_rightparticle->m_vol;
    m_children_surf     = m_leftparticle->m_surf + m_rightparticle->m_surf;
    m_children_radius   = pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
    m_children_sintering= SinteringLevel();

    //! If the coordinates of the primary particles are tracked we can directly
    //! calculate the centre-to-centre primary particle distance which is not
    //! exactly equal to the sum of the primary particle radii. Primaries are
    //! not considered to be in point contact unless there is a small overlap.
    if (m_pmodel->getTrackPrimaryCoordinates()) {
        double x1 = m_leftparticle->m_cen_bsph[0];
        double y1 = m_leftparticle->m_cen_bsph[1];
        double z1 = m_leftparticle->m_cen_bsph[2];

        double x2 = m_rightparticle->m_cen_bsph[0];
        double y2 = m_rightparticle->m_cen_bsph[1];
        double z2 = m_rightparticle->m_cen_bsph[2];

        m_distance_centreToCentre = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));    
    } else if (m_pmodel->getTrackPrimarySeparation()) {
        m_distance_centreToCentre = m_leftparticle->m_primarydiam / 2.0 + m_rightparticle->m_primarydiam / 2.0;
    }

	assert(m_distance_centreToCentre >= 0.0);

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
 * @brief       Coagulates *this* and rhs particle
 *
 * Called by coagulate, this function joins sp1 (this) and sp2 (rhs).
 * Basic properties are first calculated, followed by derived.
 *
 * @param[in] rhs   Pointer to the particle to be coagulated with this
 * @param[in] rng   Random number generator
*/
BinTreePrimary &BinTreePrimary::Fragment(const Primary &rhs, rng_type &rng)
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

//! Check for the overlap of primary particles.
/*!
 * @brief       Copy state-space and derived properties from source *
 *  @return Whether there is the overlap of primary particles.
 */
bool BinTreePrimary::checkForOverlap(BinTreePrimary &target, BinTreePrimary &bullet, int &numberOfOverlaps, double &Separation)
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
        } else {
            //! Bullet is not a leaf, call sub-nodes.
            Overlap = checkForOverlap(target, *bullet.m_leftchild, numberOfOverlaps, Separation);
            Overlap = checkForOverlap(target, *bullet.m_rightchild, numberOfOverlaps, Separation) || Overlap;

            return Overlap;
        }
    } else {
        //! Target is not a leaf.
        if (bullet.isLeaf()) {
            //! Bullet is a leaf, call target sub-nodes.
            Overlap = checkForOverlap(*target.m_leftchild, bullet, numberOfOverlaps, Separation);
            Overlap = checkForOverlap(*target.m_rightchild, bullet, numberOfOverlaps, Separation) || Overlap;

            return Overlap;
        } else {
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
bool BinTreePrimary::particlesOverlap(const Coords::Vector &p1, double r1,
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
    } else {
        Separation = sqrt(dxsqr + dysqr + dzsqr);
        return true;
    }
}

//! Calculates the radius of gyration of a particle
//! assuming primaries are point masses.
double BinTreePrimary::RadiusOfGyration() const
{
    double sum=0;
    double mass;
    double totalmass=0;
    double r2;
    double Rg;
    double rix, riy, riz, rjx, rjy, rjz;
    vector<fvector> coords;

	//! If single primary then return Rg = 0 because primaries are treated as point particles
	if(m_numprimary == 1) {
		Rg = 0.0;

	}else{

		this->GetPriCoords(coords);
    
		if (m_pmodel->getTrackPrimaryCoordinates()) {
			//! Calculation is based on Eq. (1) in R. Jullien, Transparency effects
			//! in cluster-cluster aggregation with linear trajectories, J. Phys. A
			//! 17 (1984) L771-L776. 
			for (int i = 0; i!=coords.size(); ++i) {
				for (int j = 0; j!=coords.size(); ++j) {
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
		} else {
			for (unsigned int i=0; i!=coords.size(); ++i) {
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

//! Returns a vector of primary coordinates, radius, and mass (5D).
/*!
 *  @param[in] coords The first three returned values are the cartesian x, y, z
 *                    coordinates, the final value is the radius.
 */
void BinTreePrimary::GetPriCoords(std::vector<fvector> &coords) const
{
    if (isLeaf()) {
        fvector c(5);
        c[0] = m_cen_mass[0];
        c[1] = m_cen_mass[1];
        c[2] = m_cen_mass[2];
        c[3] = m_r;
		c[4] = m_mass;
        coords.push_back(c);
    } else {
        m_leftchild->GetPriCoords(coords);
        m_rightchild->GetPriCoords(coords);
    }
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
    SetNumCarbon(source->NumCarbon());
    SetFrag(source->Frag());

    //! Set BinTreePrimary model characteristics.
    m_numprimary              = source->m_numprimary;
    m_primarydiam             = source->m_primarydiam;
    m_children_radius         = source->m_children_radius;
    m_children_vol            = source->m_children_vol;
    m_children_surf           = source->m_children_surf;
	m_free_surf				  = source->m_free_surf;
	m_sum_necks				  = source->m_sum_necks;
	m_primaryvol			  = source->m_primaryvol;
    m_distance_centreToCentre = source->m_distance_centreToCentre;
    m_cen_bsph                = source->m_cen_bsph;
    m_cen_mass                = source->m_cen_mass;
    m_r                       = source->m_r;
    m_r2                      = source->m_r2;
    m_r3                      = source->m_r3;
    m_children_sintering      = source->m_children_sintering;
    m_avg_sinter              = source->m_avg_sinter;
    m_sint_rate               = source->m_sint_rate;
    m_sint_time               = source->m_sint_time;
	m_frame_orient_z		  = source->m_frame_orient_z;
	m_frame_orient_x		  = source->m_frame_orient_x;
	m_Rg				      = source->m_Rg;
	m_tracked				  = source->m_tracked;

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
    if (target < m_leftchild->m_numprimary)
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

		//! If the centre-centre separation is tracked, the sintering level is calculated
		//! as per Lindberg et al. (J. Comp. Phys. 397, 108799, (2019)):
		//! s = R_ij / min(r_i,r_j)
		}else{
			if (m_leftparticle != NULL && m_rightparticle != NULL) {
				
				double r_i = m_leftparticle->m_primarydiam/2.0;
				double r_j = m_rightparticle->m_primarydiam/2.0;
				double d_ij =  m_distance_centreToCentre;
				double x_ij = (d_ij*d_ij - r_j*r_j + r_i*r_i)/(2.0*d_ij);
				double R_ij = sqrt(r_i*r_i - x_ij*x_ij);

				slevel = R_ij / min(r_i,r_j);
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
		//	 UpdateCache();	//this is already done at the end of merge
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

/*!
 * @brief       Checks if condition for merger is met
 *
 * @return      Boolean telling if the merger condition is met
 */
bool BinTreePrimary::MergeCondition()
{
	bool condition=false;

	if(m_leftparticle != NULL) {
		//! The condition for whether a particle has coalesced depends on whether
		//! the distance between the centres of primary particles is tracked. 
		if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()) {
			//! If coordinates are not tracked, the condition depends on whether the rounding 
			//! level exceeds an arbitrarily high threshold s > 0.95.
			condition = (m_children_sintering > 0.95);
		} else {

			double r_i = m_leftparticle->m_primarydiam/2.0;
			double r_j = m_rightparticle->m_primarydiam/2.0;
			double d_ij =  m_distance_centreToCentre;
			
			if(d_ij <= 0.0){
				//ensures that particles are merged if the sintering step overshoots
				condition = true;
			}else{
				//! If primary coordinates are tracked,
				//! primaries are merged when the neck radius is 95% of the smaller primary radius.
				//! The second condition ensures that primaries are merged even if sintering overshoots
				//! i.e. the neck crosses the centre of smaller primary
				double x_ij = (d_ij*d_ij - r_j*r_j + r_i*r_i)/(2.0*d_ij);
				double R_ij = sqrt(r_i*r_i - x_ij*x_ij);	//!neck radius
				condition = R_ij/min(r_i,r_j) >= 0.95 || ( (pow(d_ij,2.0) - pow(max(r_i,r_j),2.0) + pow(min(r_i,r_j),2.0) )/(2.0*d_ij) ) <= 0.0;
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
	//! Declare pointers for coordinate/separation tracking model
	BinTreePrimary *small_prim; //!< smaller of the merging primaries
	BinTreePrimary *big_prim;	//!< larger of the merging primaries
	BinTreePrimary *new_prim;	//!< new (merged) primary

	//initialise parameters
	double r_big,r_small; //, d_ij, x_ij;
	double r_new = 0.0;

	//! If a merging primary is tracked save the frame orientation vectors
	Coords::Vector frame_x;
	Coords::Vector frame_z;
	bool update_tracking = false;
	if (m_rightparticle->m_tracked == true){
		update_tracking = true;
		m_rightparticle->m_tracked = false;
		frame_x = m_rightparticle->m_frame_orient_x;
		frame_z = m_rightparticle->m_frame_orient_z;
	}
	else if (m_leftparticle->m_tracked == true){
		update_tracking = true;
		m_leftparticle->m_tracked = false;
		frame_x = m_leftparticle->m_frame_orient_x;
		frame_z = m_leftparticle->m_frame_orient_z;
	}

    // Make sure this primary has children to merge
    if( m_leftchild!=NULL) {

		//! Update primaries
		m_leftparticle->UpdatePrimary();
		m_rightparticle->UpdatePrimary();

		//! If the centre to centre distance is tracked we need to know which is the smaller primary of the merging pair
		if(m_leftparticle->m_primarydiam > m_rightparticle->m_primarydiam){
			small_prim = m_rightparticle;
			big_prim = m_leftparticle;
		}else{
			small_prim = m_leftparticle;
			big_prim = m_rightparticle;
		}

		r_big = big_prim->m_primarydiam/2.0;
		r_small = small_prim->m_primarydiam/2.0;

		//! If primary coordinates (or separation) are tracked then we need to estimate the radius
		//! of the new merged primary, given the new merged volume and a list of necks with neighbours.
		//  Lindberg et al., J. Comp. Phys. 397, 108799, (2019):
		if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()){

			//! Merge primary volumes and calculate new primary radius
			double V_new = small_prim->m_primaryvol + big_prim->m_primaryvol; //! Use the volume derived from the geometric properties 
			// An alternative would be to use the composition derived volume. 
			// This corrects the geometric volume to match the composition volume (see lower bound note below)
		//	double V_new = small_prim->m_vol + big_prim->m_vol; 			
			
			//! Get list of neck areas
			fvector necks;
			small_prim->GetNecks(small_prim, this, necks);
			big_prim->GetNecks(big_prim, this, necks);
			
			//! Estimate new primary diameter
			//! Use the Newton Raphson method to find new radius using the old radius as the initial guess
			double r_guess = r_big;	// Initial guess
			double r_min = 0.0;		// Lower bound
			double r_max = 10.0*r_big; // Assumed upper bound. Could also use std::numeric_limits<double>::max();
			int r_digits = static_cast<int>(0.6*std::numeric_limits<double>::digits);  // Just over half the number of digits suggested as precision in Boost manual for Newton Raphson
			const boost::uintmax_t r_maxit = 20; // Maximum number of iterations to perform
			boost::uintmax_t r_it = r_maxit;

			r_new = boost::math::tools::newton_raphson_iterate(merge_radius_functor(V_new,necks),r_guess,r_min,r_max,r_digits,r_it);

			//sanity checks
			if (r_new < r_big*0.99){ // the radius should be larger than the radius of the larger of the merging primaries
				std::cout << "BinTreePrimary::Merge: r_new < r_big \n";
				assert(r_new >= r_big*0.99);
			}
			if (r_new > r_big*2.0){ // the new radius should not too much larger than old radius
				std::cout << "BinTreePrimary::Merge: large jump in primary radius \n";
				assert(r_new <= 2.0*r_guess);
			}
			if (!std::isnormal(r_new)){
				std::cout << "BinTreePrimary::Merge: Could not calculate new radius! \n";
			}
			
			//! Get the largest neck (area)
			double max_neck = 0.0;
			if (!necks.empty())	max_neck = *std::max_element(necks.begin(), necks.end());
			//! The lower bounds on the new radius are the old radius since new volume is merged
			//! and the largest neck radius
			//  Note: if the composition dervied volume is used in estimating the new radius then
			//	the old radius should be removed as a lower bound to allow for some adjustment
			r_min = std::max(pow(max_neck / M_PI, 0.5), r_big);

			//! Impose lower bound in case the Newton method fails to find a solution
			if (r_new < r_min){
				if ((abs(r_new - r_min) / r_min) > 1e-10)
					std::cout << "BinTreePrimary::Merge: Imposing lower bound on new radius! \n";
				r_new = r_min;	
			}

			//! set diameter of larger primary to calculated diameter
			big_prim->m_primarydiam = 2.0*r_new;

		}

        if (m_leftchild==m_leftparticle && m_rightchild==m_rightparticle)
        {
            //! This node has only two primaries in its subtree, it is possible
            //! that this node is not the root node and belongs to a bigger
            //! particle.

            // Sum up the components first
            for (size_t i=0; i != m_comp.size(); i++) {
                m_comp[i] = m_leftparticle->Composition(i) +
                        m_rightparticle->Composition(i);
            }

			new_prim = this; //!< new primary

			//! Update quantities for coordinate/separation tracking model
			//! m_primarydiam of the new particle is the larger of the diameters of the merging particles
			//! If the centre to centre separation isn't tracked this will be changed to the spherical equivalent in the call to UpdatePrimary
			new_prim->m_primarydiam = big_prim->m_primarydiam;

			if(m_pmodel->getTrackPrimaryCoordinates()){
				new_prim->m_cen_bsph = big_prim->m_cen_bsph;
				new_prim->m_cen_mass = big_prim->m_cen_mass;
			}

            //! Update the pointers that pointed to the two former children
			if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()) {
		        ChangePointer(m_leftchild,this);
			    ChangePointer(m_rightchild,this);
			}else{
				//! If coordinates/separation are tracked the pointer to the larger primary is changed first
				//! so that primary properties are correctly calculated when adding neighbours
				ChangePointer(big_prim, this, this, small_prim, r_new, r_big);
				ChangePointer(small_prim, this, this, small_prim, r_new, r_small);
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
        }else{

			//! If primary coordinates or primary separations are not tracked then
			//! select subtree to keep the tree balanced
			if (!m_pmodel->getTrackPrimarySeparation() && !m_pmodel->getTrackPrimaryCoordinates()) {

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

					m_rightparticle->UpdatePrimary();

					// Set the pointers from the leftprimary to the rightprimary
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

				}else{
					// Append to right subtree
					BinTreePrimary *oldrightparticle = m_rightparticle;
					for (size_t i=0; i != m_comp.size(); i++) {
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
			}else{
			//! If the coordinates/separations are tracked then the larger primary becomes the new primary

				//! the new (merged) primary
				new_prim = big_prim;

				//! left/right flag
				bool newleft;
				if(new_prim == m_leftparticle){
					newleft = true;
				}else{
					newleft = false;
				}

				//! update composition
				BinTreePrimary *oldparticle = small_prim;
				for (size_t i=0; i != m_comp.size(); i++) {
					new_prim->m_comp[i] =
							m_leftparticle->Composition(i) +
							m_rightparticle->Composition(i);
				}

				//! update pointers to neighbours
				new_prim->ChangePointer(new_prim,new_prim,this,small_prim,r_new,r_big);
				oldparticle->ChangePointer(oldparticle,new_prim,this,small_prim,r_new,r_small);

				// Set the pointer to the parent node
				BinTreePrimary *oldchild = NULL;
				if(newleft){
					if (oldparticle->m_parent->m_leftchild==oldparticle) {
						oldparticle->m_parent->m_leftchild=m_leftchild;
					}
					else {
						oldparticle->m_parent->m_rightchild=m_leftchild;
					}
					m_leftchild->m_parent=oldparticle->m_parent;

					oldchild    = m_rightchild;
				}else{
					if (oldparticle->m_parent->m_leftchild==oldparticle) {
						oldparticle->m_parent->m_leftchild=m_rightchild;
					}
					else {
						oldparticle->m_parent->m_rightchild=m_rightchild;
					}
					m_rightchild->m_parent=oldparticle->m_parent;

					oldchild    = m_leftchild;
				}
				BinTreePrimary *oldparent = m_parent;

				// Copy the properties of the former leftchild to this node
				// so that it can be removed from the aggregate tree structure
				CopyParts(oldchild);

				// Now break the links to the tree structure in oldleftchild
				// in order to free it
				oldchild->m_leftchild   = NULL;
				oldchild->m_rightchild  = NULL;
				delete oldchild;

				m_parent = oldparent;

				if (m_leftchild!=NULL) {
					m_rightchild->m_parent=this;
					m_leftchild->m_parent=this;
				}

				delete oldparticle;
			}

        }

		//if coordinates are tracked then update the tracked radii
		if(m_pmodel->getTrackPrimaryCoordinates()){
			new_prim->setRadius(new_prim->m_primarydiam / 2.0);
		}

		//! Update tracking
		if(update_tracking == true){
			new_prim->m_tracked = true;
			new_prim->m_frame_orient_x = frame_x;
			new_prim->m_frame_orient_z = frame_z;
		}

		UpdateCache();

    }

    return *this;
}

//! Returns f(r) and f'(r) for solving the new primary radius given a new volume using the Newton method
std::pair<double, double> BinTreePrimary::merge_radius_functor::operator()(double const& r)
{
	double neck_vol(0.0);
	double neck_area(0.0);
	double fr(0.0), dr(0.0); // f(r) and f'(r)

	fvector::const_iterator ineck;
	for(ineck = a_necks.begin(); ineck!=a_necks.end(); ineck++){
		// Lower bound on primary radius is imposed by largest neck radius
		double r_min = pow(*ineck/M_PI, 0.5);

		neck_vol += 2.0*r*r*r+pow(r*r-*ineck/M_PI , 1.5)-3.0*r*r*pow(r*r-*ineck/M_PI , 0.5);
		neck_area += 2.0*r*r - r*pow(r*r-*ineck/M_PI , 0.5) - r*r*r*pow(r*r-*ineck/M_PI, -0.5);
	}

	// Calculate sums of neck volumes and areas
	fr = a_vol - 4.0*M_PI*r*r*r/3.0 + M_PI*neck_vol/3.0;
	dr = -4.0*M_PI*r*r + M_PI*neck_area;
	
	return std::make_pair(fr,dr);
}

/*!
//function to return fvector of neck radii for merged primary
//work up from primary and determine
 * @param[in] prim Pointer to merging primary
 * @param[in] node Pointer to merging node
 * @param[in] necks vector of neck areas
 */
void BinTreePrimary::GetNecks(BinTreePrimary *prim, BinTreePrimary *node, fvector &necks){

	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;
	double A_nij = 0.0;

	//! Check parent node is a neck joining prim but is not the merging neck
	if (m_parent->m_leftparticle == prim && m_parent != node) {
		//! right particle is a neighbour of prim
		r_j = m_parent->m_rightparticle->m_primarydiam/2.0;		//! neighbour radius
		x_ij = (d_ij*d_ij - r_j*r_j + r_i*r_i) / d_ij / 2.0;    //! distance from prim to neck
		A_nij = M_PI*(r_i*r_i - x_ij*x_ij);						//! neck area

		necks.push_back(A_nij);		//! Add neck area to vector

	} else if(m_parent->m_rightparticle == prim && m_parent != node) {
		//! Left primary is a neighbour of prim
		r_j = m_parent->m_leftparticle->m_primarydiam/2.0;		//! neighbour radius
		x_ij = (d_ij*d_ij - r_j*r_j + r_i*r_i) / d_ij / 2.0;    //! distance from prim to neck
		A_nij = M_PI*(r_i*r_i - x_ij*x_ij);						//! neck area

		necks.push_back(A_nij);		//! Add neck area to vector

	}

	//! Continue working up the binary tree
	if(m_parent->m_parent != NULL){
		m_parent->GetNecks(prim, node, necks);
	}
}

/*!
 * @brief       Changes pointer from source to target when centre-centre separation is tracked
 *
 * If a primary neighbours the smaller of the merging pair the centre to centre separation is
 * re-estmated as the smaller of the sum of the separation or the sum of primary radii.
 *
 * @param[in] source		Pointer to the original particle
 * @param[in] target		Pointer to the new particle
 * @param[in] node			Pointer to merging neck (non-leaf node)
 * @param[in] small_prim	Pointer to the smaller of the merging primaries
 * @param[in] r_new			New primary diameter
 * @param[in] r_old			Old primary diameter
*/
void BinTreePrimary::ChangePointer(BinTreePrimary *source, BinTreePrimary *target, BinTreePrimary *node, 
									BinTreePrimary *small_prim, double const r_new, double const r_old)
{

	if(m_rightparticle == source) {
		//! left particle is a neighbour
		if (this != node){
			//! Not the merging node

			double d_ij_old = m_distance_centreToCentre;		//!< Old centre to centre separation
			double r_j = m_leftparticle->m_primarydiam/2.0;		//!< Radius of neighbour
			double x_ji = (d_ij_old*d_ij_old - r_old*r_old + r_j*r_j)/d_ij_old/2.0;		//!< Neighbour centre to neck distance

			//! Calculate new centre to centre separation
			double d_ij_new = std::max(x_ji+ sqrt(x_ji*x_ji - r_j*r_j + r_new*r_new), x_ji - sqrt(x_ji*x_ji - r_j*r_j + r_new*r_new));

			//! Update centre to centre separation ensuring primaries are at least in point contact
			m_distance_centreToCentre = std::min(r_j + r_new, d_ij_new);
			assert(m_distance_centreToCentre >= 0.0);

			//! adjust coordinates of new neighbour and all its neighbour
			//! this translates the branch along old separation vector d_ik to new separation
			if(m_pmodel->getTrackPrimaryCoordinates()){

				Coords::Vector u_ik = UnitVector(m_leftparticle->boundSphCentre(), target->boundSphCentre());	//!< old separation unit vector
				double d_ik = Separation(m_leftparticle->boundSphCentre(), target->boundSphCentre());			//!< old separation distance
				//! Translate the neighbour
				m_leftparticle->TranslatePrimary(u_ik, d_ik - m_distance_centreToCentre);
				//! Translate all neighbours of the neighbour except the old small_prim
				m_leftparticle->TranslateNeighbours(m_leftparticle, u_ik, d_ik - m_distance_centreToCentre, small_prim);
			}

			m_rightparticle = target;

		}else{
			m_rightparticle = NULL;
		}

    }

    if(m_leftparticle == source){
		//! right particle is a neighbour
		if (this != node){
			//! Not the merging node

			double d_ij_old = m_distance_centreToCentre;			//!< Old centre to centre separation
			double r_j = m_rightparticle->m_primarydiam/2.0;		//!< Radius of neighbour
			double x_ji = (d_ij_old*d_ij_old - r_old*r_old + r_j*r_j)/d_ij_old/2.0;		//!< Neighbour centre to neck distance

			//! Calculate new centre to centre separation
			double d_ij_new = std::max(x_ji+ sqrt(x_ji*x_ji - r_j*r_j + r_new*r_new), x_ji - sqrt(x_ji*x_ji - r_j*r_j + r_new*r_new));

			//! Update centre to centre separation ensuring primaries are at least in point contact
			m_distance_centreToCentre = std::min(r_j + r_new, d_ij_new);
			
			assert(m_distance_centreToCentre >= 0.0);

			//! adjust coordinates of new neighbour and all its neighbour
			//! this translates the branch along old separation vector d_ik to new separation
			if(m_pmodel->getTrackPrimaryCoordinates()){

				Coords::Vector u_ik = UnitVector(m_rightparticle->boundSphCentre(), target->boundSphCentre());	//!< old separation unit vector
				double d_ik = Separation(m_rightparticle->boundSphCentre(), target->boundSphCentre());			//!< old separation distance
				//! Translate the neighbour
				m_rightparticle->TranslatePrimary(u_ik, d_ik - m_distance_centreToCentre);
				//! Translate all neighbours of the neighbour except the old small_prim
				m_rightparticle->TranslateNeighbours(m_rightparticle, u_ik, d_ik - m_distance_centreToCentre, small_prim);
			}

			m_leftparticle = target;

		}else{
			m_leftparticle = NULL;
		}
    }

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->ChangePointer(source, target, node, small_prim, r_new, r_old);
    }

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
		// Impose upper bound on the common surface area 
		// equal to the sum of spherical primary surfaces.
		if((m_leftparticle->m_surf + m_rightparticle->m_surf) < m_children_surf){
				m_children_surf = m_leftparticle->m_surf + m_rightparticle->m_surf;
		}
    }
    if(m_leftparticle == source){
        m_leftparticle = target;
        double sphericalsurface =
            4 * PI * pow(3 * (m_leftparticle->Volume() +
                    m_rightparticle->Volume()) /(4 * PI), TWO_THIRDS);
        m_children_surf = sphericalsurface / (m_children_sintering *
                (1.0 - TWO_ONE_THIRD) + TWO_ONE_THIRD);
		// Impose upper bound on the common surface area
		// equal to the sum of spherical primary surfaces.
		if((m_leftparticle->m_surf + m_rightparticle->m_surf) < m_children_surf){
				m_children_surf = m_leftparticle->m_surf + m_rightparticle->m_surf;
		}
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
 *  @brief Updates the properties of a primary.
 */
void BinTreePrimary::UpdatePrimary(void)
{
    //! Call parent class UpdateCache.
    Primary::UpdateCache();

    //! Set specific features of this node.
    //! If the primary centre-centre separation isn't tracked, the primary
    //! coordinates or this is a single primary m_primarydiam and m_diam are
    //! both the spherical equivalent diameter and the free surface area is the
    //! sperhical surface area.
    if(!(m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()) || m_parent == NULL){
		//To Do (csl37): should not reset properties if this is a single primary due to a merge event
		m_primarydiam = m_diam;
		m_free_surf = m_surf;
		m_primaryvol = m_vol;
		m_sum_necks = 0.0;
	}else{
		//! Update overlapping primary model
		UpdateOverlappingPrimary();
	}

    m_numprimary  = 1;

    //! Initialisation of the radius of bounding sphere which is only relevant
    //! if the primary coordinates are tracked.
    if (m_pmodel->getTrackPrimaryCoordinates()) {
        setRadius(m_primarydiam / 2.0);
    }
}

/*!
 * @brief       Updates the BinTreePrimary cache from the root node
 *
 * Works up the tree to the root node and then calls UpdateCache
 *
*/
void BinTreePrimary::UpdateCacheRoot(void){
	if(m_parent != NULL){
		m_parent->UpdateCacheRoot();
	}else{
		UpdateCache(this);
	}
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
		m_Rg = 0.0;
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
		m_primaryvol	= m_leftchild->m_primaryvol + m_rightchild->m_primaryvol;
		//titania phase transformation term
		m_phaseterm		= m_leftchild->m_phaseterm + m_rightchild->m_phaseterm;

		//calculate bounding sphere
		calcBoundSph();

		// updates m_children_radius
		if ((m_leftparticle!=NULL) && (m_rightparticle!=NULL)){
			m_children_radius = pow((3.0/(4.0*PI))*(m_leftparticle->Volume() + m_rightparticle->Volume()),(ONE_THIRD));	
		}

		//! Particle tracking
		if (m_leftchild->m_tracked == true || m_rightchild->m_tracked == true) m_tracked = true;

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
		if (this->m_parent == NULL){	//if this does not have a parent this is the root node 
             // Get spherical equivalent radius and diameter
            double spherical_radius = pow(3 * m_vol / (4*PI), ONE_THIRD);
            m_diam = 2 * spherical_radius;

            // There are m_numprimary-1 connections between the primary
            // particles
            if (m_numprimary > 1)
                m_avg_sinter = m_avg_sinter / (m_numprimary - 1);

			//! if the centre to centre distance is tracked then
			//! the surface area is the particle free surface area
			//! and the collision diameter is calculated using the
			//! primary coordinates
			if (m_pmodel->getTrackPrimaryCoordinates()) {
				m_surf = m_free_surf;
				m_dcol = CollisionDiameter();
			}else{
				
				//! If primary separations are tracked then
				//! the surface area is the particle free surface area
				if (m_pmodel->getTrackPrimarySeparation()){
					m_surf = m_free_surf;
				}else{
					//! Approxmiate the surface of the particle
					// (same as in ChangePointer)
					const double numprim_1_3 = pow(m_numprimary, -1.0 * ONE_THIRD);
					m_surf = 4 * PI * spherical_radius * spherical_radius /
						(m_avg_sinter * (1 - numprim_1_3) + numprim_1_3);
				}

				//! Calculate dcol based-on formula given in Lavvas et al. (2011)
				const double aggcolldiam = (6 * m_vol / m_surf) *
					pow(pow(m_surf, 3) / (36 * PI * m_vol * m_vol),
					(1.0 / m_pmodel->GetFractDim()));
				m_dcol = aggcolldiam;

			}

			//! Radius of gyration
			m_Rg = RadiusOfGyration();

			//! Mobility diameter
            m_dmob = MobDiameter();

        } else {
            m_diam=0;
            m_dmob=0;
			m_Rg = 0.0;
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


//! Calculates the collision diameter based on the radius of gyration
//! if primary coordinates are tracked
//	Lindberg et al., J. Comp. Phys. 397, 108799, (2019)
double BinTreePrimary::CollisionDiameter()
{
    double sum=0.0;
    double dcol=0.0;
    vector<fvector> coords;
		
	//! Calculate centre of mass
	this->calcCOM();
	//! Save centre of mass coordinates
	double COM_x = m_cen_mass[0];
	double COM_y = m_cen_mass[1];
	double COM_z = m_cen_mass[2];
	
	//! Get a list of primary coordinates
	this->GetPriCoords(coords);

	//! Calculate Rg (mass weighted)
	//! This is based on Eq. (2) in Lapuerta et al., A method to determine 
	//! the fractal dimension of diesel soot agglomerates.
	//! Journal of Colloid Interface Science, 303:149-158. 2006.
	//! A modification is made to the radius of gyration of a single primary
	//! replacing r_gp = sqrt(3/5)*r_p with r_gp = r_p, as per Eq. (4) in 
	//! Filippov et al., Fractal-like aggregates: Relation between morphology
	//! and physical properties. Journal of Colloid Interface Science, 
	//! 229:261-273, 2000.
	for (int i = 0; i!=coords.size(); ++i) {
		
		//! Add square of distance from the CoM weighted by mass
		sum += coords[i][4] * (pow((coords[i][0] - COM_x),2.0) + pow((coords[i][1] - COM_y),2.0) + pow((coords[i][2] - COM_z),2.0));
		
		//! Add square of primary radius weighted by mass
		sum += coords[i][4] * coords[i][3] * coords[i][3];

	}

	dcol = 2*sqrt(sum/m_mass);

	return dcol;
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
		double r_old = m_primarydiam / 2.0;

        // Call to Primary to adjust the state space
        n = Primary::Adjust(dcomp, dvalues, rng, n);

        // Stop doing the adjustment if n is 0.
        if (n > 0) {
            // Update only the primary
            UpdatePrimary();

			//! If the distance between the centres of primary particles or the
            //! primary coordinates are tracked, the rate of change in the
            //! primary diameter is affected by its neighbours.
			if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()) {
				
				//! Particle with more than one primary.
				if (m_parent != NULL) {	
					
					while (volOld <= m_vol){

						//! Initialisation of variables to adjust the primary diameter if the
						//! distance between the centres of primary particles is tracked.
						double r_i = m_primarydiam / 2.0;	//!< Radius of primary particle i.
						double dr_max = 0.01*r_i;			//!< Maximum change in primary radius during internal step (1% of primary radius)
						double dr = 0.0;					//!< Change in radius of i
					
						//Calculate change in volume
						dV = dr_max * m_free_surf;

						//Calculate change in radius
						if (volOld + dV > m_vol){
							dr = (m_vol - volOld)*dr_max / dV;
						}else{
							dr = dr_max;
						}

						//Update primary diameter
						m_primarydiam = 2.0* (r_i + dr);

						//! Update free surface area
						this->UpdateOverlappingPrimary();
						
						volOld += dV;
					}

					//if coordinates are tracked then update coordinate tracking properties
					if(m_pmodel->getTrackPrimaryCoordinates()){
						setRadius(m_primarydiam / 2.0);
						this->calcBoundSph();			//(csl37) are these necessary here?
						this->calcCOM();
					}

					//adjust composition of neighbours due to change in neck position
					this->AdjustNeighbours(this, m_primarydiam / 2.0 - r_old, dcomp, dvalues, rng);

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

			//! Update the cache from the root node
			UpdateCacheRoot();
        }
    }
    // Else this a non-leaf node (not a primary)
    else
    {
		
		if (m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates()) {
			
			//  Note (csl37): Adjust could be improved for deferred growth processes by distributing the new mass
			//	over all or multiple primaries rather than applying all of the deferred growth to a single primary.
			//
			//! Primary selection based on the relative free surface area
			//! Work down tree selecting left/right child based on sum of primary free surface areas under node
			//! generate random number, use bernoulli with p=free_surf(leftchild)/free_surf(this)
			boost::bernoulli_distribution<> leftRightChooser(m_leftchild->m_free_surf/m_free_surf);
			if(leftRightChooser(rng)){
				return m_leftchild->Adjust(dcomp, dvalues, rng, n);
			}else{
				return m_rightchild->Adjust(dcomp, dvalues, rng, n);
			}
			
		}else{

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

		}
    }

    // Update property cache.
    // UpdateCache(this); // (csl37) No need to do this here: UpdateCache will be called again in Particle::Adjust

    return n;

}


/*!
 * @brief       Adjusts the particle after a phase transformation event
 *
 * Analogous to the implementation in Primary. The function will
 * however descend the tree to find a primary adjust. 
 *
 * @param[in]   dcomp   Vector storing changes in particle composition
 * @param[in]   dvalues Vector storing changes in gas-phase comp
 * @param[in]   rng     Random number generator
 * @param[in]   n       Number of times for adjustment
 */
unsigned int BinTreePrimary::AdjustPhase(const fvector &dcomp,
        const fvector &dvalues, rng_type &rng, unsigned int n)
{

	if (m_leftchild == NULL && m_rightchild == NULL) {

        // Call to Primary to adjust the state space
        n = Primary::Adjust(dcomp, dvalues, rng, n);
		
        // Stop doing the adjustment if n is 0.
        if (n > 0) {
            // Update only the primary
            UpdatePrimary();
        }
    }
    // Else this a non-leaf node (not a primary)
	// Select primary to adjust
    else
    {	
        return SelectRandomSubparticle(rng)->AdjustPhase(dcomp, dvalues, rng, n);
    }

    // Update property cache.
    UpdateCache(this);

    return n;
}

//Melting point dependent phase change
void BinTreePrimary::Melt(rng_type &rng, Cell &sys){

	// Update the children
	if (m_leftchild != NULL) {
		m_leftchild->Melt(rng, sys);
		m_rightchild->Melt(rng, sys);
	}
	else{ // this is a primary

		Primary::Melt(rng, sys);
	}
}

/*! 
* @brief		Adjust composition of neighbours
*
* Called by Adjust() during a growth event to rebalance the 
* composition of a primary and its neighbours. Rebalancing is necessary
* if the position of the neck has moved due to a change in radius.
*
* @param[in]   prim		Pointer to primary adjusted
* @param[in]   delta_r  Change in radius
* @param[in]   dcomp	Vector storing changes in particle composition
* @param[in]   dvalues	Vector storing changes in gas-phase comp
* @param[in]   rng		Random number generator
*/
void BinTreePrimary::AdjustNeighbours(BinTreePrimary *prim, const double delta_r, const fvector &dcomp,const fvector &dvalues, rng_type &rng){
	
	BinTreePrimary *neighbour = NULL;

	//! Check if parent node contains a neighbour of prim
	if (m_parent->m_leftparticle == prim) {
		//! right particle is a neighbour
		neighbour = m_parent->m_rightparticle;
	} else if (m_parent->m_rightparticle == prim) {
		//! left particle is a neighbour
		neighbour = m_parent->m_leftparticle;
	} 

	//! Adjust neighbour's composition
	if (neighbour != NULL){

		double r_i = prim->m_primarydiam / 2.0;						//!< primary radius
		double r_j = neighbour->m_primarydiam / 2.0;				//!< neighbouring primary radius
		double d_ij = m_parent->m_distance_centreToCentre;			//!< centre to centre separation
		double x_ij = (d_ij*d_ij - r_j*r_j + r_i*r_i) / d_ij / 2.0;	//!< distance neck to centre of primary p_i
		double A_nij = M_PI*(r_i*r_i - x_ij*x_ij);					//!< neck area

		unsigned int max_n = 0;

		//! unit change in volume
		// this is the volume given by dcomp
		double uvol = 0.0; //!< unit change in volume
		double m = 0.0;
		for (int i = 0; i != dcomp.size(); ++i) {
			m = m_pmodel->Components(i)->MolWt() * dcomp[i] / NA;
			//! max possible change to the composition leaving 1 component
			if (i > 0) {
				max_n = min(max_n, static_cast<unsigned int>((neighbour->Composition(i) - 1.0) / dcomp[i]));
			}else{
				max_n = static_cast<unsigned int>((neighbour->Composition(i) - 1.0) / dcomp[i]);	//initial value for max_n
			}
			if (m_pmodel->Components(i)->Density() > 0.0)
				uvol += m / m_pmodel->Components(i)->Density();
		}

		//! change in volume of prim
		double dvol = A_nij * (r_i - delta_r) * delta_r / d_ij; //use the old radius (r_i - delta_r) here
		//! change in composition
		unsigned int dn = min(max_n, static_cast<unsigned int>(dvol / uvol));

		//! adjust primary compositions if change is large enough
		if (dvol > 0.0 && dn > 0){
			//! composition change vectors for the neighbour
			fvector dcomp_neighbour(dcomp.size());
			fvector dvalues_neighbour(dvalues.size());
			for (int i = 0; i != dcomp.size(); ++i) {
				dcomp_neighbour[i] = -dcomp[i];
			}
			for (int i = 0; i != dvalues.size(); ++i) {
				dvalues_neighbour[i] = -dvalues[i];
			}
			//! Adjust neighbour's composition (decrease)
			unsigned int n = Primary::Adjust(dcomp_neighbour, dvalues_neighbour, rng, dn);
			//! Adjust primary's composition (increase)
			unsigned int m = Primary::Adjust(dcomp, dvalues, rng, n);
			assert(m == n); //decrease in neighbour's composition must be equal to increase in primary's composition
		}
	}

	//! Continue working up the binary tree
	if (m_parent->m_parent != NULL){
		m_parent->AdjustNeighbours(prim,delta_r,dcomp,dvalues,rng);
	}
}


/*! Updates primary free surface area and volume
*
* @param[in]   this		Primary to update
*/
void BinTreePrimary::UpdateOverlappingPrimary(){

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
	m_primaryvol = max( M_PI*pow(m_primarydiam,3.0)/6.0 - CapVolumes , 0.0);

	//! Update sum of necks 
	m_sum_necks = SumNecks;

	//sanity checks for debugging
	assert(m_sum_necks >= 0.0);
	assert(m_free_surf >= 0.0);

}

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
void BinTreePrimary::SumCaps(BinTreePrimary *prim, double &CapAreas, double &CapVolumes, double &SumNecks){
	
	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;

	//! Check if a neighbour of prim
	if (m_parent->m_leftparticle == prim) {
		
		//! Right primary is a neighbour
		r_j = m_parent->m_rightparticle->m_primarydiam / 2.0;
		x_ij = min( ( pow(d_ij,2.0) - pow(r_j,2.0) + pow(r_i,2.0) ) / ( 2.0*d_ij ), r_i); //ensure -r_i <= x_ij <= r_i
		x_ij = max( x_ij, -r_i );

		//! Calculate cap area and add to sum
		CapAreas += 2*M_PI*(r_i*r_i - r_i*x_ij);

		//! Calculate cap volume and add to sum
		CapVolumes += M_PI * (2*pow(r_i,3.0) + pow(x_ij,3.0) - 3.0*pow(r_i,2.0)*x_ij ) /3.0;

		//! Neck area * r_i / x_ij
		SumNecks +=  abs(M_PI*(r_i*r_i - x_ij*x_ij) * r_i / x_ij);

	} else if(m_parent->m_rightparticle == prim) {
		
		//! Left primary is a neighbour
		r_j = m_parent->m_leftparticle->m_primarydiam / 2.0;
		x_ij = min( ( pow(d_ij,2.0) - pow(r_j,2.0) + pow(r_i,2.0) ) / ( 2.0*d_ij ), r_i);	//ensure -r_i <= x_ij <= r_i
		x_ij = max( x_ij, -r_i );

		//! Calculate cap area and add to sum
		CapAreas += 2*M_PI*(r_i*r_i - r_i*x_ij);

		//! Calculate cap volume and add to sum
		CapVolumes += M_PI * (2*pow(r_i,3.0) + pow(x_ij,3.0) - 3.0*pow(r_i,2.0)*x_ij ) /3.0;
	
		//! Neck area * r_i / x_ij
		SumNecks += abs(M_PI*(r_i*r_i - x_ij*x_ij) * r_i / x_ij);

	}

	//sanity checks for debugging
	assert(r_i >= 0.0);
	assert(r_j >= 0.0);
	assert(SumNecks >= 0.0);
	assert(2*M_PI*(r_i*r_i - r_i*x_ij) >= 0.0); //cap area

	//! Continue working up the binary tree
	if(m_parent->m_parent != NULL){
		m_parent->SumCaps(prim, CapAreas, CapVolumes, SumNecks);
	}
}

/*!
 * @brief       Identify neighbours and update centre to centre separation and 
 *				coordinates except for specified neighbour.
 *
 * Works up the binary tree identifying the neighbours of the adjusted primary 
 * and updates the centre to centre separation. Also sums the contribution from
 * neighbours to the free surface area.
 *
 * @param[in]   prim			Pointer to the primary being adjusted
 * @param[in]   delta_r			Change in radius of prim
 * @param[in]   prim_ignore		Pointer to the primary to be ignored
 */
void BinTreePrimary::UpdateConnectivity(BinTreePrimary *prim, double delta_r, BinTreePrimary *prim_ignore){
	
	double d_ij = m_parent->m_distance_centreToCentre;
	double r_i = prim->m_primarydiam / 2.0;
	double r_j = 0.0;
	double x_ij = 0.0;
	BinTreePrimary *neighbour = NULL;

	//! check if a neighbour of prim
	if (m_parent->m_leftparticle == prim && m_parent->m_rightparticle != prim_ignore ) {
		//! right particle is a neighbour
		neighbour = m_parent->m_rightparticle;
		r_j = neighbour->m_primarydiam / 2.0;
	} else if(m_parent->m_rightparticle == prim &&  m_parent->m_leftparticle != prim_ignore ) {
		//! left particle is a neighbour
		neighbour = m_parent->m_leftparticle;
		r_j = neighbour->m_primarydiam / 2.0;
	} else {
		//! not a neighbour
		r_j = 0.0;
	}

	if(r_j > 0.0){

		double d_ij_old = d_ij;
		x_ij = ( pow(d_ij,2.0) - pow(r_j,2.0) + pow(r_i,2.0) ) / ( 2.0*d_ij );
		//! update centre to centre separation
		//making sure centre to centre separation remains smaller than the sum of the radii
		d_ij = min(d_ij + r_i * delta_r / x_ij, r_i+r_j+delta_r);
		double d_ij_min = max(r_i - r_j, r_j - r_i); //!< minimum separation 
		//and larger than the minimum possible separation (where one primary enevelopes the other)
		d_ij = max(d_ij, d_ij_min);

		m_parent->m_distance_centreToCentre = d_ij;
		
		//! if primary coordinates are tracked then we need to update the coordinates of the neighbour 
		if (m_pmodel->getTrackPrimaryCoordinates()) {
			//! get (unit) vector separating prim and neighbour
			Coords::Vector u = UnitVector(prim->boundSphCentre(), neighbour->boundSphCentre()); 
			double delta_d = d_ij-d_ij_old; //!< change in separation (magnitude)
			//! translate the neighbour 
			neighbour->TranslatePrimary(u, delta_d);
			//! translate all neighbours of the neighbour except prim
			neighbour->TranslateNeighbours(neighbour,u, delta_d, prim);
		}
	}

	//continue working up the binary tree
	if(m_parent->m_parent != NULL){
		m_parent->UpdateConnectivity(prim, delta_r, prim_ignore);
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
        //m_parent->UpdateCache();	//updating here may cause a seg fault if a merger event occurs
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
        return SelectRandomSubparticle(rng)->Adjust(dcomp, dvalues, rng, n);
        // Generate random numbers
        //boost::bernoulli_distribution<> leftRightChooser;

        // Select particle
        //if(leftRightChooser(rng))
      //return m_leftparticle->AdjustIntPar(dcomp, dvalues, rng, n);
      //else
      //return m_rightparticle->AdjustIntPar(dcomp, dvalues, rng, n);
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

    //! The sintering model depends on whether the distance between the centres
    //! of primary particles or the coordinates of the primary particles are
    //! tracked. If tracked, sintering results in a decrease in the distance
    //! between the primaries and an increase in their diameters. If not,
    //! sintering results in a decrease in the common surface between the
    //! primaries.
    if (!(m_pmodel->getTrackPrimarySeparation() || m_pmodel->getTrackPrimaryCoordinates())) {

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

		///////////////////////////////////////////////////////
		/// References to equations in Langmuir 27:6358 (2011).
		///////////////////////////////////////////////////////

		//make sure particles are up to date
		m_leftparticle->UpdatePrimary();
		m_rightparticle->UpdatePrimary();
		
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
            if (!MergeCondition()) {

				//! Due to rounding, x_i and x_j are sometimes calculated to be larger 
				//! than the respective primary radii resulting in a negative neck area.
				//! Therefore we take the smaller of x_i and r_i.
				double x_i = min((d_ij2 - r_j2 + r_i2) / (2.0 * d_ij),r_i); //!< Eq. (3b).
				double x_j = min((d_ij2 - r_i2 + r_j2) / (2.0 * d_ij),r_j); //!< Eq. (3b).
				double A_n = M_PI * (r_i2 - pow(x_i, 2.0));        //!< Eq. (4).
			
				//declare more variables
				double dd_ij_dt=0.0;
				double R_n = 0.0;
				double r4_tau = 0.0;
				double tau = 0.0;
				double gamma_eta = 0.0;

				//! Sintering model dependent part
				//! Viscous flow model
				if(model.Type() == Processes::SinteringModel::ViscousFlow){
						
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

					//! Eq. (14a).
					dd_ij_dt = 4.0 * r_i * r_j * d_ij2 * (r_i + r_j) * gamma_eta /
								((r_i + r_j + d_ij) * (r_i4  + r_j4 - 2.0 * r_i2 * r_j2 + 4.0 * d_ij * r_i * r_j *(r_i + r_j) - d_ij2 * (r_i2 + r_j2)));

				//! Grain boundary diffusion model 
				}else if(model.Type() == Processes::SinteringModel::GBD){ 

					//! If the particles are in point contact set an initial neck radius of 1% 
					//! of the smaller primary radius, otherwise dd_ij_dt would be undefined
					if(A_n <= 0.0){
						R_n = 0.01*min(r_i,r_j);
						A_n = M_PI * R_n * R_n;
						x_i = sqrt(r_i2 - R_n * R_n);
						x_j = sqrt(r_j2 - R_n * R_n);
					}else{
						R_n = sqrt(A_n / M_PI);
					}

					//! The primary radius in the numerator cancels with the diameter dependence of tau
					//! so we can calculate this for only one of the primaries.
					//! Use smaller primary in case a minimum diameter is imposed for sintering
					//  In SintTime the diameter is calculated as 6.0 * m_vol / m_surf
					//  so r = 3.0 * m_vol / m_surf
					BinTreePrimary * small_prim;
					if (r_i <= r_j){
						small_prim = m_leftparticle;
					}else{
						small_prim = m_rightparticle;
					}
					double r4 = pow(3.0 * small_prim->m_vol / small_prim->m_surf, 4.0);
					r4_tau = r4 / model.SintTime(sys, *small_prim);

					//! J Aerosol Sci 46:7-19 (2012) Eq. (A6)
					//! dx_i_dt + dx_j_dt
					//! (this is missing a minus sign, which is accounted for below)
					dd_ij_dt = r4_tau * ( 1/(r_i - x_i) + 1/(r_j - x_j) - 2/R_n ) / A_n;

				//Other model are not coded
				}else{
					std::cout<<"Sintering model not coded"<<endl;
					break;
				}

				//! Get surface area and subtract mutual contribution
				double A_i = std::max(0.0,m_leftparticle->m_free_surf + m_leftparticle->m_sum_necks - M_PI*(r_i*r_i - x_i*x_i)*r_i/x_i);
				double A_j = std::max(0.0,m_rightparticle->m_free_surf + m_rightparticle->m_sum_necks - M_PI*(r_j*r_j - x_j*x_j)*r_j / x_j);

				//! @todo Remove derivation and replace with reference to preprint
				//!       or paper if results do get published.
				double B_i = (-r_j*A_n*A_n - x_j*A_j*A_n)/(A_i*A_j*d_ij + r_i*A_j*A_n + r_j*A_i*A_n);
				double B_j = (-r_i*A_n*A_n - x_i*A_i*A_n)/(A_j*A_i*d_ij + r_j*A_i*A_n + r_i*A_j*A_n);
				
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

				double d_ij_min = max(r_i - r_j, r_j - r_i); //!< minimum possible separation (where one primary enevelopes the other)
				double delta_dij = -(double)n * scale * dd_ij_Max; //!< Change in separation (sintering decreases d_ij hence the negative sign)
				//! Make sure that sintering doesn't overshoot
				if (d_ij + delta_dij < d_ij_min){
					delta_dij = d_ij_min - d_ij;
				}

				//! adjust separation
				m_distance_centreToCentre += delta_dij; 

				//! if coordinates are tracked then we will translate one side of the particle by the change in separation
				//! this is faster than translating both sides by half the change
				if (m_pmodel->getTrackPrimaryCoordinates()) {
					//! get direction of translation (left particle to right particle)
					Coords::Vector vector_change = UnitVector(m_leftparticle->boundSphCentre(), m_rightparticle->boundSphCentre());
					//! translate the leftparticle
					// -delta_dij because delta_dij is negative (the direction of translation is determined by the vector) 
					m_leftparticle->TranslatePrimary(vector_change, -delta_dij);
					//! translate all neighbours of the left particle except the right particle
					m_leftparticle->TranslateNeighbours(m_leftparticle,vector_change,-delta_dij,m_rightparticle);
				}
				
				//! Change in primary radii
				double delta_r_i = delta_dij * B_i;  //!< Eq. (8).
				double delta_r_j = delta_dij * B_j;  //!< Eq. (8).

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

				// Should return some sintering rate (units of m2/s expected)
				// r = dd_ij_dt;

            } else {
                break; //!do not continue to sinter.
            }
        }

		//! If coordinates are tracked update tracking radius 
		if (m_pmodel->getTrackPrimaryCoordinates()) {		
			m_leftparticle->setRadius(m_leftparticle->m_primarydiam / 2.0);
			m_rightparticle->setRadius(m_rightparticle->m_primarydiam / 2.0);
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
		diams.push_back(m_primarydiam);		//use m_primarydiam for consistency

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

		const unsigned int trueval  = 1;
		const unsigned int falseval = 0;

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

		val = m_sum_necks;
        out.write((char*)&val, sizeof(val));

		val = m_primaryvol;
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

		val = m_Rg;
		out.write((char*)&val, sizeof(val));

        val = m_children_sintering;
        out.write((char*)&val, sizeof(val));

        val = m_avg_sinter;
        out.write((char*)&val, sizeof(val));

        val = m_sint_rate;
        out.write((char*)&val, sizeof(val));

        val = m_sint_time;
        out.write((char*)&val, sizeof(val));

		// frame orientation vectors
		// Note: serialization is not actually needed because the output
		// files are not currently produced in the postpocessing step
		val = m_frame_orient_z[0];
        out.write((char*)&val, sizeof(val));

        val = m_frame_orient_z[1];
        out.write((char*)&val, sizeof(val));

        val = m_frame_orient_z[2];
        out.write((char*)&val, sizeof(val));

		val = m_frame_orient_x[0];
        out.write((char*)&val, sizeof(val));

        val = m_frame_orient_x[1];
        out.write((char*)&val, sizeof(val));

        val = m_frame_orient_x[2];
        out.write((char*)&val, sizeof(val));

		// Output if primary is tracked
		if (m_tracked) {
			out.write((char*)&trueval, sizeof(trueval));
		}
		else {
			out.write((char*)&falseval, sizeof(falseval));
		}

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
		unsigned int val_unsigned(0);

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
        m_sum_necks = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_primaryvol = val;

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
		m_Rg = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_children_sintering = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_avg_sinter = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sint_rate = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sint_time = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_frame_orient_z[0] = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_frame_orient_z[1] = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_frame_orient_z[2] = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_frame_orient_x[0] = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_frame_orient_x[1] = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_frame_orient_x[2] = val;

		 // Read if primary is tracked.
        in.read(reinterpret_cast<char*>(&val_unsigned), sizeof(val_unsigned));
        if (val_int==1) {
            m_tracked = true;
        } else {
            m_tracked = false;
        }

        // Input base class.
        Primary::Deserialize(in, model);

    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, BinTreePrimary::DeserializePrimary).");
    }
}

//! Returns true if this node is a leaf (has no children).
bool BinTreePrimary::isLeaf(void) const
{
    return (m_leftchild == NULL) && (m_rightchild == NULL);
}

//! Returns the bounding-sphere centre.
const Coords::Vector &BinTreePrimary::boundSphCentre(void) const
{
    return m_cen_bsph;
}

//! Estimates the bounding sphere position and radius using
//! Ritter's method. ~5% larger than minimum bounding sphere
//! Ritter, J. (1990). An efficient bounding sphere, Graphics Gems 
//! (Andrew S. Glassner ed.), pp. 301-303. Academic Press, Boston
void BinTreePrimary::calcBoundSph(void)
{
	if ((m_leftchild != NULL) && (m_rightchild != NULL)) {

		//! Get list of primary coordinates
		vector<fvector> coords;
		this->GetPriCoords(coords);

		//! Find 3 pairs of points with the min and max x,y,z values
		fvector min_x = coords[1];	//! initialise with first point
		fvector max_x = min_x;
		fvector min_y = min_x;
		fvector max_y = min_x;
		fvector min_z = min_x;
		fvector max_z = min_x;
		for (int i = 1; i != coords.size(); ++i) {
			if (coords[i][0] < min_x[0]) min_x = coords[i];
			if (coords[i][0] > max_x[0]) max_x = coords[i];
			if (coords[i][1] < min_y[1]) min_y = coords[i];
			if (coords[i][1] > max_y[1]) max_y = coords[i];
			if (coords[i][2] < min_z[2]) min_z = coords[i];
			if (coords[i][2] > max_z[2]) max_z = coords[i];
		}

		//! Calculate separation (squared) between min and max 
		double dx = (max_x[0] - min_x[0]);
		double dy = (max_x[1] - min_x[1]);
		double dz = (max_x[2] - min_x[2]);
		double x_sep = dx*dx + dy*dy + dz*dz;
		dx = (max_y[0] - min_y[0]);
		dy = (max_y[1] - min_y[1]);
		dz = (max_y[2] - min_y[2]);
		double y_sep = dx*dx + dy*dy + dz*dz;
		dx = (max_z[0] - min_z[0]);
		dy = (max_z[1] - min_z[1]);
		dz = (max_z[2] - min_z[2]);
		double z_sep = dx*dx + dy*dy + dz*dz;

		//! find maximum separation
		//! points p_1 and p_2 are the points with maximum separation
		double max_sep = x_sep;
		fvector p_1 = min_x;
		fvector p_2 = max_x;
		if (y_sep > max_sep){
			max_sep = y_sep;
			p_1 = min_y;
			p_2 = max_y;
		}
		if (z_sep > max_sep){
			max_sep = z_sep;
			p_1 = min_z;
			p_2 = max_z;
		}

		//! Use p_1 and p_2 for initial guess at bounding sphere
		//! bounding sphere centre
		m_cen_bsph[0] = (p_1[0] + p_2[0]) / 2.0;
		m_cen_bsph[1] = (p_1[1] + p_2[1]) / 2.0;
		m_cen_bsph[2] = (p_1[2] + p_2[2]) / 2.0;
		//! radius of bounding sphere 
		//! including the primary radii
		setRadius((sqrt(max_sep)/2.0)+p_1[3]+p_2[3]);

		//! Make a second pass through list of primaries updating the sphere
		for (int i = 1; i != coords.size(); ++i) {
			
			//! calculate distance from bounding sphere centre and add primary radius
			dx = coords[i][0] - m_cen_bsph[0];
			dy = coords[i][1] - m_cen_bsph[1];
			dz = coords[i][2] - m_cen_bsph[2];
			double r_cen = sqrt(dx*dx + dy*dy + dz*dz); //!< Distance to primary centre
			double r_out = r_cen + coords[i][3]; //!< Distance to outer edge of primary

			//! If distance from the bounding sphere centre exceeds the 
			//! bounding sphere radius then update the bounding sphere:
			//! move centre by half the difference and increase the radius by half the difference
			if (r_out > m_r){
				
				//! half the distance from current bounding sphere centre to outer edge of primary
				double delta_r = (r_out - m_r) / 2.0; 

				//! Update the bounding sphere centre:
				//! the centre is translated along the vector joining 
				//! the old centre to the primary centre by half the difference
				m_cen_bsph[0] += delta_r * dx / r_cen;
				m_cen_bsph[1] += delta_r * dy / r_cen;
				m_cen_bsph[2] += delta_r * dz / r_cen;

				//! Set new radius
				setRadius(m_r + delta_r);
			}
		}
	}
}

//! Calculates the centre-of-mass using the left and right child node values.
void BinTreePrimary::calcCOM(void)
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

//! Put the bounding-sphere at the origin.
void BinTreePrimary::centreBoundSph(void)
{
    Translate(-m_cen_bsph[0], -m_cen_bsph[1], -m_cen_bsph[2]);
}

//! Put the centre-of-mass at the origin.
void BinTreePrimary::centreCOM(void)
{
    Translate(-m_cen_mass[0], -m_cen_mass[1], -m_cen_mass[2]);
}

/*!
 *  Randomly rotates the aggregate node and child structure about its centre of
 *  mass.
 *
 *  @param[in]    theta    Rotation about the pole. 
 *  @param[in]    V        Vector for performing the reflection.
 */
void BinTreePrimary::rotateCOM(double theta, fvector V)
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

/*!
 *  @brief Sets the radius of the bounding sphere.
 *
 *  @param[in]    r    Radius of bounding sphere.
 */
void BinTreePrimary::setRadius(double r)
{
    m_r  = r;
    m_r2 = r * r;
    m_r3 = m_r2 * m_r;
}

//! Returns the bounding sphere radius.
double BinTreePrimary::Radius(void) const
{
    return m_r;
}

/*!
 *  Transform the primary particle coordinates using the transformation matrix
 *  so as to rotate it.
 *:
 *  @param[in]    mat               Transformation matrix.
 *  @param[in]    PAHTracerMatch    Flag used to indicate whether the particle
 *                                  contains the PAH to be traced.
 */
void BinTreePrimary::transform(const Coords::Matrix &mat)
{
	//! Descend binary tree to the leaf nodes, i.e. single primary particles.
	if (m_leftchild != NULL)
		m_leftchild->transform(mat);
	if (m_rightchild != NULL)
		m_rightchild->transform(mat);

	//! Rotate centre-of-mass and bounding sphere coordinates.
	m_cen_mass = mat.Mult(m_cen_mass);
	m_cen_bsph = mat.Mult(m_cen_bsph);

	//! If the primary is tracked then rotate the frame orientation
	if (m_tracked == true) {
		m_frame_orient_z = mat.Mult(m_frame_orient_z);
		m_frame_orient_x = mat.Mult(m_frame_orient_x);
	}
}

/*!
 *  Translates (moves) the aggregate node and child structure by the given
 *  amounts along the cartesian axes.
 *
 *  @param[in]    dx    Distance to translate in the x-axis.
 *  @param[in]    dy    Distance to translate in the y-axis.
 *  @param[in]    dz    Distance to translate in the z-axis.
 */
void BinTreePrimary::Translate(double dx, double dy, double dz)
{
    //! Translate child branches.
    if (m_leftchild != NULL) m_leftchild->Translate(dx, dy, dz);
    if (m_rightchild != NULL) m_rightchild->Translate(dx, dy, dz);

    //! Translate bounding sphere centre.
    m_cen_bsph.Translate(dx, dy, dz);

    //! Translate centre-of-mass.
    m_cen_mass.Translate(dx, dy, dz);
}

//! Write the coordinates of the primaries in the particle pointed to by the
//! this pointer. Units of nm.
void BinTreePrimary::writePrimaryCoordinatesRadius(void)
{
    if (m_leftchild!=NULL)
        m_leftchild->writePrimaryCoordinatesRadius();

    if (m_rightchild!=NULL)
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

/*!
 *  Calculates unit vector between two set coordinates
 *
 *  @param[in]    x_i	Coordinates
 *  @param[in]    x_j	Coordinates
 *  @param[out]   unit vector
 */
Coords::Vector BinTreePrimary::UnitVector(Coords::Vector x_i, Coords::Vector x_j)
{
	Coords::Vector delta_x;
	double len_delta_x;

	//! Calculate difference x_j - x_i
	delta_x[0] = x_j[0] - x_i[0];
	delta_x[1] = x_j[1] - x_i[1];
	delta_x[2] = x_j[2] - x_i[2];

	//! Calculate the length of the vector
	len_delta_x = sqrt(delta_x[0]*delta_x[0]+delta_x[1]*delta_x[1]+delta_x[2]*delta_x[2]);

	//! Create unit vector
	delta_x[0] /=  len_delta_x;
	delta_x[1] /=  len_delta_x;
	delta_x[2] /=  len_delta_x;

	return delta_x;
}

/*!
 *  Calculates distance between two points
 *
 *  @param[in]    x_i	Coordinates
 *  @param[in]    x_j	Coordinates
 *  @param[out]   Separation
 */
double BinTreePrimary::Separation(Coords::Vector x_i, Coords::Vector x_j)
{
	Coords::Vector delta_x;
	double len_delta_x;

	//calculate difference x_j - x_i
	delta_x[0] = x_j[0] - x_i[0];
	delta_x[1] = x_j[1] - x_i[1];
	delta_x[2] = x_j[2] - x_i[2];

	//calculate the length of the vector
	len_delta_x = sqrt(delta_x[0]*delta_x[0]+delta_x[1]*delta_x[1]+delta_x[2]*delta_x[2]);

	return len_delta_x;
}


/*!
 *  Translates a primary particle along a unit vector
 *
 *  @param[in]    u			Unit vector (direction)
 *  @param[in]    delta_d	Distance to translate
 */
void BinTreePrimary::TranslatePrimary(Coords::Vector u, double delta_d)
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

/*!
 *  Translates neighbours of a primary particle along a unit vector
 *
 *  @param[in]    prim			Primary
 *  @param[in]    u				Unit vector direction of translation
 *  @param[in]    delta_d		Magnitude of translation
 *  @param[in]    prim_ignore	Primary to ignore
 */
void BinTreePrimary::TranslateNeighbours(BinTreePrimary *prim, Coords::Vector u, double delta_d, BinTreePrimary *prim_ignore)
{
	BinTreePrimary *neighbour = NULL;
	//! Check if a neighbour of prim but not prim_ignore
	if (m_parent->m_leftparticle == prim && m_parent->m_rightparticle != prim_ignore ) {
		//! right particle is a neighbour
		neighbour = m_parent->m_rightparticle;
		//! adjust its coordinates
		neighbour->TranslatePrimary(u,delta_d);
		//! adjust its neighbours except for prim
		neighbour->TranslateNeighbours(neighbour, u, delta_d, prim);
	} else if(m_parent->m_rightparticle == prim &&  m_parent->m_leftparticle != prim_ignore ) {
		//! left particle is a neighbour
		neighbour = m_parent->m_leftparticle;
		//! adjust its coordinates
		neighbour->TranslatePrimary(u,delta_d);
		//! adjust its neighbours except for prim
		neighbour->TranslateNeighbours(neighbour, u, delta_d, prim);
	}

	//! continue working up the binary tree
	if(m_parent->m_parent != NULL){
		m_parent->TranslateNeighbours(prim, u, delta_d, prim_ignore);
	}
}

// PARTICLE TRACKING FOR VIDEOS

//! Remove primary tracking flag
void BinTreePrimary::removeTracking()
{
	//set tracking flag to false
	m_tracked = false;

	//work down bintree structure
	if (m_leftchild != NULL) m_leftchild->removeTracking();
	if (m_rightchild != NULL) m_rightchild->removeTracking();
}

//! Set tracking flag
//  A single primary in the aggregate is tracked to centre the image frame
void BinTreePrimary::setTracking()
 { 
	 m_tracked = true;

	// work done binary tree to a primary
	if (m_leftchild != NULL) m_leftchild->setTracking();
}

//! Returns the frame position and orientation, and primary coordinates
//! Used by particle tracking for videos
void BinTreePrimary::GetFrameCoords(std::vector<fvector> &coords) const
{
	if (isLeaf()) { //this is a primary
		fvector c(10);
		c[0] = m_cen_mass[0];	//primary coordinates
		c[1] = m_cen_mass[1];
		c[2] = m_cen_mass[2];
		c[3] = m_r;				//primary radius
		if (m_tracked == true){	
			//the frame orientation vectors are output
			//for the single tracked primary
			//the primary coordinates define the centre of the frame
			c[4] = m_frame_orient_x[0];	
			c[5] = m_frame_orient_x[1];
			c[6] = m_frame_orient_x[2];
			c[7] = m_frame_orient_z[0];
			c[8] = m_frame_orient_z[1];
			c[9] = m_frame_orient_z[2];
		}
		else{
			c[4] = 0.0;
			c[5] = 0.0;
			c[6] = 0.0;
			c[7] = 0.0;
			c[8] = 0.0;
			c[9] = 0.0;
		}
		//get primary composition
		fvector comp = this->Composition();
		c.insert(c.end(), comp.begin(), comp.end());

		//add to coords vector
		coords.push_back(c);
	}
	else {	//this is a non-leaf node
		m_leftchild->GetFrameCoords(coords);
		m_rightchild->GetFrameCoords(coords);
	}
}

/*!
 *  Print primary particle details and connectivity
 *
 *  @param[in]    surface			Primary connectivity
 *  @param[in]    primary_diameter	Primary details
 *  @param[in]    k					Particle counter
 */
void BinTreePrimary::PrintPrimary(vector<fvector> &surface, vector<fvector> &primary_diameter, int k) const
{
	fvector node(10);	
	fvector primary;
	fvector comp;

	if ((m_leftchild==NULL) && (m_rightchild==NULL)){
		//if leaf then print diameter
		primary.push_back((double)k+1);
		primary.push_back(m_primarydiam);
		primary.push_back(m_diam);
		primary.push_back(m_primaryvol);
		primary.push_back(m_vol);
		primary.push_back(m_free_surf);

		//primary coordinates
		vector<fvector> coords;
		this->GetPriCoords(coords); //get primary coodinates
		primary.push_back(coords[0][0]);
		primary.push_back(coords[0][1]);
		primary.push_back(coords[0][2]);
		primary.push_back(coords[0][3]);

		//primary composition
		fvector comp = Primary::Composition();
		primary.insert(primary.end(),comp.begin(),comp.end());

		primary_diameter.push_back(primary);
		
		if (m_parent == NULL){	//connectivity information for single primary case
			node[0] = k+1;
			node[1] = m_numprimary;
			node[2] = 0.0;
			node[3] = 1.0;
			node[4] = 0.0;
			node[5] = 0.0;
			node[6] = m_primarydiam/2.0;
			node[7] = 0.0;
			//print pointer to this primary as integer (assume this a unique id to the primary)
			node[8] = reinterpret_cast<uintptr_t>(this);	
			node[9] = 0.0;

			surface.push_back(node);
		}
	} else {
		
		double r_i = m_leftparticle->m_primarydiam/2.0;
		double r_j = m_rightparticle->m_primarydiam/2.0;	
		double d_ij = m_distance_centreToCentre;

		double x_ij = (d_ij*d_ij - r_j*r_j + r_i*r_i)/(2.0*d_ij);
		double R_ij = sqrt(r_i*r_i - x_ij*x_ij);	//!neck radius

		//if non-leaf node then print node and continue down the tree
		node[0] = k+1;
		node[1] = m_numprimary;
		node[2] = m_children_surf;
		node[3] = m_children_sintering;
		node[4] = d_ij;
		node[5] = R_ij;
		node[6] = r_i;
		node[7] = r_j;
		//print pointer to left primary as integer (assume this a unique id to the primary)
		node[8] = reinterpret_cast<uintptr_t>(m_leftparticle);	
		//print pointer to right primary as integer (assume this a unique id to the primary)
		node[9] = reinterpret_cast<uintptr_t>(m_rightparticle);	
		
		surface.push_back(node);
		
		//continue down binary tree
		m_leftchild->PrintPrimary(surface, primary_diameter, k);
		m_rightchild->PrintPrimary(surface, primary_diameter, k);
	}
}
