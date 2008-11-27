/*
  Author(s):      Matthew Celnik (msc37) and Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik and Markus Sander.

  File purpose:
    Implementation of the EnsembleImage class declared in the
    swp_ensemble_image.h header file.

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

#include "swp_ensemble_image.h"
#include "swp_particle_image.h"
#include "swp_imgnode.h"
#include "rng.h"
#include "string_functions.h"
#include <fstream>
#include <vector>
#include <stdexcept>
#include "swp_mechanism.h"
using namespace Sweep;
using namespace std;
using namespace Strings;
using namespace Sweep::Imaging;



// Default contructor.
EnsembleImage::EnsembleImage(void)
{
}


// Destructor.
EnsembleImage::~EnsembleImage(void)
{
}
/*
void EnsembleImage::PrintEnsemble(Cell &sys,std::ofstream &file)
{	Ensemble::iterator i;
    ParticleCache::PropID id=ParticleCache::iD;
	m_root.Clear();
	for (i=sys.Particles().begin(); i!=sys.Particles().end(); ++i) {
		if ((*(*i)).Property(id)>5e-9)
		m_root.Insert((*(*i)).Property(id),(*i));						//add 1e-8 to the radius to let the particles not collide
	}
    calc_FM_ensemble(m_root);
    m_root.CentreCOM();
	Write3dout(file,&m_root);
	
}


void EnsembleImage::Write3dout(std::ofstream &file, EnsembleImgNode *curr_node)
{
    if (file.good()) {
        string line;
		if (curr_node->IsLeaf() && curr_node->sp!=NULL) 
		{ 		
			Sweep::Imaging::ParticleImage *img= new Sweep::Imaging::ParticleImage;
			img->Construct(*curr_node->sp);
			img->Write3dout(file,1e9*curr_node->m_cen_mass[0],1e9*curr_node->m_cen_mass[1],1e9*curr_node->m_cen_mass[2]);
			delete img;
		}
		else 
		{
			Write3dout(file, curr_node->m_left);
			Write3dout(file, curr_node->m_right);
		}


    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, EnsembleImage::Write3dout).");
    }
}
*/



void EnsembleImage::PrintEnsemble(Cell &sys,std::ofstream &file, real shiftz)
{
	
	bool treedimage=false;
	Ensemble::iterator i;
    ParticleCache::PropID id=ParticleCache::iV;
	double x,y,z;
	double maxd=0;
	double boxlength;
	int numsubpart=0;
	x=0;y=0;z=0;
	m_root.Clear();
	for (i=sys.Particles().begin(); i!=sys.Particles().end(); ++i) {
		maxd=max((*(*i)).Property(id),maxd);
		numsubpart++;
		if (numsubpart>100) break;
	}
	if (treedimage==true)
	{
		maxd=1e9*6*pow(maxd,0.3333333);
		boxlength=pow(numsubpart,0.33333333);
		boxlength*=maxd;
	}
	else 
	{
		maxd=1e9*6*pow(maxd,0.3333333);
		maxd=200;
		boxlength=pow(numsubpart,0.5);
		boxlength*=maxd;
	}

	shiftz=shiftz*boxlength;
	numsubpart=0;
	for (i=sys.Particles().begin(); i!=sys.Particles().end(); ++i) {
		numsubpart++;
		if(numsubpart<100)
		{ 
			if (x<boxlength && y<boxlength)
			{
				x=x+maxd; 		
			}
			if (x>boxlength) 
			{
				x=0;
				y=y+maxd;
			}
			if (y>boxlength) 
			{
				x=0;
				y=0;
				z=z+maxd;
			}
			Sweep::Imaging::ParticleImage img;
			img.Construct(*(*i));
			//img.Write3dout(file,x+(rnd()-0.5)*(boxlength-(*(*i)).Property(id)*1e9),y+(rnd()-0.5)*(boxlength-(*(*i)).Property(id)*1e9),z+(rnd()-0.5)*(boxlength-(*(*i)).Property(id)*1e9));
			img.Write3dout(file,(x-boxlength/2)+(rnd()-0.5)*boxlength,y-boxlength/2+(rnd()-0.5)*boxlength,z+shiftz);
		
		}
	}
	
}




// Calculates the aggregate structure down from the given
// node.  Assumes that the tree leaves have been initialised
// with the correct radii, and recalculates their positions.
void EnsembleImage::calc_FM_ensemble(EnsembleImgNode &node)
{
    EnsembleImgNode *target = node.m_left;
    EnsembleImgNode *bullet = node.m_right;

    if ((target != NULL) && (bullet != NULL)) {
        // Pass calculation down binary tree left & right branches.
        calc_FM_ensemble(*target);
        calc_FM_ensemble(*bullet);

        // The first part of the collision algorithm is to
        // randomly orientate both left and right aggregates.
        // They are both then placed so that their bounding
        // spheres are at the origin.

        // Rotate left node randomly about CoM.
        real phi1   = rnd() * 2.0 * PI;
        real theta1 = ((2.0*rnd())-1.0) * PI;
        target->RotateCOM(theta1, phi1);
        
        // Rotate right node randomly about CoM.
        real phi2   = rnd() * 2.0 * PI;
        real theta2 = ((2.0*rnd())-1.0) * PI;
        bullet->RotateCOM(theta2, phi2);

        // Move both spheres so that the bounding spheres
        // sit at the origin.
        target->CentreBoundSph();
        bullet->CentreBoundSph();

        // Perform the collision of the left and right nodes.
        // This ma require several iterations if the chosen
        // x-y displacement means that the aggregates cannot
        // collide in the z-direction.
        Coords::Vector D;
        real sumr=0.0;
        bool hit = false;
        while (!hit) {
            // Need to reset target and bullet here, in case
            // the have been changed by the tree traversal
            // code below.
            target = node.m_left;
            bullet = node.m_right;
            sumr = target->Radius() + bullet->Radius();

            // Create a random displacement of the bullet node
            // in the x-y plane.  The displacement is never
            // greater than the sum of the radii, therefore they
            // should always touch.
            D[0] = ((2.0 * rnd()) - 1.0) * sumr;
            D[1] = ((2.0 * rnd()) - 1.0) * sumr;

            //// Calculate the z-position for the collision of the
            //// target and bullet.  We do this in case both the
            //// target and bullet are leaf nodes, and hence the
            //// binary tree traversal won't happen.
           // hit = calcCollZ(target->m_cen_mass, target->m_r,
           //                 bullet->m_cen_mass, bullet->m_r,
           //                 D[0], D[1], dz1);
           // if (!hit) continue; // Should never happen.
           // D[2] = target->m_cen_bsph[2] + dz1;

            // The next code determines the displacement along the z-axis
            // required for the target and bullet aggregates to touch.  This
            // requires falling down the tree progressively recalculating
            // the nearest nodes at each level, until the leaf nodes
            // are reached.
            // This next code calculates the minimum distance between the
            // target's children and bullet's children, or the target or
            // bullet if they have no children.  The two children with
            // the smallest separation are chosen as the next target
            // and bullet.
            hit = minCollZ_ensemble(*target, *bullet, D[0], D[1], D[2]);
        }

        // We have a new location for the bullet (right node), so move it.
        node.m_right->Translate(D[0], D[1], D[2]);

        // Calculate properties of this node.
        node.CalcBoundSph();
        node.CalcCOM();
        node.CentreBoundSph();
    }
}

// Calculates the minimum collision distance between
// a target and a bullet node by moving down the
// binary tree.  If the nodes collide then returns
// true, otherwise returns false.
bool EnsembleImage::minCollZ_ensemble(const EnsembleImgNode &target, 
                             const EnsembleImgNode &bullet, 
                             real dx, real dy, real &dz)
{
    bool hit=false, hit1=false;
    real dz2=0.0, dz3=0.0, dz4=0.0;

    if (target.IsLeaf()) {
        // Target is a leaf
        if (bullet.IsLeaf()) {
            // Bullet is a leaf (both leaves).
           return calcCollZ_ensemble(target.BoundSphCentre(), target.Radius(),
                            bullet.BoundSphCentre(), bullet.Radius(), 
                            dx, dy, dz);
        } else {
            // Bullet is not a leaf, call sub-nodes.
            // Calculate minimum dz for the target and the bullet left subnode.
            hit1 = calcCollZ_ensemble(target.BoundSphCentre(), target.Radius(),
                             bullet.m_left->BoundSphCentre(), bullet.m_left->Radius(), 
                             dx, dy, dz);
            if (hit1) hit = minCollZ_ensemble(target, *bullet.m_left, dx, dy, dz);
            // Calculate minimum dz for the target and the bullet right subnode.
            hit1 = calcCollZ_ensemble(target.BoundSphCentre(), target.Radius(),
                             bullet.m_right->BoundSphCentre(), bullet.m_right->Radius(), 
                             dx, dy, dz2);
            if (hit1) hit = minCollZ_ensemble(target, *bullet.m_right, dx, dy, dz2) || hit;
            // Return minimum dz.
            dz = min(dz, dz2);
            return hit;
        }
    } else {
        // Target is not a leaf.
        if (bullet.IsLeaf()) {
            // Bullet is a leaf, call target sub-nodes..
            // Calculate minimum dz for the target left subnode and the bullet.
            hit1 = calcCollZ_ensemble(target.m_left->BoundSphCentre(), target.m_left->Radius(),
                             bullet.BoundSphCentre(), bullet.Radius(), 
                             dx, dy, dz);
            if (hit1) hit = minCollZ_ensemble(*target.m_left, bullet, dx, dy, dz);
            // Calculate minimum dz for the target right subnode and the bullet.
            hit1 = calcCollZ_ensemble(target.m_right->BoundSphCentre(), target.m_right->Radius(),
                             bullet.BoundSphCentre(), bullet.Radius(), 
                             dx, dy, dz2);
            if (hit1) hit = minCollZ_ensemble(*target.m_right, bullet, dx, dy, dz2) || hit;
            // Return minimum dz.
            dz = min(dz, dz2);
            return hit;
        } else {
            // Bullet is not a leaf (neither is a leaf), check all left/right
            // collision combinations.
            // Target left and bullet left.
            hit1 = calcCollZ_ensemble(target.m_left->BoundSphCentre(), target.m_left->Radius(),
                             bullet.m_left->BoundSphCentre(), bullet.m_left->Radius(), 
                             dx, dy, dz);
            if (hit1) hit = minCollZ_ensemble(*target.m_left, *bullet.m_left, dx, dy, dz);
            // Target left and bullet right.
            hit1 = calcCollZ_ensemble(target.m_left->BoundSphCentre(), target.m_left->Radius(),
                             bullet.m_right->BoundSphCentre(), bullet.m_right->Radius(), 
                             dx, dy, dz2);
            if (hit1) hit = minCollZ_ensemble(*target.m_left, *bullet.m_right, dx, dy, dz2) || hit;
            // Target right and bullet left.
            hit1 = calcCollZ_ensemble(target.m_right->BoundSphCentre(), target.m_right->Radius(),
                             bullet.m_left->BoundSphCentre(), bullet.m_left->Radius(), 
                             dx, dy, dz3);
            if (hit1) hit = minCollZ_ensemble(*target.m_right, *bullet.m_left, dx, dy, dz3) || hit;
            // Target right and bullet right.
            hit1 = calcCollZ_ensemble(target.m_right->BoundSphCentre(), target.m_right->Radius(),
                             bullet.m_right->BoundSphCentre(), bullet.m_right->Radius(), 
                             dx, dy, dz4);
            if (hit1) hit = minCollZ_ensemble(*target.m_right, *bullet.m_right, dx, dy, dz4) || hit;
            // Returns minimum dz.
            dz = min(min(dz, dz2), min(dz3, dz4));
            return hit;
        }
    }
}

// Calculates the z-displacement of a bullet sphere for a +ve
// collision with a target sphere.  Returns true if the
// spheres collide, otherwise false.
bool EnsembleImage::calcCollZ_ensemble(const Coords::Vector &p1, real r1, 
                              const Coords::Vector &p2, real r2, 
                              real dx, real dy, real &dz)
{
    // Calculate the square of the sum of the radii.
    real sumrsqr = r1 + r2; sumrsqr *= sumrsqr;

    // Calculate dx, dy and dz.  Remember to include
    // argument contributions.
    real xdev = p2[0] - p1[0] + dx;
    real ydev = p2[1] - p1[1] + dy;
    real zdev = p2[2] - p1[2];

    // Calculate dx, dy and dz squared.
    real dxsqr = xdev * xdev;
    real dysqr = ydev * ydev;
    real dzsqr = zdev * zdev;

    // Calculate quadratic terms.
    real b = 2.0 * zdev;
    real c = dxsqr + dysqr + dzsqr - sumrsqr;

    // Calculate determinant.
    real det = (b*b) - (4.0*c);

    if (det >= 0.0) {
        // Spheres intersect.
        dz = - 0.5 * (b + sqrt(det));
        return true;
    } else {
        // Spheres do not intersect.
        dz = 1.0e10; // A large number.
        return false;
    }
}
