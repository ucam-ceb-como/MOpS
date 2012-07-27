/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ParticleImage class declared in the
    swp_particle_image.h header file.

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

#include "swp_particle_image.h"
#include "swp_bintree_primary.h"
#include "swp_silica_primary.h"
#include "swp_surfvol_primary.h"
#include "swp_PAH_primary.h"
#include "string_functions.h"
#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <fstream>
#include <vector>
#include <stdexcept>
using namespace Sweep;
using namespace Sweep::Imaging;
using namespace std;
using namespace Strings;

const real ParticleImage::m_necking = 1.000;

// CONSTRUCTORS AND DESTRUCTORS.

// Default contructor.
ParticleImage::ParticleImage(void)
{
}

// Initialising constructor.
/*ParticleImage::ParticleImage(const Particle &sp, const ParticleModel &model)
{
    Construct(sp, model);
}*/

// Destructor.
ParticleImage::~ParticleImage(void)
{
}

// PARTICLE IMAGE DATA CONSTRUCTION.

/*!
 * @brief           Construct a particle image
 *
 * @param sp        Particle object
 * @param model     Particle model describing object
 */
void ParticleImage::Construct(const Particle &sp, const ParticleModel &model)
{
    // Initialise new RNG here as the calling functions do not have access to
    // the original one
    rng_type rng(size_t(100));

    // Clear the current image data structure.
    m_root.Clear();

    if (model.AggModel() == AggModels::Spherical_ID) {
        // Spherical particle model, just draw
        // a single sphere.

        m_root.Insert(sp.SphDiameter() * 0.5e9); // Convert to nm.

    } else if (model.AggModel() == AggModels::SurfVol_ID) {
        // Surface-volume model, construct a tree
        // of identical primaries, estimating the primary
        // count and diameter.

        const AggModels::SurfVolPrimary *svp = NULL;
        svp = dynamic_cast<const AggModels::SurfVolPrimary*>(sp.Primary());

        m_root.Clear();
        for (unsigned int i=0; i!=svp->PP_Count(); ++i) {
            // Convert to radius in nm
            m_root.Insert(svp->PP_Diameter() * 0.5e9);
        }
        calc_FM(m_root, *(&rng));

    } else if (model.AggModel() == AggModels::Bintree_ID) {
        // Binary Tree generic model

        const AggModels::BintreePrimary *p;
        p = dynamic_cast<const AggModels::BintreePrimary*>(sp.Primary());
        ConstructTree(p, *(&rng));
    } else if (model.AggModel() == AggModels::Silica_ID) {
        // Silica (binary tree like) model

        const AggModels::SilicaPrimary *p;
        p = dynamic_cast<const AggModels::SilicaPrimary*>(sp.Primary());
        ConstructTree(p, *(&rng));
    } else if (model.AggModel() == AggModels::PAH_KMC_ID) {
        // PAHPP (binary tree like) model

        const AggModels::PAHPrimary *p;
        p = dynamic_cast<const AggModels::PAHPrimary*>(sp.Primary());
        ConstructTree(p, *(&rng));
    } else {
        throw std::runtime_error("Unknown particle model. (ParticleImage::Construct)");
    }

    // Set to TEM-style projection
    //m_root.Project();
}


//! Generates a projection on the zx plane (set all y to 0)
void ParticleImage::Project()
{
    // Call to ImgNode
    m_root.Project();
}

/*!
 * @brief       Writes a 3dout format file
 *
 * This made use of some custom-made OpenGL renderer by Markus Sander
 * which has since been lost. (wjm34 27/07/2012)
 *
 * @param file  File output stream
 * @param x     ?
 * @param y     ?
 * @param z     ?
 */
void ParticleImage::Write3dout(std::ofstream &file, double x, double y, double z)
{
    if (file.good()) {
        string line;
        real val = 0.0;

        // vector of arrays to store primary coordinates.  First
        // 3 values are the cartesian coordinates, final value
        // is the primary radius.
        vector<fvector> coords;

        // Get the primary coordinates from the aggregate tree.
        m_root.GetPriCoords(coords);
		file.write(line.c_str(), line.length());
        //write the radius of gyration
        double Rg=RadiusofGyration();
        line = cstr(0) + " " + cstr(0) +
               " " + cstr(0)+"\n ";
        file.write(line.c_str(), line.length());
        line = cstr(Rg)+"\n";							//write radius of gyration
        file.write(line.c_str(), line.length());
        // Write the primaries to the 3dout file.
        for (unsigned int i=0; i!=coords.size(); ++i) {
            val  = coords[i][3] * m_necking;
            line = cstr(coords[i][0]+x) + " " + cstr(coords[i][1]+y) +
                   " " + cstr(coords[i][2]+z)+"\n ";
            file.write(line.c_str(), line.length());
            line = cstr(val)+"\n";							//write radius
            file.write(line.c_str(), line.length());
        }


    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleImage::Write3dout).");
    }
}

/*
void ParticleImage::LengthWidth(double &L, double &W)
{
    //align the particle along the z axis
	double xmax=0;
	double ymax=0;
	double xmin=0;
	double ymin=0;
    double zmax=0;
	double zmin=0;
    double dx=0,dy=0,dz=0;
    //the radius of outher spheres
	double yminrad=0;
	double xminrad=0;
	double ymaxrad=0;
	double xmaxrad=0;
    double zmaxrad=0;
    double zminrad=0;
	vector<fvector> coords;
    m_root.GetPriCoords(coords);
    double dist=0, distmax=0;
	for (unsigned int i=0; i!=coords.size(); ++i) {
        dist=sqrt(coords[i][0]*coords[i][0]+coords[i][1]*coords[i][1]+coords[i][2]*coords[i][2])+coords[i][3];
		if (dist>distmax) {
            xmax=coords[i][0];ymax=coords[i][1];zmax=coords[i][2];
            distmax=dist;
        }
	}
    //added to debug ms785
    if (m_root.m_left!=NULL)
    {
        double phi=atan2(ymax,xmax);
        double theta=(PI/2)-atan(zmax/sqrt((xmax*xmax)+(ymax*ymax)));

      //  m_root.Translate(-dx,-dy,-dz);
        m_root.RotateOrigin(0,-phi-PI/2);
        ofstream out;
       // out.open("particlebeforeshift.3d");
       // this->Write3dout(out,0,0,0);
       // out.close();
        m_root.RotateOrigin(-theta,0);

      //  out.open("particleaftershift.3d");
      //  this->Write3dout(out,0,0,0);
      //  out.close();
    }
    coords.clear();
    m_root.GetPriCoords(coords);
	xmax=coords[0][0];
	xmin=xmax;
	ymax=coords[0][1];
	ymin=ymax;
    zmax=coords[0][2];
	zmin=zmax;
	for (unsigned int i=0; i!=coords.size(); ++i) {
		if (xmax<=coords[i][0]+coords[i][3]) {xmax=coords[i][0]+coords[i][3];xmaxrad=coords[i][3];}
		if (ymax<=coords[i][1]+coords[i][3]) {ymax=coords[i][1]+coords[i][3];ymaxrad=coords[i][3];}
        if (zmax<=coords[i][2]+coords[i][3]) {zmax=coords[i][2]+coords[i][3];zmaxrad=coords[i][3];}
		if (xmin>=coords[i][0]-coords[i][3]) {xmin=coords[i][0]-coords[i][3];xminrad=coords[i][3];}
		if (ymin>=coords[i][1]-coords[i][3]) {ymin=coords[i][1]-coords[i][3];yminrad=coords[i][3];}
        if (zmin>=coords[i][2]-coords[i][3]) {zmin=coords[i][2]-coords[i][3];zminrad=coords[i][3];}
	}
    dx=abs(xmax-xmin);
    dy=abs(ymax-ymin);
    dz=abs(zmax-zmin);
	L=dz;

   // W=max(dy,dx);
    W=dx;

}*/

//! Calculates the radius of gyration of a particle
double ParticleImage::RadiusofGyration()
{
    double sum=0;
    double mass;
    double totalmass=0;
    double r2;
    double Rg;
    vector<fvector> coords;
    m_root.GetPriCoords(coords);
	for (unsigned int i=0; i!=coords.size(); ++i) {
        //mass is proportional to the cube of the radius
        mass=coords[i][3]*coords[i][3]*coords[i][3];
        r2=coords[i][0]*coords[i][0]+coords[i][1]*coords[i][1]+coords[i][2]*coords[i][2];
        sum+=mass*r2;
        totalmass+=mass;
	}
    Rg=sqrt(sum/totalmass);
    return Rg;
}

// RENDERING FUNCTIONS.

// Draws the particle image to a POVRAY file.
void ParticleImage::WritePOVRAY(std::ofstream &file)
{
    if (file.good()) {
        string line;
        real val = 0.0;

        // vector of arrays to store primary coordinates.  First
        // 3 values are the cartesian coordinates, final value
        // is the primary radius.
        vector<fvector> coords;

        // Write ParticleDiameter argument to POV file.
        val  = m_root.Radius() * 2.0;
        line = "#declare ParticleDiameter = " + cstr(val) + ";\n";
        file.write(line.c_str(), line.length());

        // Write aggregate opening declaration.
        line = "#declare MyParticle = blob {\n";
        file.write(line.c_str(), line.length());

        // Write threshold radius based on necking parameter.
        val  = max(pow(1.0 - (1.0/(m_necking*m_necking)), 2.0), 1.0e-4);
        line = "  threshold " + cstr(val) + "\n";
        file.write(line.c_str(), line.length());

        // Get the primary coordinates from the aggregate tree.
        m_root.GetPriCoords(coords);

        // Write the primaries to the POV-RAY file.
        for (unsigned int i=0; i!=coords.size(); ++i) {
            val  = coords[i][3] * m_necking;
            line = "sphere {<" + cstr(coords[i][0]) + ", " + cstr(coords[i][1]) +
                   ", " + cstr(coords[i][2]) + ">, " + cstr(val) + ", 1.0}\n";
            file.write(line.c_str(), line.length());
        }

        // Write closing brace for MyParticle declaration.
        line = "}\n";
        file.write(line.c_str(), line.length());

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleImage::WritePOVRAY).");
    }
}


// AGGREGATE SPHERE-TREE CONSTRUCTORS (FREE-MOLECULAR).


/*!
 * @brief       Generate the free-molecular structure of a particle
 *
 * Calculates the aggregate structure down from the given
 * node.  Assumes that the tree leaves have been initialised
 * with the correct radii, and recalculates their positions.
 *
 * @param node  Pointer to node of ImgNode tree
 * @param rng   Random number generator
 */
void ParticleImage::calc_FM(ImgNode &node, Sweep::rng_type &rng)
{
    ImgNode *target = node.m_left;
    ImgNode *bullet = node.m_right;

    if ((target != NULL) && (bullet != NULL)) {
        // Pass calculation down binary tree left & right branches.
        calc_FM(*target,  rng);
        calc_FM(*bullet,  rng);

        // The first part of the collision algorithm is to
        // randomly orientate both left and right aggregates.
        // They are both then placed so that their bounding
        // spheres are at the origin.

        // Rotate left node randomly about CoM.
        // Generate a random number on [0,1)-real-interval
        boost::uniform_01<rng_type&, real> uniformGenerator(rng);
        real phi1   = uniformGenerator() * 2.0 * PI;
        real theta1 = ((2.0*uniformGenerator())-1.0) * PI;
        target->RotateCOM(theta1, phi1);

        // Rotate right node randomly about CoM.
        real phi2   = uniformGenerator() * 2.0 * PI;
        real theta2 = ((2.0*uniformGenerator())-1.0) * PI;
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
            D[0] = ((2.0 * uniformGenerator()) - 1.0) * sumr;
            D[1] = ((2.0 * uniformGenerator()) - 1.0) * sumr;

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
            hit = minCollZ(*target, *bullet, D[0], D[1], D[2]);
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
/*!
 * @brief           Calculates the minimum collision distance
 *
 * Calculates the minimum collision distance between
 * a target and a bullet node by moving down the
 * binary tree.  If the nodes collide then returns
 * true, otherwise returns false.
 *
 * @param target    Target node
 * @param bullet    Bullet node
 * @param dx        ?
 * @param dy        ?
 * @param dz        ?
 * @return          Have the nodes collided?
 */
bool ParticleImage::minCollZ(const ImgNode &target,
                             const ImgNode &bullet,
                             real dx, real dy, real &dz)
{
    bool hit=false, hit1=false;
    real dz2=0.0, dz3=0.0, dz4=0.0;

    if (target.IsLeaf()) {
        // Target is a leaf
        if (bullet.IsLeaf()) {
            // Bullet is a leaf (both leaves).
           return calcCollZ(target.BoundSphCentre(), target.Radius(),
                            bullet.BoundSphCentre(), bullet.Radius(),
                            dx, dy, dz);
        } else {
            // Bullet is not a leaf, call sub-nodes.
            // Calculate minimum dz for the target and the bullet left subnode.
            hit1 = calcCollZ(target.BoundSphCentre(), target.Radius(),
                             bullet.m_left->BoundSphCentre(), bullet.m_left->Radius(),
                             dx, dy, dz);
            if (hit1) hit = minCollZ(target, *bullet.m_left, dx, dy, dz);
            // Calculate minimum dz for the target and the bullet right subnode.
            hit1 = calcCollZ(target.BoundSphCentre(), target.Radius(),
                             bullet.m_right->BoundSphCentre(), bullet.m_right->Radius(),
                             dx, dy, dz2);
            if (hit1) hit = minCollZ(target, *bullet.m_right, dx, dy, dz2) || hit;
            // Return minimum dz.
            dz = min(dz, dz2);
            return hit;
        }
    } else {
        // Target is not a leaf.
        if (bullet.IsLeaf()) {
            // Bullet is a leaf, call target sub-nodes..
            // Calculate minimum dz for the target left subnode and the bullet.
            hit1 = calcCollZ(target.m_left->BoundSphCentre(), target.m_left->Radius(),
                             bullet.BoundSphCentre(), bullet.Radius(),
                             dx, dy, dz);
            if (hit1) hit = minCollZ(*target.m_left, bullet, dx, dy, dz);
            // Calculate minimum dz for the target right subnode and the bullet.
            hit1 = calcCollZ(target.m_right->BoundSphCentre(), target.m_right->Radius(),
                             bullet.BoundSphCentre(), bullet.Radius(),
                             dx, dy, dz2);
            if (hit1) hit = minCollZ(*target.m_right, bullet, dx, dy, dz2) || hit;
            // Return minimum dz.
            dz = min(dz, dz2);
            return hit;
        } else {
            // Bullet is not a leaf (neither is a leaf), check all left/right
            // collision combinations.
            // Target left and bullet left.
            hit1 = calcCollZ(target.m_left->BoundSphCentre(), target.m_left->Radius(),
                             bullet.m_left->BoundSphCentre(), bullet.m_left->Radius(),
                             dx, dy, dz);
            if (hit1) hit = minCollZ(*target.m_left, *bullet.m_left, dx, dy, dz);
            // Target left and bullet right.
            hit1 = calcCollZ(target.m_left->BoundSphCentre(), target.m_left->Radius(),
                             bullet.m_right->BoundSphCentre(), bullet.m_right->Radius(),
                             dx, dy, dz2);
            if (hit1) hit = minCollZ(*target.m_left, *bullet.m_right, dx, dy, dz2) || hit;
            // Target right and bullet left.
            hit1 = calcCollZ(target.m_right->BoundSphCentre(), target.m_right->Radius(),
                             bullet.m_left->BoundSphCentre(), bullet.m_left->Radius(),
                             dx, dy, dz3);
            if (hit1) hit = minCollZ(*target.m_right, *bullet.m_left, dx, dy, dz3) || hit;
            // Target right and bullet right.
            hit1 = calcCollZ(target.m_right->BoundSphCentre(), target.m_right->Radius(),
                             bullet.m_right->BoundSphCentre(), bullet.m_right->Radius(),
                             dx, dy, dz4);
            if (hit1) hit = minCollZ(*target.m_right, *bullet.m_right, dx, dy, dz4) || hit;
            // Returns minimum dz.
            dz = min(min(dz, dz2), min(dz3, dz4));
            return hit;
        }
    }
}

/*!
 * @brief       Calculates the z-displacement of a sphere
 *
 * Calculates the z-displacement of a bullet sphere for a +ve
 * collision with a target sphere.  Returns true if the
 * spheres collide, otherwise false.
 *
 * @param p1    Coordinates of sphere 1
 * @param r1    Radius of sphere 1
 * @param p2    Coordinates of sphere 2
 * @param r2    Radius of sphere 2
 * @param dx    ?
 * @param dy    ?
 * @param dz    ?
 * @return      Have the nodes collided?
 */
bool ParticleImage::calcCollZ(const Coords::Vector &p1, real r1,
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
