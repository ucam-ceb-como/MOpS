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
#include "rng.h"
#include "string_functions.h"
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
ParticleImage::ParticleImage(const Particle &sp)
{
    Construct(sp);
}

// Destructor.
ParticleImage::~ParticleImage(void)
{
}

// PARTICLE IMAGE DATA CONSTRUCTION.

// Constructs the particle image from the given particle.
void ParticleImage::Construct(const Particle &sp)
{
    switch(m_creg) {
        case FreeMol:
            // Free-molecular is the default regime.
        default:
            constructAgg_FM(sp);
    }
}

// Constructs a random particle image.
void ParticleImage::ConstructRandom(real minrad, real maxrad, unsigned int n)
{
    // Clear the current image data structure.
    m_root.Clear();

    // Generate the random primaries.
    for (unsigned int i=0; i!=n; ++i) {
        real r = minrad + rnd()*(maxrad-minrad);
        m_root.Insert(r);
    }

    // Calculate aggregate structure
    switch(m_creg) {
        case FreeMol:
            // Free-molecular is the default regime.
        default:
            calc_FM(m_root);
    }
}




void ParticleImage::Write3dout(std::ofstream &file)
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

        // Write the primaries to the 3dout file.
        for (unsigned int i=0; i!=coords.size(); ++i) {
            val  = coords[i][3] * m_necking;
            line = cstr(coords[i][0]) + " " + cstr(coords[i][1]) + 
                   " " + cstr(coords[i][2])+"\n";
            file.write(line.c_str(), line.length());
            line = cstr(val)+"\n";							//write diameter
            file.write(line.c_str(), line.length());
        }

       
        line = "\n";
        file.write(line.c_str(), line.length());

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleImage::Write3dout).");
    }
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

        // Write include command for tem_single.pov, which defines TEM style
        // for a single particle.
        line = "#include \"particle.pov\"\n";
        file.write(line.c_str(), line.length());
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleImage::WritePOVRAY).");
    }
}


// AGGREGATE SPHERE-TREE CONSTRUCTORS (FREE-MOLECULAR).

// Constructs a PNode sphere-tree aggregate from the given 
// particle using free-molecular collision dynamics.
void ParticleImage::constructAgg_FM(const Particle &sp)
{
    // Clear the current image data structure.
    m_root.Clear();

    // Need to determine the type of particle this is.
    if (sp.ParticleModel()->UseSubPartTree()) {
        // This particle uses the sub-particle tree, so 
        // we can usethat to construct the sphere-tree for
        // image output.
        // TODO:  Complete sub-particle tree TEM output.
		constructSubParttree(&sp);

		cout <<"use subparttree to create the image tree";

    } else {
        const AggModels::SurfVolPrimary *svp = NULL;
        const AggModels::PriPartPrimary *ppp = NULL;

        switch(sp.ParticleModel()->AggModel()) {
            case AggModels::Spherical_ID:
                // Spherical particle model, just draw
                // a single sphere.
                m_root.Insert(sp.SphDiameter() * 0.5e9); // Convert to nm.
                return;
            case AggModels::SurfVol_ID:
                // Surface-volume model, construct a tree
                // of identical primaries, estimating the primary
                // count and diameter.
                svp = dynamic_cast<const AggModels::SurfVolPrimary*>(sp.Primary());
                if (svp != NULL) uniformAgg_FM(svp->PP_Count(), svp->PP_Diameter());
                break;
            case AggModels::PriPartList_ID:
                // Primary-particle list is known.  Use that
                // to construct the output sphere-tree.
                ppp = dynamic_cast<const AggModels::PriPartPrimary*>(sp.Primary());
                if (ppp != NULL) constructAgg_FM(*ppp);
                break;
        }
    }
}
void ParticleImage::constructSubParttree(const SubParticle *sp)
{
	//Copy the subparticle tree into the img tree
    m_root.CopySPT(sp);
    // Use the free-molecular regime to calculate the
    // aggregate structure.
    calc_FM(m_root);
    m_root.CentreCOM();
}

// Constructs a PNode sphere-tree aggregate from the given 
// pri-part list primary using free-molecular collision dynamics.
void ParticleImage::constructAgg_FM(const AggModels::PriPartPrimary &pri)
{
    // Get number of primary particles.
    unsigned int n = pri.PriCount();

    // Create a vector of primary radii.
    fvector radii;
    for (unsigned int i=0; i!=n; ++i) {
        radii.push_back(pri.PriDiameter(i)*0.5e9); // Convert to nm.
    }

    // Randomly add the primaries to the image aggregate tree.
    m_root.Clear();
    while (radii.size() > 0) {
        int j = irnd(0, radii.size()-1);
        m_root.Insert(radii[j]);
        radii.erase(radii.begin()+j);
    }

    // Use the free-molecular regime to calculate the
    // aggregate structure.
    calc_FM(m_root);
    m_root.CentreCOM();
}

// Constructs a PNode sphere-tree aggregate with uniform 
// primaries (equal diameter).  The diameter and primary
// count are passed as arguments.
void ParticleImage::uniformAgg_FM(unsigned int n, real d)
{
    real r = d * 0.5e9; // Convert to nm.

    // Add n identical primaries to the image data structure.
    m_root.Clear();
    for (unsigned int i=0; i!=n; ++i) {
        m_root.Insert(r);
    }

    // Use the free-molecular regime to calculate the
    // aggregate structure.
    calc_FM(m_root);
}

// Calculates the aggregate structure down from the given
// node.  Assumes that the tree leaves have been initialised
// with the correct radii, and recalculates their positions.
void ParticleImage::calc_FM(ImgNode &node)
{
    ImgNode *target = node.m_left;
    ImgNode *bullet = node.m_right;

    if ((target != NULL) && (bullet != NULL)) {
        // Pass calculation down binary tree left & right branches.
        calc_FM(*target);
        calc_FM(*bullet);

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

// Calculates the z-displacement of a bullet sphere for a +ve
// collision with a target sphere.  Returns true if the
// spheres collide, otherwise false.
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
