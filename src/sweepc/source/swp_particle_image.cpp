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

const double ParticleImage::m_necking = 1.000;

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

        //! At the moment tracking of the distance between the centres of
        //! primary particles does not apply to the surface-volume model.
        calc_FM(m_root, *(&rng), false);

    } else if (model.AggModel() == AggModels::BinTree_ID ||
            model.AggModel() == AggModels::BinTreeSilica_ID) {
        // Binary Tree generic model

        const AggModels::BinTreePrimary *p;
        p = dynamic_cast<const AggModels::BinTreePrimary*>(sp.Primary());
        ConstructTree(p, *(&rng), model.getTrackPrimaryCoordinates());

    } else if (model.AggModel() == AggModels::PAH_KMC_ID) {
        // PAHPP (binary tree like) model

        const AggModels::PAHPrimary *p;
        p = dynamic_cast<const AggModels::PAHPrimary*>(sp.Primary());
        ConstructTree(p, *(&rng), model.getTrackPrimaryCoordinates());
    } else {
        throw std::runtime_error("Unknown particle model. (ParticleImage::Construct)");
    }

    // Set to TEM-style projection
    //m_root.Project();
}

/*!
 *  This function is like a limited assignment operator, except that the
 *  children are not copied and the pointers to the particles may need
 *  adjusting after this method has finished.
 *
 *  @param[in,out] node   Pointer to node of ImgNode tree.
 *  @param[in]     source Pointer to the primary to be copied.
 */
template <class ParticleClass>
void ParticleImage::CopyParts(ImgNode &node, const ParticleClass *source)
{
	//! Since m_cen_bshp and m_cen_mass are vectors do unit conversion in
    //! function.
    node.setBoundSph(source->m_cen_bsph);
    node.setCOM(source->m_cen_mass);

	//! Units of nm.
    node.setRadius(source->m_primarydiam*0.5e9);
    node.setDistance(source->m_distance_centreToCentre*1.0e9);
}

/*!
 *  Recursively copy the tree for non-leaf nodes.
 *
 *  @param[in,out] node   Pointer to node of ImgNode tree.
 *  @param[in]     source Pointer to the primary to be copied.
 */
template <class ParticleClass>
void ParticleImage::CopyTree(ImgNode &node, const ParticleClass *source)
{
    //! Create the new left and right children with nothing in them.
    node.m_leftchild = new ImgNode();
    node.m_rightchild = new ImgNode();

    //! Copy the properties such as the volume, surface area and list of
    //! constituent PAH molecules.
    CopyParts(*node.m_leftchild, source->m_leftchild);
    CopyParts(*node.m_rightchild, source->m_rightchild);

    // Set the pointers to the parent.
    node.m_leftchild->m_parent=&node;
    node.m_rightchild->m_parent=&node;

    //! The left and right particle are set further down in UpdateAllPointers.
    //! These are the pointers that specify which primary particles touch each
    //! other in the aggregate structure.
    node.m_leftparticle=NULL;
    node.m_rightparticle=NULL;

    //! Recurse to copy the subtrees.
    if (source->m_leftchild->m_leftchild!=NULL)
        CopyTree(*node.m_leftchild, source->m_leftchild);

    if (source->m_rightchild->m_leftchild!=NULL)
        CopyTree(*node.m_rightchild, source->m_rightchild);

    //! Set the leftparticle and rightparticle.
    UpdateAllPointers(node, source);
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
        double val = 0.0;

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

/*!
 *  Each node contains two pointers (m_leftparticle and m_rightparticle)
 *  to primary particles that are connected by this node.
 *  This function is used when the entire particle tree is duplicated.
 *  It sets the pointers in the copied node, that the connectivity
 *  of the primary particles in this node is the same as in the original node.
 *
 *  @todo give this method a more accurate name.
 *
 *  @param[in,out] node     Pointer to node of ImgNode tree.
 *  @param[in]     original Pointer to the primary to be copied.
 */
template <class ParticleClass>
void ParticleImage::UpdateAllPointers(ImgNode &node, const ParticleClass *original)
{
    //! The primary has no children => there are no left and right particles.
    if (original->m_leftchild == NULL) {
        //! Since this is not a connecting node it does not have left and
        //!right particles.
        node.m_leftparticle=NULL;
        node.m_rightparticle=NULL;
    } else {
        //! Find the route to m_leftparticle in the original tree.
        std::stack<bool> route = recordPath(original->m_leftparticle, original);

        //! Now follow the same route down the new tree to find the new left particle.
        node.m_leftparticle = descendPath(&node, route);

        //! Find the route to m_rightparticle in the original tree.
        route = recordPath(original->m_rightparticle, original);

        //! Now follow the same route down the new tree to find the new right particle.
        node.m_rightparticle = descendPath(&node, route);
    }
}

/*!
 *  This is a helper function for UpdateAllPointers.
 *  It climbs up the tree from bottom to top recording a route
 *  suitable for use in call to @see descendPath.
 *
 *  @param[in] bottom Tree node from which to start climbing.
 *  @param[in] top    Tree node at which to stop climbing.
 *
 *  @pre top must be above bottom in a tree.
 *
 *  @return Stack that can be used to descend the same path by moving to the
 *          left child each time the top of the stack is true.
 */
template <class ParticleClass>
std::stack<bool> ParticleImage::recordPath(const ParticleClass* bottom, const ParticleClass* const top) {
    std::stack<bool> wasLeftChild;

    while(bottom != top) {
        //! check whether bottom was a left child of its parent.
        wasLeftChild.push(bottom == bottom->m_parent->m_leftchild);

        //! Climb one level up the tree.
        bottom = bottom->m_parent;
    }
    return wasLeftChild;
}

/*!
 *  @param[in]     here           Point in tree from which to start descent.
 *  @param[in,out] takeLeftBranch Instructions for which child to move to at each level.
 *
 *  @return The node at the bottom of the path.
 *
 *  @pre  here must be a node of tree in which takeLeftBranch is a valid path.
 *  @post takeLeftBranch.empty() == true.
 */
template <class ParticleClass>
ParticleClass* ParticleImage::descendPath(ParticleClass *here, std::stack<bool> &takeLeftBranch) {
    while(!takeLeftBranch.empty()) {
        //! Move one step down the tree in the instructed direction.
        if(takeLeftBranch.top())
            here = here->m_leftchild;
        else
            here = here->m_rightchild;

        //! This instuction has now been processed.
        takeLeftBranch.pop();
    }
    return here;
}

// RENDERING FUNCTIONS.

// Draws the particle image to a POVRAY file.
void ParticleImage::WritePOVRAY(std::ofstream &file)
{
    if (file.good()) {
        string line;
        double val = 0.0;

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
 *  @brief Generate the free-molecular structure of a particle.
 *
 *  Calculates the aggregate structure down from the given
 *  node. Assumes that the tree leaves have been initialised
 *  with the correct radii, and recalculates their positions.
 *
 *  @param[in,out] node                   Pointer to node of ImgNode tree.
 *  @param[in]     rng                    Random number generator.
 *  @param[in]     trackPrimaryCoordinates Flag used to indicate whether to track primary coordinates.
 */
void ParticleImage::calc_FM(ImgNode &node, Sweep::rng_type &rng, const bool trackPrimaryCoordinates)
{
    ImgNode *target = node.m_leftchild;
    ImgNode *bullet = node.m_rightchild;

    if ((target != NULL) && (bullet != NULL)) {
        //! Pass calculation down binary tree left & right branches.
        calc_FM(*target, rng, trackPrimaryCoordinates);
        calc_FM(*bullet, rng, trackPrimaryCoordinates);

        //! The first part of the collision algorithm is to
        //! randomly orientate both left and right aggregates.
        //! They are both then placed so that their bounding
        //! spheres are at the origin.

        //! Rotate left node randomly about CoM.
        //! Generate a random number on [0,1)-double-interval.
        boost::uniform_01<rng_type&, double> uniformGenerator(rng);
        double phi1   = uniformGenerator() * 2.0 * PI;
        double theta1 = ((2.0*uniformGenerator())-1.0) * PI;
        target->RotateCOM(theta1, phi1);

        //! Rotate right node randomly about CoM.
        double phi2   = uniformGenerator() * 2.0 * PI;
        double theta2 = ((2.0*uniformGenerator())-1.0) * PI;
        bullet->RotateCOM(theta2, phi2);

        //! Move both spheres so that the bounding spheres
        //! sit at the origin.
        target->CentreBoundSph();
        bullet->CentreBoundSph();

        //! Perform the collision of the left and right nodes.
        //! This may require several iterations if the chosen
        //! x-y displacement means that the aggregates cannot
        //! collide in the z-direction.
        Coords::Vector D;
        double sumr=0.0;
        bool hit = false;
        while (!hit) {
            //! Need to reset target and bullet here, in case
            //! they have been changed by the tree traversal
            //! code below.
            target = node.m_leftchild;
            bullet = node.m_rightchild;

            if (!trackPrimaryCoordinates) {
                sumr = target->Radius() + bullet->Radius();
            } else {
                sumr = node.m_distance_centreToCentre;
            }

            //! Create a random displacement of the bullet node
            //! in the x-y plane.  The displacement is never
            //! greater than the sum of the radii, therefore they
            //! should always touch.
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

            if (!trackPrimaryCoordinates) {
                //! The next code determines the displacement along the z-axis
                //! required for the target and bullet aggregates to touch.  This
                //! requires falling down the tree progressively recalculating
                //! the nearest nodes at each level, until the leaf nodes
                //! are reached.
                //! This next code calculates the minimum distance between the
                //! target's children and bullet's children, or the target or
                //! bullet if they have no children.  The two children with
                //! the smallest separation are chosen as the next target
                //! and bullet.
                hit = minCollZ(*target, *bullet, D[0], D[1], D[2]);
            } else {
                //! By including pointers to a node's left and right particles,
                //! we can directly calculate the displacement in the z
                //! direction for these two particles which maintains the
                //! particle's connectivity thus negating the need to call
                //! calcCollZ function through minCollZ.
                hit = calcCollZ(node.m_leftparticle->BoundSphCentre(), node.m_leftparticle->Radius(),
                                node.m_rightparticle->BoundSphCentre(), node.m_rightparticle->Radius(),
                                D[0], D[1], D[2], sumr, trackPrimaryCoordinates);
            }
        }

        //! We have a new location for the bullet (right node), so move it.
        node.m_rightchild->Translate(D[0], D[1], D[2]);

        //! Calculate properties of this node.
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
                             double dx, double dy, double &dz)
{
    bool hit=false, hit1=false;
    double dz2=0.0, dz3=0.0, dz4=0.0;

    if (target.IsLeaf()) {
        // Target is a leaf
        if (bullet.IsLeaf()) {
            // Bullet is a leaf (both leaves).
           return calcCollZ(target.BoundSphCentre(), target.Radius(),
                            bullet.BoundSphCentre(), bullet.Radius(),
                            dx, dy, dz, 0.0, false);
        } else {
            // Bullet is not a leaf, call sub-nodes.
            // Calculate minimum dz for the target and the bullet left subnode.
            hit1 = calcCollZ(target.BoundSphCentre(), target.Radius(),
                             bullet.m_leftchild->BoundSphCentre(), bullet.m_leftchild->Radius(),
                             dx, dy, dz, 0.0, false);
            if (hit1) hit = minCollZ(target, *bullet.m_leftchild, dx, dy, dz);
            // Calculate minimum dz for the target and the bullet right subnode.
            hit1 = calcCollZ(target.BoundSphCentre(), target.Radius(),
                             bullet.m_rightchild->BoundSphCentre(), bullet.m_rightchild->Radius(),
                             dx, dy, dz2, 0.0, false);
            if (hit1) hit = minCollZ(target, *bullet.m_rightchild, dx, dy, dz2) || hit;
            // Return minimum dz.
            dz = min(dz, dz2);
            return hit;
        }
    } else {
        // Target is not a leaf.
        if (bullet.IsLeaf()) {
            // Bullet is a leaf, call target sub-nodes..
            // Calculate minimum dz for the target left subnode and the bullet.
            hit1 = calcCollZ(target.m_leftchild->BoundSphCentre(), target.m_leftchild->Radius(),
                             bullet.BoundSphCentre(), bullet.Radius(),
                             dx, dy, dz, 0.0, false);
            if (hit1) hit = minCollZ(*target.m_leftchild, bullet, dx, dy, dz);
            // Calculate minimum dz for the target right subnode and the bullet.
            hit1 = calcCollZ(target.m_rightchild->BoundSphCentre(), target.m_rightchild->Radius(),
                             bullet.BoundSphCentre(), bullet.Radius(),
                             dx, dy, dz2, 0.0, false);
            if (hit1) hit = minCollZ(*target.m_rightchild, bullet, dx, dy, dz2) || hit;
            // Return minimum dz.
            dz = min(dz, dz2);
            return hit;
        } else {
            // Bullet is not a leaf (neither is a leaf), check all left/right
            // collision combinations.
            // Target left and bullet left.
            hit1 = calcCollZ(target.m_leftchild->BoundSphCentre(), target.m_leftchild->Radius(),
                             bullet.m_leftchild->BoundSphCentre(), bullet.m_leftchild->Radius(),
                             dx, dy, dz, 0.0, false);
            if (hit1) hit = minCollZ(*target.m_leftchild, *bullet.m_leftchild, dx, dy, dz);
            // Target left and bullet right.
            hit1 = calcCollZ(target.m_leftchild->BoundSphCentre(), target.m_leftchild->Radius(),
                             bullet.m_rightchild->BoundSphCentre(), bullet.m_rightchild->Radius(),
                             dx, dy, dz2, 0.0, false);
            if (hit1) hit = minCollZ(*target.m_leftchild, *bullet.m_rightchild, dx, dy, dz2) || hit;
            // Target right and bullet left.
            hit1 = calcCollZ(target.m_rightchild->BoundSphCentre(), target.m_rightchild->Radius(),
                             bullet.m_leftchild->BoundSphCentre(), bullet.m_leftchild->Radius(),
                             dx, dy, dz3, 0.0, false);
            if (hit1) hit = minCollZ(*target.m_rightchild, *bullet.m_leftchild, dx, dy, dz3) || hit;
            // Target right and bullet right.
            hit1 = calcCollZ(target.m_rightchild->BoundSphCentre(), target.m_rightchild->Radius(),
                             bullet.m_rightchild->BoundSphCentre(), bullet.m_rightchild->Radius(),
                             dx, dy, dz4, 0.0, false);
            if (hit1) hit = minCollZ(*target.m_rightchild, *bullet.m_rightchild, dx, dy, dz4) || hit;
            // Returns minimum dz.
            dz = min(min(dz, dz2), min(dz3, dz4));
            return hit;
        }
    }
}

/*!
 *  @brief Calculates the z-displacement of a sphere.
 *
 *  Calculates the z-displacement of a bullet sphere for a +ve
 *  collision with a target sphere. Returns true if the
 *  spheres collide, otherwise false.
 *
 *  @param[in]  p1                     Coordinates of sphere 1.
 *  @param[in]  r1                     Radius of sphere 1.
 *  @param[in]  p2                     Coordinates of sphere 2.
 *  @param[in]  r2                     Radius of sphere 2
 *  @param[in]  dx                     Bullet x displacement.
 *  @param[in]  dy                     Bullet y displacement.
 *  @param[out] dz                     Bullet z displacement.
 *  @param[in]  distanceCentreToCentre Distance between the centres of neighbouring primary particles.
 *  @param[in]  trackPrimaryCoordinates Flag used to indicate whether to track primary coordinates.
 *
 *  @return Have the nodes collided?
 */
bool ParticleImage::calcCollZ(const Coords::Vector &p1, double r1,
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