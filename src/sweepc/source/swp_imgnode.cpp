/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ImgNode class declared in the
    swp_imgnode.h header file.

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

#include "swp_imgnode.h"
#include "swp_particle.h"
#include "swp_primary.h"


using namespace Sweep;
using namespace Sweep::Imaging;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ImgNode::ImgNode(void)
: m_distance_centreToCentre(0.0), m_r(0.0), m_r2(0.0), m_r3(0.0), m_mindepth(0), m_parent(NULL), 
  m_leftchild(NULL), m_rightchild(NULL), m_leftparticle(NULL), m_rightparticle(NULL)
{
}

// Destructor.
ImgNode::~ImgNode(void)
{
    delete m_leftchild;
    delete m_rightchild;
}


// TREE BUILDING (PRIMARY ADDITION).

// Clears the node back to initial state.  Also clears
// children.
void ImgNode::Clear(void)
{
    // Set size to zero.
    m_distance_centreToCentre = 0.0;
    m_r  = 0.0;
    m_r2 = 0.0;
    m_r3 = 0.0;

    // Centre coordinates.
    m_cen_bsph[0] = 0.0;
    m_cen_bsph[1] = 0.0;
    m_cen_bsph[2] = 0.0;
    m_cen_mass[0] = 0.0;
    m_cen_mass[1] = 0.0;
    m_cen_mass[2] = 0.0;

    // Clear children.
    if (m_leftchild != NULL) {
        m_leftchild->Clear();
        delete m_leftchild; m_leftchild = NULL;
    }
    if (m_rightchild != NULL) {
        m_rightchild->Clear();
        delete m_rightchild; m_rightchild = NULL;
    }
}



void ImgNode::CopySPT(const Particle *sp)
{
    setRadius( (sp->Primary()->Property(Sweep::iDsph))*0.5e9);       //convert to nm
}

// Adds a primary into the tree.  Algorithm maintains
// well-balanced binary tree structure by observing
// the minimum branch depth below node.
void ImgNode::Insert(double radius)
{
    // Check for an empty tree.
    if ((m_leftchild == NULL) && (m_rightchild == NULL)) {
        if (m_r == 0.0) {
            // This node is empty, just set its radius.
            setRadius(radius);
            m_mindepth = 0;
        } else {
            // This node is not empty, but also has
            // no children.

            // Set left child to have the properties
            // of this node.
            m_leftchild = new ImgNode();
            m_leftchild->setRadius(m_r);
            m_leftchild->m_parent = this;

            // Add a new right node for the new primary.
            m_rightchild = new ImgNode();
            m_rightchild->setRadius(radius);
            m_rightchild->m_parent = this;

            // Mark this node as not calculated.
           setRadius(0.0);

           // Set the tree depth.
           m_mindepth = 1;
        }
    } else {
        // This node has children.
        if (m_leftchild->m_mindepth <= m_rightchild->m_mindepth) {
            // Add primary to the left branch if branch
            // depths are equal, or left branch is shallower.
            m_leftchild->Insert(radius);
        } else {
            // Add primary to the right branch if it is
            // shallower.
            m_rightchild->Insert(radius);
        }
        // Recalulate the minimum tree depth of this node.
        m_mindepth = min(m_leftchild->m_mindepth, m_rightchild->m_mindepth)+1;
    }
}

// Returns true if this node is a leaf (has no children).
bool ImgNode::IsLeaf(void) const
{
    return (m_leftchild==NULL) && (m_rightchild==NULL);
}

// Returns a vector of primary coordinates (4D).  The first
// three returned values are the cartesian x,y,z coordinates,
// the final value is the radius.
void ImgNode::GetPriCoords(std::vector<fvector> &coords) const
{
    if (IsLeaf()) {
        fvector c(4);
        c[0] = m_cen_mass[0];
        c[1] = m_cen_mass[1];
        c[2] = m_cen_mass[2];
        c[3] = m_r;
        coords.push_back(c);
    } else {
        m_leftchild->GetPriCoords(coords);
        m_rightchild->GetPriCoords(coords);
    }
}


// COORDINATE MANIPULATION.

//!Projects the particle to the x-z plane (sets all y to zero)
void ImgNode::Project()
{
    m_cen_mass[1]=0;
    if (!IsLeaf()) { 
        m_leftchild->Project();
        m_rightchild->Project();
    }
}

// Translates (moves) the aggregate node and child structure
// by the given amounts along the cartesian axes.
void ImgNode::Translate(double dx, double dy, double dz)
{
    // Translate child branches.
    if (m_leftchild != NULL) m_leftchild->Translate(dx, dy, dz);
    if (m_rightchild != NULL) m_rightchild->Translate(dx, dy, dz);
    // Translate bounding sphere and centre-of-mass centres.
    m_cen_bsph.Translate(dx, dy, dz);
    m_cen_mass.Translate(dx, dy, dz);
}

// Put the centre-of-mass at the origin.
void ImgNode::CentreCOM(void)
{
    Translate(-m_cen_mass[0], -m_cen_mass[1], -m_cen_mass[2]);
}

// Put the bounding-sphere at the origin.
void ImgNode::CentreBoundSph(void)
{
    Translate(-m_cen_bsph[0], -m_cen_bsph[1], -m_cen_bsph[2]);
}

// Rotates the aggregate node and child structure about its centre
// of mass by the given angles (spherical coordinates).
void ImgNode::RotateCOM(double dtheta, double dphi)
{
    // Move the aggregate so that its centre-of-mass
    // is at the origin.  Store the coordinates, so that
    // they can be restored afterwards.
    Coords::Vector D(m_cen_mass);
    Translate(-D.X(), -D.Y(), -D.Z());

    // Create transformation matrix.
    Coords::Matrix M;
    M.SetIdentity();
    M.Rotate(dtheta, dphi);

    // Rotate child nodes.
    if (m_leftchild != NULL) m_leftchild->transform(M);
    if (m_rightchild != NULL) m_rightchild->transform(M);

    // Rotate bounding-sphere coordinates.
    m_cen_bsph = M.Mult(m_cen_bsph);

    // Restore centre-of-mass coordinates.
    Translate(D.X(), D.Y(), D.Z());
}

// Rotates the aggregate node and child structure about its bounding
// sphere centre by the given angles (spherical coordinates).
void ImgNode::RotateBoundSph(double dtheta, double dphi)
{
    // Move the aggregate so that its bounding sphere
    // is at the origin.  Store the coordinates, so that
    // they can be restored afterwards.
    Coords::Vector D(m_cen_bsph);
    Translate(-D.X(), -D.Y(), -D.Z());

    // Create transformation matrix.
    Coords::Matrix M;
    M.SetIdentity();
    M.Rotate(dtheta, dphi);

    // Rotate child nodes.
    if (m_leftchild != NULL) m_leftchild->transform(M);
    if (m_rightchild != NULL) m_rightchild->transform(M);

    // Rotate centre-of-mass coordinates.
    m_cen_mass = M.Mult(m_cen_mass);

    // Restore bounding-sphere coordinates.
    Translate(D.X(), D.Y(), D.Z());
}

// Rotates the aggregate node and child structure about the
// coordinate system origin by the given angles (spherical coordinates).
void ImgNode::RotateOrigin(double dtheta, double dphi)
{
    // Create transformation matrix.
    Coords::Matrix M;
    M.SetIdentity();
    M.Rotate(dtheta, dphi);
    // Perform transformation.
    transform(M);
}

// Transforms the node coordinates using the given
// transformation matrix.
void ImgNode::transform(const Coords::Matrix &mat)
{
    // Rotate child nodes.
    if (m_leftchild != NULL) 
        m_leftchild->transform(mat);
    if (m_rightchild != NULL) 
        m_rightchild->transform(mat);
    // Rotate centre-of-mass and bounding sphere coords.
    m_cen_mass = mat.Mult(m_cen_mass);
    m_cen_bsph = mat.Mult(m_cen_bsph);
}


// BOUNDING SPHERE.

// Returns the bounding-sphere centre.
const Coords::Vector &ImgNode::BoundSphCentre(void) const
{
    return m_cen_bsph;
}

// Returns the bounding sphere radius.
double ImgNode::Radius(void) const
{
    return m_r;
}

// Sets the radius of the bounding sphere.
void ImgNode::setRadius(double r)
{
    m_r  = r;
    m_r2 = r * r;
    m_r3 = m_r2 * m_r;
}

//! Set the bounding sphere. Units of nm.
void ImgNode::setBoundSph(Coords::Vector bsph)
{
    m_cen_bsph[0] = bsph[0] * 1.0e9;
    m_cen_bsph[1] = bsph[1] * 1.0e9;
    m_cen_bsph[2] = bsph[2] * 1.0e9;
}

//! Set the centre-of-mass. Units of nm.
void ImgNode::setCOM(Coords::Vector mass)
{
    m_cen_mass[0] = mass[0] * 1.0e9;
    m_cen_mass[1] = mass[1] * 1.0e9;
    m_cen_mass[2] = mass[2] * 1.0e9;
}

//! Set the distance between the centres of primary particles.
void ImgNode::setDistance(double distance)
{
    m_distance_centreToCentre = distance;
}

// Calculates the bounding sphere position and radius using
// the left and right child node values.
void ImgNode::CalcBoundSph(void)
{
    if ((m_leftchild != NULL) && (m_rightchild != NULL)) {
        // Calculate bounding spheres of children.
        m_leftchild->CalcBoundSph();
        m_rightchild->CalcBoundSph();

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


// CENTRE-OF-MASS.

// Returns the centre-of-mass.
const Coords::Vector &ImgNode::CentreOfMass(void) const
{
    return m_cen_mass;
}

// Calculates the centre-of-mass using the left and right child
// node values.
void ImgNode::CalcCOM(void)
{
    if ((m_leftchild != NULL) && (m_rightchild != NULL)) {
        // Calculate centres-of-mass of left and
        // right children.
        m_leftchild->CalcCOM();
        m_rightchild->CalcCOM();

        // Mass is proportional to r^3.  Calculated total
        // mass (inverse) of left and right children.
    //    double invtotmass = 1.0 / (m_left->m_r3 + m_right->m_r3);
        m_mass=m_leftchild->m_mass + m_rightchild->m_mass;
        double invtotmass = 1.0 / m_mass;
        // Now calculate centre of mass.
        for (unsigned int i=0; i!=3; ++i) {
            //m_cen_mass[i]  = m_left->m_cen_mass[i] * m_left->m_r3;
          //  m_cen_mass[i] += m_right->m_cen_mass[i] * m_right->m_r3;
            m_cen_mass[i]  = m_leftchild->m_cen_mass[i] * m_leftchild->m_mass;
            m_cen_mass[i] += m_rightchild->m_cen_mass[i] * m_rightchild->m_mass;
            m_cen_mass[i] *= invtotmass;
        }
    } else {
        // If there are no children, then the centre-of-mass and
        // bounding sphere centre are the same.
        m_cen_mass[0] = m_cen_bsph[0];
        m_cen_mass[1] = m_cen_bsph[1];
        m_cen_mass[2] = m_cen_bsph[2];
        m_mass=m_r3;
    }
}
