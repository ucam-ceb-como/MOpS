/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ImgNode class is used to construct a binary tree of sphere which
    describe an aggregate particle structure consisting of spherical
    primary particles.  It is a critical component of the ParticleImage
    class.
*/

#ifndef SWEEP_IMGNODE_H
#define SWEEP_IMGNODE_H

#include "swp_params.h"
#include "swp_coords.h"

namespace Sweep
{
namespace Imaging
{
class ParticleImage;

class ImgNode
{
friend class ParticleImage;
public:
    // Constructors.
    ImgNode(void); // Default constructor.

    // Destructors.
    ~ImgNode(void);


    // TREE BUILDING (PRIMARY ADDITION).

    // Clears the data from the node.  Also clears the
    // child nodes, but does not delete them.
    void Clear(void);

    // Adds a primary into the tree.  Algorithm maintains
    // well-balanced binary tree structure by observing
    // the minimum branch depth below node.
    void Insert(real radius);

    // Returns true if this node is a leaf (has no children).
    bool IsLeaf(void) const;

    // Returns a vector of primary coordinates (4D).  The first
    // three returned values are the cartesian x,y,z coordinates,
    // the final value is the radius.
    void GetPriCoords(std::vector<fvector> &coords) const;


    // COORDINATE MANIPULATION.

    // Translates (moves) the aggregate node and child structure
    // by the given amounts along the cartesian axes.
    void Translate(real dx, real dy, real dz);

    // Put the centre-of-mass at the origin.
    void CentreCOM(void);

    // Put the bounding-sphere at the origin.
    void CentreBoundSph(void);

    // Rotates the aggregate node and child structure about its centre
    // of mass by the given angles (spherical coordinates).
    void RotateCOM(real dtheta, real dphi);

    // Rotates the aggregate node and child structure about its bounding
    // sphere centre by the given angles (spherical coordinates).
    void RotateBoundSph(real dtheta, real dphi);

    // Rotates the aggregate node and child structure about the
    // coordinate system origin by the given angles (spherical coordinates).
    void RotateOrigin(real dtheta, real dphi);
    

    // BOUNDING SPHERE.

    // Returns the bounding-sphere centre.
    const Coords::Vector &BoundSphCentre(void) const;

    // Returns the bounding sphere radius.
    real Radius(void) const;

    // Calculates the bounding sphere position and radius using
    // the left and right child node values.
    void CalcBoundSph(void);

    // CENTRE-OF-MASS.

    // Returns the centre-of-mass.
    const Coords::Vector &CentreOfMass(void) const;

    // Calculates the centre-of-mass using the left and right child
    // node values.
    void CalcCOM(void);

private:
    // Size.
    real m_r;  // Bounding sphere radius of aggregate/primary.
    real m_r2; // r squared (useful for efficient collision detection computation).
    real m_r3; // r cubed (useful for calculating centre-of-mass).

    // Position.
    Coords::Vector m_cen_bsph; // Bounding-sphere centre.
    Coords::Vector m_cen_mass; // Centre-of-mass coordinates.

    // Minimum branch depth beneath this node.
    unsigned int m_mindepth;

    // Binary tree (divide & conquer) structure.
    ImgNode *m_parent;         // Parent node.
    ImgNode *m_left, *m_right; // Left and right child nodes.

    // Sets the radius of the bounding sphere.
    void setRadius(real r);

    // Transforms the node coordinates using the given
    // transformation matrix.
    void transform(const Coords::Matrix &mat);
};
};
};

#endif
