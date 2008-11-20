/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Imager class draws particle images and saves them to file.  Currently
    only POVRAY output files are generated, though it is expected that
    other formats will be added later.

    Particles are drawn as fractal aggregates using the primary-particle
    information.  Different particle models generate different primary-particle
    data structures, which are all catered for in this class.  The Draw
    functions automatically detect what form of data is available and draw
    the aggregate images accordingly.

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

#ifndef SWEEP_IMAGER_H
#define SWEEP_IMAGER_H

#include "swp_params.h"
#include "swp_particle.h"
#include "swp_surfvol_primary.h"
#include "swp_pripart_primary.h"
#include "swp_imgnode.h"
#include <string>
#include <fstream>
#include <vector>

namespace Sweep
{
namespace Imaging
{
class ParticleImage
{
public:
    // Enumeration of possible collision regimes.  Different regimes
    // use different collision algorithms, so may produce differently
    // shaped aggregates.
    enum CollRegime{FreeMol=0};

    // Constructors.
    ParticleImage(void);               // Default constructor.
    ParticleImage(const Particle &sp); // Initialising constructor.

    // Destructors.
    ~ParticleImage(void);


    // COLLISION REGIME.

    // Returns the collision regime used to calculate the
    // particle image.
    CollRegime CollisionRegime(void) const;

    // Sets the collision regime used to calculate the
    // particle image.  If the collision regime is being 
    // changed then this function also forces a recalculation
    // of the image aggregate structure, assuming that it 
    // has already been initialised.
    void SetCollRegime(CollRegime creg);


    // PARTICLE IMAGE DATA CONSTRUCTION.

    // Constructs the particle image from the given particle.
    void Construct(const Particle &sp);

    // Constructs a random particle image.
    void ConstructRandom(
        real minrad,   // Minimum primary radius.
        real maxrad,   // Maximum primary radius.
        unsigned int n // Number of primaries to generate.
        );


	void constructSubParttree(const SubParticle *sp);
	void copysptinsert(const SubParticle *sp);
	void Write3dout(std::ofstream &file, double x, double y, double z);


    // RENDERING FUNCTIONS.

    // Draws the particle image to a POVRAY file.
    void WritePOVRAY(
        std::ofstream &file // Output file stream.
        );

private:
    // The collision regime used to calculate the image
    // aggregate structure.  The default is free-molecular.
    CollRegime m_creg;

    // The image aggregate structure root node.
    ImgNode m_root;

    // Amount of necking between particles in output.
    static const real m_necking;

    // AGGREGATE SPHERE-TREE CONSTRUCTORS (FREE-MOLECULAR).

    // Constructs a PNode sphere-tree aggregate from the given 
    // particle using free-molecular collision dynamics.
    void constructAgg_FM(const Particle &sp);

    // Constructs a PNode sphere-tree aggregate from the given 
    // pri-part list primary using free-molecular collision dynamics.
    void constructAgg_FM(const AggModels::PriPartPrimary &pri);

    // Constructs a PNode sphere-tree aggregate with uniform 
    // primaries (equal diameter).  The diameter and primary
    // count are passed as arguments.
    void uniformAgg_FM(
        unsigned int n, // Primary particle count.
        real d          // Primary particle diameter.
        );

    // Calculates the aggregate structure down from the given
    // node.  Assumes that the tree leaves have been initialised
    // with the correct radii, and recalculates their positions.
    static void calc_FM(ImgNode &node);

    // Calculates the z-displacement of a bullet sphere for a +ve
    // collision with a target sphere.  Returns true if the
    // spheres collide, otherwise false.
    static bool calcCollZ(
        const Coords::Vector &p1, // Positional vector of sphere 1.
        real r1,                  // Radius of sphere 1.
        const Coords::Vector &p2, // Positional vector of sphere 2.
        real r2,                  // Radius of sphere 2.
        real dx, real dy,         // Sphere 2 x and y displacements.
        real &dz                  // The output z-axis displacement of the bullet (+ve).
        );

    // Calculates the minimum collision distance between
    // a target and a bullet node by moving down the
    // binary tree.  If the nodes collide then returns
    // true, otherwise returns false.
    static bool minCollZ(
        const ImgNode &target, // Target node.
        const ImgNode &bullet, // Bullet node.
        real dx, real dy,      // Bullet x-y displacements.
        real &dz               // Return minimum distance.
        );

    // OUTPUT FUNCTIONS.

};
};
};

#endif
