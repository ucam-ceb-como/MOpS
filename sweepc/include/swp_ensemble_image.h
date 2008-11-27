/*
  Author(s):      Matthew Celnik (msc37) and Markus Sander (ms785)
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

#ifndef SWEEP_ENSEMBLE_IMAGER_H
#define SWEEP_ENSEMBLE_IMAGER_H

#include "swp_params.h"
#include "swp_particle.h"
#include "swp_surfvol_primary.h"
#include "swp_pripart_primary.h"
#include "swp_ensemble_imgnode.h"
#include <string>
#include <fstream>
#include <vector>


namespace Sweep
{
namespace Imaging
{
class EnsembleImage
{
public:

    // Constructors.
    EnsembleImage(void);               // Default constructor.
    EnsembleImage(const Particle &sp); // Initialising constructor.

    // Destructors.
    ~EnsembleImage(void);


    // COLLISION REGIME.

    // Returns the collision regime used to calculate the
    // particle image.
    

    // Sets the collision regime used to calculate the
    // particle image.  If the collision regime is being 
    // changed then this function also forces a recalculation
    // of the image aggregate structure, assuming that it 
    // has already been initialised.
	void PrintEnsemble(Cell &sys,std::ofstream &file, real shiftz);
	void Write3dout(std::ofstream &file, EnsembleImgNode *curr_node);





private:

    // The collision regime used to calculate the image
    // aggregate structure.  The default is free-molecular.
    // The image aggregate structure root node.
    EnsembleImgNode m_root;


    // Calculates the aggregate structure down from the given
    // node.  Assumes that the tree leaves have been initialised
    // with the correct radii, and recalculates their positions.
    static void calc_FM_ensemble(EnsembleImgNode &node);

    // Calculates the z-displacement of a bullet sphere for a +ve
    // collision with a target sphere.  Returns true if the
    // spheres collide, otherwise false.
    static bool calcCollZ_ensemble(
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
    static bool minCollZ_ensemble(
        const EnsembleImgNode &target, // Target node.
        const EnsembleImgNode &bullet, // Bullet node.
        real dx, real dy,      // Bullet x-y displacements.
        real &dz               // Return minimum distance.
        );
	

};
};
};

#endif
