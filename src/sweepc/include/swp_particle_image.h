/*!
 * @file    swp_particle_image.h
 * @author  Matthew S Celnik.
 * @brief   Draws TEM-style images of particles
 *
 *  Copyright (C) 2008 Matthew S Celnik.
 *  re-integrated into Sweep by William J Menz 2012
 *
 *  File purpose:
 *    The Imager class draws particle images and saves them to file.  Currently
 *    only POVRAY output files are generated, though it is expected that
 *    other formats will be added later.
 *
 *    Particles are drawn as fractal aggregates using the primary-particle
 *    information.  Different particle models generate different primary-particle
 *    data structures, which are all catered for in this class.  The Draw
 *    functions automatically detect what form of data is available and draw
 *    the aggregate images accordingly.
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

#ifndef SWEEP_IMAGER_H
#define SWEEP_IMAGER_H

#include "swp_params.h"
#include "swp_particle.h"
#include "swp_imgnode.h"
#include "swp_particle_model.h"
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
    // Constructors.
    ParticleImage(void);               // Default constructor.

    // Destructors.
    ~ParticleImage(void);


    // PARTICLE IMAGE DATA CONSTRUCTION.

    // Constructs the particle image from the given particle.
    void Construct(const Particle &sp, const ParticleModel &model);

    // RENDERING FUNCTIONS.
    //! Writes 3dout file (Markus Sander), deprecated
    void Write3dout(std::ofstream &file, double x, double y, double z);

    //! Draws the particle image to a POVRAY file.
    void WritePOVRAY(
        std::ofstream &file // Output file stream.
        );

private:

    // The image aggregate structure root node.
    ImgNode m_root;

    // Amount of necking between particles in output.
    static const real m_necking;

    // AGGREGATE SPHERE-TREE CONSTRUCTORS (FREE-MOLECULAR).

    //! Generate the free-molecular structure of a particle
    static void calc_FM(ImgNode &node, Sweep::rng_type &rng);

    /*!
     * @brief       Loop helper function to construct particle images
     *
     * Recursively loops though the binary tree structure to generate
     * the required image
     *
     * @param p     Pointer to particle for which to construct tree
     */
    template <class ParticleClass>
    void ConstructTreeLoop(const ParticleClass *p) {
        // Loop through tree structure
        if (p->m_leftchild != NULL && p->m_rightchild != NULL) {
            ConstructTreeLoop(p->m_leftchild);
            ConstructTreeLoop(p->m_rightchild);
        } else {
            //convert to nm, store the radius not the diameter
            m_root.Insert(p->SphDiameter()*0.5e9);
        }
    };

    /*!
     * @brief      Constructs a binary tree image
     *
     * Used for generating images of particles with a binary tree
     * structure, e.g. PAHPrimary or SilicaPrimary
     *
     * @param p    Pointer to particle for which to construct tree
     * @param rng  Random number generator
     */
    template <class ParticleClass>
    void ConstructTree(const ParticleClass *p, Sweep::rng_type &rng) {

        m_root.Clear();

        // Call the helper function to generate structure
        ConstructTreeLoop(p);

        // Use the free-molecular regime to calculate the
        // aggregate structure.
        calc_FM(m_root, rng);
        m_root.CentreCOM();
    };

    //void LengthWidth(double &L, double &W);

    //! Calculates the radius of gyration
    double RadiusofGyration();

    //! Generates a projection on the zx plane (set all y to 0)
    void Project();


    //! Calculates the z-displacement of a sphere
    static bool calcCollZ(
        const Coords::Vector &p1, // Positional vector of sphere 1.
        real r1,                  // Radius of sphere 1.
        const Coords::Vector &p2, // Positional vector of sphere 2.
        real r2,                  // Radius of sphere 2.
        real dx, real dy,         // Sphere 2 x and y displacements.
        real &dz                  // The output z-axis displacement of the bullet (+ve).
        );

    //! Calculates the minimum collision distance
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
