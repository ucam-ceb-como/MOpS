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
#include <stack> //!< Used by descendPath and recordPath.

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

protected:
    //! Sets the pointers to the primary particles correct after a copy event.
    template <class ParticleClass>
    void UpdateAllPointers(ImgNode &node, const ParticleClass *original);

private:
    // The image aggregate structure root node.
    ImgNode m_root;

    // Amount of necking between particles in output.
    static const double m_necking;

    // AGGREGATE SPHERE-TREE CONSTRUCTORS (FREE-MOLECULAR).

    //! Generate the free-molecular structure of a particle.
    static void calc_FM(ImgNode &node, Sweep::rng_type &rng, const bool trackPrimarySeparation);

    /*!
     *  @brief Constructs a binary tree image.
     *
     *  Used for generating images of particles with a binary tree
     *  structure, e.g. PAHPrimary or SilicaPrimary.
     *
     *  @param[in] p                      Pointer to particle for which to construct tree.
     *  @param[in] rng                    Random number generator.
     *  @param[in] trackPrimarySeparation Flag used to indicate whether to track primary separation.
     */
    template <class ParticleClass>
    void ConstructTree(const ParticleClass *p, Sweep::rng_type &rng, const bool trackPrimaryCoordinates) {
        m_root.Clear();

		//! Copy properties of the original particle to the root node.
        CopyParts(m_root, p);
        
		//! Copy the rest of the tree if particle is made up of more than one
        //! primary.
        if (p->m_leftchild != NULL)
            CopyTree(m_root, p);

        //! Generate structure of aggregate if coordinates of the primary
        //! particles are not tracked on-the-fly.
        if (!trackPrimaryCoordinates) {
            calc_FM(m_root, rng, trackPrimaryCoordinates);
            m_root.CentreCOM();
		}
    };
    
    //! Copies the node without the children.
    template <class ParticleClass>
    void CopyParts(ImgNode &node, const ParticleClass *source);

    //! Copies the subtree of a node.
    template <class ParticleClass>
    void CopyTree(ImgNode &node, const ParticleClass *source);

    //! Find the path through the tree from node top to node bottom.
    template <class ParticleClass>
    static std::stack<bool> recordPath(const ParticleClass* bottom, const ParticleClass* const top);

    //! Follow a path down the tree.
    template <class ParticleClass>
    static ParticleClass* descendPath(ParticleClass *here, std::stack<bool> &takeLeftBranch);

    //void LengthWidth(double &L, double &W);

    //! Calculates the radius of gyration
    double RadiusofGyration();

    //! Generates a projection on the zx plane (set all y to 0)
    void Project();

    //! Calculates the z-displacement of a sphere.
    static bool calcCollZ(
        const Coords::Vector &p1,         //!< Positional vector of sphere 1.
        double r1,                        //!< Radius of sphere 1.
        const Coords::Vector &p2,         //!< Positional vector of sphere 2.
        double r2,                        //!< Radius of sphere 2.
        double dx, double dy,             //!< Sphere 2 x and y displacements.
        double &dz,                       //!< The output z-axis displacement of the bullet (+ve).
        double distanceCentreToCentre,    //!< Distance between the centres of primary particles.
        const bool trackPrimaryCoordinates //!< Flag used to indicate whether to track primary coordinates.
        );

    //! Calculates the minimum collision distance.
    static bool minCollZ(
        const ImgNode &target, //!< Target node.
        const ImgNode &bullet, //!< Bullet node.
        double dx, double dy,  //!< Bullet x-y displacements.
        double &dz             //!< Return minimum distance.
        );

    // OUTPUT FUNCTIONS.

};
};
};

#endif
