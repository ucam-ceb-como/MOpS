/*!
 * \file   swp_tree_cache.h
 * \author Robert I A Patterson
 *  Copyright (C) 2010 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Vector like structure for summing double properties in binary tree
 *
 Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Prof Markus Kraft
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

#ifndef SWEEP_TREE_CACHE_H
#define SWEEP_TREE_CACHE_H

#include "swp_params.h"
#include "swp_property_indices.h"

namespace Sweep
{
    //forward declaration
    class Particle;

/*!
 * \brief   Summary properties that can be summed and used to calculate overall population rates
 *
 * The intended use is to store and facilitate the summation of numerically useful
 * particle properties.  At the time of writing the expected use was in the leaves
 * of binary trees which are used to cache partial sums and use them for efficient
 * inversion of arbitrary discrete distributions defined by those sums.  See the
 * binary tree documentation or pages 46-49 of "Numerical Modelling of Soot Formation",
 * R I A Patterson, 2007, which is the PhD thesis of the original author.  It should
 * be available electronically from http://como.cheng.cam.ac.uk/...
 *
 * Note that this class is a cache, it contains no original data.  The names of the
 * data members and accessors are suggestive of particular physical interpretations,
 * because these interpretations were the intention of the original developers.
 * However a cache holds what is put into it.  For example, consider a summation
 * tree based on spherical particles:  For every leaf node (a leaf node corresponds
 * to one particle SphDiameter == (6*Volume/pi)^(1/3), however, once one starts
 * summing this relation ceases to hold, even though the cache correctly holds
 * the sums of SphDiameter and Volume, because the relationship is non-linear.
 *
 * Currently has 13 data members
 */
class TreeCache
{
public:
    //! Initialise with zeros
    TreeCache();

    //! Initialise with values from a particle
    TreeCache(const Sweep::Particle &part);

    //! Increment elementwise
    TreeCache &operator+=(const TreeCache &rhs);

    //! Add together the values from two caches element by element
    const TreeCache operator+(const TreeCache &rhs) const;

    //! Resets the cache to its "empty" condition.
    void Clear(void);


    // BASIC PROPERTIES.

    //! Returns the particle equivalent sphere diameter.
    double SphDiameter(void) const;

    //! Returns the collision diameter.
    double CollDiameter(void) const;

    //! Returns the mobility diameter.
    double MobDiameter(void) const;

    //! Returns the surface area.
    double SurfaceArea(void) const;

    //! Cached for surface area of volume equivalent spheres
    double SphSurfaceArea(void) const;

    //! Returns the cached volume.
    double Volume(void) const;

    //! Returns the cached mass.
    double Mass(void) const;

    //! Returns the number of carbon atoms.
    int NumCarbon(void) const;

    //! Returns fragmentation flag.
    int Frag(void) const;
    // Returns the property with the given ID.
    double Property(PropID id) const;

	//double Sintered();

    // BASIC DERIVED PROPERTY OVERWRITES.

    /*// Sets the spherical particle diameter
    void SetSphDiameter(double diam);

    // Sets the collision diameter of the particle.
    void SetCollDiameter(double dcol);

    // Sets the mobility diameter.
    void SetMobDiameter(double dmob);

    // Sets the surface area, subject to minimum spherical area condition.
    void SetSurfaceArea(double surf);

    // Sets the volume.
    void SetVolume(double vol);

    // Sets the mass.
    void SetMass(double m);



    // COLLISION RATE CALCULATION PARTICLE PROPERTIES.

    // Collision diameter squared (cm2).
    double CollDiamSquared() const;

    // Inverse collision diameter (cm-1).
    double InvCollDiam() const;

    // Inverse squared collision diameter (cm-2).
    double InvCollDiamSquared() const;

    // Inverse of square root of mass (g-1/2).
    double InvSqrtMass() const;

    // Collision diameter squared times the inverse square root of mass.
    double CollDiamSqrdInvSqrtMass() const;*/


    // READ/WRITE/COPY.

    // Writes the object to a binary stream.
    //virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    /*virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );*/

private:

    // Basic properties
    //! Equivalent spherical diameter.
    double m_sphdiam;

    //! Collision diameter.
    double m_dcol;

    //! Mobility diameter.
    double m_dmob;

    //! Surface area.
    double m_surf;

    //! Volume.
    double m_vol;

    //! Mass.
    double m_mass;

	//! Number of Carbon atoms
	unsigned int m_numcarbon;

    //! Fragmentation flag.
	unsigned int m_frag;

    // Collision rate calculation particle properties.
    double m_dcolsqr;      // Collision diameter squared.
    double m_inv_dcol;     // Inverse collision diameter.
    double m_inv_dcolsqr;  // Inverse of the diameter squared.
    double m_inv_sqrtmass; // Inverse of the square-root of the mass.
    double m_d2_m_1_2;     // D^2 * M^-1/2.

    //! Statistical weight
    double m_weight;

	// The free surface available for other particles to sinter
	double m_freesurface;
};
}

#endif
