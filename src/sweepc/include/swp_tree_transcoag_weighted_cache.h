/*!
 * \file   swp_tree_transcoag_weighted_cache.h
 * \author William J Menz
 *  Copyright (C) 2011 William J Menz.
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

#ifndef SWEEP_TREE_TRANSCOAG_WEIGHTED_CACHE_H
#define SWEEP_TREE_TRANSCOAG_WEIGHTED_CACHE_H

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
 * The terms mass and weight have very different meanings.  Mass refers to the
 * ordinary mass of a particle (SI unit kg); statistical weight specifies how
 * many physical particles are represented by one simulation particle
 * (SI unit m^-3).
 *
 * Currently has 21 data members
 */
class TreeTransCoagWeightedCache
{
public:
    //! Initialise with zeros
    TreeTransCoagWeightedCache();

    //! Initialise with values from a particle
    TreeTransCoagWeightedCache(const Sweep::Particle &part);

    //! Increment elementwise
    TreeTransCoagWeightedCache &operator+=(const TreeTransCoagWeightedCache &rhs);

    //! Add together the values from two caches element by element
    const TreeTransCoagWeightedCache operator+(const TreeTransCoagWeightedCache &rhs) const;

    //! Resets the cache to its "empty" condition.
    void Clear(void);

    //!Returns the property with the given ID.
    double Property(PropID id) const;

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

    //! Number of carbons.
    int m_numcarbon;

	//! Fragmentation flag.
	int m_frag;
    // Collision rate calculation particle properties.
    double m_dcolsqr;      //!< Collision diameter squared.
    double m_inv_dcol;     //!< Inverse collision diameter.
    double m_inv_dcolsqr;  //!< Inverse of the diameter squared.
    double m_inv_sqrtmass; //!< Inverse of the square-root of the mass.
    double m_d2_m_1_2;     //!< D^2 * M^-1/2.

    // Particle properties required for the weighted transition kernel (dcol used)
    double m_weight;            //!< Statistical weight
    double m_weight_mass;       //!< W * M
    double m_d_w;               //!< D * M
    double m_d2_w;              //!< D^2 * W
    double m_d_1_w;             //!< D^-1 * W
    double m_d_2_w;             //!< D^_2 * W
    double m_m_1_2_w;           //!< M^-1/2 * W
    double m_d2m_1_2_w;         //!< D^2 * M^-1/2 * W

    // SilicaPrimary properties
    
    //! Number of sites to which reaction is proportional
    double m_sites;

    //! Sintering rate of a particle
    double m_sinterrate;

    // Silicon model properties
    //! Coverage fraction of active sites (ratio of components)
    double m_coverage;

	int m_select;

	//Titania model properties
	double m_phaseterm;	//Anatase fraction ^ (2/3) * total composition
};
}

#endif
