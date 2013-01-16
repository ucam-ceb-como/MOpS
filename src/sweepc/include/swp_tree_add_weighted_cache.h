/*!
 * \file   swp_tree_add_weighted_cache.h
 * \author Robert I A Patterson
 *  Copyright (C) 2011 Robert I A Patterson.
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

#ifndef SWEEP_TREE_ADD_WEIGHTED_CACHE_H
#define SWEEP_TREE_ADD_WEIGHTED_CACHE_H

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
 * Currently has 3 data members
 *
 * This class is (was originally a slimmed down version of TreeWeightedCache for
 * use with the additive and constant coagulation kernels.
 */
class TreeAddWeightedCache
{
public:
    //! Initialise with zeros
    TreeAddWeightedCache();

    //! Initialise with values from a particle
    TreeAddWeightedCache(const Sweep::Particle &part);

    //! Increment elementwise
    TreeAddWeightedCache &operator+=(const TreeAddWeightedCache &rhs);

    //! Add together the values from two caches element by element
    const TreeAddWeightedCache operator+(const TreeAddWeightedCache &rhs) const;

    //! Resets the cache to its "empty" condition.
    void Clear(void);

    // Returns the property with the given ID.
    double Property(PropID id) const;

private:
    //! Mass.
    double m_mass;

    //! Statistical weight
    double m_weight;

    //! Product of statistical weight and mass
    double m_weight_mass;
};
}

#endif
