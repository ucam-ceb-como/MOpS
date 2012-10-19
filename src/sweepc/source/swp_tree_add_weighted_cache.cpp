/*!
 * \file   swp_tree_add_weighted_cache.cpp
 * \author Robert I A Patterson
 *  Copyright (C) 2010 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
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

#include "swp_tree_add_weighted_cache.h"

#include "swp_particle.h"

#include <stdexcept>

// CONTRUCTORS AND DESTRUCTORS.

/*!
 * The cache is meant to be used to sum up quantities, so the natural
 * initial state is zero.
 */
Sweep::TreeAddWeightedCache::TreeAddWeightedCache()
: m_mass(0.0)
, m_weight(0.0)
, m_weight_mass(0.0)
{}

/*!
 * @param[in]	part	Particle supplying physical data with which to initialise
 *
 * This method queries the particle for basic quantities and then computes
 * derived quantities so that different caches can be written for different
 * kernels or in situations where different distributions are required.
 */
Sweep::TreeAddWeightedCache::TreeAddWeightedCache(const Sweep::Particle &part)
{
    // Quantities that must be provided by the particle
    // Effectively this code defines an interface that the Particle
    // class must provide.
    m_mass    = part.Mass();

    m_weight = part.getStatisticalWeight();
    m_weight_mass = m_weight * m_mass;
}

// OPERATOR OVERLOADS.


// Addition-assignment operator (TreeAddWeightedCache RHS).
//   This function is not used to coagulate particles, rather
//   it is used to sum up particle properties in the ensemble binary tree.
//
Sweep::TreeAddWeightedCache &Sweep::TreeAddWeightedCache::operator+=(const TreeAddWeightedCache &rhs)
{
    // Sum cache variables.
    m_mass += rhs.m_mass;
    m_weight       += rhs.m_weight;
    m_weight_mass  += rhs.m_weight_mass;

    return *this;
}

// Addition operator (TreeAddWeightedCache RHS).
const Sweep::TreeAddWeightedCache Sweep::TreeAddWeightedCache::operator+(const TreeAddWeightedCache &rhs) const {
    // Use copy constructor and += operator to define.
    return TreeAddWeightedCache(*this) += rhs;
}

// CLEAR THE PARTICLE CACHE.

// Resets the particle cache to its "empty" condition.
void Sweep::TreeAddWeightedCache::Clear(void)
{
    // Clear derived properties.
    m_mass = 0.0;
    m_weight       = 0.0;
    m_weight_mass  = 0.0;
}

/*!
 *  Returns one of the values stored in the cache
 */
double Sweep::TreeAddWeightedCache::Property(PropID id) const
{
    switch (id) {
        case iM:      // Mass.
            return m_mass;
        case iW:
        	return m_weight;
        case iWM:
        	return m_weight_mass;
        case -1:
            // Special case property, used to select particles
            // uniformly.
            return 1.0;
        default:
            return 0.0;
    }
}

