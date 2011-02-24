/*!
 * \file   swp_secondary_primary_coag.cpp
 * \author Robert I A Patterson
 *  Copyright (C) 2010 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Process for free molecular coagulation of small (secondary) particles with particles from the main population

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

#include "swp_tree_cache.h"

#include "swp_particle.h"

// CONTRUCTORS AND DESTRUCTORS.

/*!
 * The cache is meant to be used to sum up quantities, so the natural
 * initial state is zero.
 */
Sweep::TreeCache::TreeCache()
: m_sphdiam(0.0)
, m_dcol(0.0)
, m_dmob(0.0)
, m_surf(0.0)
, m_vol(0.0)
, m_mass(0.0)
, m_dcolsqr(0.0)
, m_inv_dcol(0.0)
, m_inv_dcolsqr(0.0)
, m_inv_sqrtmass(0.0)
, m_d2_m_1_2(0.0)
, m_freesurface(0.0)
, m_numcarbon(0)
{}

/*!
 * @param[in]	part	Particle supplying physical data with which to initialise
 *
 * This method queries the particle for basic quantities and then computes
 * derived quantities so that different caches can be written for different
 * kernels or in situations where different distributions are required.
 */
Sweep::TreeCache::TreeCache(const Sweep::Particle &part)
{
    // Quantities that must be provided by the particle
    // Effectively this code defines an interface that the Particle
    // class must provide.
    m_sphdiam = part.SphDiameter();
    m_dcol    = part.CollDiameter();
    m_dmob    = part.MobDiameter();
    m_surf    = part.SurfaceArea();
    m_vol     = part.Volume();
    m_mass    = part.Mass();

    // The particle does not currently provide this data (although it stores it)
    m_freesurface = 0.0;
    m_numcarbon = 0;


    // Derived quantites that are needed to the typical transition
    // regime coagulation kernel.
    m_dcolsqr      = m_dcol * m_dcol;
    m_inv_dcol     = 1.0 / m_dcol;
    m_inv_dcolsqr  = 1.0 / m_dcolsqr;
    m_inv_sqrtmass = 1.0 / std::sqrt(m_mass);
    m_d2_m_1_2     = m_dcolsqr * m_inv_sqrtmass;
}

// OPERATOR OVERLOADS.


// Addition-assignment operator (TreeCache RHS).
//   This function is not used to coagulate particles, rather
//   it is used to sum up particle properties in the ensemble binary tree.
//
Sweep::TreeCache &Sweep::TreeCache::operator+=(const TreeCache &rhs)
{
    // Sum cache variables.
    m_sphdiam += rhs.m_sphdiam;
    m_dcol += rhs.m_dcol;
    m_dmob += rhs.m_dmob;
    m_surf += rhs.m_surf;
    m_vol  += rhs.m_vol;
    m_mass += rhs.m_mass;
    m_dcolsqr      += rhs.m_dcolsqr;
    m_inv_dcol     += rhs.m_inv_dcol;
    m_inv_dcolsqr  += rhs.m_inv_dcolsqr;
    m_inv_sqrtmass += rhs.m_inv_sqrtmass;
    m_d2_m_1_2     += rhs.m_d2_m_1_2;

    m_freesurface += rhs.m_freesurface;

    m_numcarbon += rhs.m_numcarbon;

    return *this;
}

// Addition operator (TreeCache RHS).
const Sweep::TreeCache Sweep::TreeCache::operator+(const TreeCache &rhs) const {
    // Use copy constructor and += operator to define.
    return TreeCache(*this) += rhs;
}

// CLEAR THE PARTICLE CACHE.

// Resets the particle cache to its "empty" condition.
void Sweep::TreeCache::Clear(void)
{
    // Clear derived properties.
    m_sphdiam = 0.0;
    m_dcol = 0.0;
    m_dmob = 0.0;
    m_surf = 0.0;
    m_vol  = 0.0;
    m_mass = 0.0;
    m_dcolsqr      = 0.0;
    m_inv_dcol     = 0.0;
    m_inv_dcolsqr  = 0.0;
    m_inv_sqrtmass = 0.0;
    m_d2_m_1_2     = 0.0;

    m_freesurface = 0.0;
    m_numcarbon = 0;
}



// BASIC DERIVED PARTICLE PROPERTIES.

// Returns the particle equivalent sphere diameter.
Sweep::real Sweep::TreeCache::SphDiameter(void) const {return m_sphdiam;}

// Returns the collision diameter.
Sweep::real Sweep::TreeCache::CollDiameter(void) const {return m_dcol;}

// Rethrns the mobility diameter.
Sweep::real Sweep::TreeCache::MobDiameter(void) const {return m_dmob;}

// Returns the surface area.
Sweep::real Sweep::TreeCache::SurfaceArea(void) const {return m_surf;}

// Returns the equivalent sphere surface area, based
// on the volume.
Sweep::real Sweep::TreeCache::SphSurfaceArea(void) const
{
    return PI * pow(m_vol * 6.0 / PI, TWO_THIRDS);
}

// Returns the volume.
Sweep::real Sweep::TreeCache::Volume(void) const {return m_vol;}

// Returns the mass.
Sweep::real Sweep::TreeCache::Mass(void) const {return m_mass;}


// Returns the property with the given ID.
Sweep::real Sweep::TreeCache::Property(PropID id) const
{
    switch (id) {
        case iDsph:      // Equivalent sphere diameter.
            return m_sphdiam;
        case iDcol:   // Collision diameter.
            return m_dcol;
        case iDmob:   // Mobility diameter.
            return m_dmob;
        case iS:      // Surface area.
            return m_surf;
        case iV:      // Volume.
            return m_vol;
        case iM:      // Mass.
            return m_mass;
        // Collision rate properties:
        case iD2:
            return m_dcolsqr;
        case iD_1:
            return m_inv_dcol;
        case iD_2:
            return m_inv_dcolsqr;
        case iM_1_2:
            return m_inv_sqrtmass;
        case iD2_M_1_2:
            return m_d2_m_1_2;
		case iFS:
			return m_freesurface;
		case iNumCarbon:
		    return m_numcarbon;
        case -1:
            // Special case property, used to select particles
            // uniformly.
            return 1.0;
        default:
            return 0.0;
    }
}

