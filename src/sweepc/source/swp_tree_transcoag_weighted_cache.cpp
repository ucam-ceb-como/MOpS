/*!
 * \file   swp_tree_weighted_cache.cpp
 * \author William J Menz
 *  Copyright (C) 2010 William J Menz
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief Implementation of the binary tree structure for the storage
 *        of particle properties for the weighted transition kernel.
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

 * This cache is one of a number of options that can be used in the binary tree
 * for particle properties.  At the time of writing these is a version that does
 * not include weights and a slimmed down cache which is only suitable for the
 * additive and constant coagulation kernels, and thus not useful for engineering
 * applications.  The memory consumed by the instances of the cache in the binary
 * tree makes up a significant proportion of the total memory used by the program
 * in the case of simple particle models, it is probably less significant when
 * individual primary particles are stored.
 */

#include "swp_tree_transcoag_weighted_cache.h"
#include "swp_PAH_primary.h"

#include "swp_particle.h"

#include <stdexcept>

// CONTRUCTORS AND DESTRUCTORS.

/*!
 * The cache is meant to be used to sum up quantities, so the natural
 * initial state is zero.
 */
Sweep::TreeTransCoagWeightedCache::TreeTransCoagWeightedCache()
: m_sphdiam(0.0)
, m_dcol(0.0)
, m_dmob(0.0)
, m_surf(0.0)
, m_vol(0.0)
, m_mass(0.0)
, m_numcarbon(0)
, m_frag(0), m_dcolsqr(0.0)
, m_inv_dcol(0.0)
, m_inv_dcolsqr(0.0)
, m_inv_sqrtmass(0.0)
, m_d2_m_1_2(0.0)
, m_weight(0.0)
, m_weight_mass(0.0)
, m_d_w(0.0)
, m_d2_w(0.0)
, m_d_1_w(0.0)
, m_d_2_w(0.0)
, m_m_1_2_w(0.0)
, m_d2m_1_2_w(0.0)
, m_sites(0.0)
, m_sinterrate(0.0)
, m_coverage(0.0)
, m_select(0.0)
, m_phaseterm(0.0)
{}

/*!
 * Calculation of basic and derived particle properties.
 * 
 * \param[in]   part    Particle supplying physical data with which to initialise
 *
 * This method queries the particle for basic quantities and then computes
 * derived quantities so that different caches can be written for different
 * kernels or in situations where different distributions are required.
 */
Sweep::TreeTransCoagWeightedCache::TreeTransCoagWeightedCache(const Sweep::Particle &part)
{
    /**
     * Quantities that must be provided by the particle.
     * Effectively this code defines an interface that the Particle class must
     * provide.
     */
	int incepRing;
	bool cut;

	// If the particle represents a single PAH that has fewer than the number of rings
	// specified by the inceptionThreshold, do not include this particle's parameters
	// in the binary tree summation. As the efficiency for the coagulation of this particle
	// with any other particle is zero. This check increases the efficiency of the majorant kernel
	if (part.Primary()->AggID() == AggModels::PAH_KMC_ID){
		incepRing = part.Primary()->ParticleModel()->inceptionThreshold();
		if (part.NumRings() < incepRing){
			cut = true;
		}
		else{
			cut = false;
		}
	}
	else{
		cut = false;
	}

    m_sphdiam   = part.SphDiameter();
	if (cut){
		m_dcol = 0.0;
	}
	else{
		m_dcol = part.CollDiameter();
	}
    m_dmob      = part.MobDiameter();
    m_surf      = part.SurfaceArea();
    m_vol       = part.Volume();
    m_mass      = part.Mass();
    m_numcarbon = part.NumCarbon();
    m_frag      = part.Frag();      //!< Fragmentation flag.

    // Derived quantites that are needed to the typical transition
    // regime coagulation kernel.
	m_dcolsqr = m_dcol * m_dcol;
	if (m_dcol != 0){
		m_inv_dcol = 1.0 / m_dcol;
		m_inv_dcolsqr = 1.0 / m_dcolsqr;
		m_inv_sqrtmass = 1.0 / std::sqrt(m_mass);
		m_d2_m_1_2 = m_dcolsqr * m_inv_sqrtmass;
	}
	else{
		m_inv_dcol = 0.0;
		m_inv_dcolsqr = 0.0;
		m_inv_sqrtmass = 0.0;
		m_d2_m_1_2 = 0.0;
	}

    // Quantities associated with statistical weighting
	if (cut){
		m_weight = 0.0;
	}
	else{
		m_weight = part.getStatisticalWeight();
	}
	m_weight_mass = m_weight * m_mass;
	m_d_w = m_dcol * m_weight;
	m_d2_w = m_dcolsqr * m_weight;
	m_d_1_w = m_inv_dcol * m_weight;
	m_d_2_w = m_inv_dcolsqr * m_weight;
	m_m_1_2_w = m_inv_sqrtmass * m_weight;
	m_d2m_1_2_w = m_d2_m_1_2 * m_weight;

    // Silica parameters
    m_sites =       part.GetSites();
    m_sinterrate =  part.GetSintRate();

    // Silicon parameters
    m_coverage = part.GetCoverageFraction();

	if (cut){
		m_select = 0;
	}
	else{
		m_select = m_weight;
	}

	//Titania parameters
	m_phaseterm = part.GetPhaseTerm();
}

// OPERATOR OVERLOADS.

/**
 * Addition-assignment operator (TreeWeightedCache RHS).
 * This function is not used to coagulate particles, rather
 * it is used to sum up particle properties in the ensemble binary tree.
 */
Sweep::TreeTransCoagWeightedCache &Sweep::TreeTransCoagWeightedCache::operator+=(const TreeTransCoagWeightedCache &rhs)
{
    //! Sum cache variables.
    m_sphdiam      += rhs.m_sphdiam;
    m_dcol         += rhs.m_dcol;
    m_dmob         += rhs.m_dmob;
    m_surf         += rhs.m_surf;
    m_vol          += rhs.m_vol;
    m_mass         += rhs.m_mass;
    m_numcarbon    += rhs.m_numcarbon;
    m_frag         += rhs.m_frag;
    m_numcarbon += rhs.m_numcarbon;
    m_dcolsqr      += rhs.m_dcolsqr;
    m_inv_dcol     += rhs.m_inv_dcol;
    m_inv_dcolsqr  += rhs.m_inv_dcolsqr;
    m_inv_sqrtmass += rhs.m_inv_sqrtmass;
    m_d2_m_1_2     += rhs.m_d2_m_1_2;
    m_weight       += rhs.m_weight;
    m_weight_mass  += rhs.m_weight_mass;
    m_d_w          += rhs.m_d_w;
    m_d2_w         += rhs.m_d2_w;
    m_d_1_w        += rhs.m_d_1_w;
    m_d_2_w        += rhs.m_d_2_w;
    m_m_1_2_w      += rhs.m_m_1_2_w;
    m_d2m_1_2_w    += rhs.m_d2m_1_2_w;
    m_sites        += rhs.m_sites;
    m_sinterrate   += rhs.m_sinterrate;
    m_coverage     += rhs.m_coverage;
	m_select       += rhs.m_select;
	m_phaseterm    += rhs.m_phaseterm;

    return *this;
}

//! Aaddition operator (TreeWeightedCache RHS).
const Sweep::TreeTransCoagWeightedCache Sweep::TreeTransCoagWeightedCache::operator+(const TreeTransCoagWeightedCache &rhs) const {
    // Use copy constructor and += operator to define.
    return TreeTransCoagWeightedCache(*this) += rhs;
}

// CLEAR THE PARTICLE CACHE.

//! Resets the particle cache to its "empty" condition.
void Sweep::TreeTransCoagWeightedCache::Clear(void)
{
    //! Clear derived properties.
    m_sphdiam      = 0.0;
    m_dcol         = 0.0;
    m_dmob         = 0.0;
    m_surf         = 0.0;
    m_vol          = 0.0;
    m_mass         = 0.0;
    m_numcarbon    = 0;
    m_frag         = 0;
    m_numcarbon = 0;
    m_dcolsqr      = 0.0;
    m_inv_dcol     = 0.0;
    m_inv_dcolsqr  = 0.0;
    m_inv_sqrtmass = 0.0;
    m_d2_m_1_2     = 0.0;
    m_weight       = 0.0;
    m_weight_mass  = 0.0;
    m_d_w          = 0.0;
    m_d2_w         = 0.0;
    m_d_1_w        = 0.0;
    m_d_2_w        = 0.0;
    m_m_1_2_w      = 0.0;
    m_d2m_1_2_w    = 0.0;
    m_sites        = 0.0,
    m_sinterrate   = 0.0;
    m_coverage     = 0.0;
	m_select = 0.0;
	m_phaseterm = 0.0;
}

/**
 * Returns one of the values stored in the cache.
 * 
 * \param[in]    id    PropID corresponding to particle property of interest.
 * \return       value of the property requested.
 */
double Sweep::TreeTransCoagWeightedCache::Property(PropID id) const
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

        //! Number of carbon atoms.
        case iNumCarbon:
            return m_numcarbon;

        //! Fragmentation flag.
        case iFrag:
            return m_frag;

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
        case iW:
            return m_weight;
        case iWM:
            return m_weight_mass;
        case iDW:
            return m_d_w;
        case iD2W:
            return m_d2_w;
        case iD_1W:
            return m_d_1_w;
        case iD_2W:
            return m_d_2_w;
        case iM_1_2W:
            return m_m_1_2_w;
        case iD2_M_1_2W:
            return m_d2m_1_2_w;
        case iASN:
            return m_sites;
        case iSintRate:
            return m_sinterrate;
        case iCoverage:
            return m_coverage;
        case iFS:
            throw std::logic_error("Free surface no longer cached (TreeWeightedCache::Property)");
            return 0.0;   
		case iUniform1:
			return m_select;
		case -1:
            // Special case property, used to select particles
            // uniformly.
            return 1.0;
		case iAn_2_3_comp:
			return m_phaseterm;
        default:
            return 0.0;
    }
}

