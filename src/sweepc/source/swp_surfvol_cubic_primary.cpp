/*!
 * @file swp_surfvol_cubic_primary.cpp
 * @author Robert Patterson
 *
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2013 Robert I A Patterson

  File purpose:
    Implementation of a particle the approximates cuboidal crystals by
    storing their surface area and volume
    @brief Implementation of cuboidal crystals described by volume and surface area.

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

#include "swp_surfvol_cubic_primary.h"

#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"

#include <stdexcept>
#include <boost/random/poisson_distribution.hpp>

using namespace Sweep;
using namespace Sweep::AggModels;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
SurfVolCubicPrimary::SurfVolCubicPrimary(void)
{
}

// Initialising constructor.
SurfVolCubicPrimary::SurfVolCubicPrimary(double time, const Sweep::ParticleModel &model)
: Primary(time, model)
{
}

// Copy constructor.
SurfVolCubicPrimary::SurfVolCubicPrimary(const SurfVolCubicPrimary &copy)
{
    *this = copy;
}

// Stream-reading constructor.
SurfVolCubicPrimary::SurfVolCubicPrimary(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
SurfVolCubicPrimary::~SurfVolCubicPrimary()
{
    releaseMem();
}


// OPERATOR OVERLOADS.

// Assignment operator (Primary RHS).
SurfVolCubicPrimary &SurfVolCubicPrimary::operator=(const Primary &rhs)
{
    operator=(dynamic_cast<const SurfVolCubicPrimary&>(rhs));
    return *this;
}

// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType SurfVolCubicPrimary::AggID(void) const {return AggModels::SurfVolCubic_ID;}

// BASIC DERIVED PARTICLE PROPERTIES.

// Calculates the derived properties from the unique properties.  This
// function is broadly similar to the version in the spherical Primary
// class except that it uses the surface-volume model.  Therefore the
// surface area is not altered and the collision diameter is calculated
// using the arithmetic mean function.
//
// Mobility diameter:
// This calculation is based on the work of Rogak et al., 1993 Aer. Sci.
// Tech. 18:25-47, who give the calculation of dmob in the FM, SF and
// transition regime.
void SurfVolCubicPrimary::UpdateCache(void)
{
    // Store the old surface area, this value may be 0 if the particle is in
    // the process of initialisation
    double s = m_surf;

    // Pretend that the primary is spherical and set the cache
    // accordingly.
    Primary::UpdateCache();

    m_surf = std::max(s, m_surf);

   
    // Calculate diameters.
    m_dcol = (m_diam + sqrt(m_surf / PI)) * 0.5;

    // Calculate mobility diameter
    m_dmob = PP_Diameter();
    if (false) {
        // SF regime mobility diameter
        m_dmob *= 0.9 * sqrt(m_pmodel->GetFractDim() / (m_pmodel->GetFractDim() + 2));
        m_dmob *= pow(PP_Count(), (1.0/m_pmodel->GetFractDim()));
    } else {
        // FM regime mobility diameter
        m_dmob *= sqrt(0.802*(PP_Count()-1) + 1);
    }

    if (m_dmob < m_diam) m_dmob = m_diam;
}

// Returns the number of primary particles if the aggregate is assumed
// to consist of mono-sized primaries.
unsigned int SurfVolCubicPrimary::PP_Count(void) const
{
    // Note the minimum number of primary particles must be 1.
    return std::max(1u, (unsigned int)((m_surf * m_surf * m_surf) /
                                       (36.0 * PI * m_vol * m_vol)));
}

// Returns the primary particle diameter if the aggregate is assumed
// to consist of mono-sized primaries.
double SurfVolCubicPrimary::PP_Diameter(void) const
{
    // This should always be <= equiv. sphere diameter.
    return 6.0 * m_vol / m_surf;
}


// OPERATIONS.

// Adjusts the primary with the given composition and 
// tracker values changes n times.  If the particle cannot be adjust
// n times, then this function returns the number of times
// it was adjusted.
unsigned int SurfVolCubicPrimary::Adjust(const fvector &dcomp, const fvector &dvalues, rng_type &rng,
                                    unsigned int n)
{
    // Store initial surface and volume
    double surfOld = m_surf;
    double volOld  = m_vol;

    // Adjust the particle assuming that it is spherical.
    n = Primary::Adjust(dcomp, dvalues, rng, n);


    // Calculate change in volume.
    double dvol = 0.0;
    for (unsigned int i=0; i!=dcomp.size(); ++i) {
        dvol += dcomp[i] * m_pmodel->Components(i)->MolWt() / 
                m_pmodel->Components(i)->Density();
    }
    dvol *= (double)n / NA;

    // Calculate change in surface area.
    double invRadius = 0.0;
    if (dvol > 0.0) {
        // Inverse growth radius.
        invRadius = sqrt(4.0 * PI / surfOld);
    } else {
        // Inverse oxidation radius.    
        invRadius = surfOld / (3.0 * volOld);
    }

    // Save new surface area.
    double s = surfOld + (2.0 * dvol * invRadius);

    // Set correct surface area, which was incorrectly set by
    // Primary::Adjust.
    m_surf    = s;

    // This has a knock-on affect of changing the collision diameter.
    // Note, we can avoid recalling UpdateCache() here, because only
    // a couple of cached values will have changed.
    m_dcol = (m_diam + sqrt(m_surf / PI)) * 0.5;
    m_dmob = m_dcol;

    return n;
}

/*!
 * Combines this primary with another.
 *
 * \param[in]       rhs         Particle to add to current instance
 * \param[in,out]   rng         Random number generator
 *
 * \return      Reference to the current instance after rhs has been added
 */
SurfVolCubicPrimary &SurfVolCubicPrimary::Coagulate(const Primary &rhs, rng_type &rng)

{
    // Store the resultant surface area.
    double s = m_surf + rhs.SurfaceArea();

    // Perform the coagulation.
    Primary::Coagulate(rhs, rng);

    // One sixth of the surface area is assumed to be lost as a result of two faces, one
    // from each incoming particle covering each other.  More sophisticated formulae
    // could be implemented here.
    m_surf = s * 0.8333333333333333333;

    // This has a knock-on affect of changing the collision diameter.
    // Note, we can avoid recalling UpdateCache() here, because only
    // a couple of cached values will have changed.
    m_dcol = (m_diam + sqrt(m_surf / PI)) * 0.5;
    m_dmob = m_dcol;

    return *this;
}

/*!
 * Combines this primary with another.
 *
 * \param[in]       rhs         Particle to add to current instance
 * \param[in,out]   rng         Random number generator
 *
 * \return      Reference to the current instance after rhs has been added
 */
SurfVolCubicPrimary &SurfVolCubicPrimary::Fragment(const Primary &rhs, rng_type &rng)

{
    // Store the resultant surface area.
    double s = m_surf + rhs.SurfaceArea();

    // Perform the coagulation.
    Primary::Fragment(rhs, rng);

    // One sixth of the surface area is assumed to be lost as a result of two faces, one
    // from each incoming particle covering each other.  More sophisticated formulae
    // could be implemented here.
    m_surf = s * 0.8333333333333333333;

    // This has a knock-on affect of changing the collision diameter.
    // Note, we can avoid recalling UpdateCache() here, because only
    // a couple of cached values will have changed.
    m_dcol = (m_diam + sqrt(m_surf / PI)) * 0.5;
    m_dmob = m_dcol;

    return *this;
}

// READ/WRITE/COPY.

//! Returns a copy of the model data.
SurfVolCubicPrimary *const SurfVolCubicPrimary::Clone(void) const
{
    return new SurfVolCubicPrimary(*this);
}

// Writes the object to a binary stream.
void SurfVolCubicPrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output base class.
        Primary::Serialize(out);
    } else {
        throw std::invalid_argument("Output stream not ready "
                                    "(Sweep, SurfVolCubicPrimary::Serialize).");
    }
}

// Reads the object from a binary stream.
void SurfVolCubicPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    if (in.good()) {
        // Read base class.
        Primary::Deserialize(in, model);
    } else {
        throw std::invalid_argument("Input stream not ready "
                                    "(Sweep, SurfVolCubicPrimary::Deserialize).");
    }
}
