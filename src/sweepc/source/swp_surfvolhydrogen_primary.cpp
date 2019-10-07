/*!
 * @file swp_surfvolhydrogen_primary.cpp
 * @author Robert I A Patterson Robert.Patterson@wias-berlin.de
 *
  
  Copyright (C) 2012 Robert I A Patterson

  File purpose:
    Implementation of the SurfVolHydrogenPrimary class declared in the
    swp_surfvol_primary.h header file.
    @brief Implementation of particle described by its volume and surface area.

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

#include "swp_primary.h"
#include "swp_surfvolhydrogen_primary.h"
#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"

#include <stdexcept>
#include <boost/random/poisson_distribution.hpp>

using namespace Sweep;
using namespace Sweep::AggModels;


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
SurfVolHydrogenPrimary::SurfVolHydrogenPrimary(void)
: iH(1u)
{
}

// Initialising constructor.
SurfVolHydrogenPrimary::SurfVolHydrogenPrimary(double time, const Sweep::ParticleModel &model)
: Primary(time, model)
, iH(1u)
{
}

// Copy constructor.
SurfVolHydrogenPrimary::SurfVolHydrogenPrimary(const SurfVolHydrogenPrimary &copy)
: iH(1u)
{
    *this = copy;
}

// Stream-reading constructor.
SurfVolHydrogenPrimary::SurfVolHydrogenPrimary(std::istream &in, const Sweep::ParticleModel &model)
: iH(1u)
{
    Deserialize(in, model);
}

// Default destructor.
SurfVolHydrogenPrimary::~SurfVolHydrogenPrimary()
{
    releaseMem();
}


// OPERATOR OVERLOADS.

// Assignment operator (Primary RHS).
SurfVolHydrogenPrimary &SurfVolHydrogenPrimary::operator=(const Primary &rhs)
{
    // Attempt to cast the primary to a SurfVolHydrogenPrimary.  This will
    // throw an exception if the cast fails.
    operator=(dynamic_cast<const SurfVolHydrogenPrimary&>(rhs));

    return *this;
}

// Assignment operator (SurfVolHydrogenPrimary RHS).
SurfVolHydrogenPrimary &SurfVolHydrogenPrimary::operator=(const SurfVolHydrogenPrimary &rhs)
{
    // First copy everything for spherical primaries.
    Primary::operator=(rhs);

    return *this;
}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType SurfVolHydrogenPrimary::AggID(void) const {return AggModels::SurfVolHydrogen_ID;}


// BASIC DERIVED PARTICLE PROPERTIES.

// Calculates the derived properties from the unique properties.  This
// function is broadly similar to the version in the spherical Primary
// class except that it uses the surface-volume model.  Therefore the
// surface area is not altered and the collision diameter is calculated
// using the arithmetic mean function.
void SurfVolHydrogenPrimary::UpdateCache(void)
{
    // Store the correct surface area.
    const double s = m_surf;

    // Pretend that the primary is spherical and set the cache
    // accordingly.  This will set m_surf the the surface are
    // of an equivalent volume sphere.
    Primary::UpdateCache();

    // The surface area is now set to that of a sphere. Set 
    // correct surface area (minimum is spherical surface area).
    m_surf = std::max(s, m_surf);

   
    // Calculate diameters.
    m_dcol = (m_diam + sqrt(m_surf / PI)) * 0.5;
    m_dmob = m_dcol; // TODO:  Correct expression for Dmob.
}

// Returns the equivalent spherical particle surface area.
double SurfVolHydrogenPrimary::SphSurfaceArea(void) const {return PI * std::pow(6.0 * m_vol / PI, TWO_THIRDS);}

// Returns the number of primary particles if the aggregate is assumed
// to consist of mono-sized primaries.
unsigned int SurfVolHydrogenPrimary::PP_Count(void) const
{
    // Note the minimum number of primary particles must be 1.
    return std::max(1u, (unsigned int)((m_surf * m_surf * m_surf) /
                                       (36.0 * PI * m_vol * m_vol)));
}

// Returns the primary particle diameter if the aggregate is assumed
// to consist of mono-sized primaries.
double SurfVolHydrogenPrimary::PP_Diameter(void) const
{
    // This should always be <= equiv. sphere diameter.
    return 6.0 * m_vol / m_surf;
}


// OPERATIONS.

// Adjusts the primary with the given composition and 
// tracker values changes n times.  If the particle cannot be adjust
// n times, then this function returns the number of times
// it was adjusted.
unsigned int SurfVolHydrogenPrimary::Adjust(const fvector &dcomp, const fvector &dvalues, rng_type &rng,
                                    unsigned int n)
{
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
        invRadius = sqrt(4.0 * PI / m_surf);
    } else {
        // Inverse oxidation radius.    
        invRadius = m_surf / (3.0 * m_vol);
    }

    // Save new and old surface areas
    const double oldS = m_surf;
    const double s = m_surf + (2.0 * dvol * invRadius);

    // Adjust the particle assuming that it is spherical.
    Primary::Adjust(dcomp, dvalues, rng, n);

    // Maintain the concentration of H per unit surface area during oxidation
    if(dvol < 0.0) {
        fvector comp = Composition();
        //comp[iH] = static_cast<int>(s * comp[iH] / oldS + 0.5);
        comp[iH] = s * comp[iH] / oldS;
        SetComposition(comp);
    }


    // Set correct surface area, which was incorrectly set by
    // Primary::Adjust.
    m_surf    = std::max(s, m_surf);

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
SurfVolHydrogenPrimary &SurfVolHydrogenPrimary::Coagulate(const Primary &rhs, rng_type &rng)

{
    // Store the resultant surface area.
    const double s = m_surf + rhs.SurfaceArea();

    // Perform the coagulation.
    Primary::Coagulate(rhs, rng);

    // The spherical particle Coagulate() routine has set the
    // surface area incorrectly.  We now replace the surface area
    // to the correct point-contact value.
    m_surf    = std::max(m_surf, s);

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
SurfVolHydrogenPrimary &SurfVolHydrogenPrimary::Fragment(const Primary &rhs, rng_type &rng)

{
    // Store the resultant surface area.
    const double s = m_surf + rhs.SurfaceArea();

    // Perform the coagulation.
    Primary::Fragment(rhs, rng);

    // The spherical particle Coagulate() routine has set the
    // surface area incorrectly.  We now replace the surface area
    // to the correct point-contact value.
    m_surf    = std::max(m_surf, s);

    // This has a knock-on affect of changing the collision diameter.
    // Note, we can avoid recalling UpdateCache() here, because only
    // a couple of cached values will have changed.
    m_dcol = (m_diam + sqrt(m_surf / PI)) * 0.5;
    m_dmob = m_dcol;

    return *this;
}

// This routine sinters the Primary for the given length of
// time using the provided sintering model.
void SurfVolHydrogenPrimary::Sinter(double dt, Cell &sys,
                            const Processes::SinteringModel &model,
                            rng_type &rng,
                            double wt)
{
  throw std::runtime_error("sintering not implemented at present for SurfVolHydrogenPrimary");
}


// READ/WRITE/COPY.

// Returns a copy of the model data.
SurfVolHydrogenPrimary *const SurfVolHydrogenPrimary::Clone(void) const
{
    return new SurfVolHydrogenPrimary(*this);
}

/*!
 * \return      Number of surface Hydrogens (active sites)
 */
double SurfVolHydrogenPrimary::GetSites() const
{
    return Composition()[iH];
}


// Writes the object to a binary stream.
void SurfVolHydrogenPrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output base class.
        Primary::Serialize(out);
    } else {
        throw std::invalid_argument("Output stream not ready (Sweep, SurfVolHydrogenPrimary::Serialize).");
    }
}

// Reads the object from a binary stream.
void SurfVolHydrogenPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                // Read base class.
                Primary::Deserialize(in, model);
                break;
            default:
                throw std::runtime_error("Serialized version number is invalid (Sweep, SurfVolHydrogenPrimary::Deserialize).");
        }
    } else {
        throw std::invalid_argument("Input stream not ready (Sweep, SurfVolHydrogenPrimary::Deserialize).");
    }
}
