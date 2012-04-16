/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Primary class declared in the
    swp_primary.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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
#include "swp_model_factory.h"
#include "swp_cell.h"

#include <stdexcept>
#include <memory.h>

using namespace Sweep;

using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Primary::Primary(void)
: m_pmodel(NULL), m_createt(0.0), m_time(0.0), m_diam(0.0), m_dcol(0.0), 
  m_dmob(0.0), m_surf(0.0), m_vol(0.0), m_mass(0.0)
{
}

// Initialising constructor.
Primary::Primary(real time, const Sweep::ParticleModel &model)
{
    init();
    m_pmodel  = &model;
    
    // Set particle create and update times.
    m_createt = time;
    m_time    = time;

    // Resize component and tracker vectors to match number of
    // each in the particle model.
    m_comp.resize(model.ComponentCount(), 0.0);
    m_values.resize(model.TrackerCount(), 0.0);
}

// Copy constructor.
Primary::Primary(const Primary &copy) 
{
    init();
    *this = copy;
}

// Stream-reading constructor.
Primary::Primary(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
Primary::~Primary()
{
    releaseMem();
}


// OPERATOR OVERLOADS.

// Assignment operator.
Primary &Primary::operator=(const Primary &rhs)
{
    if (this != &rhs) {
        // Check if the RHS uses the same particle model before copying
        // components and tracker values.
        if (rhs.m_pmodel == m_pmodel) {
            // Copy composition.
            memcpy(&m_comp[0], &rhs.m_comp[0], sizeof(real)*m_comp.size());        
            // Copy tracker variables.
            memcpy(&m_values[0], &rhs.m_values[0], sizeof(real)*m_values.size());
        } else {
            // Set the particle model.
            m_pmodel = rhs.m_pmodel;
            // Copy components.
            m_comp.assign(rhs.m_comp.begin(), rhs.m_comp.end());
            // Copy tracker variables.
            m_values.assign(rhs.m_values.begin(), rhs.m_values.end());
        }

        // Copy primary time.
        m_createt = rhs.m_createt;
        m_time    = rhs.m_time;

        // Copy the derived properties.
        m_diam = rhs.m_diam;
        m_dcol = rhs.m_dcol;
        m_dmob = rhs.m_dmob;
        m_surf = rhs.m_surf;
        m_vol  = rhs.m_vol;
        m_mass = rhs.m_mass;
    }
    return *this;
}

// DEFINING PARTICLE MODEL.

// Returns the particle model used to create this primary.
const Sweep::ParticleModel *const Primary::ParticleModel(void) const
{
    return m_pmodel;
}


// PARTICLE COMPOSITION.

// Returns the composition vector.
const fvector &Primary::Composition() const 
{
	return m_comp;
}

// Returns the ith component value.  Returns 0.0 if i is invalid.
real Primary::Composition(unsigned int i) const
{
	if (i < m_comp.size()) {
		return m_comp[i];
	} else {
		return 0.0;
	}
}

// Sets the composition vector.
void Primary::SetComposition(const Sweep::fvector &comp)
{
	m_comp.assign(comp.begin(), comp.end());
}


// TRACKER VARIABLE VALUES.

// Returns the tracker value vector.
const fvector &Primary::Values() const
{
	return m_values;
}

// Returns the ith tracker variable value.  Returns 0.0 if i is invalid.
real Primary::Values(unsigned int i) const
{
	if (i < m_values.size()) {
		return m_values[i];
	} else {
		return 0.0;
	}
}

// Sets the values vector.
void Primary::SetValues(const fvector &vals)
{
    m_values.assign(vals.begin(), vals.end());
}

// Sets the ith trackervalue.
void Primary::SetValue(unsigned int i, real val)
{
	if (i < m_values.size()) {
		m_values[i] = val;
	}
}


// PRIMARY CREATE TIME.

// Returns the primary create time.
real Primary::CreateTime() const {return m_createt;}


// LAST UPDATE TIME.

// Returns the primary last update time.
real Primary::LastUpdateTime() const {return m_time;}

// Sets the last update time of the primary.
void Primary::SetTime(real t) {m_time = t;}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType Primary::AggID(void) const {return AggModels::Spherical_ID;}

// BASIC DERIVED PARTICLE PROPERTIES.

// Calculates the derived properties from the unique properties.
void Primary::UpdateCache(void)
{
    real m = 0.0;

    // Loop over composition and calculate mass and volume.
    m_mass = m_vol = 0.0;
    for (unsigned int i=0; i!=m_pmodel->ComponentCount(); ++i) {
        m = m_pmodel->Components(i)->MolWt() * m_comp[i] / NA;
        m_mass += m;
        if (m_pmodel->Components(i)->Density() > 0.0)
            m_vol  += m / m_pmodel->Components(i)->Density();
    }
    
    // Calculate other properties (of sphere).
    m_diam = pow(6.0 * m_vol / PI, ONE_THIRD);
    m_dcol = m_diam;
    m_dmob = m_diam;
    m_surf = PI * m_diam * m_diam;
}

// Returns the particle equivalent sphere diameter.
real Primary::SphDiameter(void) const {return m_diam;}

// Returns the collision diameter.
real Primary::CollDiameter(void) const {return m_dcol;}

// Rethrns the mobility diameter.
real Primary::MobDiameter(void) const {return m_dmob;}

// Returns the surface area.
real Primary::SurfaceArea(void) const {return m_surf;}

// Returns the equivalent sphere surface area, based
// on the volume.
real Primary::SphSurfaceArea(void) const
{
    return PI * m_diam * m_diam;
}

// Returns the volume.
real Primary::Volume(void) const {return m_vol;}

// Returns the mass.
real Primary::Mass(void) const {return m_mass;}

// Returns the property with the given ID.
real Primary::Property(const Sweep::PropID id) const
{
    switch (id) {
        case iDsph:      // Equivalent sphere diameter.
            return m_diam;
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
        default:
            return 0.0;
    }
}


// BASIC DERIVED PROPERTY OVERWRITES.

// Sets the spherical particle diameter
void Primary::SetSphDiameter(real diam) {m_diam = diam;}

// Sets the collision diameter of the particle.
void Primary::SetCollDiameter(real dcol) {m_dcol = dcol;}

// Sets the mobility diameter.
void Primary::SetMobDiameter(real dmob) {m_dmob = dmob;}

// Sets the surface area, subject to minimum spherical area condition.
void Primary::SetSurfaceArea(real surf) {m_surf = surf;}

// Sets the volume.
void Primary::SetVolume(real vol) {m_vol = vol;}

// Sets the mass.
void Primary::SetMass(real m) {m_mass = m;}

/*!
 * Check that this primary is a physically valid particle.  This currently
 * means checking that the amount of each component is within a range specifed
 * for that component.
 *
 *@return   True iff primary particle is valid as a physical particle
 */
bool Primary::IsValid() const {
    for(unsigned int i = 0; i < m_comp.size(); ++i) {
        // Check each component value, but stop as soon as an invalid value
        // is found
        if(!m_pmodel->Components(i)->IsValidValue(m_comp[i]))
            return false;
    }

    // Check one of the cached values, all the m_??? could be added here
    // if necessary.
    return (m_mass > 0);
}

// OPERATIONS.

// Adjusts the primary with the given composition and 
// tracker values changes n times.  If the particle cannot be adjust
// n times, then this function returns the number of times
// it was adjusted.
unsigned int Primary::Adjust(const fvector &dcomp, const fvector &dvalues, rng_type &rng, unsigned int n)
{
	unsigned int i = 0;

	// Add the components.
	for (i=0; i!=min(m_comp.size(),dcomp.size()); ++i) {
		m_comp[i] += dcomp[i] * (real)n;
	}

	// Add the tracker values.
	for (i=0; i!=min(m_values.size(),dvalues.size()); ++i) {
		m_values[i] += dvalues[i] * (real)n;
	}

    // Update property cache.
    Primary::UpdateCache();

    return n;
}

// Adjusts the particle n times for IntParticle reaction
unsigned int Primary::AdjustIntPar(const fvector &dcomp, const fvector &dvalues, rng_type &rng, unsigned int n)
{
	unsigned int i = 0;

	// Add the components.
	for (i=0; i!=min(m_comp.size(),dcomp.size()); ++i) {
		m_comp[i] += dcomp[i] * (real)n;
	}

	// Add the tracker values.
	for (i=0; i!=min(m_values.size(),dvalues.size()); ++i) {
		m_values[i] += dvalues[i] * (real)n;
	}

    // Update property cache.
    Primary::UpdateCache();

    return n;
}

/*!
 *  Combines this primary with another.
 *
 *  Note the very strange behaviour when the primaries do not
 *  have equal particle model pointers.  Users should also
 *  be very careful about not mixing different types that inherit
 *  from PrimaryParticle.
 *
 * \param[in]       rhs         Primary particle to add to current instance
 * \param[in,out]   rng         Random number generator
 *
 * \return      Reference to the current instance after rhs has been added
 */
Primary &Primary::Coagulate(const Primary &rhs, rng_type &rng)
{
    // Check if the RHS uses the same particle model.  If not, then
    // just use the assignment operator because you can't add apples 
    // and bananas!
    if (rhs.m_pmodel == m_pmodel) {
        // Add the components.
        for (unsigned int i=0; i!=min(m_comp.size(),rhs.m_comp.size()); ++i) {
            m_comp[i] += rhs.m_comp[i];
        }

        // Add the tracker values.
        for (unsigned int i=0; i!=min(m_values.size(),rhs.m_values.size()); ++i) {
            m_values[i] += rhs.m_values[i];
        }

        // Create time is the earliest time.
        m_createt = min(m_createt, rhs.m_createt);
        m_time    = min(m_time, rhs.m_time);
    } else {
        // Different particle models!
        *this = rhs;
        std::cerr << "Sweep::Primary::Coagulate called for particles with models " << m_pmodel
                  << " and " << rhs.m_pmodel << std::endl;
    }
    Primary::UpdateCache();
    return *this;
}

// This routine sinters the Primary for the given length of
// time using the provided sintering model.
void Primary::Sinter(real dt, Cell &sys,
                     const Processes::SinteringModel &model,
                     rng_type &rng,
                     real wt)
{
    // Spherical primaries don't sinter.
	
    return;
}

/*!
 * @brief       Gets the ratio of second to first component number
 *
 * A coarse estimation of the coverage fraction; as the amount of the second
 * component relative to the amount of the first. Returns 1.0 for single comp.
 *
 * @return      Ratio of number of second to first component
 */
real Primary::GetCoverageFraction() const
{
    real val = 1.0;

    if (m_pmodel->ComponentCount() > 1 && m_comp[0] > 0.0) {
        val = m_comp[1]/m_comp[0];
    }

    return val;
}

// READ/WRITE/COPY.

// Returns a copy of the model data.
Primary *const Primary::Clone(void) const
{
    return new Primary(*this);
}

// Writes the object to a binary stream.
void Primary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write number of components.
        unsigned int n = (unsigned int)m_comp.size();
        out.write((char*)&n, sizeof(n));

        // Write components.
        double val = 0.0;
        for (unsigned int i=0; i!=n; ++i) {
            val = (double)m_comp[i];
            out.write((char*)&val, sizeof(val));
        }

        // Write number of tracker values.
        n = (unsigned int)m_values.size();
        out.write((char*)&n, sizeof(n));

        // Write values.
        for (unsigned int i=0; i!=n; ++i) {
            val = (double)m_values[i];
            out.write((char*)&val, sizeof(val));
        }

        // Write create time.
        val = (double)m_createt;
        out.write((char*)&val, sizeof(val));

        // Write last update time.
        val = (double)m_time;
        out.write((char*)&val, sizeof(val));

        // Write equivalent sphere diameter.
        val = (double)m_diam;
        out.write((char*)&val, sizeof(val));

        // Write collision diameter.
        val = (double)m_dcol;
        out.write((char*)&val, sizeof(val));

        // Write mobility diameter.
        val = (double)m_dmob;
        out.write((char*)&val, sizeof(val));

        // Write surface area.
        val = (double)m_surf;
        out.write((char*)&val, sizeof(val));

        // Write volume.
        val = (double)m_vol;
        out.write((char*)&val, sizeof(val));

        // Write mass.
        val = (double)m_mass;
        out.write((char*)&val, sizeof(val));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Primary::Serialize).");
    }
}

// Reads the object from a binary stream.
void Primary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    releaseMem();
    m_pmodel = &model;

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double val = 0.0;

        switch (version) {
            case 0:
                // Read number of components.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read components.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_comp.push_back((real)val);
                }

                // Read number of tracker values.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read values.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_values.push_back((real)val);
                }

                // Read create time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_createt = (real)val;

                // Read last update time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_time = (real)val;

                // Read equivalent sphere diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_diam = (real)val;

                // Read collision diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_dcol = (real)val;

                // Read mobility diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_dmob = (real)val;

                // Read surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_surf = (real)val;

                // Read volume.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_vol = (real)val;

                // Read mass.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_mass = (real)val;
    
                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Primary::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Primary::Deserialize).");
    }
}


// DATA MANAGEMENT.

// Release all memory associated with object.
void Primary::releaseMem(void)
{
    m_comp.clear();
    m_values.clear();
}

// Initialisation routine.
void Primary::init(void)
{
    m_pmodel = NULL;
    m_createt = 0.0;
    m_time = 0.0;
    m_diam = 0.0; // Equivalent spherical diameter.
    m_dcol = 0.0; // Collision diameter.
    m_dmob = 0.0; // Mobility diameter.
    m_surf = 0.0; // Surface area.
    m_vol = 0.0;  // Volume.
    m_mass = 0.0; // Mass.
    releaseMem();
}
