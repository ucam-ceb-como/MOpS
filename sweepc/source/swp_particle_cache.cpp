/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ParticleCache class declared in the
    swp_particle_cache.h header file.

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

#include "swp_particle_cache.h"
#include "swp_particle_model.h"
#include "swp_model_factory.h"
#include "swp_primary.h"

#include <vector>
#include <stdexcept>
#include <memory.h>

using namespace Sweep;

using namespace std;

// CONTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
ParticleCache::ParticleCache()
{
	init();
}

// Initialising constructor (using Primary particle).
ParticleCache::ParticleCache(const Primary &pri)
  : m_diam(pri.SphDiameter()),
    m_dcol(pri.CollDiameter()),
    m_dmob(pri.MobDiameter()),
    m_surf(pri.SurfaceArea()),
    m_vol (pri.Volume()),
    m_mass(pri.Mass()),
    m_dcolsqr(m_dcol * m_dcol),
    m_inv_dcol(1.0 / m_dcol),
    m_inv_dcolsqr(1.0 / m_dcolsqr),
    m_inv_sqrtmass(1.0 / sqrt(m_mass)),
    m_d2_m_1_2(m_dcolsqr * m_inv_sqrtmass),
    m_numsubpart(0),
    m_numcarbon(static_cast<unsigned int> (pri.Composition(iNumCarbon)))
{}

// Stream-reading constructor.
ParticleCache::ParticleCache(std::istream &in, const Sweep::ParticleModel &model)
{
    init();
    Deserialize(in, model);
}

// OPERATOR OVERLOADS.

// Assignment operator (Primary RHS).
ParticleCache &ParticleCache::operator=(const Sweep::Primary &rhs)
{
    // Copy the derived properties.
    m_diam = rhs.SphDiameter();
    m_dcol = rhs.CollDiameter();
    m_dmob = rhs.MobDiameter();
    m_surf = rhs.SurfaceArea();
    m_vol  = rhs.Volume();
    m_mass = rhs.Mass();

    // Calculate collision-rate properties.
    m_dcolsqr      = m_dcol * m_dcol;
    m_inv_dcol     = 1.0 / m_dcol;
    m_inv_dcolsqr  = 1.0 / m_dcolsqr;
    m_inv_sqrtmass = 1.0 / sqrt(m_mass);
    m_d2_m_1_2     = m_dcolsqr * m_inv_sqrtmass;

    m_numsubpart = 0;
    m_numcarbon = static_cast<unsigned int> (rhs.Composition(iNumCarbon));

    return *this;
}

// Addition-assignment operator (ParticleCache RHS).
//   This function is not used to coagulate particles, rather
//   it is used to sum up particle properties in the ensemble binary tree.
//
// It is the responsibility of the caller to make sure that both particle
// cache objects use the same particle model.  It is not meaningful to add
// caches for different particle models.
ParticleCache &ParticleCache::operator+=(const Sweep::ParticleCache &rhs)
{
    // Sum cache variables.
    m_diam += rhs.m_diam;
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
    m_numsubpart += rhs.m_numsubpart;

    m_numcarbon += rhs.m_numcarbon;

    return *this;
}

// Addition-assignment operator (Primary RHS).
//   This function is not used to coagulate particles, rather
//   it is used to sum up particle properties in the ensemble binary tree.
ParticleCache &ParticleCache::operator+=(const Sweep::Primary &rhs)
{
        // Get cache properties from primary.
        real diam = rhs.SphDiameter();
        real dcol = rhs.CollDiameter();
        real dmob = rhs.MobDiameter();
        real surf = rhs.SurfaceArea();
        real vol  = rhs.Volume();
        real mass = rhs.Mass();

        // Sum cache variables.
        m_diam += diam;
        m_dcol += dcol;
        m_dmob += dmob;
        m_surf += surf;
        m_vol  += vol;
        m_mass += mass;

        // Sum collision-rate properties.
        m_dcolsqr      += dcol * dcol;
        m_inv_dcol     += 1.0 / dcol;
        m_inv_dcolsqr  += 1.0 / (dcol * dcol);
        m_inv_sqrtmass += 1.0 / sqrt(mass);
        m_d2_m_1_2     += (dcol * dcol) / sqrt(mass);


        m_numcarbon += static_cast<unsigned int>(rhs.Composition(iNumCarbon));

    return *this;
}

// Addition operator (ParticleCache RHS).
const ParticleCache ParticleCache::operator+(const Sweep::ParticleCache &rhs) const {
	// Use copy constructor and += operator to define.
	return ParticleCache(*this) += rhs;
}

// Addition operator (Primary RHS).
const ParticleCache ParticleCache::operator+(const Sweep::Primary &rhs) const {
	// Use copy constructor and += operator to define.
	return ParticleCache(*this) += rhs;
}


// CLEAR THE PARTICLE CACHE.

// Resets the particle cache to its "empty" condition.
void ParticleCache::Clear(void)
{
    // Clear derived properties.
    m_diam = 0.0;
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
    m_numsubpart = 0;
    m_numcarbon = 0;
}

// BASIC DERIVED PARTICLE PROPERTIES.

// Returns the particle equivalent sphere diameter.
real ParticleCache::SphDiameter(void) const {return m_diam;}

// Returns the collision diameter.
real ParticleCache::CollDiameter(void) const {return m_dcol;}

// Rethrns the mobility diameter.
real ParticleCache::MobDiameter(void) const {return m_dmob;}

// Returns the surface area.
real ParticleCache::SurfaceArea(void) const {return m_surf;}

// Returns the equivalent sphere surface area, based
// on the volume.
real ParticleCache::SphSurfaceArea(void) const
{
    return PI * pow(m_vol * 6.0 / PI, TWO_THIRDS);
}

// Returns the volume.
real ParticleCache::Volume(void) const {return m_vol;}

// Returns the mass.
real ParticleCache::Mass(void) const {return m_mass;}


// COLLISION RATE EXPRESSION PROPERTIES.

// Collision diameter squared (cm2).
real ParticleCache::CollDiamSquared() const {return m_dcolsqr;}

// Inverse collision diameter (cm-1).
real ParticleCache::InvCollDiam() const {return m_inv_dcol;}

// Inverse squared collision diameter (cm-2).
real ParticleCache::InvCollDiamSquared() const {return m_inv_dcolsqr;}

// Inverse of square root of mass (g-1/2).
real ParticleCache::InvSqrtMass() const {return m_inv_sqrtmass;}

// Collision diameter squared times the inverse square root of mass.
real ParticleCache::CollDiamSqrdInvSqrtMass() const {return m_d2_m_1_2;}


// Returns the property with the given ID.
real ParticleCache::Property(PropID id) const
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


// BASIC DERIVED PROPERTY OVERWRITES.




// Sets the spherical particle diameter
void ParticleCache::SetSphDiameter(real diam) {m_diam = diam;}

// Sets the collision diameter of the particle.
void ParticleCache::SetCollDiameter(real dcol)
{
    m_dcol         = dcol;
    m_dcolsqr      = m_dcol * m_dcol;
    m_inv_dcol     = 1.0 / m_dcol;
    m_inv_dcolsqr  = 1.0 / m_dcolsqr;
    m_inv_sqrtmass = 1.0 / sqrt(m_mass);
    m_d2_m_1_2     = m_dcolsqr * m_inv_sqrtmass;
}

// Sets the mobility diameter.
void ParticleCache::SetMobDiameter(real dmob) {m_dmob = dmob;}

// Sets the surface area, subject to minimum spherical area condition.
void ParticleCache::SetSurfaceArea(real surf) {m_surf = surf;}

// Sets the volume.
void ParticleCache::SetVolume(real vol) {m_vol = vol;}

// Sets the mass.
void ParticleCache::SetMass(real m)
{
    m_mass         = m;
    m_inv_sqrtmass = 1.0 / sqrt(m_mass);
    m_d2_m_1_2     = m_dcolsqr * m_inv_sqrtmass;
}


// READ/WRITE/COPY.

// Writes the object to a binary stream.
void ParticleCache::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        double val = 0.0;

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

        // Output the collision diameter squared.
        double v = (double)m_dcolsqr;
        out.write((char*)&v, sizeof(v));

        // Output the inverse collision diameter.
        v = (double)m_inv_dcol;
        out.write((char*)&v, sizeof(v));

        // Output the inverse collision diameter squared.
        v = (double)m_inv_dcolsqr;
        out.write((char*)&v, sizeof(v));

        // Output the inverse square root of mass.
        v = (double)m_inv_sqrtmass;
        out.write((char*)&v, sizeof(v));

        // Output the D^2 / sqrt(M).
        v = (double)m_d2_m_1_2;
        out.write((char*)&v, sizeof(v));

        out.write(reinterpret_cast<const char*>(&m_freesurface), sizeof(m_freesurface));
        out.write(reinterpret_cast<const char*>(&m_numsubpart), sizeof(m_numsubpart));
        out.write(reinterpret_cast<const char*>(&m_numcarbon), sizeof(m_numcarbon));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleCache::Serialize).");
    }
}

// Reads the object from a binary stream.
void ParticleCache::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    releaseMem();
    //m_pmodel = &model;

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double       val = 0.0;

        switch (version) {
            case 0:
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

                // Read the collision diameter squared.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_dcolsqr = (real)val;

                // Read the inverse collision diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_inv_dcol = (real)val;

                // Read the inverse collision diameter squared.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_inv_dcolsqr = (real)val;

                // Read the inverse square root of mass.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_inv_sqrtmass = (real)val;

                // Read the D^2 / sqrt(M).
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_d2_m_1_2 = (real)val;

                // Read the free surface.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_freesurface = (real)val;

                // Read the number of sub paricles.
                in.read(reinterpret_cast<char*>(&m_numsubpart), sizeof(m_numsubpart));

                // Read number of C atoms.
                in.read(reinterpret_cast<char*>(&m_numcarbon), sizeof(m_numcarbon));

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ParticleCache::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ParticleCache::Deserialize).");
    }
}


// MEMORY MANAGEMENT.

// Release all memory associated with the ParticleData object.
void ParticleCache::releaseMem(void)
{
    //nothing to do
}


// Initialisation routine.
void ParticleCache::init(void)
{

	m_diam       = 0.0;
	m_dcol       = 0.0;
	m_dmob       = 0.0;
	m_surf       = 0.0;
	m_vol        = 0.0;
	m_mass       = 0.0;
    m_dcolsqr      = 0.0;
    m_inv_dcol     = 0.0;
    m_inv_dcolsqr  = 0.0;
    m_inv_sqrtmass = 0.0;
    m_d2_m_1_2     = 0.0;


    m_freesurface = 0.0;
    m_numsubpart = 0;

    m_numcarbon = 0;
}
