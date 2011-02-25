/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ParticleCache class defines the data structure of a sub-particle.  It contains
    the particle composition, basic properties and the additional properties
    required for the various particle models.  This class is the interface between
    the particle class and the ensemble, as both those classes contain data objects
    of type ParticleCache.

    The particle cache contains only a sub-set of information about the sub-particle
    structure.  It is designed as a data store only which is applicable to the
    numerics of the ensemble binary tree (in terms of property sums), but the data
    stored within the ParticleCache class may not be physical.  For example the diameter
    could be the sum of two sub-particle diameters, rather than the true aggregate
    collision diameter.  The SubParticle class inherits the ParticleCache and imposes
    strict physical data requirements, that is it represents the true particle
    structure, rather than a numerical construct.

    The is necessarily some overlap between the ParticleCache class and the Primary
    class.  Both contain basic data such as mass, volume and diameter.  However, there
    are two distinctions:  first, the ParticleCache contains only a subset of the
    Primary sub-model data.  This prevents complex models with large memory requirements
    being stored at every branch in the sub-particle and ensemble trees.  Secondly the
    ParticleCache contains an aggregation model cache, which holds cached data from the
    primary particle model, for example the surface-volume model.

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

#ifndef SWEEP_PARTICLE_CACHE_H
#define SWEEP_PARTICLE_CACHE_H

#include "swp_params.h"
#include "swp_property_indices.h"

#include <iostream>

namespace Sweep
{
// forward declarations
    class Primary;
    class ParticleModel;

/*!
 * \brief   Summary properties that can be summed and used to calculate overall population rates
 */
class ParticleCache
{
public:
    //! Create cache of zeros
    ParticleCache();

    //! Initialise from primary particle
    ParticleCache(const Primary &pri);

    ParticleCache(                            // Stream-reading constructor.
        std::istream &in,                     //  - Input stream.
        const Sweep::ParticleModel &model     //  - Model to which this cache subscribes.
        );

    // Operators.
    ParticleCache &operator=(const Primary &rhs);

    // Resets the particle cache to its "empty" condition.
    void Clear(void);


    // BASIC PROPERTIES.

    // Returns the particle equivalent sphere diameter.
    real SphDiameter(void) const;

    // Returns the collision diameter.
    real CollDiameter(void) const;

    // Returns the mobility diameter.
    real MobDiameter(void) const;

    // Returns the surface area.
    real SurfaceArea(void) const;

    // Returns the equivalent sphere surface area, based
    // on the volume.
    real SphSurfaceArea(void) const;

    // Returns the volume.
    real Volume(void) const;

    // Returns the mass.
    real Mass(void) const;

    // Returns the property with the given ID.
    real Property(PropID id) const;

	real Sintered();

    // BASIC DERIVED PROPERTY OVERWRITES.

    // Sets the spherical particle diameter
    void SetSphDiameter(real diam);

    // Sets the collision diameter of the particle.
    void SetCollDiameter(real dcol);

    // Sets the mobility diameter.
    void SetMobDiameter(real dmob);

    // Sets the surface area, subject to minimum spherical area condition.
    void SetSurfaceArea(real surf);

    // Sets the volume.
    void SetVolume(real vol);

    // Sets the mass.
    void SetMass(real m);



    // COLLISION RATE CALCULATION PARTICLE PROPERTIES.

    // Collision diameter squared (cm2).
    real CollDiamSquared() const;

    // Inverse collision diameter (cm-1).
    real InvCollDiam() const;

    // Inverse squared collision diameter (cm-2).
    real InvCollDiamSquared() const;

    // Inverse of square root of mass (g-1/2).
    real InvSqrtMass() const;

    // Collision diameter squared times the inverse square root of mass.
    real CollDiamSqrdInvSqrtMass() const;


    // READ/WRITE/COPY.

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

protected:
    // Operators.
    ParticleCache &operator+=(const ParticleCache &rhs);
    ParticleCache &operator+=(const Primary &rhs);
    const ParticleCache operator+(const ParticleCache &rhs) const;
    const ParticleCache operator+(const Primary &rhs) const;

    // Basic properties (calculated from above properties).
    real m_diam; // Equivalent spherical diameter.
    real m_dcol; // Collision diameter.
    real m_dmob; // Mobility diameter.
    real m_surf; // Surface area.
    real m_vol;  // Volume.
    real m_mass; // Mass.

    // Collision rate calculation particle properties.
    real m_dcolsqr;      // Collision diameter squared.
    real m_inv_dcol;     // Inverse collision diameter.
    real m_inv_dcolsqr;  // Inverse of the diameter squared.
    real m_inv_sqrtmass; // Inverse of the square-root of the mass.
    real m_d2_m_1_2;     // D^2 * M^-1/2.



    // The free surface available for other particles to sinter
    real m_freesurface;
    //number of subparticles below this node (no longer used, should be removed)
    unsigned int m_numsubpart;

    //! Number of Carbon atoms
    unsigned int m_numcarbon;

    // Release all memory associated with the object.
    void releaseMem(void);

    // Initialisation routine.
    void init(void);
};
}

#endif
