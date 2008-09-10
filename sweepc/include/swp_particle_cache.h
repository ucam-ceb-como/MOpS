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
#include "swp_particle_model.h"
#include "swp_primary.h"
#include "swp_submodel_type.h"
#include "swp_submodel_cache.h"
#include "swp_aggmodel_cache.h"
#include <iostream>

namespace Sweep
{
// Forward declare the Ensemble and TreeNode classes as containers.
class Ensemble;
class TreeNode;

class ParticleCache
{
friend class Ensemble;
friend class TreeNode;

public:
    // Enumeration of ParticleCache properties which can be accessed using
    // the Property() function.
    enum PropID {
        iUniform=-1, // Special Case:  Always returns 1.0.  Used to select particles uniformly.
        iCTime,  // Create time.
        iLUTime, // Last update time.
        iD,      // Equivalent sphere diameter.
        iDcol,   // Collision diameter.
        iDmob,   // Mobility diameter.
        iS,      // Surface area.
        iV,      // Volume.
        iM,      // Mass.
        // The next properties are provided for calculation of
        // collision rates.
        iD2,      // Collision diameter squared.
        iD_1,     // Inverse collision diameter.
        iD_2,     // Inverse of the diameter squared.
        iM_1_2,   // Inverse of the square-root of the mass.
        iD2_M_1_2, // D^2 * M^-1/2.


    };

    // Constructors.
    //   The default constructor is protected to prevent sub-particles
    //   being created without knowledge of the defining particle model.
    ParticleCache( // Initialising constructor.
        real time,                        // Create time.
        const Sweep::ParticleModel &model // Defining particle model.
        );
    ParticleCache(const Primary &pri);        // Initialising constructor (from primary particle).
    ParticleCache(const ParticleCache &copy); // Copy constructor.
    ParticleCache(                            // Stream-reading constructor.
        std::istream &in,                     //  - Input stream.
        const Sweep::ParticleModel &model     //  - Model to which this cache subscribes.
        );

    // Destructors.
    virtual ~ParticleCache(void);

    // Operators.
    ParticleCache &operator=(const ParticleCache &rhs);
    ParticleCache &operator=(const Primary &rhs);
    ParticleCache &operator+=(const ParticleCache &rhs);
    ParticleCache &operator+=(const Primary &rhs);
    const ParticleCache operator+(const ParticleCache &rhs) const;
    const ParticleCache operator+(const Primary &rhs) const;

    // Resets the particle cache to its "empty" condition.
    void Clear(void);


    // DEFINING PARTICLE MODEL.
    
    // Returns the particle model used to define this particle cache.
    const Sweep::ParticleModel *const ParticleModel(void) const;


    // COMPOSITION.

    // Returns the composition vector.
    const fvector &Composition(void) const;

    // Returns the ith component value.  Returns 0.0 if i invalid.
    real Composition(unsigned int i) const;

    // Sets the composition vector.
    void SetComposition(const fvector &comp);


    // TRACKER VALUES.

    // Returns the values vector.
    const fvector &Values(void) const;

    // Returns the ith value.  Returns 0.0 if i invalid.
    real Values(unsigned int i) const;

    // Sets the values vector.
    void SetValues(const fvector &vals);


    // PARTICLE CREATE TIME.

    // Returns the particle create time.
    real CreateTime(void) const;


    // LAST UPDATE TIME.

    // Returns the last update time of the particle.
    real LastUpdateTime(void) const;

    // Sets the last update time of the particle.
    void SetTime(real t);


    // AGGREGATION-MODEL CACHE.

    // Returns the aggregation model data.
    const AggModels::AggModelCache *const AggCache(void) const;


    // SUB-MODEL CACHE.

    // Returns the model data.
    const SubModels::SubModelCacheMap &SubModelCache(void) const;

    // Returns the data for the idth model.  Returns NULL if id is invalid.
    const SubModels::SubModelCache *const SubModel(SubModels::SubModelType id) const;


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

    // Creates a copy of the particle data object.
    virtual ParticleCache *const Clone(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

protected:
    // Defining particle model.
    const Sweep::ParticleModel *m_pmodel;

    // Unique properties.
    fvector m_comp;   // Composition.
    fvector m_values; // Other particle values (defined at run time).
    real m_createt;   // Time at which particle was created (earliest part).
    real m_time;      // Last time particle was updated.  Required for LPDA.



    // Additional particle models.
    AggModels::AggModelCache *m_aggcache;    // The aggregation-model used by the particle.
    SubModels::SubModelCacheMap m_submodels; // Cache of data required by different sub-models.

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

    // ParticleCache class cannot be created without knowledge of the
    // defining particle model.
    ParticleCache(void);

    // Release all memory associated with the object.
    void releaseMem(void);

    // Initialisation routine.
    void init(void);
};
};

#endif
