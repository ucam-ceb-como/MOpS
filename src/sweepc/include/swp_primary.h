/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Primary class defines the smallest sub-particle unit.  It uses 
    the spherical particle model, though is designed so that derived
    classes can relax that assumption.  It stores the fundamental
    data structures required to describe different particle models.

    In the context of the sub-particle tree, primaries can be thought of
    as the leaves of the tree.  The SubParticle class describes the branches
    and the Particle class describes the trunk.

    A Primary particle subscribes to a particular ParticleModel object.  The
    particle model defines what constitutes a primary and what sub-models are
    enabled.  A Primary is constructed with knowledge of the particle model,
    and the particle model cannot be changed.  It is expected that the
    data in the ParticleModel object does not change once Primaries have
    been created (e.g. by adding/removing components), and this class is
    defined uner this assumption.

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

#ifndef SWEEP_PRIMARY_H
#define SWEEP_PRIMARY_H

#include "swp_params.h"
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_particle_model.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_sintering_model.h"
#include <iostream>

namespace Sweep
{
// Forward declaration
class Cell;

class Primary
{
public:
    // Enumeration of properties which can be accessed using
    // the Property() function.  These are fundamental properties
    // which are not altered by sub-models.
    enum PropID {
        iCTime,  // Create time.
        iLUTime, // Last update time.
        iD,      // Equivalent sphere diameter.
        iDcol,   // Collision diameter.
        iDmob,   // Mobility diameter.
        iS,      // Surface area.
        iV,      // Volume.
        iM       // Mass.
    };

    // Constructors.
    //   Note:  Default constructor is protected to prevent a
    //          Primary being created without knowledge of the
    //          defining particle model.
    Primary(                              // Initialising constructor.
        real time,                        //  - Create time.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );
    Primary(const Primary &copy); // Copy constructor.
    Primary(                              // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );

    // Destructors.
    virtual ~Primary(void);

    // Operators.
    virtual Primary &operator=(const Primary &rhs);

    // DEFINING PARTICLE MODEL.
    
    // Returns the particle model used to create this primary.
    const Sweep::ParticleModel *const ParticleModel(void) const;


    // PRIMARY COMPOSITION.

    // Returns the composition vector.
    const fvector &Composition(void) const;

    // Returns the ith component value.  Returns 0.0 if i invalid.
    real Composition(unsigned int i) const;

    // Sets the composition vector.
    void SetComposition(const fvector &comp);


    // PRIMARY TRACKER VALUES.

    // Returns the values vector.
    const fvector &Values(void) const;

    // Returns the ith value.  Returns 0.0 if i invalid.
    real Values(unsigned int i) const;

    // Sets the values vector.
    void SetValues(const fvector &vals);

    // Sets the ith trackervalue.
    void SetValue(unsigned int i, real val);


    // PRIMARY CREATE TIME.

    // Returns the particle create time.
    real CreateTime(void) const;


    // LAST UPDATE TIME.
    //   This is provided to assist numerics only.  It is not available
    //   to other programs/libraries outside sweep.

    // Returns the last update time of the particle.
    real LastUpdateTime(void) const;

    // Sets the last update time of the particle.
    void SetTime(real t);


    // AGGREGATION MODEL.

    // Returns the aggregation model which this primary describes.
    virtual AggModels::AggModelType AggID(void) const;

    // Creates an aggregation data cache for this primary type.
    virtual AggModels::AggModelCache *const CreateAggCache() const;


    // BASIC DERIVED PROPERTIES.

    // Calculates the derived properties from the unique properties.
    virtual void UpdateCache(void);

    // Returns the equivalent-sphere diameter.
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

    //! Check particle still meets physical conditions for being a particle.
    bool IsValid() const;

    // OPERATIONS.

    // Adjusts the primary with the given composition and 
    // tracker values changes n times.  If the particle cannot be adjust
    // n times, then this function returns the number of times
    // it was adjusted.
    virtual unsigned int Adjust(
        const fvector &dcomp,   // Composition changes.
        const fvector &dvalues, // Tracker variable changes.
        unsigned int n=1        // Number of times to perform adjustment.
        );

    // Combines this primary with another.  This is also the
    // implementation of the + and += operators.
    virtual Primary &Coagulate(const Primary &rhs,
                               int (*rand_int)(int, int),
                               real(*rand_u01)());

    // This routine sinters the Primary for the given length of
    // time using the provided sintering model.
    virtual void Sinter(
        real dt, // Delta-t for sintering to occur.
        const Cell &sys, // System which defines primary's environment.
        const Processes::SinteringModel &model, // Sintering model to use.
        real (*rand_u01)() // Uniform [0,1] sample generator
        );

    // READ/WRITE/COPY.

    // Returns a copy of the primary.
    virtual Primary *const Clone(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

protected:
    // Particle model used to define the Primary.
    const Sweep::ParticleModel *m_pmodel;

    // Unique properties.
    fvector m_comp;   // Primary composition.
    fvector m_values; // Other primary values (defined at run time).
    real m_createt;   // Time at which primary was created.
    real m_time;      // Last time primary was updated.  Required for LPDA.

    // Basic derived properties (calculated from above properties).
    real m_diam; // Equivalent spherical diameter.
    real m_dcol; // Collision diameter.
    real m_dmob; // Mobility diameter.
    real m_surf; // Surface area.
    real m_vol;  // Volume.
    real m_mass; // Mass.

    // Primary class cannot be created without knowledge of the
    // particle model, therefore default constructor is protected.
    Primary(void);


    // DERIVED PROPERTIES.

    // Sets the particle cache to that of a spherical particle.
    void setSphereCache(void);


    // MEMORY MANAGEMENT.

    // Release all memory associated with object.
    void releaseMem(void);

    // Initialisation routine.
    void init(void);


};
};

#endif
