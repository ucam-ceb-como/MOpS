/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The SubParticle class is integral to the sub-particle tree paradigm.  It is
    the tree node in a binary tree structure in which each sub-particle consists
    of either two further sub-particles or a primary particle.  In this way
    advanced particle structure and processes (such as sintering) may be modelled.

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

#ifndef SWEEP_SUBPARTICLE_H
#define SWEEP_SUBPARTICLE_H

#include "swp_params.h"
#include "swp_sintering_model.h"
#include "swp_particle_model.h"
#include "swp_particle_cache.h"
#include "swp_primary.h"
#include <iostream>

namespace Sweep
{
class Cell;

class SubParticle : public ParticleCache
{
public:
    // Constructors.
    SubParticle(                          // Initialising constructor.
        real time,                        //   - Create time.
        const Sweep::ParticleModel &model //   - Defining particle model.
        );
    SubParticle(Sweep::Primary &pri);     // Initialising constructor (from primary particle).
    SubParticle(const SubParticle &copy); // Copy constructor.
    SubParticle(                          // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Model to which this particle subscribes.
        );

    // Destructors.
    virtual ~SubParticle(void);

    // Operators.
    SubParticle &operator=(const SubParticle &rhs);
    SubParticle &operator=(const Sweep::Primary &rhs);
    SubParticle &operator+=(const SubParticle &rhs);
    SubParticle &operator+=(const Sweep::Primary &rhs);
    const SubParticle operator+(const SubParticle &rhs) const;
    const SubParticle operator+(const Sweep::Primary &rhs) const;


    // PARENT SUB-PARTICLE.

    // Returns the sub-particle's parent.
    SubParticle *const Parent(void);
    const SubParticle *const Parent(void) const;


    // PRIMARY PARTICLE CHILD.

    // Returns a pointer to the child primary particle, NULL if 
    // this sub-particle has no primary.
    Sweep::Primary *const Primary(void);
    const Sweep::Primary *const Primary(void) const;


    // CHILD SUB-PARTICLES IN SUB-PARTICLE TREE.

    // Returns a pointer to the left child.
    SubParticle *const Left(void);
    const SubParticle *const Left(void) const;

    // Returns a pointer to the right child.
    SubParticle *const Right(void);
    const SubParticle *const Right(void) const;

    // Selects a leaf from the sub-particle tree below this
    // sub-particle using the given property ID to weight the
    // primaries.
    Sweep::SubParticle *const SelectLeaf(
        SubModels::SubModelType model_id, // Sub-model ID for weighting parameter.
        int id                            // Weighting property ID.
        );


    // PROPERTY SETTING OVERRIDES (FROM PARTICLE CACHE).

    // Sets the composition vector.
    void SetComposition(const fvector &comp);

    // Sets the values vector.
    void SetValues(const fvector &vals);

    // Sets the last update time of the particle.
    void SetTime(real t);

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


    // PARTICLE OPERATIONS.

    // Adjusts the particle with the given composition and 
    // values changes n times.  If the particle cannot be adjust
    // n times, then this function returns the number of times
    // it was adjusted. The function proceeds down the sub-particle tree
    // until it reaches a primary particle leaf.  At this point
    // that primary particle is adjusted using the given values.
    unsigned int Adjust(
        const fvector &dcomp,             // Composition changes.
        const fvector &dvalues,           // Tracker variable changes.
        SubModels::SubModelType model_id, // Sub-model ID for weighting parameter.
        int id,                           // Weighting property ID.
        unsigned int n                    // Number of times to perform adjustment.
        );

    // Combines this particle with another.
    virtual SubParticle &Coagulate(const SubParticle &sp);

    // Sinters the sub-particle for the given time using the given
    // sintering model.
    virtual void Sinter(
        real dt,         // Delta-t for sintering.
        const Cell &sys, // System which defines particle's environment.
        const Processes::SinteringModel &model // Sintering model to use.
        );


    // PARTICLE UPDATE AND CHECKING.

    // Recalculates the derived properties from the 
    // unique properties.  This function moves down the tree
    // from the root to the leaves.
    void UpdateCache(void);

    // Tells the parent sub-particle to update its cache of
    // derived properties.  This operation is passed up the
    // sub-particle tree to the root particle.
    void UpdateTree(void);

    // Check the that the particle is valid by querying the
    // validity conditions of the models and ensuring that it 
    // contains any components.
    bool IsValid() const;


    // READ/WRITE/COPY.

    // Creates a copy of the particle data object.
    virtual SubParticle *const Clone(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

protected:
    // Parent sub-particle.
    SubParticle *m_parent;

    // Left and right children in the sub-particle tree structure.
    SubParticle *m_leftchild, *m_rightchild;

    // Primary particle of which this sub-particle comprises.  A sub-particle
    // may contain one primary, or exactly two sub-particles (see above).
    Sweep::Primary *m_primary;

    // SubParticle class cannot be created without knowledge of the
    // components and the values.
    SubParticle(void);


    // MEMORY MANAGEMENT.

    // Release all memory associated with the SubParticle object.
    void releaseMem(void);

    // Initialisation routine.
    void init(void);


    // PARENT SUB-PARTICLE.

    // Sets the sub-particle's parent.
    void setParent(SubParticle *const p);


    // PRIMARY PARTICLE CHILD.

    // Sets the pointer to the child primary particle.
    void setPrimaryPtr(Sweep::Primary *const pri);

    // Selects a leaf from the sub-particle tree below this
    // sub-particle using the given property ID to weight the
    // primaries.  This particular function is used internally by the SubParticle class
    // to reuse the random variable used to select the correct leaf which
    // is generated by the root sub-particle.  
    Sweep::SubParticle *const selLeafLoop(
        SubModels::SubModelType model_id, // Sub-model ID for weighting parameter.
        int id,                           // Weighting property ID.
        real r                            // Random leaf-selection number.
       );


    // CHILD SUB-PARTICLES IN SUB-PARTICLE TREE.

    // Sets the pointer to the left child.
    void setLeftPtr(SubParticle *const sub);

    // Sets the pointer to the right child.
    void setRightPtr(SubParticle *const sub);


    // OPERATIONS.

    // Adjusts the particle with the given composition and 
    // values changes n times.  If the particle cannot be adjust
    // n times, then this function returns the number of times
    // it was adjusted. The function proceeds down the sub-particle tree
    // until it reaches a primary particle leaf.  At this point
    // that primary particle is adjusted using the given values.  This
    // particular function is used internally by the SubParticle class
    // to reuse the random variable used to select the correct leaf which
    // is generated by the root sub-particle.
    unsigned int adjustLoop(
        const fvector &dcomp,             // Composition changes.
        const fvector &dvalues,           // Tracker variable changes.
        SubModels::SubModelType model_id, // Sub-model ID for weighting parameter.
        int id,                           // Weighting property ID.
        real r,                           // Random leaf-selection number.
        unsigned int n                    // Number of times to perform adjustment.
        );


};
};

#endif
