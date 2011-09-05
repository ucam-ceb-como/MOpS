/*
  Author(s):      Matthew Celnik (msc37) and Markus Sander (ms785)
  Project:        sweep (population balance solver)

  File purpose:
    The SubParticle class is integral to the sub-particle tree paradigm.  It is
    the tree node in a binary tree structure in which each sub-particle consists
    of either two further sub-particles or a primary particle.  In this way
    advanced particle structure and processes (such as sintering) may be modelled.
*/

#ifndef SWEEP_SUBPARTICLE_H
#define SWEEP_SUBPARTICLE_H

#include "swp_params.h"
#include "swp_sintering_model.h"
#include "swp_particle_model.h"
#include "swp_primary.h"
#include "swp_property_indices.h"

#include <iostream>


namespace Sweep
{
class Cell;

/*!
 * \brief Class that implements some of the Particle interface.
 *
 * This class is now obsolete.  It was originally used to store the
 * connectivity structure of primary particle aggregates, but this
 * connectivity is now handled in the primary particle classes.
 *
 * This class will eventually removed so that Particle inherits
 * directly from ParticleCache, without having SubParticle in the
 * inheritance hierachy any more.
 *
 * DO NOT add methods or data to this class, since the whole
 * class is due to be removed. (riap2 28 jan 2011)
 *
 */
class SubParticle
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

    // PRIMARY PARTICLE CHILD.

    // Returns a pointer to the child primary particle, NULL if
    // this sub-particle has no primary.
    Sweep::Primary *const Primary(void);
    const Sweep::Primary *const Primary(void) const;

    // BASIC PROPERTIES.

    //! Returns the particle equivalent sphere diameter.
    real SphDiameter(void) const;

    //! Returns the collision diameter.
    real CollDiameter(void) const;

    //! Returns the mobility diameter.
    real MobDiameter(void) const;

    //! Returns the surface area.
    real SurfaceArea(void) const;

    //! Returns the equivalent sphere surface area, based on the volume.
    real SphSurfaceArea(void) const;

    //! Returns the volume.
    real Volume(void) const;

    //! Returns the mass.
    real Mass(void) const;

    //! Returns the property with the given ID.
    real Property(Sweep::PropID id) const;



    // PROPERTY SETTING OVERRIDES (FROM PARTICLE CACHE).

    // Sets the last update time of the particle.
    void SetTime(real t);

	// Sets the sintering level of the children
    void SetSintering(real);

	// Returns the sintering level of the childen
    real GetSintering();

    // COMPOSITION.

    //! Returns the composition vector.
    const fvector &Composition(void) const;

    //! Returns the ith component value.  Returns 0.0 if i invalid.
    real Composition(unsigned int i) const;


    // TRACKER VALUES.

    //! Returns the tracked values vector.
    const fvector &Values(void) const;

    //! Returns the ith tracked value.  Returns 0.0 if i invalid.
    real Values(unsigned int i) const;

    // PARTICLE CREATE TIME.

    // Returns the particle create time.
    real CreateTime(void) const;


    // LAST UPDATE TIME.

    // Returns the last update time of the particle.
    real LastUpdateTime(void) const;

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
        unsigned int n                    // Number of times to perform adjustment.
        );

    // Combines this particle with another.
    virtual SubParticle &Coagulate(const SubParticle &sp,
                                   int (*rand_int)(int, int),
                                   real(*rand_u01)());


    // Sinters the sub-particle for the given time using the given
    // sintering model.
    virtual void Sinter(
        real dt,         // Delta-t for sintering.
        Cell &sys, // System which defines particle's environment.
        const Processes::SinteringModel &model, // Sintering model to use.
        real (*rand_u01)() // Generate U[0,1] samples
        );

    // PARTICLE UPDATE AND CHECKING.

    // Recalculates the derived properties from the
    // unique properties.  This function moves down the tree
    // from the root to the leaves.
    void UpdateCache(void);

    // Check the that the particle is valid by querying the
    // validity conditions of the models and ensuring that it
    // contains any components.
    bool IsValid() const;

    //! Aggregation model data
    const AggModels::AggModelCache& AggCache(void) const {return *m_aggcache;}


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

	real avgeomdiam(double oneovernumsubpart);

protected:

    // Primary particle of which this sub-particle comprises (since sub-particles
    // must contain exactly one primary, it might be possible to make this a member)..
    Sweep::Primary *m_primary;

    // SubParticle class cannot be created without knowledge of the
    // components and the values.
    SubParticle(void);

    // MEMORY MANAGEMENT.

    // Release all memory associated with the SubParticle object.
    void releaseMem(void);

    // Initialisation routine.
    void init(void);

    // PRIMARY PARTICLE CHILD.

    // Sets the pointer to the child primary particle.
    void setPrimaryPtr(Sweep::Primary *const pri);

private:
    //! Time at which particle was created (earliest part).
    real m_createt;

    //! Last time particle was updated.  Required for LPDA.
    real m_time;

    //! The aggregation-model used by the particle.
    AggModels::AggModelCache *m_aggcache;
};
};

#endif
