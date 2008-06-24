/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The AggModelCache class is a base model cache for different primary
    particle aggregation models.
*/

#ifndef SWEEP_AGGMODEL_CACHE_H
#define SWEEP_AGGMODEL_CACHE_H

#include "swp_params.h"
#include "swp_aggmodel_type.h"
#include <iostream>

namespace Sweep
{
// Forward declare parent particle cache class.
class ParticleCache;
class Primary;

namespace AggModels
{
class AggModelCache
{
public:
    // Destructors.
    virtual ~AggModelCache(void);

    // Operators.
    virtual AggModelCache& operator=(const AggModelCache &rhs);
    virtual AggModelCache& operator=(const Primary &rhs) = 0;
    virtual AggModelCache& operator+=(const AggModelCache &rhs);
    virtual AggModelCache& operator+=(const Primary &rhs) = 0;

    // DATA MANAGEMENT.

    // Resets the cache to its "empty" condition.
    virtual void Clear(void) = 0;


    // PARENT OBJECT.

    // Returns a pointer to the parent particle data.
    ParticleCache *const Parent(void) const;

    // Sets the parent particle data.
    void SetParent(ParticleCache &parent);


    // READ/WRITE/COPY.

    // Returns a correctly typed reference to the model cache object.
//    virtual AggModelCache &TypedRef(void) = 0;
//    virtual const AggModelCache &TypedRef(void) const = 0;

    // Creates a copy of the particle data object.
    virtual AggModelCache *const Clone(void) const = 0;

    // Returns the aggregation model type ID.  Required for serialization.
    virtual AggModelType ID(void) const = 0;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const = 0;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,     // Input stream.
        ParticleCache &parent // Parent object.
        ) = 0;

//    virtual AggModelCache &CopyTo(AggModelCache &lhs) const = 0;

protected:
    // Constructors (protected because this class should not be created,
    //               instead create derived classes).
    //    Note default constructor is protected to prevent an AggModelCache
    //    object being created without knowledge of parent object.
    AggModelCache(void);                      // Default constructor.
    AggModelCache(ParticleCache &parent);     // Initialising constructor.
    AggModelCache(const AggModelCache &copy); // Copy constructor.
    AggModelCache(            // Stream-reading constructor.
        std::istream &in,     //  - Input stream.
        ParticleCache &parent //  - Parent object.
        );

    // Parent object.
    ParticleCache *m_parent;
};
};
};

#endif
