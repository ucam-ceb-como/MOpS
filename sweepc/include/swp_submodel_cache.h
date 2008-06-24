/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ModelData class defines the additional data which is
    added to a particle to enable a given model.
*/

#ifndef SWEEP_SUBMODEL_CACHE_H
#define SWEEP_SUBMODEL_CACHE_H

#include "swp_params.h"
#include "swp_submodel_type.h"
#include "swp_submodel.h"
#include <vector>
#include <map>
#include <iostream>

namespace Sweep
{
// Forward declare parent object.
class ParticleCache;

namespace SubModels
{
// Declare base class for storing sub-model data caches.
class SubModelCache
{
public:
    // Constructors.
    SubModelCache(ParticleCache &parent);     // Default constructor.
    SubModelCache(const SubModelCache &copy); // Copy constructor.

    // Destructors.
    virtual ~SubModelCache(void);

    // Operators.
    virtual SubModelCache &operator=(const SubModelCache &rhs) = 0;
    virtual SubModelCache &operator=(const SubModel &rhs) = 0;
    virtual SubModelCache &operator+=(const SubModelCache &rhs) = 0;
    virtual SubModelCache &operator+=(const SubModel &rhs) = 0;
//    virtual const SubModelCache operator+(const SubModelCache &rhs) const = 0;
//    virtual const SubModelCache operator+(const SubModel &rhs) const = 0;

    // Resets the model data to the default state.
    virtual void Clear() = 0;


    // PARENT OBJECT.

    // Returns a pointer to the parent particle data.
    ParticleCache *const Parent(void) const;

    // Sets the parent particle data.
    void SetParent(ParticleCache &parent);


    // PROPERTIES.

    // Returns the property with the given ID.
    virtual real Property(unsigned int id) const = 0;


    // READ/WRITE/COPY.

    // Returns a copy of the model data.
    virtual SubModelCache *const Clone(void) const = 0;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    virtual SubModelType ID(void) const = 0;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const = 0;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,     // Input stream.
        ParticleCache &parent // Parent object.
        ) = 0;

protected:
    // Can't create a SubModelCache object independently of a
    // parent ParticleCache.
    SubModelCache(void);

private:
    // Parent object.
    ParticleCache *m_parent;
};

typedef std::vector<SubModelCache*> SubModelCachePtrVector;
typedef std::map<SubModelType, SubModelCache*> SubModelCacheMap;
};
};

#endif
