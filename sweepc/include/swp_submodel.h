/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The SubModel class defines the interface for a particle sub-model (layer).
*/

#ifndef SWEEP_SUBMODEL_H
#define SWEEP_SUBMODEL_H

#include "swp_params.h"
#include "swp_submodel_type.h"
#include <vector>
#include <iostream>

namespace Sweep
{
// Forward declare Primary, ParticleCache classes.
class Primary;
class ParticleCache;

namespace SubModels
{
// Forward declare SubModelCache class.
class SubModelCache;

class SubModel
{
public:
    // Constructors.
    //   Default constructor is protected to prevent a sub-model
    //   being created without knowledge of the parent primary.
    SubModel(Primary &parent);       // Initialising constructor.
    SubModel(const SubModel &copy);  // Copy constructor.
    SubModel(                        // Stream-reading constructor.
        const std::istream &in,      //  - Input stream.
        Primary &parent              //  - Parent primary particle.
        );

    // Destructor.
    virtual ~SubModel(void);

    // Operators.
    virtual SubModel &operator=(const SubModel &rhs);
//    virtual SubModel &operator+=(const SubModel &rhs)=0;
//    virtual const SubModel operator+(const SubModel &rhs) const=0;


    // PARENT PRIMARY PARTICLE.

    // Returns the parent primary particle.
    const Primary *const Parent(void) const;

    // Sets the parent primary particle.
    void SetParent(Primary &parent);


    // PARTICLE UPDATES.

    // Updates the model data according to the rules of the model
    // given the changes to particle composition and values.
    // Performs the update n times.
    virtual void Update(
        const fvector &dcomp, // Changes to the composition.
        const fvector &dval,  // Changes to the tracker variables.
        unsigned int n=1      // Number of times to perform update.
        ) = 0;

    // Updates the data according to the model rules for
    // coagulation.
    virtual SubModel &Coagulate(const SubModel &rhs) = 0;


    // CACHE UPDATES.

    // Updates the sub-model cached variables.
    virtual void UpdateCache(void) = 0;

    // Creates a model cache which subscribes to this particle model.
    virtual SubModelCache *const CreateCache(ParticleCache &pcache) const = 0;


    // READ/WRITE/COPY.

    // Returns a copy of the model data.
    virtual SubModel *const Clone(void) const = 0;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    virtual SubModelType ID(void) const = 0;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const = 0;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in, // Input stream.
        Primary &parent   // Parent object.
        ) = 0;

protected:
    // Parent primary particle.
    Primary *m_parent;

    // Sub-models cannot be created without knowledge of the parent
    // primary particle.
    SubModel(void);
};

typedef std::vector<SubModel*> SubModelPtrVector;
typedef std::map<SubModelType, SubModel*> SubModelMap;
};
};

#endif
