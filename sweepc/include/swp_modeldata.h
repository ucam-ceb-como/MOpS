/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ModelData class defines the additional data which is
    added to a particle to enable a given model.
*/

#ifndef SWEEP_MODELDATA_H
#define SWEEP_MODELDATA_H

#include "swp_params.h"
#include "swp_modeltype.h"
#include "swp_model.h"
#include <vector>
#include <map>
#include <iostream>

namespace Sweep
{
class ParticleData; // Forward declare parent object.

class IModelData
{
public:
    // Constructors.
    IModelData(ParticleData &parent);   // Default constructor.
    IModelData(const IModelData &copy); // Copy constructor.

    // Destructors.
    virtual ~IModelData(void);

    // Operators.
    virtual IModelData &operator=(const IModelData &rhs);
    virtual IModelData &operator+=(const IModelData &rhs);
    virtual const IModelData operator+(const IModelData &rhs) const;


    // PARENT ParticleData OBJECT.

    // Returns a pointer to the parent particle data.
    ParticleData *const Parent(void) const;

    // Sets the parent particle data.
    void SetParent(ParticleData &parent);


    // PROPERTIES.

    // Returns the property with the given ID.
    virtual real Property(unsigned int id) const = 0;


    // MODEL WHICH USES THIS DATA.

    virtual const IModel &Model(void) const = 0;


    // READ/WRITE/COPY.

    // Returns a copy of the model data.
    virtual IModelData *const Clone(void) const = 0;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    virtual ModelType ID(void) const = 0;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const = 0;

    // Reads the object from a binary stream.
    virtual void Deserialize(std::istream &in) = 0;

protected:
    // Can't create a ModelData object independently of a
    // parent ParticleData.
    IModelData(void);

private:
    ParticleData *m_parent;
};

typedef std::vector<IModelData*> ModelDataPtrVector;
typedef std::map<ModelType,IModelData*> ModelMap;
};

#endif
