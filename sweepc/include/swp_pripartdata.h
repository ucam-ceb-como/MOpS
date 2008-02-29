/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ModelData class defines the additional data which is
    added to a particle to enable a given model.
*/

#ifndef SWEEP_PRIPARTDATA_H
#define SWEEP_PRIPARTDATA_H

#include "swp_primary.h"
#include "swp_modeltype.h"
#include "swp_modeldata.h"
#include "swp_pripartmodel.h"
#include <vector>

namespace Sweep
{
class PriPartData : public IModelData
{
friend class PriPartModel;

public:
    // Constructors.
    PriPartData(ParticleData &parent);  // Default constructor.
    PriPartData(const PriPartData &copy); // Copy constructor.
    PriPartData(             // Stream-reading constructor.
        std::istream &in,    //  - Input stream.
        ParticleData &parent //  - Parent ParticleData object.
        );

    // Destructors.
    ~PriPartData(void);

    // Operators.
    PriPartData &operator=(const PriPartData &rhs);
    PriPartData &operator+=(const PriPartData &rhs);
    const PriPartData operator+(const PriPartData &rhs) const;

    // Resets the model data to the default state.
    void Clear();


    // PROPERTIES.

    // Returns the vector of primary particles.
    std::vector<Primary> &Primaries(void);
    const std::vector<Primary> &Primaries(void) const;

    // Returns the property with the given ID.
    real Property(unsigned int id) const;


    // MODEL WHICH USES THIS DATA.

    // Returns the PriPartModel object which operator on this data.
    const PriPartModel &Model(void) const;


    // READ/WRITE/COPY.

    // Returns a copy of the model data.
    PriPartData *const Clone(void) const;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    ModelType ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

protected:
    // Can't create a ModelData object independently of a
    // parent ParticleData.
    PriPartData(void);

private:
    std::vector<Primary> m_primaries;
};

};

#endif
