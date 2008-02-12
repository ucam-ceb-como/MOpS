/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ModelData class defines the additional data which is
    added to a particle to enable a given model.
*/

#ifndef SWEEP_PRIPARTMODELDATA_H
#define SWEEP_PRIPARTMODELDATA_H

#include "swp_primary.h"
#include "swp_modeltype.h"
#include "swp_modeldata.h"
#include "swp_pripartmodel.h"
#include <vector>

namespace Sweep
{
class PriPartModelData : public IModelData
{
friend PriPartModel;

public:
    // Constructors.
    PriPartModelData(ParticleData &parent);  // Default constructor.
    PriPartModelData(const PriPartModelData &copy); // Copy constructor.
    PriPartModelData(        // Stream-reading constructor.
        std::istream &in,    //  - Input stream.
        ParticleData &parent //  - Parent ParticleData object.
        );

    // Destructors.
    ~PriPartModelData(void);

    // Operators.
    PriPartModelData &operator=(const PriPartModelData &rhs);
    PriPartModelData &operator+=(const PriPartModelData &rhs);
    const PriPartModelData operator+(const PriPartModelData &rhs) const;

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
    PriPartModelData *const Clone(void) const;

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
    PriPartModelData(void);

private:
    std::vector<Primary> m_primaries;
};

};

#endif
