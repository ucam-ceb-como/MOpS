/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The PriPartData class holds particle data for the simple primary-particle
    model.
*/

#ifndef SWEEP_PRIPART_DATA_H
#define SWEEP_PRIPART_DATA_H

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

    // Returns the number of primary particles.
    unsigned int Count(void) const;

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
    // Vector of primary particles.
    std::vector<Primary> m_primaries;

    // The index of the component in the parent ParticleData object which
    // is used to define the mass and surface area of the primary particles.
    // The default value is 0 (the first component).
    unsigned int m_icomp;
};

};

#endif
