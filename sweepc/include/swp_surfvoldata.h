/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The SurfVolData class is a specialisation of the CoagModelData
    class which holds additional data for the surface-volume particle
    model.
*/

#ifndef SWEEP_SURFVOL_DATA_H
#define SWEEP_SURFVOL_DATA_H

#include "swp_params.h"
#include "swp_coagmodeldata.h"
#include "swp_surfvolmodel.h"
#include <iostream>

namespace Sweep
{
class SurfVolData : public CoagModelData
{
friend class SurfVolModel;

public:
    // Constructors.
    SurfVolData(ParticleData &parent);    // Default constructor.
    SurfVolData(const SurfVolData &copy); // Copy constructor.
    SurfVolData(             // Stream-reading constructor.
        std::istream &in,    //  - Input stream.
        ParticleData &parent //  - Parent ParticleData object.
        );

    // Destructors.
    ~SurfVolData(void);

    // Operators.
    SurfVolData &operator=(const SurfVolData &rhs);
    SurfVolData &operator+=(const SurfVolData &rhs);
    const SurfVolData operator+(const SurfVolData &rhs) const;

    // Resets the model data to the default state.
    void Clear();


    // COAGULATION MODEL PARTICLE PROPERTIES.

    // Returns the equivalent spherical surface area.
    real SphSurfaceArea(void) const;

    // Returns the actual surface area.
    real SurfaceArea(void) const;


    // MODEL WHICH USES THIS DATA.

    // Returns the SurfVolModel which operates on this data.
    const SurfVolModel &Model(void) const;


    // READ/WRITE/COPY.

    // Returns a copy of the data.
    SurfVolData *const Clone(void) const;

    // Returns the model ID.
    ModelType ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

private:
    // Point contact model particle properties.
    real m_sphsurf; // Equivalent sphere surface area.
    real m_surf;    // Actual surface area.

    // Can't create a SurfVolData without knowledge
    // of the parent ParticleData.
    SurfVolData(void);
};
};

#endif
