/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The PointContactData class is a specialisation of the CoagModelData
    class which holds additional data for the surface-volume particle
    model.
*/

#ifndef SWEEP_POINTCONTACTDATA_H
#define SWEEP_POINTCONTACTDATA_H

#include "swp_params.h"
#include "swp_coagmodeldata.h"
#include "swp_pointcontactmodel.h"
#include <iostream>

namespace Sweep
{
class PointContactData : public CoagModelData
{
friend class PointContactModel;

public:
    // Constructors.
    PointContactData(ParticleData &parent); // Default constructor.
    PointContactData(const PointContactData &copy); // Copy constructor.
    PointContactData(        // Stream-reading constructor.
        std::istream &in,    //  - Input stream.
        ParticleData &parent //  - Parent ParticleData object.
        );

    // Destructors.
    ~PointContactData(void);

    // Operators.
    PointContactData &operator=(const PointContactData &rhs);
    PointContactData &operator+=(const PointContactData &rhs);
    const PointContactData operator+(const PointContactData &rhs) const;

    // Resets the model data to the default state.
    void Clear();


    // COAGULATION MODEL PARTICLE PROPERTIES.

    // Returns the equivalent spherical surface area.
    real SphSurfaceArea(void) const;

    // Returns the actual surface area.
    real SurfaceArea(void) const;


    // MODEL WHICH USES THIS DATA.

    // Returns the PointContactModel which operates on this data.
    const PointContactModel &Model(void) const;


    // READ/WRITE/COPY.

    // Returns a copy of the data.
    PointContactData *const Clone(void) const;

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

    // Can't create a PointContactData without knowledge
    // of the parent ParticleData.
    PointContactData(void);
};
};

#endif
