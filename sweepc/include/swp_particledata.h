/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ParticleData class defines the data structure of a particle.  It contains
    the particle composition, basic properties and the additional properties
    required for the various particle models.  This class is the interface between
    the particle class and the ensemble, as both those classes contain data objects
    of type ParticleData.
*/

#ifndef SWEEP_PARTICLEDATA_H
#define SWEEP_PARTICLEDATA_H

#include "swp_params.h"
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_modeltype.h"
#include "swp_modeldata.h"
#include "swp_coagmodeldata.h"
#include <iostream>

namespace Sweep
{
class Ensemble;
class TreeNode;
class Mechanism;

class ParticleData
{
friend Ensemble;
friend TreeNode;

public:
    // Enumeration of ParticleData properties which can be accessed using
    // the Property() function.
    enum PropertyID {
        iCTime,  // Create time.
        iLUTime, // Last update time.
        iD,      // Equivalent sphere diameter.
        iDcol,   // Collision diameter.
        iDmob,   // Mobility diameter.
        iS,      // Surface area.
        iV,      // Volume.
        iM       // Mass.
    };

    // Constructors.
    ParticleData(                            // Default constructor.
        const CompPtrVector &components,
        const TrackPtrVector &trackers
        ); 
    ParticleData(const ParticleData &copy);  // Copy constructor.
    ParticleData(
        std::istream &in,     // Input stream.
        const Mechanism &mech // Mechanism which defines components and trackers.
        );

    // Destructors.
    virtual ~ParticleData(void);

    // Operators.
    ParticleData &operator=(const ParticleData &rhs);
    ParticleData &operator+=(const ParticleData &rhs);
    const ParticleData operator+(const ParticleData &rhs) const;

    // Resets the particle data to its "empty" condition.
    void Clear(void);


    // DEFINING COMPONENTS.
    
    // Returns the components vector.
    const CompPtrVector *const Components(void) const;

    // Sets the components vector.
    void SetComponents(const CompPtrVector &comps);


    // DEFINING VALUES.

    // Returns the tracker values vector.
    const TrackPtrVector *const Trackers(void) const;

    // Sets the values/trackers vector.
    void SetTrackers(const TrackPtrVector &track);


    // PARTICLE COMPOSITION.

    // Returns the composition vector.
    const fvector &Composition(void) const;

    // Returns the ith component value.  Returns 0.0 if i invalid.
    real Composition(unsigned int i) const;

    // Sets the composition vector.
    void SetComposition(const fvector &comp);


    // PARTICLE VALUES.

    // Returns the values vector.
    const fvector &Values(void) const;

    // Returns the ith value.  Returns 0.0 if i invalid.
    real Values(unsigned int i) const;

    // Sets the values vector.
    void SetValues(const fvector &vals);


    // PARTICLE CREATE TIME.

    // Returns the particle create time.
    real CreateTime(void) const;

    // Sets the particle create time.
    void SetCreateTime(real t);


    // LAST UPDATE TIME.

    // Returns the last update time of the particle.
    real LastUpdateTime(void) const;

    // Sets the last update time of the particle.
    void SetTime(real t);


    // COAGULATION MODEL CACHE.

    // Returns the coagulation model data.
    CoagModelData *const CoagModelCache(void);
    const CoagModelData *const CoagModelCache(void) const;

    // Sets the coagulation model data.  This ParticleData object
    // then takes control of the data for destruction purposes.
    void SetCoagModelCache(CoagModelData &data);


    // MODEL CACHE.

    // Returns the model data.
    ModelMap &ModelCache(void);

    // Returns the data for the ith model.  Returns NULL if i is invalid.
    IModelData *const ModelCache(ModelType id);
    const IModelData *const ModelCache(ModelType id) const;

    // Add a model to the particle definition.  This ParticleData object
    // then takes control of the model for destruction purposes.
    void AddModel(IModelData &model);

    // Removes the ith model from the particle definition.
//    void RemoveModel(unsigned int i);


    // BASIC PROPERTIES.

    // Returns the particle equivalent sphere diameter.
    real SphDiameter(void) const;

    // Returns the collision diameter.
    real CollDiameter(void) const;

    // Returns the mobility diameter.
    real MobDiameter(void) const;

    // Returns the surface area.
    real SurfaceArea(void) const;

    // Returns the volume.
    real Volume(void) const;

    // Returns the mass.
    real Mass(void) const;

    // Returns the property with the given ID.
    real Property(PropertyID id) const;


    // BASIC DERIVED PROPERTY OVERWRITES.

    // Sets the spherical particle diameter
    void SetSphDiameter(real diam);

    // Sets the collision diameter of the particle.
    void SetCollDiameter(real dcol);

    // Sets the mobility diameter.
    void MobDiameter(real dmob);

    // Sets the surface area, subject to minimum spherical area condition.
    void SetSurfaceArea(real surf);

    // Sets the volume.
    void SetVolume(real vol);

    // Sets the mass.
    void SetMass(real m);


    // READ/WRITE/COPY.

    // Creates a copy of the particle data object.
    virtual ParticleData *const Clone(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,
        const Mechanism &mech
        );

protected:
    const CompPtrVector *m_components; // Components used to define properties.
    const TrackPtrVector *m_trackers;  // Tracker values used to define particle data.

    // Unique properties.
    fvector m_comp;   // Particle composition.
    fvector m_values; // Other particle values (defined at run time).
    real m_createt;   // Time at which particle was created.
    real m_time;      // Last time particle was updated.  Required for LPDA.

    // Additional particle models.
    CoagModelData *m_coag; // The coagulation model used by the particle.
    ModelMap m_models;     // Cache of data required by different particle models.

    // Basic properties (calculated from above properties).
    real m_diam; // Equivalent spherical diameter.
    real m_dcol; // Collision diameter.
    real m_dmob; // Mobility diameter.
    real m_surf; // Surface area.
    real m_vol;  // Volume.
    real m_mass; // Mass.

    // ParticleData class cannot be created without knowledge of the
    // components and the values.
    ParticleData(void);

    // Release all memory associated with the ParticleData object.
    void releaseMem(void);

    // Initialisation routine.
    void init(void);
};
};

#endif
