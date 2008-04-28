/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The SurfVolModel class implements the surface-volume model
    for coagulation and surface growth.
*/

#ifndef SWEEP_SURFVOL_MODEL_H
#define SWEEP_SURFVOL_MODEL_H

#include "swp_params.h"
#include "swp_coagmodel.h"

namespace Sweep
{
class SurfVolModel : public CoagModel
{
public:
    // Performs a coagulation event on the two given particles.
    // The resultant particle is stored in the location of p1.
    void Perform(Particle &p1, Particle &p2) const;


    // PARTICLE UPDATES.

    // Updates the particle according to the rules of the model
    // given the changes to particle composition and values.
    void UpdateParticle(
        Particle &p,          // The particle which is being updated.
        const fvector &dcomp, // Changes to the composition.
        const fvector &dval   // Changes to the tracker variables.
        ) const;

    // Updates the particle according to the rules of the model
    // given the changes to particle composition and values.
    // Performs the update n times.
    void UpdateParticle(
        Particle &p,          // The particl which is being updated.
        const fvector &dcomp, // Changes to the composition.
        const fvector &dval,  // Changes to the tracker variables.
        unsigned int n        // Number of times to perform update.
        ) const;

    // Updates the first particle according to the model rules for
    // coagulation.  The first particle is the recipient of the second.
    void CoagParticles(
        Particle &p1, 
        const Particle &p2
        ) const;

    // Recalculates those model properties which are function of
    // other particle properties.
    void UpdateCache(ParticleData &p) const;


    // SINGLETON IMPLEMENTATION.

    // Returns the one and only SurfVolModel object.
    static SurfVolModel &Instance(void);

protected:
    // Class implements the Singleton idiom, therefore
    // the default constructor, copy constructor and assignment
    // operator have to be made inaccessible.
    SurfVolModel(void);
    SurfVolModel(const SurfVolModel &copy);
    SurfVolModel &operator=(const SurfVolModel &rhs);
    ~SurfVolModel(void);

private:
    // Number of stats in the point contact model.
    static const unsigned int STAT_COUNT = 1;
};
};

#endif
