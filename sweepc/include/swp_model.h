/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The IModel class defines the interface for a particle model.
*/

#ifndef SWEEP_IMODEL_H
#define SWEEP_IMODEL_H

#include "swp_params.h"

namespace Sweep
{
class Particle;
class ParticleData;

class IModel
{
public:
    // Updates the particle according to the rules of the model
    // given the changes to particle composition and values.
    virtual void UpdateParticle(
        Particle &p,          // The particle which is being updated.
        const fvector &dcomp, // Changes to the composition.
        const fvector &dval   // Changes to the tracker variables.
        ) const = 0;

    // Updates the particle according to the rules of the model
    // given the changes to particle composition and values.
    // Performs the update n times.
    virtual void UpdateParticle(
        Particle &p,          // The particl which is being updated.
        const fvector &dcomp, // Changes to the composition.
        const fvector &dval,  // Changes to the tracker variables.
        unsigned int n        // Number of times to perform update.
        ) const = 0;

    // Updates the first particle according to the model rules for
    // coagulation.  The first particle is the recipient of the second.
    virtual void CoagParticles(
        Particle &p1,
        const Particle &p2
        ) const = 0;

    // Recalculates those model properties which are functions of
    // other particle properties.
    virtual void UpdateCache(ParticleData &data) const = 0;
protected:
    // IModels should implement the Singleton idiom, therefore
    // the default constructor, copy constructor and assignment
    // operator have to be made inaccessible.
    IModel(void);
    IModel(const IModel &copy);
    const IModel &operator=(const IModel &rhs);
    virtual ~IModel(void);
};
};

#endif
