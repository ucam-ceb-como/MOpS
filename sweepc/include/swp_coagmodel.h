/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The CoagModel class is a specialisation of the Model class for
    coagulation processes.
*/

#ifndef SWEEP_COAGMODEL_H
#define SWEEP_COAGMODEL_H

#include "swp_params.h"
#include "swp_model.h"

namespace Sweep
{
class CoagModel : public IModel
{
public:
    // Performs a coagulation event on the two given particles.
    // The resultant particle is stored in the location of p1.
    virtual void Perform(Particle &p1, Particle &p2) const;

    // Updates the particle according to the rules of the model
    // given the changes to particle composition and values.
    virtual void UpdateParticle(
        Particle &p,          // The particle which is being updated.
        const fvector &dcomp, // Changes to the composition.
        const fvector &dval   // Changes to the tracker variables.
        ) const;

    // Updates the particle according to the rules of the model
    // given the changes to particle composition and values.
    // Performs the update n times.
    virtual void UpdateParticle(
        Particle &p,          // The particl which is being updated.
        const fvector &dcomp, // Changes to the composition.
        const fvector &dval,  // Changes to the tracker variables.
        unsigned int n        // Number of times to perform update.
        ) const;

    // Updates the first particle according to the model rules for
    // coagulation.  The first particle is the recipient of the second.
    virtual void CoagParticles(
        Particle &p1, 
        const Particle &p2
        ) const;

    // Recalculates those model properties which are function of
    // other particle properties.
    virtual void UpdateCache(ParticleData &p) const;

    // Returns the one and only CoagModel object.
    static CoagModel &Instance(void);
protected:
    // CoagModel class implements the Singleton idiom, therefore
    // the default constructor, copy constructor and assignment
    // operator have to be made inaccessible.
    CoagModel(void);
    CoagModel(const CoagModel &copy);
    CoagModel &operator=(const CoagModel &rhs);
    ~CoagModel(void);
};
};

#endif
