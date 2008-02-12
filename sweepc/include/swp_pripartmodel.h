/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The PriPartModel class defines the spherical primary particle model
    which describes aggregate constituent primary particles of variable
    size.
*/

#ifndef SWEEP_PRIPARTMODEL_H
#define SWEEP_PRIPARTMODEL_H

#include "swp_params.h"
#include "swp_model.h"
#include "swp_primary.h"
#include "swp_component.h"
#include <vector>

namespace Sweep
{
class PriPartModel : public IModel
{
public:
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

    // Recalculates those model properties which are functions of
    // other particle properties.
    virtual void UpdateCache(ParticleData &data) const;

    // Returns the one and only PriPartModel object.
    static PriPartModel &Instance(void);
protected:
    // PriPartModel implements the Singleton idiom, therefore
    // the default constructor, copy constructor and assignment
    // operator have to be made inaccessible.
    PriPartModel(void);
    PriPartModel(const PriPartModel &copy);
    const PriPartModel &operator=(const PriPartModel &rhs);
    virtual ~PriPartModel(void);

private:
    // Distributes mass over a vector of primaries.
    static void distMass(
        std::vector<Primary> &pri,  // Vector of primary particles to receive mass.
        real dmass,                 // Mass change.
        const Component *const comp // Component which defines the mass.
        );
};
};

#endif
