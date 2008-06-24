/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The Particle class represents the top node in the sub-particle tree.  This
    object is placed in the Ensemble and exhibits an interface which allows
    bulk particle properties to be summed and stored in the ensemble binary tree.
*/

#ifndef SWEEP_PARTICLE_H
#define SWEEP_PARTICLE_H

#include "swp_params.h"
#include "swp_particle_model.h"
#include "swp_subparticle.h"
#include <vector>
#include <iostream>

namespace Sweep
{
class Ensemble;
class Mechanism;

class Particle : public SubParticle
{
public:
	// Constructors.
    Particle(                             // Initialising constructor.
        real time,                        // Create time.
        const Sweep::ParticleModel &model // Defining particle model.
        );
    Particle(Sweep::Primary &pri);        // Initialising constructor (from primary).
    Particle(const Particle &copy);       // Copy constructor.
    Particle(                             // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Model to which this particle subscribes.
        );

	// Destructor.
    ~Particle(void);

    // Operators.
    Particle &operator=(const Particle &rhs);
    Particle &operator+=(const Particle &rhs);
    const Particle operator+(const Particle &rhs) const;


    // PARENT ENSEMBLE.

    // Returns the parent ensemble.
    const Sweep::Ensemble *const Ensemble(void) const;

    // Sets the parent ensemble.
    void SetEnsemble(Sweep::Ensemble &ens);


    // READ/WRITE/COPY.

    // Creates a clone of the particle.
    Particle *const Clone() const;

private:
    // Parent ensemble.
    Sweep::Ensemble *m_ensemble;

    // Can't create a particle without knowledge of the components
    // and the tracker variables.
    Particle(void);
};

typedef std::vector<Particle*> PartPtrVector;
};

#endif
