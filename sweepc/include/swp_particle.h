/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a default particle used by sweep.  Particles
    are assumed spherical and are only described by their volume.

    This class serves as the base class for all other stochastic particle types
    defined in sweep.
*/

#ifndef SWEEP_PARTICLE_H
#define SWEEP_PARTICLE_H

#include "swp_params.h"
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_particledata.h"
#include <vector>
#include <iostream>

namespace Sweep
{
class Ensemble;
class Mechanism;

class Particle : public ParticleData
{
public:
	// Constructors.
    Particle( // Default constructor.
        const CompPtrVector &components, 
        const TrackPtrVector &trackers
        );  
    Particle(const Particle &copy); // Copy constructor.
    Particle(                 // Stream-reading constructor.
        std::istream &in,     //  - Input stream.
        const Mechanism &mech //  - Mechanism which defines components and trackers.
        );

	// Destructor.
    ~Particle(void);

    // Operators.
    Particle &operator=(const Particle &rhs);
    Particle &operator+=(const Particle &rhs);
    const Particle operator+(const Particle &rhs) const;


    // PARTICLE ADJUSTMENT.

    // Adjusts the particle with the given composition and value
    // changes.
    void Adjust(
        const fvector &dcomp, 
        const fvector &dvalues
        );

    // Adjusts the particle with the given composition and 
    // values changes n times.  If the particle cannot be adjust
    // n times, then this function returns the number of times
    // it was adjusted.
    unsigned int Adjust(
        const fvector &dcomp, 
        const fvector &dvalues, 
        unsigned int n
        );

    // Combines this particle with another.
    Particle &Coagulate(const Particle &sp);


    // PARTICLE UPDATE AND CHECKING.

    // Recalculates the derived properties from the 
    // unique properties.
    void UpdateCache(void);

    // Check the that the particle is valid by querying the
    // validity conditions of the models and ensuring that it 
    // contains any components.
    bool IsValid() const;


    // PARENT ENSEMBLE.

    // Returns the parent ensemble.
    const Sweep::Ensemble *const Ensemble(void) const;

    // Sets the parent ensemble.
    void SetEnsemble(Sweep::Ensemble &ens);


    // READ/WRITE/COPY.

    // Creates a clone of the particle.
    Particle *const Clone();

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
