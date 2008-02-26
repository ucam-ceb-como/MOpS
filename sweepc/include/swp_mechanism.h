/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Stochastic mechanism for Sweep.  The mechanism holds all the processes
    which can be enacted on a system with a particle ensemble.  It also holds
    auxilliary info which defines how those processes work.

    The mechanism class defines routines for calculating process rates and perform
    a process.  It also contains the Linear Process Deferment Algorithm (LPDA) for
    updating linear processes that have been removed from the stochastic jump process.

    The mechanism also defines particle components and additional particle values.  Particle
    values are particle properties which can be changed by processes, but which are not 
    components, therefore they do not contribute to particle mass.
*/

#ifndef SWEEP_MECHANISM_H
#define SWEEP_MECHANISM_H

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_process.h"
#include "swp_particleprocess.h"
#include "swp_inception.h"
#include "swp_coagulation.h"
#include "swp_modeltype.h"
#include "sprog.h"
#include <vector>
#include <string>

namespace Sweep
{
class Mechanism
{
public:
	// Constructors.
	Mechanism(void);                  // Default Constructor.
	Mechanism(const Mechanism &copy); // Copy-Constructor.

    // Destructor.
    ~Mechanism(void);

    // Operators.
    Mechanism &operator=(const Mechanism &rhs);


	// CHEMICAL SPECIES.

    // Returns the chemical species vector.
    const Sprog::SpeciesPtrVector *const Species(void) const;

    // Sets the chemical species vector.
    void SetSpecies(const Sprog::SpeciesPtrVector &sp);


	// COMPONENT DEFINITIONS.

    // Returns the number of components in the mechanism.
    unsigned int ComponentCount(void) const;

    // Returns the vector of particle components.
    const CompPtrVector &Components() const;

    // Returns the component with the given index.
    const Component *const Components(unsigned int i) const;

    // Returns the index of the component with the 
	// given name in the mechanism if found, otherwise 
	// return negative.
    int ComponentIndex(const std::string &name) const;

    // Adds a component to the mechanism and returns the index
    // of the component.
    unsigned int AddComponent(Component &comp);

    // Overwrites the ith component with that given.  Previous
    // component is deleted from memory.
    void ReplaceComponent(
        unsigned int i, // Index of component to overwrite.
        Component &comp // New component.
        );

    // Sets the particle components vector.
    void SetComponents(const CompPtrVector &comps);


	// TRACKER VARIABLES.

    // Returns the number of tracker variables.
    unsigned int TrackerCount(void) const;

    // Returns the vector of tracker variables.
    const TrackPtrVector &Trackers(void) const;

    // Returns the ith tracker variable.
    const Tracker *const Trackers(unsigned int i) const;

    // Returns the index of the tracker variable with the given name 
	// on success, otherwise returns negative.
    int GetTrackerIndex(const std::string &name) const;

    // Adds a tracker variable to the mechanism.
    void AddTracker(Tracker &track);

    // Replaces the tracker at the given index with the given tracker
    // object.
    void ReplaceTracker(
        unsigned int i, // Index of tracker variable to overwrite.
        Tracker &track  // New track variable.
        );

    // Sets the vector of tracker variables.
    void SetTrackers(const TrackPtrVector &track);


	// INCEPTIONS.

    // Returns the vector of inceptions.
    const IcnPtrVector &Inceptions(void) const;

    // Returns the inception with the given index.
    const Inception *const Inceptions(unsigned int i) const;

    // Adds an inception to the mechanism.
    void AddInception(Inception &icn);


    // PARTICLE PROCESSES.

    // Returns the vector of particle processes.
    const PartProcPtrVector &Processes(void) const;

    // Returns the process with the given index.
    const ParticleProcess *const Processes(unsigned int i) const;

    // Adds a process to the mechanism.
    void AddProcess(ParticleProcess &p);


    // COAGULATIONS.

    // Adds a coagulation process to the mechanism.
    void AddCoagulation();


    // PROCESS INFORMATION.

    // Returns the number of processes (including 
    // inceptions) in the mechanism.
    unsigned int ProcessCount(void) const;

    // Returns the number of terms in all process rate expressions.
    unsigned int TermCount(void) const;

    // Returns true if the mechanism contains deferred (LPDA) 
    // processes otherwise false.
    bool AnyDeferred(void) const;

    // Checks all processes to see if any are deferred.
    void CheckDeferred(void) const;


    // PARTICLE MODELS.

    // Returns the set of particle model ID used by this mechanism
    const ModelTypeSet &Models(void) const;

    // Returns true if the mechanism include the given model.
    bool ContainsModel(ModelType id) const;

    // Adds a model to the mechanism.  Any subsequent particles
    // created with this mechanism will use this model.
    void AddModel(ModelType id);


	// RATE CALCULATION.

    // Get total rates of all processes.  Returns the sum of
    // all rates.
    real CalcRates(
        real t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        fvector &rates   // Return vector for process rates.
        ) const;

    // Get total rates of all processes.  Returns the sum of
    // all rates.  Uses supplied gas-phase conditions rather than
    // those in the given system.
    real CalcRates(
        real t,          // Time at which to get rates.
        const Sprog::Thermo::IdealGas &gas, // Use these gas-phase conditions instead.
        const Cell &sys, // System cell for which to get rates.
        fvector &rates   // Return vector for process rates.
        ) const;


    // Get total rates of non-deferred processes.  Returns the sum
    // of all rates.
    real CalcJumpRates(
        real t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        fvector &rates   // Return vector for process rates.
        ) const;

    // Get total rates of non-deferred processes.  Returns the sum
    // of all rates.  Uses supplied gas-phase conditions rather than
    // those in the given system.
    real CalcJumpRates(
        real t,          // Time at which to get rates.
        const Sprog::Thermo::IdealGas &gas, // Use these gas-phase conditions instead.
        const Cell &sys, // System cell for which to get rates.
        fvector &rates   // Return vector for process rates.
        ) const;


	// PERFORMING THE PROCESSES.

    // Performs the Process specified.  Process index could be
    // an inception, particle process or a coagulation event.
	void DoProcess(
        unsigned int i, // Index of process to perform.
        real t,         // Current time (s).
        Cell &sys       // System to update (includes ensemble).
        ) const;


    // LINEAR PROCESS DEFERMENT ALGORITHM.

    // Performs linear update algorithm on the 
    // given system up to given time.
    void LPDA(
        real t,   // Time up to which to integrate.
        Cell &sys // System to update.
        ) const;

    // Performs linear update algorithm on the given system up to given time,
    // with the given chemical conditions rather than those in the 
    // given system.
    void LPDA(
        real t,    // Time up to which to integrate.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        Cell &sys // System to update.
        ) const;

    // Performs linear process updates on a particle in the given system.
    void UpdateParticle(
        Particle &sp, // Particle to update.
        Cell &sys,    // System to which the particle belongs.
        real t        // Time up to which to integrate.
        ) const;

    // Performs linear process updates on a particle in the given system,
    // with the current chemical conditions precalculated.
    void UpdateParticle(
        Particle &sp, // Particle to update.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        Cell &sys,    // System to which the particle belongs.
        real t        // Time up to which to integrate.
        ) const;


    // PARTICLE FUNCTIONS.
    
    // Creates a new particle and sets it up with all the models
    // required by the mechanism.
    Particle *const CreateParticle(void) const;

    // READ/WRITE/COPY.

    // Creates a copy of the mechanism.
    Mechanism *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

private:
    // True if the mechanism contains deferred (LPDA)
    // processes, otherwise false. 
    mutable bool m_anydeferred;

    // The species used to define the processes and the particles.
    const Sprog::SpeciesPtrVector *m_species;

    // Particle component and tracker variable definitions.
    CompPtrVector m_components; // The components used to build particles.
    TrackPtrVector m_trackers;  // Tracker variables.

    // Processes in mechanism.
    IcnPtrVector m_inceptions;     // Inception process list.
    PartProcPtrVector m_processes; // Particle process list.
    Coagulation *m_coag;           // Coagulation process.

    // Auxilliary information about the processes.
    int m_icoag;              // Index of first coagulation process in mechanism.
    unsigned int m_termcount; // the rate term count of all processes.

    // Set of models which all particles produced by this mechanism
    // must use.
    ModelTypeSet m_models;

    // Process counters.
    mutable std::vector<unsigned int> m_proccount, m_fictcount; 

    // Clears the current mechanism from memory.
    void releaseMem(void);
};
};
#endif
