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

#include <vector>
#include <string>
#include <set>
#include "swpparams.h"
#include "swpprocess.h"
#include "swpsystem.h"
#include "swpcomponent.h"
#include "gpcspecieslist.h"
#include "swpinception.h"
#include "swpcoagulation.h"

namespace Sweep
{
class Mechanism
{
public:
    // Enumeration of different particle types available to mechanisms.
    enum ParticleModel {SphericalParticle, SurfaceVolume};
    // enumeration of different models which can be loaded into sweep.
    enum Model {HACA, CNT, PAH};
protected:
    bool m_anydeferred;              // Are any of the mechanism's processes deferred?
    vector<Inception*> m_inceptions; // Inception process list.
    vector<Process*> m_processes;    // Process list.
    vector<Coagulation*> m_coags;    // Coagulation process list.
    vector<Component*> m_components; // The components used to build particles.
    int m_icoag;                     // Index of first coagulation process in mechanism.
    SpeciesList *m_species;          // The species used to define the processes.
    unsigned int m_termcount;        // the rate term count of all processes.
    vector<string> m_valuenames;     // Names of additional values stored per particle.
    ParticleModel m_pmodel;          // Particle model used by this mechanism.
    set<Model> m_reqmodels;          // Models required by this mechanism.
public:  // Default constructor and destructor.
    Mechanism(void);
    ~Mechanism(void);
public:
    /* Calculates the rates of all processes. */
    int GetRates(vector<real> &rates, const real t, const System &sys) const;
    /* Calculates the rates of all stochastic jump processes that are not deferred. */
    int GetStochasticRates(vector<real> &rates, const real t, const System &sys) const;
    /* Calculates the combined rate of all stochastic processes (not
       deferred) given the individual process rates. */
    real JumpRate(const vector<real> &rates) const;
public:
    /* Sets the list of species used by the mechanism to define processes. */
    void SetSpeciesList(SpeciesList &list);
    /* Returns reference to the species list used by the mechanism. */
    SpeciesList &GetSpeciesList(void) const;
public:
    /* Adds an inception to the mechanism. */
    int AddInception(Inception *icn);
    /* Returns the inception with the given index. */
    inline Inception &GetInception(const unsigned int i);
    /* Adds a process to the mechanism. */
    int AddProcess(Process *p);
    /* Returns the process with the given index. */
    Process &GetProcess(const unsigned int i);
    /* Adds the coagulation process to the mechanism. */
    int AddCoagulation();
    /* Returns the number of processes (including inceptions) in the mechanism. */
    inline unsigned int ProcessCount(void) const {
        return (unsigned int)m_inceptions.size() + (unsigned int)m_processes.size();
    };
    /* Returns the number of terms of all processes. */
    inline unsigned int TermCount(void) const {return m_termcount;};
    /* Returns true if the mechanism contains deferred (LPDA) processes
       otehrwise false. */
    inline bool AnyDeferred(void) const {return m_anydeferred;};
    /* Checks all processes to see if any are deferred. */
    void CheckDeferred(void);
public:
    /* Peforms the process with the given index on the given particle
       ensemble. */
    int DoProcess(const unsigned int i, // Index of process to perform.
                  const real t,         // Current time (s).
                  System &sys           // System to update (includes ensemble).
                 ) const;
    /* Performs linear update algorithm on the given system up to
       given time. */
    int LPDA(const real t, System &sys);
    /* Performs linear update algorithm on the given system up to given time,
       with the current chemical conditions precalculated. */
    int LPDA(const real t, const vector<real> &chem, const real T, const real P, 
             const vector<real> &sums, System &sys);
    /* Performs linear process updates on a particle in the given system. */
    int UpdateParticle(DefaultParticle &sp, System &sys, const real t);
    /* Performs linear process updates on a particle in the given system,
       with the current chemical conditions precalculated. */
    int UpdateParticle(DefaultParticle &sp, System &sys, const real t, 
                       const vector<real> &chem, const real T, const real P, 
                       const vector<real> &sums);
public: // Property get/set;
    /* Returns the number of particle components defined in the mechanism. */
    inline unsigned int ComponentCount(void) const {return (unsigned int)m_components.size();};
    /* Returns a reference to the particle component with the given index. */
    Component const &GetComponent(const unsigned int i) const;
    /* Returns a reference to the vector of all particle components defined in mechanism. */
    vector<Component*> &GetComponents() {return m_components;};
    /* Returns the index of the component with the given name in the mechanism
       if found, otherwise return negative. */
    int GetComponentIndex(const string &name) const;
    /* Adds a particle component to the mechanism. */
    unsigned int AddComponent(Component &comp);
    /* Sets the particle component with the given index. */
    void SetComponent(const unsigned int i, Component &comp);
    /* Sets the vector of all components. */
    void SetComponents(const vector<Component*> &comps);
public: // Value name get and sets.
    /* Returns the number of particle values defined in the mechanism. */
    inline unsigned int ValueCount(void) const {return (unsigned int)m_valuenames.size();};
    /* Adds a new value name to the mechanism. */
    void AddValueName(const string &name);
    /* Sets the value name at the given index. */
    void SetValueName(const unsigned int i, const string &name);
    /* Gets the value name at the given index. */
    string GetValueName(const unsigned int i) const;
    /* Returns the index of the tracker value with the given name on
       success, otherwise returns negative. */
    int GetValueIndex(const string &name) const;
public:
    /* Get the enum ID of the particle model currently used by the mechanism. */
    inline ParticleModel GetParticleModel(void) const {return m_pmodel;};
    /* Sets the ID of the particle model used by the mechanism. */
    inline void SetParticleModel(const ParticleModel model) {m_pmodel=model;};
public: // Routines to check which models are required by this mechanism.
    /* Adds a model to the list of required models. */
    void AddReqdModel(const Model m);
    /* Removes a model from the list of required models. */
    void RemoveReqdModel(const Model m);
    /* Checks the mechanism to see if it requires the given model.  Returns true
       if it does. */
    bool RequiresModel(const Model m) const;
    /* Initialises the required models for this mechanism. Returns 
       negative on error. */
    int InitReqdModels(void);
};

inline Inception &Mechanism::GetInception(const unsigned int i) 
{
    // Check the index is within range then return the inception.
    if (i < (unsigned int)m_inceptions.size()) {
        return *m_inceptions[i];
    } else {
        throw out_of_range("Index is outside the list of inceptions.");
    }
};

inline Process &Mechanism::GetProcess(const unsigned int i) 
{
    // Check the index is within range then return the process.
    if (i < (unsigned int)m_processes.size()) {
        return *m_processes[i];
    } else {
        throw out_of_range("Index is outside list of processes.");
    }
};

/* Particle value name routines. */

inline void Mechanism::AddValueName(const string &name)
{
    m_valuenames.push_back(name);
}

inline void Mechanism::SetValueName(const unsigned int i, const string &name)
{
    m_valuenames.reserve(i+1);
    m_valuenames[i] = name;
}

inline string Mechanism::GetValueName(const unsigned int i) const
{
    if (i<m_valuenames.size()) {
        return m_valuenames[i];
    } else {
        return "";
    }
}

/* Required model routines. */

inline void Mechanism::AddReqdModel(const Sweep::Mechanism::Model m)
{
    m_reqmodels.insert(m);
}

inline void Mechanism::RemoveReqdModel(const Sweep::Mechanism::Model m)
{
    m_reqmodels.erase(m);
}

inline bool Mechanism::RequiresModel(const Sweep::Mechanism::Model m) const
{
    set<Model>::const_iterator i = m_reqmodels.find(m);
    return (i!=m_reqmodels.end());
}

};

#endif