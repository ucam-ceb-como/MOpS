/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Stochastic mechanism for Sweep.
*/

#ifndef SWEEP_MECHANISM_H
#define SWEEP_MECHANISM_H

#include <vector>
#include <string>
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
    enum ParticleModel {SphericalParticle, SurfaceVolume};
protected:
    bool m_anydeferred;
    vector<Inception*> m_inceptions; // Inception process list.
    vector<Process*> m_processes;    // Process list.
    vector<Coagulation*> m_coags;    // Coagulation process list.
    vector<Component*> m_components; // The components used to build particles.
    int m_icoag;                     // Index of first coagulation process in mechanism.
    SpeciesList *m_species;          // The species used to define the processes.
    unsigned int m_termcount;        // the rate term count of all processes.
    vector<string> m_valuenames;     // Names of additional values stored per particle.
    ParticleModel m_pmodel;           // Particle model used by this mechanism.
public:
    Mechanism(void);
    ~Mechanism(void);
public: // Overrides.
    int GetRates(vector<real> &rates, const real t, const System &sys) const;
    int GetStochasticRates(vector<real> &rates, const real t, const System &sys) const;
    /* Peforms the process with the given index on the given particle
       ensemble. */
    /* Calculates the combined rate of all stochastic processes (not
       deferred) given the individual process rates. */
    real JumpRate(const vector<real> &rates) const;
public:
    void SetSpeciesList(SpeciesList &list);
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
    /* Adds a coagulation process to the mechanism. */
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
    inline unsigned int ComponentCount(void) const {return (unsigned int)m_components.size();};
    Component const &GetComponent(const unsigned int i) const;
    vector<Component*> &GetComponents() {return m_components;};
    /* Returns the index of the component with the given name in the mechanism
       if found, otherwise return negative. */
    int GetComponentIndex(const string &name) const;
    unsigned int AddComponent(Component &comp);
    void SetComponent(const unsigned int i, Component &comp);
    void SetComponents(const vector<Component*> &comps);
public: // Value name get and sets.
    inline unsigned int ValueCount(void) const {return (unsigned int)m_valuenames.size();};
    void AddValueName(const string &name);
    void SetValueName(const unsigned int i, const string &name);
    string GetValueName(const unsigned int i) const;
    /* Returns the index of the tracker value with the given name on
       success, otherwise returns negative. */
    int GetValueIndex(const string &name) const;
public:
    inline ParticleModel GetParticleModel(void) const {return m_pmodel;};
    inline void SetParticleModel(const ParticleModel model) {m_pmodel=model;};
};

inline Inception &Mechanism::GetInception(const unsigned int i) 
{
    if (i < (unsigned int)m_inceptions.size()) {
        return *m_inceptions[i];
    } else {
        throw out_of_range("Index is outside the list of inceptions.");
    }
};

inline Process &Mechanism::GetProcess(const unsigned int i) 
{
    if (i < (unsigned int)m_processes.size()) {
        return *m_processes[i];
    } else {
        throw out_of_range("Index is outside list of processes.");
    }
};

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
};

#endif