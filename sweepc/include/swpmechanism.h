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
	/* CONSTRUCTORS AND DESTRUCTOR. */
	Mechanism(void);                      // Default Constructor.
	Mechanism(const Mechanism & mech);    // Copy-Constructor.
    ~Mechanism(void);                     // Destructor.


    enum ParticleModel {SphericalParticle, SurfaceVolume};   // enumeration of particle model used.


	/* OVERRIDES. */
    int GetRates(vector<real> &rates, const real t, const System &sys) const;           // Get total rates.
    int GetStochasticRates(vector<real> &rates, const real t, const System &sys) const; /* Get total rates of 
																						non-deferred processes */
    real JumpRate(const vector<real> &rates) const;  /* Calculates the combined rate of all stochastic processes (not
                                                     deferred) given the individual process rates. */


	/* SPECIESLIST FUNCTIONS */
	void SetSpeciesList(SpeciesList &list);     // Set m_species point to the given SpeciesList.
    SpeciesList &GetSpeciesList(void) const;    // Return reference to m_species.


	/* PROCESS FUNCTIONS */
    int AddInception(Inception *icn);                      // Adds an inception to the mechanism.
    inline Inception &GetInception(const unsigned int i);  // Returns the inception with the given index.
    int AddProcess(Process *p);                            // Adds a process to the mechanism.
    Process &GetProcess(const unsigned int i);             // Returns the process with the given index.
    int AddCoagulation();                                  // Adds a coagulation process to the mechanism.

    // Returns the number of processes (including inceptions) in the mechanism.
    inline unsigned int ProcessCount(void) const {                                
        return (unsigned int)m_inceptions.size() + (unsigned int)m_processes.size();
    };
    inline unsigned int TermCount(void) const {return m_termcount;};  // Returns the number of terms of all processes.
    inline bool AnyDeferred(void) const {return m_anydeferred;};      /* Returns true if the mechanism contains 
																	  deferred (LPDA) processes otherwise false. */
    void CheckDeferred(void);                                         // Checks all processes to see if any are deferred.



	/* PERFORMING THE PROCESSES. */

    // Performs the Process specified.
	int DoProcess(const unsigned int i, // Index of process to perform.
                  const real t,         // Current time (s).
                  System &sys           // System to update (includes ensemble).
                 ) const;

    int LPDA(const real t, System &sys);  // Performs linear update algorithm on the given system up to given time.
    /* Performs linear update algorithm on the given system up to given time,
       with the current chemical conditions precalculated. */
    int LPDA(const real t, const vector<real> &chem, const real T, const real P, 
             const vector<real> &sums, System &sys);

    /* Performs linear process updates on a particle in the given system. */
    int UpdateParticle(Particle &sp, System &sys, const real t);
    /* Performs linear process updates on a particle in the given system,
       with the current chemical conditions precalculated. */
    int UpdateParticle(Particle &sp, System &sys, const real t, 
                       const vector<real> &chem, const real T, const real P, 
                       const vector<real> &sums);

	/* PROPERTY GET/SET FUNCTIONS*/
    inline unsigned int ComponentCount(void) const {            // Return size of m_components.
		return (unsigned int)m_components.size();
	};
    Component const &GetComponent(const unsigned int i) const;  /* returns m_components[i], but if i not valid index, 
	                                                            then returns m_components[0] */
    vector<Component*> &GetComponents() {return m_components;}; // returns m_components.
    int GetComponentIndex(const string &name) const;            /* Returns the index of the component with the 
																given name in the mechanism if found, otherwise 
																return negative. */
    unsigned int AddComponent(Component &comp);   // Adds a component to m_components - returns previous size of vector?
    void SetComponent(const unsigned int i, Component &comp);  // Set m_components[i] to be the given component.
    void SetComponents(const vector<Component*> &comps);       // Assign comps vector to m_components vector.


	/* VALUE NAME GET AND SETS */
    inline unsigned int ValueCount(void) const {return (unsigned int)m_valuenames.size();};  // Returns m_valuenames.size()
    void AddValueName(const string &name);                       // Add the given string as a new element of m_valuenames
    void SetValueName(const unsigned int i, const string &name); // Set ith element of m_valuenames with given string.
    string GetValueName(const unsigned int i) const;             // Get m_valuenames[i].
    
    int GetValueIndex(const string &name) const;   /* Returns the index of the tracker value with the given name 
												   on success, otherwise returns negative. */
    inline ParticleModel GetParticleModel(void) const {return m_pmodel;};       // Return m_pmodel.
    inline void SetParticleModel(const ParticleModel model) {m_pmodel=model;};  // Set particle model with given model.

protected:

	/* PROTECTED VARIABLES. */
    bool m_anydeferred;              // True if the mechanism contains deferred (LPDA) processes otherwise false. 
    vector<Inception*> m_inceptions; // Inception process list.
    vector<Process*> m_processes;    // Process list.
    vector<Coagulation*> m_coags;    // Coagulation process list.
    vector<Component*> m_components; // The components used to build particles.
    int m_icoag;                     // Index of first coagulation process in mechanism.
    SpeciesList *m_species;          // The species used to define the processes.
    unsigned int m_termcount;        // the rate term count of all processes.
    vector<string> m_valuenames;     // Names of additional values stored per particle.
    ParticleModel m_pmodel;          // Particle model used by this mechanism.
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
