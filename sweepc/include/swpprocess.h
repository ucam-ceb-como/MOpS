/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    A specialised process for Sweep, which acts on Sweep::System
    objects.
*/

#ifndef SWEEP_PROCESS_H
#define SWEEP_PROCESS_H

#include "swpparams.h"
#include "swpsystem.h"
#include "swpparticle1d.h"
#include "rng.h"
#include <map>

namespace Sweep
{
class Mechanism;

class Process
{
public:

	/* CONSTRUCTORS AND DESTRUCTOR */
    Process(void);                  // Default Constructor.
	Process(const Process & proc);  // Copy-constructor.
    ~Process(void);                 // Destructor.

	/* GET/SET MECHANISMS */
	inline Mechanism &GetMechanism() {return *m_mech;};            // Return pointer to the parent Mechanism m_mech.
    inline void SetMechanism(Mechanism &mech) {m_mech = &mech;};   // Set m_mech to given Mechanism.



	/* TOTAL RATE CALCULATIONS - PURE VIRTUAL FUNCTIONS */
    virtual real Rate(const real t, const System &sys) const = 0;  // Returns rate of the process for the given system.

	// Calculates the process rate given the chemical conditions.
    virtual real Rate(const real t, const vector<real> &chem,  const real T,
                      const real P, const vector<real> &sums, const System &sys) const = 0;

	// SINGLE PARTICLE RATE CALCULATIONS - PURE VIRTUAL FUNCTIONS.
    /* Returns the rate of the process for the given particle in
       the system. Process must be linear in particle number. */
    virtual real Rate(const real t,             // Current time (s).
                      const System &sys,        // System for which to get rate.
                      const Particle &sp        // Address of particle for which to calculate rate.
                     ) const = 0;
	/* Function used to calculate relevant rate for the above 'rate' function */
    virtual real Rate(const real t,               // Current time (s).
                      const vector<real> &chem,   // Species concentrations.
                      const real T, const real P, // Temperature and pressure.
                      const vector<real> &sums,   // Particle sums from the ensemble.
                      const System &sys,          // System for which to get rate.
                      const Particle &sp          // Address of particle for which to calculate rate.
                     ) const = 0;


	/* RATE TERM CALCULATIONS. */
    virtual unsigned int TermCount(void) const = 0;  // Returns the number of rate terms for this process.
    /* Calculates the rate terms given an iterator to a real vector. The 
       iterator is advanced to the position after the last term for this
       process. */
    virtual void RateTerms(const real t, const System &sys, vector<real>::iterator &iterm) const = 0;
    virtual void RateTerms(const real t, const vector<real> &chem, const real T, 
                           const real P, const vector<real> &sums, const System &sys, 
                           vector<real>::iterator &iterm) const = 0;


	/* PERFORMING THE PROCESS */
    /* Performs the process on the given system.  The responsible rate term is given
       by index.  Returns 0 on success, otherwise negative. */
    virtual int Perform(const real t, System &sys, const unsigned int iterm) const = 0;
    /* Performs the process on a given particle in the system.  Particle
       is given by index.  The process is performed n times. */
    virtual int Perform(const real t,         // Current time (s).
                        System &sys,          // System for which to perform process. 
                        Particle &sp,         // Address of particle for which to perform process.
                        const unsigned int n  // Number of times to perform the process.
                       ) const = 0;

	
    inline bool IsDeferred(void) const {return m_defer;}; // Returns TRUE if process should be deferred, otherwise false.
    virtual void SetDeferred(const bool defer);  // Sets m_defer to FALSE.

	/* Determines whether a rate is ficticious given the majorant and true values. */
    static inline bool Ficticious(const real majk, const real truek) {
        return !((majk*rnd()) < truek);
    };


	/* PROPERTY GETS */
    /* Returns the number of reactants. */
    inline unsigned int ReactantCount(void) const {return (unsigned int)m_reac.size();};


	/* PROPERTY SETS */
    void SetReactant(const unsigned int id, const int stoich);  // Adds a reactant to the inception.
    void SetProduct(const unsigned int id, const int stoich);   // Adds a product to the inception.

	/* Applies the changes to the gas-phase to the given chemical system. */
    void ApplyToSystem(System &sys) const;
    /* Applies the changes to the gas-phase to the given chemical system n times. */
    void ApplyToSystem(System &sys, const unsigned int n) const;

protected:
    Mechanism *m_mech;             // Pointer to the parent Mechanism.
    bool m_defer;                  // Is the process solved by LPDA?
    map<unsigned int, int> m_reac; // Reactant species stoichiometry.
    map<unsigned int, int> m_prod; // Product species stoichiometry.

};


inline void Process::ApplyToSystem(System &sys) const
{
    // Apply changes to gas-phase chemistry.
    map<unsigned int,int>::const_iterator j;
    real NAvol = NA * sys.SampleVolume();
    for (j=m_reac.begin(); j!=m_reac.end(); j++)
        sys.AdjustSpeciesConc(j->first, -(real)(j->second) / (NAvol));
    for (j=m_prod.begin(); j!=m_prod.end(); j++)
        sys.AdjustSpeciesConc(j->first, (real)(j->second) / (NAvol));
}

inline void Process::ApplyToSystem(System &sys, const unsigned int n) const
{
    // Apply changes to gas-phase chemistry.
    map<unsigned int,int>::const_iterator j;
    real NAvol = NA * sys.SampleVolume();
    for (j=m_reac.begin(); j!=m_reac.end(); j++)
        sys.AdjustSpeciesConc(j->first, -(real)(j->second * n) / (NAvol));
    for (j=m_prod.begin(); j!=m_prod.end(); j++)
        sys.AdjustSpeciesConc(j->first, (real)(j->second * n) / (NAvol));
}
};

#endif
