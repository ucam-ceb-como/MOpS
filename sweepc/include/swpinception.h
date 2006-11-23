/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of an inception process.  Inceptions are modelled as the coagulation
    of 2 gas-phase species.  The rate is given by the collision kernel of these species,
    therefore one must provide their masses and diameters.  Kernel parameters are
    calculated internally.
*/

#ifndef SWEEP_INCEPTION_H
#define SWEEP_INCEPTION_H

#include "swpparams.h"
#include "swpprocess.h"
#include "swpcomponent.h"
#include "swpparticlechanger.h"
#include <map>

using namespace std;

namespace Sweep
{
class Inception : public Process, public ParticleChanger
{
protected:
    real m_a;                   // Rate constant.
    real m_kfm, m_ksf1, m_ksf2; // Free-mol and slip-flow kernel parameters.

public:
    Inception(void);
    ~Inception(void);
    /* Initialises the inception reaction. */
    void Initialise(const map<unsigned int, int> &reac, // Gas-phase reactants.
                    const map<unsigned int, int> &prod, // Gas-phase products.
                    const real a,                       // Rate constant.
                    const real m1, const real m2,       // Masses of two incepting species (g).
                    const real d1, const real d2,       // Diameters of two incepting species (cm).
                    const vector<real> &comp,           // Component counts of new particle.
                    const vector<real> &values,         // Other values for new particle.
                    vector<Component*> &components);     // Reference to component vector used to define process.
public: // Rate calculation.
    /* Returns the rate of the process for the given system. */
    real Rate(const real t, const System &sys) const;
    /* Calculates the process rate given the chemical conditions explicitly. */
    real Rate(const real t, const vector<real> &chem, const real T, 
              const real P, const vector<real> &sums, const System &sys) const;
    /* A faster rate calculation routine for Inception events only.  Requires all the
       parameters that would otherwise be calculated by the routine to be passed as
       arguments, */
    real Rate(const vector<real> &chem, const real sqrtT, const real T_mu, 
              const real T_P, const real vol) const;
public: // Single particle rates (invalid for inceptions, so return invalid negative value).
    /* Returns the rate of the process for the given particle in
       the system. Process must be linear in particle number. */
    inline real Rate(const real t, const System &sys, const DefaultParticle &sp) const {return -1.0;};
    inline real Rate(const real t, const vector<real> &chem, const real T, const real P, 
                     const vector<real> &sums, const System &sys, const DefaultParticle &sp) const {
        return -1.0;
    };
public: // Rate term calculation.
    /* Returns the number of rate terms for this process. */
    inline unsigned int TermCount(void) const {return 1;};
    /* Calculates the rate terms give an iterator to a real vector.  The 
       iterator is advanced to the postition after the last term for this
       process. */
    void RateTerms(const real t, const System &sys, vector<real>::iterator &iterm) const;
    void RateTerms(const real t, const vector<real> &chem, const real T, 
                   const real P, const vector<real> &sums, const System &sys,
                   vector<real>::iterator &iterm) const;
public:
    /* Performs the process on the given system.  The responsible rate term is given
       by index.  Returns 0 on success, otherwise negative. */
    int Perform(const real t, System &sys, const unsigned int iterm) const;
    /* Particle Perform() routine is invalid for an inception. */
    inline int Perform(const real t, System &sys, DefaultParticle &sp, 
                       const unsigned int n) const {return -1;};
public: // Property gets.
    inline real A() const {return m_a;};
public: // Property sets.
    inline void SetDeferred(const bool defer) {m_defer=false;}
    /* Sets the rate constant. */
    inline void SetA(const real a) {m_a=0.5*a;};
    /* Sets the coagulation kernel constants given incepting species
       masses and diameters. */
    void SetInceptingSpecies(real m1, real m2, real d1, real d2);
};

inline real Inception::Rate(const vector<real> &chem, const real sqrtT, const real T_mu, 
                            const real T_P, const real vol) const
{
    // Temperature and pressure dependence.
    real fm = sqrtT * m_kfm;
    real sf = T_mu * (m_ksf1 + (T_P*m_ksf2));
    real rate = m_a * (fm*sf) / (fm+sf) * vol;

    // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        rate *= pow(NA*chem[(*i).first], (real)(*i).second);
    }

    return rate;
}

};

#endif