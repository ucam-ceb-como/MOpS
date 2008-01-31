/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Defines a condensation reaction.  Condensation reactions are treated
    differently from surface reactions as they are modelled as free molecular
    collisions.  The assumptions made in this model are:

    1.  The colliding species is much smaller than the recipient particle.

    Before the condensation process can be used it must be provided the mass and
    diameter of the condensing species in order to calculate the rate terms.  The
    SetCondensingSpecies() function is used for this.
*/

#ifndef SWEEP_CONDENSATION_H
#define SWEEP_CONDENSATION_H

#include "swpprocess.h"
#include "swpparticlechanger.h"

namespace Sweep
{
class Condensation : public Process, public ParticleChanger
{
protected:
    const static real COND_MAJ; // See def'n in cpp file.
public:
    static const unsigned int TERM_COUNT = 3;
protected:
    real m_a; // Rate constant.
    real m_kfm1, m_kfm2, m_kfm3; // Free-mol term parameters.
public: // Default constructor and destructor.
    Condensation(void);
    ~Condensation(void);
public:
    /* Initialises the condensation reaction. */
    virtual void Initialise(const map<unsigned int, int> &reac, // Gas-phase reactants.
                            const map<unsigned int, int> &prod, // Gas-phase products.
                            const real a,                       // Rate constant.
                            const real m,                       // Mass of condensating species (g).
                            const real d,                       // Diameter of condensating species (cm).
                            const vector<real> &comp,           // Component counts of new particle.
                            const vector<real> &values,         // Other values for new particle.
                            vector<Component*> &components);    // Reference to component vector used to define process.
    /* Clears memory associated with the reaction. */
    virtual void Destroy(void);
public: // Rate calculation.
    /* Returns the rate of the process for the given system. */
    virtual real Rate(const real t, const System &sys) const;
    /* Calculates the process rate given the chemical conditions. */
    virtual real Rate(const real t, const vector<real> &chem, const real T, 
                      const real P, const vector<real> &sums, const System &sys) const;
    /* Returns the rate of the process for the given particle in
       the system. Process must be linear in particle number. */
    virtual real Rate(const real t,             // Current time (s).
                      const System &sys,        // System for which to get rate.
                      const Particle &sp // Index of particle for which to calculate rate.
                     ) const;
    /* Returns the rate of the process for the given particle in
       the system with precalculated chemical conditions. Process
       must be linear in particle number. */
    virtual real Rate(const real t,             // Current time (s).
                      const vector<real> &chem,
                      const real T, 
                      const real P,
                      const vector<real> &sums,
                      const System &sys,        // System for which to get rate.
                      const Particle &sp // Index of particle for which to calculate rate.
                     ) const;
    real MajorantRate(const real t, const System &sys, const Particle &sp) const;
    real MajorantRate(const real t, const vector<real> &chem, const real T, const real P, 
                      const vector<real> &sums, const System &sys, const Particle &sp) const;
public: // Rate term calculation.
    /* Returns the number of rate terms for this process. */
    inline unsigned int TermCount(void) const {return TERM_COUNT;};
    /* Calculates the rate terms give an iterator to a real vector.  The 
       iterator is advanced to the postition after the last term for this
       process. */
    void RateTerms(const real t, const System &sys, vector<real>::iterator &iterm) const;
    /* Calculates the rate terms give an iterator to a real vector.  The 
       iterator is advanced to the postition after the last term for this
       process. Chemical conditions are precalculated. */
    void RateTerms(const real t, const vector<real> &chem, const real T, 
                   const real P, const vector<real> &sums, const System &sys, 
                   vector<real>::iterator &iterm) const;
public: // Doing the process.
    /* Performs the process on the given system.  The responsible rate term is given
       by index.  Returns 0 on success, otherwise negative. */
    virtual int Perform(const real t, System &sys, const unsigned int iterm) const;
    /* Performs the process on a given particle in the system.  Particle
       is given by index.  The process is performed n times. */
    virtual int Perform(const real t,         // Current time (s).
                        System &sys,          // System for which to perform process. 
                        Particle &sp,  // Index of particle for which to perform process.
                        const unsigned int n  // Number of times to perform the process.
                       ) const;
public: // Property gets/sets.
    /* Returns the pre-exp rate constant. */
    inline real A() const {return m_a / NA;};
    /* Sets the pre-exp rate constant. */
    inline void SetA(const real a) {m_a = a * NA;};
public: // Property sets.
    /* Sets the coagulation kernel constants given incepting species'
       mass and diameter. */
    void SetCondensingSpecies(const real m, const real d);
};
};

#endif