/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a surface reaction process.  Surface reactions are single
    particle processes, hence are deferred by default.
*/

#ifndef SWEEP_SURFRXN_H
#define SWEEP_SURFRXN_H

#include "swpparams.h"
#include "swpprocess.h"
#include "swpcomponent.h"
#include "swpparticlechanger.h"

using namespace std;

namespace Sweep
{
class SurfaceReaction : public Process, public ParticleChanger
{
protected:
    const static real SURF_MAJ; // See def'n below class.
protected:
    real m_a, m_n, m_e; // Arrhenius rate parameters.
    unsigned int m_pid;  // ID of particle property on which this reaction depends.
public:
    SurfaceReaction(void);
    virtual ~SurfaceReaction(void);
    /* Initialises the surface reaction reaction. */
    virtual void Initialise(const map<unsigned int, int> &reac, // Gas-phase reactants.
                            const map<unsigned int, int> &prod, // Gas-phase products.
                            const real a, const real n,         // Arrhenius rate parameters.
                            const real e,                       // ''.
                            const vector<real> &comp,           // Component counts of new particle.
                            const vector<real> &values,         // Other values for new particle.
                            const unsigned int pid,             // Index of particle property.
                            vector<Component*> &components);     // Reference to component vector used to define process.
    /* Clears memory associated with the reaction. */
    virtual void Destroy(void);
public:
    /* Returns the rate of the process for the given system. */
    virtual real Rate(const real t, const System &sys) const;
    /* Calculates the process rate given the chemical conditions. */
    virtual real Rate(const real t, const vector<real> &chem, const real T, 
                      const real P, const vector<real> &sums, const System &sys) const;
    /* Returns the rate of the process for the given particle in
       the system. Process must be linear in particle number. */
    virtual real Rate(const real t, const System &sys, const Particle &sp) const;
    virtual real Rate(const real t, const vector<real> &chem, const real T, 
                      const real P, const vector<real> &sums, const System &sys, 
                      const Particle &sp) const;
    real MajorantRate(const real t, const System &sys, const Particle &sp) const;
    real MajorantRate(const real t, const vector<real> &chem, const real T, 
                      const real P, const vector<real> &sums, const System &sys, 
                      const Particle &sp) const;
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
    /* Performs the process on a given particle in the system.  Particle
       is given by index.  The process is performed n times. */
    int Perform(const real t, System &sys, Particle &sp, const unsigned int n) const;
public: // Property gets.
    inline real A() const {return m_a;};
    inline real n() const {return m_n;};
    inline real E() const {return m_e;};
    inline unsigned int PropertyID(void) {return m_pid;};
public: // Property sets.
    /* Sets the rate constant. */
    inline void SetArrhenius(const real a, const real n, const real e) {
        m_a=a; m_n=n; m_e=e;
    }
    /* Sets the particle property index. */
    inline virtual void SetPropertyID(const unsigned int i) {m_pid = i;};
};
};

#endif