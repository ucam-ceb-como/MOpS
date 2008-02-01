/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a PAH surface reaction process.  This is a specialised
    process specifically written to simulate to growth of PAH structures
    in soot particles.  It requires certain elements to have been defined
    in the mechanism in order to work.
*/

#ifndef SWEEP_PAHRXN_H
#define SWEEP_PAHRXN_H

#include "swpsurfacereaction.h"

using namespace std;

namespace Sweep
{
class PAHReaction : public Sweep::SurfaceReaction
{
public:
    /* Enumeration of different possible modes for a PAH reaction.  Modes dictate
       the process rate calculation and how the process is performed. */
    enum PAHReactionMode {EdgeGrowth, ArmchairGrowth, R5Addition, R5EdgeConv, 
                          R5ArmchairConv, EdgeOxidation, None=-1};
protected:
    PAHReactionMode m_mode; // The PAH reaction mode used by this reaction.
public: // Default constructor and destructor.
    PAHReaction(void);
    ~PAHReaction(void);
public:
    /* Initialises the PAH surface reaction. */
    virtual void Initialise(const map<unsigned int, int> &reac, // Gas-phase reactants.
                            const map<unsigned int, int> &prod, // Gas-phase products.
                            const real a, const real n,         // Arrhenius rate parameters.
                            const real e,                       // ''.
                            const vector<real> &comp,           // Component counts of new particle.
                            const vector<real> &values,         // Other values for new particle.
                            const unsigned int pid,             // Index of particle property.
                            vector<Component*> &components,     // Reference to component vector used to define process.
                            const PAHReactionMode mode);        // Reaction mode for this process.
public:
    /* Sets the PAH reaction mode. */
    inline void SetMode(const PAHReactionMode mode) {m_mode = mode;};
    /* Returns the PAH reaction mode. */
    inline PAHReactionMode Mode(void) {return m_mode;};
public:
    /* Returns the rate of the process for the given system. */
    virtual real Rate(const real t, const System &sys) const;
    /* Calculates the process rate given the chemical conditions. */
    virtual real Rate(const real t, const vector<real> &chem, const real T, 
                      const real P, const vector<real> &sums, const System &sys) const;
    /* Returns the rate of the process for the given particle in
       the system. Process must be linear in particle number. */
    virtual real Rate(const real t, const System &sys, const DefaultParticle &sp) const;
    virtual real Rate(const real t, const vector<real> &chem, const real T, 
                      const real P, const vector<real> &sums, const System &sys, 
                      const DefaultParticle &sp) const;
public:
    /* Performs the process on the given system.  The responsible rate term is given
       by index.  Returns 0 on success, otherwise negative. */
    int Perform(const real t, System &sys, const unsigned int iterm) const;
    /* Performs the process on a given particle in the system.  Particle
       is given by index.  The process is performed n times. */
    int Perform(const real t, System &sys, DefaultParticle &sp, const unsigned int n) const;
};
}
#endif
