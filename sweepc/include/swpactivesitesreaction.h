/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a surface reaction which includes terms for active sites.  The
    active sites concentration is provided by a function pointer, which must be
    provided before the process is used.
*/

#ifndef SWEEP_ACTIVESITESRXN_H
#define SWEEP_ACTIVESITESRXN_H

#include "swpparams.h"
#include "swpcomponent.h"
#include "swpsurfacereaction.h"
#include "swpsystem.h"
#include <map>

using namespace std;

namespace Sweep
{
class ActiveSitesReaction : public Sweep::SurfaceReaction
{
public:
    /* Pointer to the function which returns the number of active
       sites per unit surface area of particles. */
    typedef real (*ActiveSitesFnPtr) (const real t, const System &sys, 
                                      const vector<real> &chem, const real T, 
                                      const real P, const vector<real> &sums);
protected:
    ActiveSitesFnPtr m_psitesfn;
public:
    ActiveSitesReaction(void);
    ~ActiveSitesReaction(void);
    /* Initialises the inception reaction. */
    void Initialise(const map<unsigned int, int> &reac, // Gas-phase reactants.
                    const map<unsigned int, int> &prod, // Gas-phase products.
                    const real a, const real n,         // Arrhenius rate parameters.
                    const real e,                       // ''.
                    const vector<real> &comp,           // Component counts of new particle.
                    const vector<real> &values,         // Other values for new particle.
                    const unsigned int pid,             // Index of particle property.
                    vector<Component*> &components,     // Reference to component vector used to define process.
                    ActiveSitesFnPtr pfn);              // Pointer to the active sites function.
    /* Sets the function used to calculate active sites density. */
    inline void SetActiveSitesFn(ActiveSitesFnPtr pfn) {m_psitesfn = pfn;};
    /* Returns the rate of the process for the given system. */
    real Rate(const real t, const System &sys) const;
    /* Calculates the process rate given the chemical conditions. */
    real Rate(const real t, const vector<real> &chem, const real T, 
              const real P, const vector<real> &sums, const System &sys) const;
    /* Returns the rate of the process for the given particle in
       the system. Process must be linear in particle number. */
    real Rate(const real t, const System &sys, const DefaultParticle &sp) const;
    real Rate(const real t, const vector<real> &chem, const real T, 
              const real P, const vector<real> &sums, const System &sys, 
              const DefaultParticle &sp) const;
};
};

#endif