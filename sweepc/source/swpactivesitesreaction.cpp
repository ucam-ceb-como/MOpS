#include "swpactivesitesreaction.h"

using namespace Sweep;

ActiveSitesReaction::ActiveSitesReaction(void)
{
    // Initialise the active sites function pointer to NULL.
    m_psitesfn = NULL;
}

ActiveSitesReaction::~ActiveSitesReaction(void)
{
}

void ActiveSitesReaction::Initialise(const std::map<unsigned int,int> &reac, 
                                     const std::map<unsigned int,int> &prod, 
                                     const Sweep::real a, const Sweep::real n, 
                                     const Sweep::real e, const std::vector<Sweep::real> &comp, 
                                     const std::vector<Sweep::real> &values, const unsigned int pid, 
                                     std::vector<Component*> &components, 
                                     Sweep::ActiveSitesReaction::ActiveSitesFnPtr pfn)
{
    // Call the base class initialisation routine.
    SurfaceReaction::Initialise(reac, prod, a, n, e, comp, values, pid, components);
    // Set the pointer to the active sites function.
    m_psitesfn = pfn;
}

Sweep::real ActiveSitesReaction::Rate(const real t, const System &sys) const
{
    // Get conditions from the system.
    vector<real> chem, sums;
    real T, P;
    sys.GetConditions(t, chem, T, P); sys.ConstEnsemble().GetSums(sums);

    // Return the surface reaction rate multiplied by the actives sites concentration per
    // unit surface area.
    return m_psitesfn(t, sys, chem, T, P, sums) * SurfaceReaction::Rate(t, chem, T, P, sums, sys);
}

real ActiveSitesReaction::Rate(const real t, const vector<real> &chem, const real T, 
                               const real P, const vector<real> &sums, const System &sys) const
{
    // Return the surface reaction rate multiplied by the actives sites concentration per
    // unit surface area.
    return m_psitesfn(t, sys, chem, T, P, sums) * SurfaceReaction::Rate(t, chem, T, P, sums, sys);
}

Sweep::real ActiveSitesReaction::Rate(const real t, const Sweep::System &sys, const DefaultParticle &sp) const
{
    // Get conditions from the system.
    vector<real> chem, sums;
    real T, P;
    sys.GetConditions(t, chem, T, P); sys.ConstEnsemble().GetSums(sums);

    // Return the surface reaction rate multiplied by the actives sites concentration per
    // unit surface area.
    return m_psitesfn(t, sys, chem, T, P, sums) * SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
}

Sweep::real ActiveSitesReaction::Rate(const real t, const vector<real> &chem, const real T, 
                                      const real P, const vector<real> &sums, const Sweep::System &sys, 
                                      const DefaultParticle &sp) const
{
    // Return the surface reaction rate multiplied by the actives sites concentration per
    // unit surface area.
    return m_psitesfn(t, sys, chem, T, P, sums) * SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
}