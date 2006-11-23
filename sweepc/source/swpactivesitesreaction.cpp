#include "swpactivesitesreaction.h"

using namespace Sweep;

ActiveSitesReaction::ActiveSitesReaction(void)
{
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
    SurfaceReaction::Initialise(reac, prod, a, n, e, comp, values, pid, components);
    m_psitesfn = pfn;
}

Sweep::real ActiveSitesReaction::Rate(const real t, const System &sys) const
{
    // Get conditions from the system.
    vector<real> chem, sums;
    real T, P;
    sys.GetConditions(t, chem, T, P); sys.ConstEnsemble().GetSums(sums);

    return m_psitesfn(t, sys, chem, T, P, sums) * SurfaceReaction::Rate(t, chem, T, P, sums, sys);
}

real ActiveSitesReaction::Rate(const real t, const vector<real> &chem, const real T, 
                               const real P, const vector<real> &sums, const System &sys) const
{
    return m_psitesfn(t, sys, chem, T, P, sums) * SurfaceReaction::Rate(t, chem, T, P, sums, sys);
}

Sweep::real ActiveSitesReaction::Rate(const real t, const Sweep::System &sys, const DefaultParticle &sp) const
{
    // Get conditions from the system.
    vector<real> chem, sums;
    real T, P;
    sys.GetConditions(t, chem, T, P); sys.ConstEnsemble().GetSums(sums);

    return m_psitesfn(t, sys, chem, T, P, sums) * SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
}

Sweep::real ActiveSitesReaction::Rate(const real t, const vector<real> &chem, const real T, 
                                      const real P, const vector<real> &sums, const Sweep::System &sys, 
                                      const DefaultParticle &sp) const
{
    return m_psitesfn(t, sys, chem, T, P, sums) * SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
}