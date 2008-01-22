#include "swppahreaction.h"
#include "swppahparticle.h"
#include "swpmechanism.h"
#include "swppahmodel.h"

using namespace Sweep;

PAHReaction::PAHReaction(void)
{
    m_mode = None;
}

PAHReaction::~PAHReaction(void)
{
}

void PAHReaction::Initialise(const std::map<unsigned int,int> &reac, 
                             const std::map<unsigned int,int> &prod, 
                             const Sweep::real a, const Sweep::real n, 
                             const Sweep::real e, const vector<Sweep::real> &comp, 
                             const vector<Sweep::real> &values, const unsigned int pid,
                             vector<Component*> &components, const PAHReactionMode mode)
{
    SurfaceReaction::Initialise(reac, prod, a, n, e, comp, values, pid, components);
    m_mode = mode;
}

real PAHReaction::Rate(const real t, const System &sys) const
{
    // Get conditions from the system.
    vector<real> chem, sums;
    real T, P;
    sys.GetConditions(t, chem, T, P); sys.ConstEnsemble().GetSums(sums);

    // Return the reaction rate.
    return PAHReaction::Rate(t, chem, T, P, sums, sys);
}

real PAHReaction::Rate(const real t, const vector<real> &chem, const real T, 
                       const real P, const vector<real> &sums, const System &sys) const
{
    // The reaction rate depends on the reaction mode.  This is because different
    // sites react at different rates, hence have different radical site fractions.
    switch (m_mode) {
        case EdgeGrowth :
            return PAHModel::EdgeRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys);
        case ArmchairGrowth :
            return PAHModel::ArmchairRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys);
        case R5Addition :
            return PAHModel::ZigZagRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys);
        case R5EdgeConv :
            return PAHModel::EdgeRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys);
        case R5ArmchairConv :
            return PAHModel::ArmchairRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys);
        case EdgeOxidation :
            return PAHModel::EdgeRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys);
        default:
            return SurfaceReaction::Rate(t, chem, T, P, sums, sys);
    }
}

real PAHReaction::Rate(const real t, const System &sys, const DefaultParticle &sp) const
{
    // Get conditions from the system.
    vector<real> chem, sums;
    real T, P;
    sys.GetConditions(t, chem, T, P); sys.ConstEnsemble().GetSums(sums);

    // Return the surface reaction rate multiplied by the actives sites concentration per
    // unit surface area.
    return PAHReaction::Rate(t, chem, T, P, sums, sys, sp);
}

Sweep::real PAHReaction::Rate(const real t, const vector<real> &chem, const real T, 
                              const real P, const vector<real> &sums, const Sweep::System &sys, 
                              const DefaultParticle &sp) const
{
    // The reaction rate depends on the reaction mode.  This is because different
    // sites react at different rates, hence have different radical site fractions.
    switch (m_mode) {
        case EdgeGrowth :
            return PAHModel::EdgeRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
        case ArmchairGrowth :
            return PAHModel::ArmchairRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
        case R5Addition :
            return PAHModel::ZigZagRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
        case R5EdgeConv :
            return PAHModel::EdgeRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
        case R5ArmchairConv :
            return PAHModel::ArmchairRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
        case EdgeOxidation :
            return PAHModel::EdgeRadicalFraction(chem, T, P) * 
                   SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
        default:
            return SurfaceReaction::Rate(t, chem, T, P, sums, sys, sp);
    }
}

int PAHReaction::Perform(const real t, System &sys, const unsigned int iterm) const
{
    // Select a particle from the ensemble
    int i = sys.Ensemble().SelectParticle(m_pid);

    if (i >= 0) {
        PAHParticle *sp = (PAHParticle*)sys.Ensemble().GetParticle(i);
        real majr = MajorantRate(t, sys, *sp);

        // Update particle with deferred processes.
        if (m_mech->AnyDeferred()) {
            m_mech->UpdateParticle(*sp, sys, t);
        }

        // Check that the particle is still valid.
        if (sp->IsValid()) {
            real truer = Rate(t, sys, *sp);

            if (!Ficticious(majr, truer)) {
                
                // Update the particle using the PAH model.
                unsigned int n = 1;
                switch (m_mode) {
                    case EdgeGrowth :
                        sp->PerformEdgeRingGrowth(n);
                    case ArmchairGrowth :
                        sp->PerformArmchairRingGrowth(n);
                    case R5Addition :
                        sp->PerformR5Addition(n);
                    case R5EdgeConv :
                        sp->PerformR5EdgeConversion(n);
                    case R5ArmchairConv :
                        sp->PerformR5ArmchairConversion(n);
                    case EdgeOxidation :
                        sp->PerformEdgeOxidation(n);
                }

                // n will be set to zero if the particle had insufficient required
                // sites to perform the process.  In this case we do not update the
                // particle further.
                if (n>0) {
                    // Adjust particle and update ensemble.
                    AdjustParticle(*sp, t, 1);
                    sys.Ensemble().Update(i);

                    // Apply changes to gas-phase chemistry.
                    ApplyToSystem(sys);
                }
            }
        } else {
            // If not valid then remove the particle.
            sys.Ensemble().RemoveParticle(i);
        }
    } else {
        // Failed to select a particle.
        return -1;
    }

    return 0;
}

int PAHReaction::Perform(const real t, System &sys, DefaultParticle &sp, 
                         const unsigned int n) const
{
    // Update the particle using the PAH model.
    unsigned int m = n;
    switch (m_mode) {
        case EdgeGrowth :
            ((PAHParticle&)sp).PerformEdgeRingGrowth(m);
        case ArmchairGrowth :
            ((PAHParticle&)sp).PerformArmchairRingGrowth(m);
        case R5Addition :
            ((PAHParticle&)sp).PerformR5Addition(m);
        case R5EdgeConv :
            ((PAHParticle&)sp).PerformR5EdgeConversion(m);
        case R5ArmchairConv :
            ((PAHParticle&)sp).PerformR5ArmchairConversion(m);
        case EdgeOxidation :
            ((PAHParticle&)sp).PerformEdgeOxidation(m);
    }

    // m will be set to zero if the particle had insufficient required
    // sites to perform the process.  In this case we do not update the
    // particle further.
    if (m>0) {
        // Adjust particle.
        AdjustParticle(sp, t, m);

        // Apply changes to gas-phase chemistry.
        ApplyToSystem(sys, m);
    }

    return 0;
}

