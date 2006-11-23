#include "swpcoagulation.h"
#include "swpmechanism.h"
#include "swpensemble.h"

using namespace Sweep;

Coagulation::Coagulation(void)
{
    m_defer = false;
}

Coagulation::~Coagulation(void)
{
}

Sweep::real Coagulation::Rate(const Sweep::real t, const System &sys) const 
{
    real n = (real)sys.ParticleCount();

    if (n > 1.0) {
        real T = sys.GetTemperature(t);
        vector<real> sums; sys.ConstEnsemble().GetSums(sums);
        return Rate(sums, n, sqrt(T), T/Viscosity(T), T/sys.GetPressure(t), sys.SampleVolume());
    } else {
        return 0.0;
    }
}

real Coagulation::Rate(const real t, const vector<real> &chem, const real T, 
                       const real P, const vector<real> &sums, const System &sys) const
{
    real n = (real)sys.ParticleCount();

    if (n > 1.0) {
        return Rate(sums, n, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
    } else {
        return 0.0;
    }
}

void Coagulation::RateTerms(const real t, const System &sys, vector<real>::iterator &iterm) const
{
    real n = (real)sys.ParticleCount();

    if (n > 1.0) {
        real T = sys.GetTemperature(t);
        vector<real> sums; sys.ConstEnsemble().GetSums(sums);
        RateTerms(sums, n, sqrt(T), T/Viscosity(T), T/sys.GetPressure(t), sys.SampleVolume(), iterm);
    } else {
        // No coagulation as there are too few particles.
        for(int i=0; i<6; i++,iterm++) *iterm = 0.0;
    }
}

void Coagulation::RateTerms(const real t, const vector<real> &chem, const real T, 
                            const real P, const vector<real> &sums, const System &sys, 
                            vector<real>::iterator &iterm) const
{
    real n = (real)sys.ParticleCount();

    if (n > 1.0) {
        RateTerms(sums, n, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume(), iterm);
    } else {
        // No coagulation as there are too few particles.
        for(int i=0; i<6; i++,iterm++) *iterm = 0.0;
    }
}

int Coagulation::Perform(const Sweep::real t, Sweep::System &sys, const unsigned int iterm) const
{
    // Select property by which to choose particles (-1 means
    // choose uniformly.
    int rid1=-1, rid2=-1;
    MajorantType maj;
    switch (iterm) {
        case 0 :
            rid1 = -1;
            rid2 = -1;
            maj = SlipFlow;
            break;
        case 1 :
            rid1 = DefaultParticle::iD;
            rid2 = DefaultParticle::iD_1;
            maj = SlipFlow;
            break;
        case 2 :
            rid1 = -1;
            rid2 = DefaultParticle::iD_1;
            maj = SlipFlow;
            break;
        case 3 :
            rid1 = DefaultParticle::iD;
            rid2 = DefaultParticle::iD_2;
            maj = SlipFlow;
            break;
        case 4 :
            rid1 = -1;
            rid2 = DefaultParticle::iD2M_1_2;
            maj = FreeMol;
            break;
        case 5 :
            rid1 = DefaultParticle::iD2;
            rid2 = DefaultParticle::iM_1_2;
            maj = FreeMol;
            break;
        default :
            rid1 = -1;
            rid2 = -1;
            maj = SlipFlow;
            break;
    }

    // Save current T and P.
    real T = sys.GetTemperature(t);
    real P = sys.GetPressure(t);

    // Get first particle and update it.
    int ip1 = sys.Ensemble().SelectParticle(rid1);
    DefaultParticle *sp1, *sp1old;
    if (ip1 >= 0) {
        sp1 = sys.Ensemble().GetParticle(ip1);
        sp1old = &sp1->CreateCopy();
        m_mech->UpdateParticle(*sp1, sys, t);
    } else {
        // Failed to choose a particle.
        return -1;
    }

    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
        sys.Ensemble().RemoveParticle(ip1);
        ip1 = -1;
    }

    // Get unique second particle and update it.
    int ip2 = ip1, guard=0;
    DefaultParticle *sp2, *sp2old;
    while ((ip2 == ip1) && (guard<1000)) {
        guard++;
        ip2 = sys.Ensemble().SelectParticle(rid2);
    }
    if ((ip2>=0) && (ip2!=ip1)) {
        sp2 = sys.Ensemble().GetParticle(ip2);
        sp2old = &sp2->CreateCopy();
        m_mech->UpdateParticle(*sp2, sys, t);
    } else {
        // Failed to select a unique particle.
        delete sp1old;
        return -1;
    }

    // Check validity of particles after update.
    if (!sp2->IsValid()) {
        // Must remove second particle now.
        sys.Ensemble().RemoveParticle(ip2);
        ip2 = -1;
    }

    if ((ip1>=0) && (ip2>=0)) {
        // Must check for ficticious event now.
        real majk = CoagKernel(sp1old, sp2old, T, P, maj);
        real truek = CoagKernel(sp1, sp2, T, P, None);
        if (!Ficticious(majk,truek)) {
            // We can now coagulate the particles, remember to
            // remove second particle afterwards.
            *sp1 += *sp2;
            sp1->SetTime(t);
            sys.Ensemble().RemoveParticle(ip2);
        } else {
            delete sp1old; delete sp2old;
            return 1;
        }
    } else {
        // One or both particles were invalidated on update,
        // but that's not a problem.
    }

    delete sp1old; delete sp2old;
    return 0;
}

/* Coagulation kernels. */

Sweep::real Coagulation::CoagKernel(const Sweep::DefaultParticle *sp1, const Sweep::DefaultParticle *sp2, 
                             const Sweep::real T, const Sweep::real P, const Sweep::Coagulation::MajorantType maj) const
{
    real fm, sf;
    switch (maj) {
        case None :
            fm = FreeMolKernel(sp1, sp2, T, P, false);
            sf = SlipFlowKernel(sp1, sp2, T, P, false);
            return (fm*sf)/(fm+sf);
        case FreeMol :
            return FreeMolKernel(sp1, sp2, T, P, true);
        case SlipFlow :
            return SlipFlowKernel(sp1, sp2, T, P, true);
    }
    return 0.0;
}

Sweep::real Coagulation::FreeMolKernel(const Sweep::DefaultParticle *sp1, const Sweep::DefaultParticle *sp2, 
                                const Sweep::real T, const Sweep::real P, const bool maj) const
{
    if (maj) {
        return CFMMAJ * CFM * sqrt(T) * 
              (sp1->InvSqrtMass()+sp2->InvSqrtMass()) *
              (sp1->CollDiamSquared()+sp2->CollDiamSquared());
    } else {
        return CFM * sqrt(T * 
               ((1.0/sp1->Mass())+(1.0/sp2->Mass()))) * 
               pow(sp1->CollisionDiameter()+sp2->CollisionDiameter(),2.0);
    }
}

Sweep::real Coagulation::SlipFlowKernel(const Sweep::DefaultParticle *sp1, const Sweep::DefaultParticle *sp2, 
                                 const Sweep::real T, const Sweep::real P, const bool maj) const
{
    return ((1.257 * KNUDSEN_K * T * 
           (sp1->InvCollDiamSquared()+sp2->InvCollDiamSquared()) / P) +
           (sp1->InvCollDiam() + sp2->InvCollDiam())) * CSF * T * 
           (sp1->CollisionDiameter()+sp2->CollisionDiameter()) / Viscosity(T);
}