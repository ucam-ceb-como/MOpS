#include "swpcoagulation.h"
#include "swpmechanism.h"
#include "swpensemble.h"

using namespace Sweep;

Coagulation::Coagulation(void)
{
    // The coagulation process cannot be deferred (LPDA) because it is a 2 particle
    // process.
    m_defer = false;
}

Coagulation::~Coagulation(void)
{
	// Call the abstract base class destructor.
}

Sweep::real Coagulation::Rate(const Sweep::real t, const System &sys) const 
{
    // Get the number of particles in the system.
    real n = (real)sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1.0) {
        // Get system properties required to calculate coagulation rate.
        real T = sys.GetTemperature(t);
        vector<real> sums; sys.ConstEnsemble().GetSums(sums);

        // Calculate the rate.
        return Rate(sums, n, sqrt(T), T/Viscosity(T), T/sys.GetPressure(t), sys.SampleVolume());
    } else {
        return 0.0;
    }
}

real Coagulation::Rate(const real t, const vector<real> &chem, const real T, 
                       const real P, const vector<real> &sums, const System &sys) const
{
    // Get the number of particles in the system.
    real n = (real)sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1.0) {
        // All system properties required are passed as arguments, so just
        // calculate the rate.
        return Rate(sums, n, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
    } else {
        return 0.0;
    }
}

void Coagulation::RateTerms(const real t, const System &sys, vector<real>::iterator &iterm) const
{
    // Get the number of particles in the system.
    real n = (real)sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1.0) {
        // Get system properties required to calculate coagulation rate.
        real T = sys.GetTemperature(t);
        vector<real> sums; sys.ConstEnsemble().GetSums(sums);

        // Calculate the rate terms.
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
    // Get the number of particles in the system.
    real n = (real)sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1.0) {
        // All system properties required are passed as arguments, so just
        // calculate the rate terms.
        RateTerms(sums, n, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume(), iterm);
    } else {
        // No coagulation as there are too few particles.
        for(int i=0; i<6; i++,iterm++) *iterm = 0.0;
    }
}

int Coagulation::Perform(const Sweep::real t, Sweep::System &sys, const unsigned int iterm) const
{
    // Select properties by which to choose particles (-1 means
    // choose uniformly).  Note we need to choose 2 particles.  There
    // are six possible rate terms to choose from; 4 slip-flow and 2
    // free molecular.
    int rid1=-1, rid2=-1;
    MajorantType maj;
    switch (iterm) {
        case 0 :
            rid1 = -1;
            rid2 = -1;
            maj = SlipFlow;
            break;
        case 1 :
            rid1 = Particle::iD;
            rid2 = Particle::iD_1;
            maj = SlipFlow;
            break;
        case 2 :
            rid1 = -1;
            rid2 = Particle::iD_1;
            maj = SlipFlow;
            break;
        case 3 :
            rid1 = Particle::iD;
            rid2 = Particle::iD_2;
            maj = SlipFlow;
            break;
        case 4 :
            rid1 = -1;
            rid2 = Particle::iD2M_1_2;
            maj = FreeMol;
            break;
        case 5 :
            rid1 = Particle::iD2;
            rid2 = Particle::iM_1_2;
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

    // Choose and get first particle, then update it.
    int ip1 = sys.Ensemble().SelectParticle(rid1);
    Particle *sp1, *sp1old;
    if (ip1 >= 0) {
        sp1 = sys.Ensemble().GetParticle(ip1);

        // Create a copy of the particle before updating.
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

        // Invalidating the index tells this routine not to perform coagulation.
        ip1 = -1;
    }

    // Choose and get unique second particle, then update it.  Note, we are allowed to do
    // this even if the first particle was invalidated.
    int ip2 = ip1, guard=0;
    Particle *sp2, *sp2old;
    while ((ip2 == ip1) && (guard<1000)) {
        guard++;
        ip2 = sys.Ensemble().SelectParticle(rid2);
    }
    if ((ip2>=0) && (ip2!=ip1)) {
        sp2 = sys.Ensemble().GetParticle(ip2);

        // Create a copy of the particle before updating.
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

        // Invalidating the index tells this routine not to perform coagulation.
        ip2 = -1;
    }

    // Check that both the particles are still valid.
    if ((ip1>=0) && (ip2>=0)) {
        // Must check for ficticious event now by calculate the original
        // majorant rate and the current (after updates) true rate.
        real majk = CoagKernel(sp1old, sp2old, T, P, maj);
        real truek = CoagKernel(sp1, sp2, T, P, None);
        if (!Ficticious(majk, truek)) {
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

Sweep::real Coagulation::CoagKernel(const Sweep::Particle *sp1, const Sweep::Particle *sp2, 
                             const Sweep::real T, const Sweep::real P, const Sweep::Coagulation::MajorantType maj) const
{
    // This routine calculates the coagulation kernel for two particles.  The kernel
    // type is chosen by the majorant type requested.

    real fm, sf;
    switch (maj) {
        case None :
            // No majorant, so return half harmonic mean of free molecular and
            // slip-flow kernels (non-majorant).
            fm = FreeMolKernel(sp1, sp2, T, P, false);
            sf = SlipFlowKernel(sp1, sp2, T, P, false);
            return (fm*sf)/(fm+sf);
        case FreeMol :
            // Free molecular majorant.
            return FreeMolKernel(sp1, sp2, T, P, true);
        case SlipFlow :
            // Slip-flow majorant.
            return SlipFlowKernel(sp1, sp2, T, P, true);
    }

    // Invalid majorant, return zero.
    return 0.0;
}

Sweep::real Coagulation::FreeMolKernel(const Sweep::Particle *sp1, const Sweep::Particle *sp2, 
                                const Sweep::real T, const Sweep::real P, const bool maj) const
{
    // This routine calculate the free molecular coagulation kernel for two particles.
    // There are two forms of kernel; a majorant form and a non-majorant form.

    if (maj) {
        // The majorant form is always >= the non-majorant form.
        return CFMMAJ * CFM * sqrt(T) * 
              (sp1->InvSqrtMass()+sp2->InvSqrtMass()) *
              (sp1->CollDiamSquared()+sp2->CollDiamSquared());
    } else {
        return CFM * sqrt(T * 
               ((1.0/sp1->Mass())+(1.0/sp2->Mass()))) * 
               pow(sp1->CollisionDiameter()+sp2->CollisionDiameter(),2.0);
    }
}

Sweep::real Coagulation::SlipFlowKernel(const Sweep::Particle *sp1, const Sweep::Particle *sp2, 
                                 const Sweep::real T, const Sweep::real P, const bool maj) const
{
    // For the slip-flow kernel the majorant and non-majorant forms are identical.
    return ((1.257 * KNUDSEN_K * T * 
           (sp1->InvCollDiamSquared()+sp2->InvCollDiamSquared()) / P) +
           (sp1->InvCollDiam() + sp2->InvCollDiam())) * CSF * T * 
           (sp1->CollisionDiameter()+sp2->CollisionDiameter()) / Viscosity(T);
}