#include "swpcondensation.h"
#include <cmath>
#include "swpmechanism.h"

using namespace Sweep;

const Sweep::real Condensation::COND_MAJ = 2.0;

Condensation::Condensation(void)
{
    // Clear the current data from the condensation.
    Destroy();

    // Assume the condensation is simulated as a deferred process (LPDA).
    m_defer = true;
}

Condensation::~Condensation(void)
{
    Destroy();
}

void Condensation::Initialise(const map<unsigned int, int> &reac,const map<unsigned int, int> &prod,
                              const real a, const real m, const real d, const vector<real> &comp,           
                              const vector<real> &values, vector<Component*> &components)
{
    // Destroy the current data in the condensation.
    Destroy();

    // Set up the condensation variables.
    m_reac.insert(reac.begin(),reac.end());
    m_prod.insert(reac.begin(),reac.end());
    m_a = a * NA;
    SetCondensingSpecies(m, d);
    m_comp.assign(comp.begin(), comp.end());
    m_values.assign(values.begin(), values.end());
    m_components = &components;
}

void Condensation::Destroy()
{
    // Clears all memory and resets all variables associated with
    // the condensation.
    m_reac.clear();
    m_prod.clear();
    m_a = NA;
    m_kfm1 = 0.0;
    m_kfm2 = 0.0;
    m_kfm3 = 0.0;
    m_comp.clear();
    m_values.clear();
    m_components = NULL;
}

real Condensation::Rate(const real t, const System &sys) const
{
    // Get required system conditions for rate calculation.
    vector<real> chem; real T, P;
    sys.GetConditions(t, chem, T, P);
    vector<real> sums; sys.ConstEnsemble().GetSums(sums);

    // Calculate the rate.
    return Rate(t, chem, T, P, sums, sys);
}

real Condensation::Rate(const real t, const vector<real> &chem, const real T, 
                        const real P, const vector<real> &sums, const System &sys) const
{
    // Calculate temperature terms.
    real sqrtT = sqrt(T);
    real cterm = m_a * sqrtT;

     // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        cterm *= pow(chem[i->first], (real)i->second);
    }

    // Free molecular terms.
    cterm *= (m_kfm1 * sys.ParticleCount()) + (m_kfm2 * sums[Particle::iD]) +
             (m_kfm3 * sums[Particle::iD2]);

    // If the mechanism contains any deferred processes then we must use the
    // majorant form of the rate, in order to account for any changes to
    // particles during deferred process updates.
    if (m_mech->AnyDeferred()) {
        return cterm * COND_MAJ;
    } else {
        return cterm;
    }
}

real Condensation::Rate(const real t, const System &sys, const Particle &sp) const
{
    // Calculate temperature terms.
    real sqrtT = sqrt(sys.GetTemperature(t));
    real cterm = m_a * sqrtT, trm[3];

    // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator j;
    for (j=m_reac.begin(); j!=m_reac.end(); j++) {
        cterm *= pow(sys.GetSpeciesConc((*j).first,t), (real)j->second);
    }

    // Free molecular terms.
    trm[0] = cterm * m_kfm1;
    trm[1] = cterm * (m_kfm2 * sp.CollisionDiameter());
    trm[2] = cterm * (m_kfm3 * sp.CollDiamSquared());
    cterm *= m_kfm1 + (m_kfm2 * sp.CollisionDiameter()) +
             (m_kfm3 * sp.CollDiamSquared());
    return trm[0] + trm[1] + trm[2];
}

real Condensation::Rate(const real t, const vector<real> &chem, const real T, 
                        const real P, const vector<real> &sums, const System &sys, 
                        const Particle &sp) const
{
    // Calculate temperature terms.
    real sqrtT = sqrt(sys.GetTemperature(t));
    real cterm = m_a * sqrtT;

     // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator j;
    for (j=m_reac.begin(); j!=m_reac.end(); j++) {
        cterm *= pow(sys.GetSpeciesConc((*j).first,t), (real)j->second);
    }

    // Free molecular terms.
    cterm *= m_kfm1 + (m_kfm2 * sp.CollisionDiameter()) +
             (m_kfm3 * sp.CollDiamSquared());
    return cterm;
}

Sweep::real Condensation::MajorantRate(const real t, const System &sys, 
                                       const Particle &sp) const
{
    // Return the single particle rate multiplied by the 
    // condensation majorant factor.
    return Rate(t,sys,sp) * COND_MAJ;
}

real Condensation::MajorantRate(const real t, const vector<real> &chem, const real T, const real P,
                                const vector<real> &sums, const System &sys, const Particle &sp) const
{
    // Return the single particle rate multiplied by the 
    // condensation majorant factor.
    return Rate(t, chem, T, P, sums, sys, sp) * COND_MAJ;
}

void Condensation::RateTerms(const real t, const System &sys, vector<real>::iterator &iterm) const
{
    // Calculate temperature terms.
    real sqrtT = sqrt(sys.GetTemperature(t));
    real cterm = m_a * sqrtT;

     // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        cterm *= pow(sys.GetSpeciesConc((*i).first,t), (real)i->second);
    }

    // If the mechanism contains any deferred processes then we must use the
    // majorant form of the rate, in order to account for any changes to
    // particles during deferred process updates.
    if (m_mech->AnyDeferred()) cterm *= COND_MAJ;

    // Free molecular terms.
    *iterm++ = m_kfm1 * cterm * sys.ParticleCount();
    *iterm++ = m_kfm2 * cterm * sys.ConstEnsemble().GetSum(Particle::iD);
    *iterm++ = m_kfm3 * cterm * sys.ConstEnsemble().GetSum(Particle::iD2);
}

void Condensation::RateTerms(const Sweep::real t, const std::vector<real> &chem, 
                             const Sweep::real T, const Sweep::real P, 
                             const vector<real> &sums, 
                             const Sweep::System &sys, vector<real>::iterator &iterm) const
{
    // Calculate temperature terms.
    real sqrtT = sqrt(T);
    real cterm = m_a * sqrtT;

     // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        cterm *= pow(chem[i->first], (real)i->second);
    }

    // If the mechanism contains any deferred processes then we must use the
    // majorant form of the rate, in order to account for any changes to
    // particles during deferred process updates.
    if (m_mech->AnyDeferred()) cterm *= COND_MAJ;

    // Free molecular terms.
    *iterm = m_kfm1 * cterm * sys.ParticleCount(); iterm++;
    *iterm = m_kfm2 * cterm * sums[Particle::iD]; iterm++;
    *iterm = m_kfm3 * cterm * sums[Particle::iD2]; iterm++;
}

int Condensation::Perform(const real t, System &sys, const unsigned int iterm) const
{
    // Select particle based on which term was called.
    int i;
    switch(iterm) {
        case 0 :
            i = sys.Ensemble().SelectParticle();
        case 1 :
            i = sys.Ensemble().SelectParticle(Particle::iD);
        case 2 :
            i = sys.Ensemble().SelectParticle(Particle::iD2);
        default :
            i = sys.Ensemble().SelectParticle();
    }

    // Check for a valid particle (i>=0).
    if (i >= 0) {
        Particle *sp = sys.Ensemble().GetParticle(i);
        real majr = 0.0;
        if (m_mech->AnyDeferred()) {
            majr = MajorantRate(t, sys, *sp);
            // Update particle with deferred processes.
            m_mech->UpdateParticle(*sp, sys, t);
        } else {
            majr = Rate(t, sys, *sp);
        }

        // Check that the particle is still valid.
        if (sp->IsValid()) {
            // Get the true process rate (after updates).
            real truer = Rate(t, sys, *sp);

            // Check that the event is not ficticious by comparing the
            // majorant rate with the true rate.
            if (!Ficticious(majr, truer)) {
                // Adjust particle.
                AdjustParticle(*sp, t, 1);
                sys.Ensemble().Update(i);

                // Apply changes to gas-phase chemistry.
                map<unsigned int,int>::const_iterator j;
                real vol = sys.SampleVolume();
                for (j=m_reac.begin(); j!=m_reac.end(); j++)
                    sys.AdjustSpeciesConc(j->first, -(real)(j->second) / (NA * vol));
                for (j=m_prod.begin(); j!=m_prod.end(); j++)
                    sys.AdjustSpeciesConc(j->first, (real)(j->second) / (NA * vol));
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

int Condensation::Perform(const real t, System &sys, Particle &sp,
                          const unsigned int n) const
{
    // Adjust particle.
    AdjustParticle(sp, t, n);

    // Apply changes to gas-phase chemistry.
    ApplyToSystem(sys, n);
    map<unsigned int,int>::const_iterator j;

    return 0;
}

/* Property sets. */

void Condensation::SetCondensingSpecies(const real m, const real d)
{
    // Calculate the free-mol terms for condensation.  This must be done
    // before the condensation process is used.
    m_kfm3 = CFM / sqrt(m);
    m_kfm2 = d * m_kfm3 * 2.0;
    m_kfm1 = d * m_kfm2 / 2.0;
}