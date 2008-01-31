#include "swpsurfacereaction.h"
#include "swpparams.h"
#include "swpparticle1d.h"
#include "swpmechanism.h"
#include "swpensemble.h"
#include <cmath>

using namespace Sweep;
using namespace std;

const real SurfaceReaction::SURF_MAJ = 2.0;

SurfaceReaction::SurfaceReaction(void)
{
    Destroy();
    m_defer = true;
}

SurfaceReaction::~SurfaceReaction(void)
{
    Destroy();
}

void SurfaceReaction::Initialise(const std::map<unsigned int,int> &reac, 
                                 const std::map<unsigned int,int> &prod, 
                                 const Sweep::real a, const Sweep::real n, 
                                 const Sweep::real e, const vector<Sweep::real> &comp, 
                                 const vector<Sweep::real> &values, const unsigned int pid,
                                 vector<Component*> &components)
{
    Destroy();
    m_reac.insert(reac.begin(),reac.end());
    m_prod.insert(reac.begin(),reac.end());
    m_a = a;
    m_n = n;
    m_e = e;
    m_comp.assign(comp.begin(), comp.end());
    m_values.assign(values.begin(), values.end());
    m_pid = pid;
    m_components = &components;
}

void SurfaceReaction::Destroy()
{
    m_reac.clear();
    m_prod.clear();
    m_a = 1.0;
    m_n = 0.0;
    m_e = 0.0;
    m_pid = -1;
}

Sweep::real SurfaceReaction::Rate(const real t, const System &sys) const
{
    // Rate constant.
    real rate = m_a;

    // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        rate *= pow(sys.GetSpeciesConc((*i).first,t), (real)(*i).second);
    }

    // Tempearature dependance.
    real T = sys.GetTemperature(t);
    rate *= pow(T, m_n) * exp(-m_e / (R * T));

    // Paticle dependence.
    vector<real> sums; sys.ConstEnsemble().GetSums(sums);
    rate *= sums[m_pid];

    if (m_mech->AnyDeferred()) {
        return rate * SURF_MAJ;
    } else {
        return rate;
    }
}

real SurfaceReaction::Rate(const real t, const vector<real> &chem, const real T, 
                           const real P, const vector<real> &sums, const System &sys) const
{
    // Rate constant.
    real rate = m_a;

    // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        rate *= pow(chem[(*i).first], (real)(*i).second);
    }

    // Tempearature dependance.
    rate *= pow(T, m_n) * exp(-m_e / (R * T));

    // Paticle dependence.
    rate *= sums[m_pid];

    if (m_mech->AnyDeferred()) {
        return rate * SURF_MAJ;
    } else {
        return rate;
    }
}

Sweep::real SurfaceReaction::Rate(const real t, const System &sys, const Particle &sp) const
{
    // Rate constant.
    real rate = m_a;

    // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator j;
    for (j=m_reac.begin(); j!=m_reac.end(); j++) {
        rate *= pow(sys.GetSpeciesConc((*j).first,t), (real)(*j).second);
    }

    // Tempearature dependance.
    real T = sys.GetTemperature(t);
    rate *= pow(T, m_n) * exp(-m_e / (R * T));

    // Paticle dependence.
    rate *= sp.GetProperty(m_pid);

    return rate;
}

Sweep::real SurfaceReaction::Rate(const real t, const vector<real> &chem, const real T, 
                                  const real P, const vector<real> &sums, const System &sys, 
                                  const Particle &sp) const
{
    // Rate constant.
    real rate = m_a;

    // Chemical species concentration dependence.
    map<unsigned int,int>::const_iterator j;
    for (j=m_reac.begin(); j!=m_reac.end(); j++) {
        rate *= pow(sys.GetSpeciesConc((*j).first,t), (real)(*j).second);
    }

    // Tempearature dependance.
    rate *= pow(T, m_n) * exp(-m_e / (R * T));

    // Paticle dependence.
    rate *= sp.GetProperty(m_pid);

    return rate;
}

Sweep::real SurfaceReaction::MajorantRate(const real t, const System &sys, 
                                          const Particle &sp) const
{
    return Rate(t,sys,sp) * SURF_MAJ;
}

Sweep::real SurfaceReaction::MajorantRate(const real t, const vector<real> &chem, const real T, 
                                          const real P, const vector<real> &sums, const System &sys, 
                                          const Particle &sp) const
{
    return Rate(t,chem,T,P,sums,sys,sp) * SURF_MAJ;
}

void SurfaceReaction::RateTerms(const real t, const System &sys, vector<real>::iterator &iterm) const
{
    *iterm = Rate(t, sys);
    iterm++;
}

void SurfaceReaction::RateTerms(const real t, const vector<real> &chem, const real T, 
                                const real P, const vector<real> &sums, const System &sys, 
                                vector<real>::iterator &iterm) const
{
    *iterm = Rate(t, chem, T, P, sums, sys);
    iterm++;
}

int SurfaceReaction::Perform(const Sweep::real t, Sweep::System &sys, const unsigned int iterm) const
{
    int i = sys.Ensemble().SelectParticle(m_pid);

    if (i >= 0) {
        Particle *sp = sys.Ensemble().GetParticle(i);
        real majr = MajorantRate(t, sys, *sp);

        // Update particle with deferred processes.
        if (m_mech->AnyDeferred()) {
            m_mech->UpdateParticle(*sp, sys, t);
        }

        // Check that the particle is still valid.
        if (sp->IsValid()) {
            real truer = Rate(t, sys, *sp);

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

int SurfaceReaction::Perform(const real t, System &sys, 
                             Particle &sp, const unsigned int n) const
{
    // Adjust particle.
    AdjustParticle(sp, t, n);

    // Apply changes to gas-phase chemistry.
    ApplyToSystem(sys, n);

    return 0;
}