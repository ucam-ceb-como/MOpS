#include "swpinception.h"
#include "swpsystem.h"
#include "swpparams.h"
#include "swpparticle1d.h"
#include "swpsvparticle.h"
#include "swpmechanism.h"

using namespace Sweep;

Inception::Inception(void)
{
    // give default values.  Note the incepting species have
    // not been defined so the rate kernels are zero.
    m_a = 0.5;
    m_kfm = 0.0;
    m_ksf1 = 0.0;
    m_ksf2 = 0.0;
    m_defer = false;
}

Inception::~Inception(void)
{
    // Clear memory associated with inception.
    m_comp.clear();
    m_values.clear();
    m_components = NULL;
}

void Inception::Initialise(const std::map<unsigned int,int> &reac, const std::map<unsigned int,int> &prod, 
                           const Sweep::real a, const Sweep::real m1, const Sweep::real m2, const Sweep::real d1, 
                           const Sweep::real d2, const std::vector<Sweep::real> &comp, const std::vector<Sweep::real> &values, 
                           std::vector<Component*> &components)
{
    // Fill the inception with all the variables required to define it.
    m_reac.clear();
    m_reac.insert(reac.begin(), reac.end());
    m_prod.clear();
    m_prod.insert(prod.begin(), prod.end());
    m_a = 0.5 * a;
    SetInceptingSpecies(m1, m2, d1, d2);
    m_comp.assign(comp.begin(), comp.end());
    m_values.assign(values.begin(), values.end());
    m_components = &components;
}

Sweep::real Inception::Rate(const Sweep::real t, const System &sys) const 
{
    // Get the current chemical conditions.
    vector<real> chem;
    real T, P;
    sys.GetConditions(t, chem, T, P);

    // Calculate the inception rate.
    return Rate(chem, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
}

real Inception::Rate(const real t, const vector<real> &chem, const real T, 
                     const real P, const vector<real> &sums, const System &sys) const
{
    // Chemical conditions have been precalculated, so just calculate the rate.
    return Rate(chem, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
}

void Inception::RateTerms(const real t, const System &sys, vector<real>::iterator &iterm) const
{
    // Get the current chemical conditions.
    vector<real> chem;
    real T, P;
    sys.GetConditions(t, chem, T, P);

    // Calculate the single rate term and advance iterator.
    *iterm = Rate(chem, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
    iterm++;
}

void Inception::RateTerms(const real t, const vector<real> &chem, const real T, 
                          const real P, const vector<real> &sums, const System &sys, 
                          vector<real>::iterator &iterm) const
{
    // Chemical conditions have been precalculated, so just calculate rate term and
    // advance iterator.
    *iterm = Rate(chem, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
    iterm++;
}

int Inception::Perform(const Sweep::real t, System &sys, const unsigned int iterm) const 
{
    // This routine performs the inception on the given chemical system..

    // Create a new particle of the type specified by the mechanism.
    Particle *sp;
    switch (m_mech->GetParticleModel())
    {
        case Mechanism::SphericalParticle:
            sp = new Particle();
            break;
        case Mechanism::SurfaceVolume:
            sp = new SVParticle();
            break;
    }
    
    // Initialise the new particle.
    sp->Initialise(*m_components, (int)m_values.size());
    sp->SetCreateTime(t);
    SetParticle(*sp, t);

    // Add particle to system's ensemble.
    sys.Ensemble().AddParticle(*sp);

    // Update gas-phase chemistry of system.
    map<unsigned int,int>::const_iterator i;
    real vol = sys.SampleVolume();
    for (i=m_reac.begin(); i!=m_reac.end(); i++)
        sys.AdjustSpeciesConc(i->first, -(real)(i->second) / (NA * vol));
    for (i=m_prod.begin(); i!=m_prod.end(); i++)
        sys.AdjustSpeciesConc(i->first, (real)(i->second) / (NA * vol));
    
    return 0;
}

/* Property sets. */

void Inception::SetInceptingSpecies(real m1, real m2, real d1, real d2)
{
    // This routine sets the free-mol and slip flow kernel parameters given
    // the mass and diameter of the incepting species.
    real invd1=1.0/d1, invd2=1.0/d2;
    m_kfm = CFM * sqrt((1.0/m1) + (1.0/m2)) * pow(d1+d2, 2.0);
    m_ksf1 = CSF * (d1+d2);
    m_ksf2 = 1.257 * KNUDSEN_K * m_ksf1 * ((invd1*invd1) + (invd2*invd2));
    m_ksf1 = m_ksf1 * (invd1+invd2);
}