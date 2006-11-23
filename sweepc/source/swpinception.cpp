#include "swpinception.h"
#include "swpsystem.h"
#include "swpparams.h"
#include "swpparticle1d.h"
#include "swpsvparticle.h"
#include "swpmechanism.h"

using namespace Sweep;

Inception::Inception(void)
{
    m_a = 0.5;
    m_kfm = 0.0;
    m_ksf1 = 0.0;
    m_ksf2 = 0.0;
    m_defer = false;
}

Inception::~Inception(void)
{
    m_comp.clear();
    m_values.clear();
    m_components = NULL;
}

void Inception::Initialise(const std::map<unsigned int,int> &reac, const std::map<unsigned int,int> &prod, 
                           const Sweep::real a, const Sweep::real m1, const Sweep::real m2, const Sweep::real d1, 
                           const Sweep::real d2, const std::vector<Sweep::real> &comp, const std::vector<Sweep::real> &values, 
                           std::vector<Component*> &components)
{

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
    vector<real> chem;
    real T, P;
    sys.GetConditions(t, chem, T, P);
    return Rate(chem, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
}

real Inception::Rate(const real t, const vector<real> &chem, const real T, 
                     const real P, const vector<real> &sums, const System &sys) const
{
    return Rate(chem, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
}

void Inception::RateTerms(const real t, const System &sys, vector<real>::iterator &iterm) const
{
    vector<real> chem;
    real T, P;
    sys.GetConditions(t, chem, T, P);
    *iterm = Rate(chem, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
    iterm++;
}

void Inception::RateTerms(const real t, const vector<real> &chem, const real T, 
                          const real P, const vector<real> &sums, const System &sys, 
                          vector<real>::iterator &iterm) const
{
    *iterm = Rate(chem, sqrt(T), T/Viscosity(T), T/P, sys.SampleVolume());
    iterm++;
}

int Inception::Perform(const Sweep::real t, System &sys, const unsigned int iterm) const 
{
    // Create a new particle.
    DefaultParticle *sp;
    switch (m_mech->GetParticleModel())
    {
        case Mechanism::SphericalParticle:
            sp = new DefaultParticle();
            break;
        case Mechanism::SurfaceVolume:
            sp = new SVParticle();
            break;
    }
        
    sp->Initialise(*m_components, (int)m_values.size());
    sp->SetCreateTime(t);
    SetParticle(*sp, t);

    // Add particle to ensemble.
    sys.Ensemble().AddParticle(*sp);

    // Update gas-phase chemistry.
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
    real invd1=1.0/d1, invd2=1.0/d2;
    m_kfm = CFM * sqrt((1.0/m1) + (1.0/m2)) * pow(d1+d2, 2.0);
    m_ksf1 = CSF * (d1+d2);
    m_ksf2 = 1.257 * KNUDSEN_K * m_ksf1 * ((invd1*invd1) + (invd2*invd2));
    m_ksf1 = m_ksf1 * (invd1+invd2);
}