#include "swpsystem.h"
#include "swpparticle1d.h"
#include "swpensemble.h"

#include <fstream>

using namespace Sweep;

System::System(void)
{
    m_smpvol = 1.0;
    m_species = NULL;
}

System::~System(void)
{
    m_ensemble.Destroy();
}

Sweep::Ensemble &System::Ensemble(void) 
{
    return m_ensemble;
}

const Sweep::Ensemble &System::ConstEnsemble(void) const 
{
    return m_ensemble;
}

unsigned int System::ParticleCount(void) const 
{
    return m_ensemble.Count();
}

void System::SetSpeciesList(SpeciesList &list)
{
    m_species = &list;
}

SpeciesList &System::GetSpeciesList()
{
    return *m_species;
}

real System::SampleVolume() const
{
    return m_smpvol * m_ensemble.Scaling();
}

int System::SetM0(const real m0)
{
    if ((m_ensemble.Count() > 0) && (m0 > 0.0)) {
        m_smpvol = (real)m_ensemble.Count() / m0;
        return 0;
    } else {
        // The ensemble contains no particles, so assume this
        // is the maximum M0.        
        return SetMaxM0(m0);
    }
}

int System::SetMaxM0(const real m0)
{
    if ((m_ensemble.Capacity() > 0) && (m0 > 0.0)) {
        m_smpvol = (real)m_ensemble.Capacity() / m0;
        return 0;
    } else {
        // The ensemble has not yet been initialised, hence guess
        // unit volume and report an error.
        m_smpvol = 1.0;
        return -1;
    }
}

void System::Reset(const real m0)
{
    m_ensemble.Clear();
    SetMaxM0(m0);
}
