#include "swp_system.h"

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
Cell::Cell(void)
: m_smpvol(1.0)
{
}

// Default constructor (public).
Cell::Cell(const Sprog::SpeciesPtrVector &sp)
: Sprog::Thermo::IdealGas(sp)
{
}

// Copy constructor.
Cell::Cell(const Cell &copy)
{
    *this = copy;
}

// Default destructor.
Cell::~Cell(void)
{
}


// OPERATOR OVERLOADS.

// Assignment operator.
Cell &Cell::operator=(const Sweep::Cell &rhs)
{
    if (this != &rhs) {
        Sprog::Thermo::IdealGas::operator=(rhs);
        m_ensemble = rhs.m_ensemble;
        m_smpvol   = rhs.m_smpvol;
    }
    return *this;
}


// THE GAS-PHASE INTERFACE.

// Returns the description of the gas-phase mixture.
const Sprog::Thermo::IdealGas &Cell::GasPhase(void) const
{
    return *this;
}

// Sets the gas-phase mixture.
void Cell::SetGasPhase(const Sprog::Thermo::IdealGas &gas)
{
    Sprog::Thermo::IdealGas::operator=(gas);
}

// Adjusts the concentration of the ith species.
void Cell::AdjustConc(unsigned int i, real dc)
{
    unsigned int k;

    // Precalculate DC / density.
    real dc_rho = dc / *m_pdens;

    // Calculate change to all mole fractions k < i.
    for (k=0; k<i; ++k) {
        m_data[k] -= dc_rho * m_data[k];
    }

    // Calculate change for ith species.
    m_data[i] += dc_rho * (1.0 - m_data[i]);

    // Calculate change for all mole fractions k > i.
    for (k=i+1; k<m_species->size(); ++k) {
        m_data[k] -= dc_rho * m_data[k];
    }
}

// Adjusts the concentration of all species.
void Cell::AdjustConcs(const fvector &dc)
{
    // Calculate total change in density.
    real drho = 0.0;
    unsigned int k;
    for (k=0; k!=m_species->size(); ++k) {
        drho += dc[k];
    }

    // Calculate changes to the mole fractions.
    real invrho = 1.0 / *m_pdens;
    for (k=0; k!=m_species->size(); ++k) {
        m_data[k] += (invrho * dc[k]) - (invrho * m_data[k] * drho);
    }
}


// THE PARTICLE ENSEMBLE.

// Returns the particle ensemble.
Ensemble &Cell::Particles(void) 
{
    return m_ensemble;
}

// Returns the particle count.
unsigned int Cell::ParticleCount(void) const 
{
    return m_ensemble.Count();
}


// SCALING ROUTINES INCL. SAMPLE VOLUME.

// Returns the real system to stochastic system scaling factor.
real Cell::SampleVolume() const
{
    return m_smpvol * m_ensemble.Scaling();
}

// Sets the number density which the ensemble represents.
int Cell::SetM0(const real m0)
{
    if ((m_ensemble.Count() > 0) && (m0 > 0.0)) {
        m_smpvol = (real)m_ensemble.Count() / m0;
        m_ensemble.ResetScaling();
        return 0;
    } else {
        // The ensemble contains no particles, so assume this
        // is the maximum M0.        
        return SetMaxM0(m0);
    }
}

// Sets the number density which the full 
// ensemble would represent.
int Cell::SetMaxM0(const real m0)
{
    if ((m_ensemble.Capacity() > 0) && (m0 > 0.0)) {
        m_smpvol = m_ensemble.Scaling() * (real)m_ensemble.Capacity() / m0;
        m_ensemble.ResetScaling();
        return 0;
    } else {
        // The ensemble has not yet been initialised, hence guess
        // unit volume and report an error.
        m_smpvol = 1.0;
        return -1;
    }
}

void Cell::Reset(const real m0)
{
    m_ensemble.Clear();
    SetMaxM0(m0);
}
