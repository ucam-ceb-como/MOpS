#include "swp_inception.h"
#include "swp_mechanism.h"
#include <cmath>

using namespace Sweep;

// Free-molecular enhancement factor.
const real Inception::m_efm = 2.2; // 2.2 is for soot.

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Inception::Inception(void)
: m_a(0.5), m_kfm(0.0), m_ksf1(0.0), m_ksf2(0.0)
{
}

// Copy constructor.
Inception::Inception(const Sweep::Inception &copy)
{
    *this = copy;
}

// Default destructor.
Inception::~Inception(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
Inception &Inception::operator =(const Inception &rhs)
{
    if (this != &rhs) {
        Process::operator =(rhs);
        m_a    = rhs.m_a;
        m_kfm  = rhs.m_kfm;
        m_ksf1 = rhs.m_ksf1;
        m_ksf2 = rhs.m_ksf2;
    }
    return *this;
}


// RATE CONSTANT.

// Returns the rate constant.
real Inception::A(void) const {return m_a;}

// Sets the rate constant.
void Inception::SetA(real a) {m_a = a;}


// INCEPTION KERNEL.

// Sets the coagulation kernel constants given incepting species
// masses and diameters.
void Inception::SetInceptingSpecies(real m1, real m2, real d1, real d2)
{
    // This routine sets the free-mol and slip flow kernel parameters given
    // the mass and diameter of the incepting species.
    real invd1=1.0/d1, invd2=1.0/d2;
    m_kfm  = m_efm * CFM * sqrt((1.0/m1) + (1.0/m2)) * pow(d1+d2, 2.0);
    m_ksf1 = CSF * (d1+d2);
    m_ksf2 = 1.257 * m_ksf1 * ((invd1*invd1) + (invd2*invd2));
    m_ksf1 = m_ksf1 * (invd1+invd2);
}


// PROPERTIES OF INCEPTED PARTICLES.

// Returns the composition vector of the new particle.
const fvector &Inception::ParticleComp(void) const {return m_newcomp;}

// Returns the amount of the ith component of the new particle.
real Inception::ParticleComp(unsigned int i) const
{
    if (i < m_mech->ComponentCount()) {
        return m_newcomp[i];
    } else {
        return 0.0;
    }
}

// Sets the particle composition vector.
void Inception::SetParticleComp(const fvector &comp)
{
    m_newcomp.assign(comp.begin(), comp.end());
}

// Sets the amount of the ith component in the new particle.
void Inception::SetParticleComp(unsigned int i, real comp)
{
    if (i < m_mech->ComponentCount()) {
        m_newcomp[i] = comp;
    }
}

// Returns the tracker variable vector of the new particle.
const fvector &Inception::ParticleTrackers(void) const
{
    return m_newvals;
}

// Returns the value of the ith tracker variable of the
// new particle.
real Inception::ParticleTrackers(unsigned int i) const
{
    if (i < m_mech->TrackerCount()) {
        return m_newvals[i];
    } else {
        return 0.0;
    }
}

// Sets the new particle tracker variable vector.
void Inception::SetParticleTrackers(const fvector &track)
{
    m_newvals.assign(track.begin(), track.end());
}

// Sets the value of the ith tracker variable in the
// new particle.
void Inception::SetParticleTracker(unsigned int i, real track)
{
    if (i < m_mech->TrackerCount()) {
        m_newvals[i] = track;
    }
}


// TOTAL RATE CALCULATIONS.

// Returns rate of the process for the given system.
real Inception::Rate(real t, const Cell &sys) const 
{
    // Get the current chemical conditions.
    real T = sys.Temperature();
    real P = sys.Pressure();

    // Calculate the rate.
    return Rate(sys.MoleFractions(), sys.Density(), sqrt(T), 
                T/ViscosityAir(T), MeanFreePathAir(T,P), 
                sys.SampleVolume());
}

// Calculates the process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real Inception::Rate(const real t, const Sprog::Thermo::IdealGas &gas, 
                     const Cell &sys) const
{
    // Get the current chemical conditions.
    real T = gas.Temperature();
    real P = gas.Pressure();

    // Calculate the rate.
    return Rate(gas.MoleFractions(), gas.Density(), sqrt(T), 
                T/ViscosityAir(T), MeanFreePathAir(T,P), 
                sys.SampleVolume());
}

// Calculates the rate of multiple inceptions given a
// vector of inceptions and an iterator to a vector of
// reals for output.
real Inception::CalcRates(real t, const Cell &sys, 
                          const IcnPtrVector &icns, 
                          fvector &rates, unsigned int start)
{
    // Precalculate some values.
    real T     = sys.Temperature();
    real P     = sys.Pressure();
    real sqrtT = sqrt(T);
    real T_mu  = T / ViscosityAir(T);
    real MFP   = MeanFreePathAir(T,P);
    real vol   = sys.SampleVolume();

    IcnPtrVector::const_iterator p;
    fvector::iterator i = (rates.begin()+start);
    real sum = 0.0;
    for (p=icns.begin(); p!=icns.end(); ++p,++i) {
        *i = (*p)->Rate(sys.MoleFractions(), sys.Density(), T, 
                        T_mu, MFP, vol);
        sum += *i;
    }
    return sum;
}

// A faster rate calculation routine for Inception events only.  Requires all the
// parameters that would otherwise be calculated by the routine to be passed as
// arguments.
real Inception::Rate(const fvector &fracs, real density, real sqrtT, 
                     real T_mu, real MVP, real vol) const
{
    // Temperature and pressure dependence.
    real fm   = sqrtT * m_kfm;
    real sf   = T_mu  * (m_ksf1 + (MVP*m_ksf2));
    real rate = m_a   * (fm*sf) / (fm+sf) * vol;

    // Chemical species concentration dependence.
    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        for (unsigned int j=0; j!=i->second; ++j) {
            rate *= NA * density * fracs[(*i).first];
        }
    }

    return rate;
}


// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int Inception::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
real Inception::RateTerms(const real t, const Cell &sys, 
                          fvector::iterator &iterm) const
{
    // Get the current chemical conditions.
    real T = sys.Temperature();
    real P = sys.Pressure();

    // Calculate the single rate term and advance iterator.
    *iterm = Rate(sys.MoleFractions(), sys.Density(), sqrt(T), 
                  T/ViscosityAir(T), MeanFreePathAir(T,P), 
                  sys.SampleVolume());
    return *(iterm++);
}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  The given chemical conditions are used instead of those
// in the given system object.
real Inception::RateTerms(const real t, const Sprog::Thermo::IdealGas &gas,
                          const Cell &sys, fvector::iterator &iterm) const
{
    // Get the current chemical conditions.
    real T = gas.Temperature();
    real P = gas.Pressure();

    // Calculate rate term and advance iterator.
    *iterm = Rate(gas.MoleFractions(), gas.Density(), sqrt(T), 
                  T/ViscosityAir(T), MeanFreePathAir(T,P), 
                  sys.SampleVolume());
    return *(iterm++);
}


// PERFORMING THE PROCESS.

// Performs the process on the given system.  The responsible rate term is given
// by index.  Returns 0 on success, otherwise negative.
int Inception::Perform(real t, Cell &sys, unsigned int iterm) const 
{
    // This routine performs the inception on the given chemical system.

    // Create a new particle of the type specified
    // by the system ensemble.
    Particle *sp = m_mech->CreateParticle();
    
    // Initialise the new particle.
    sp->SetCreateTime(t);
    sp->SetTime(t);
    sp->SetComposition(m_newcomp);
    sp->SetValues(m_newvals);
    sp->UpdateCache();

    // Add particle to system's ensemble.
    sys.Particles().Add(*sp);

    // Update gas-phase chemistry of system.
    adjustGas(sys);

    return 0;
}
