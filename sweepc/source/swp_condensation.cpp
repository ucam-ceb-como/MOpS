#include "swp_condensation.h"
#include "swp_mechanism.h"
#include <cmath>

using namespace Sweep;

const unsigned int Condensation::TERM_COUNT = 3;
const real Condensation::m_majfactor       = 2.0;
const real Condensation::m_efm             = 2.2;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Condensation::Condensation(void)
: m_a(1.0), m_kfm1(0.0), m_kfm2(0.0), m_kfm3(0.0)
{
    // Assume the condensation is simulated as a deferred process (LPDA).
    m_defer = true;
}

// Copy constructor.
Condensation::Condensation(const Condensation &copy)
{
    *this = copy;
}

// Default destructor.
Condensation::~Condensation(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
Condensation &Condensation::operator=(const Condensation &rhs)
{
    if (this != &rhs) {
        ParticleProcess::operator =(rhs);
        m_a    = rhs.m_a;
        m_kfm1 = rhs.m_kfm1;
        m_kfm2 = rhs.m_kfm2;
        m_kfm3 = rhs.m_kfm3;
    }
    return *this;
}


// RATE CONSTANT AND PARAMETERS.

// Returns the fixed rate constant.
real Condensation::A() const {return m_a;}

// Sets the fixed rate constant.
void Condensation::SetA(real a) {m_a = a;}

// Sets the coagulation kernel parameters given the mass and
// collision diameter of the condensing species.
void Condensation::SetCondensingSpecies(const real m, const real d)
{
    // Calculate the free-mol terms for condensation.  This must be done
    // before the condensation process is used.
    m_kfm3 = m_efm * CFM / sqrt(m);
    m_kfm2 = d * m_kfm3 * 2.0;
    m_kfm1 = d * m_kfm2 / 2.0;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

// Returns rate of the process for the given system.
real Condensation::Rate(real t, const Cell &sys) const
{
    return Rate(t, sys, sys);
}

// Calculates the process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real Condensation::Rate(real t, const Sprog::Thermo::IdealGas &gas, const Cell &sys) const
{
    // Calculate temperature terms.
    real cterm = m_a * sqrt(gas.Temperature());

     // Chemical species concentration dependence.
    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        for (unsigned int j=0; j<i->second; ++j) {
            cterm *= gas.MolarConc(i->first);
        }
    }

    // Free molecular terms.
    cterm *= (m_kfm1 * sys.ParticleCount()) + 
             (m_kfm2 * sys.Particles().GetSum(ParticleData::iDcol)) +
             (m_kfm3 * sys.Particles().GetSum(CoagModel_ID, CoagModelData::iD2));

    // If the mechanism contains any deferred processes then we must use the
    // majorant form of the rate, in order to account for any changes to
    // particles during deferred process updates.
    if (m_mech->AnyDeferred()) {
        return cterm * m_majfactor;
    } else {
        return cterm;
    }
}


// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
real Condensation::Rate(real t, const Cell &sys, const Particle &sp) const
{
    return Rate(t, sys, sys, sp);
}

// Returns rate of the process for the given particle using the
// given chemical conditions rather than those conditions in the
// the given system.
real Condensation::Rate(real t, const Sprog::Thermo::IdealGas &gas, 
                        const Cell &sys, const Particle &sp) const
{
    // Calculate temperature terms.
    real cterm = m_a * sqrt(gas.Temperature());
//    real trm[3];

    // Chemical species concentration dependence.
    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        for (unsigned int j=0; j<i->second; ++j) {
            cterm *= gas.MolarConc(i->first);
        }
    }

    // Free molecular terms.
//    trm[0] = cterm * m_kfm1;
//    trm[1] = cterm * (m_kfm2 * sp.CollDiameter());
//    trm[2] = cterm * (m_kfm3 * sp.CoagModelCache->CollDiamSquared());
    cterm *= m_kfm1 + 
             (m_kfm2 * sp.CollDiameter()) +
             (m_kfm3 * sp.CoagModelCache()->CollDiamSquared());
    return cterm; //trm[0] + trm[1] + trm[2];
}

// Returns majorant rate of the process for the given system.
real Condensation::MajorantRate(real t, const Cell &sys, const Particle &sp) const
{
    // Return the single particle rate multiplied by the 
    // condensation majorant factor.
    return Rate(t, sys, sys, sp) * m_majfactor;
}

// Calculates the majorant process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real Condensation::MajorantRate(real t, const Sprog::Thermo::IdealGas &gas, 
                                const Cell &sys, const Particle &sp) const
{
    // Return the single particle rate multiplied by the 
    // condensation majorant factor.
    return Rate(t, gas, sys, sp) * m_majfactor;
}


// RATE TERM CALCULATIONS.
//   These routines return the individual rate terms for a 
//   process, which may have multiple terms (e.g. condensation).

// Returns the number of rate terms for this process.
unsigned int Condensation::TermCount(void) const {return TERM_COUNT;}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.
real Condensation::RateTerms(real t, const Cell &sys, 
                             fvector::iterator &iterm) const
{
    return RateTerms(t, sys, sys, iterm);
}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  The given chemical conditions are used instead of those
// in the given system object.
real Condensation::RateTerms(real t, const Sprog::Thermo::IdealGas &gas, 
                             const Cell &sys, fvector::iterator &iterm) const
{
    // Calculate temperature terms.
    real cterm = m_a * sqrt(gas.Temperature());

     // Chemical species concentration dependence.
    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        for (unsigned int j=0; j<i->second; ++j) {
            cterm *= gas.MolarConc(i->first);
        }
    }

    // If the mechanism contains any deferred processes then we must use the
    // majorant form of the rate, in order to account for any changes to
    // particles during deferred process updates.
    if (m_mech->AnyDeferred()) cterm *= m_majfactor;

    // Free molecular terms.
    real sum = 0.0;
    sum += *(iterm++) = m_kfm1 * cterm * sys.ParticleCount();
    sum += *(iterm++) = m_kfm2 * cterm * sys.Particles().GetSum(ParticleData::iDcol);
    sum += *(iterm++) = m_kfm3 * cterm * 
                        sys.Particles().GetSum(CoagModel_ID, CoagModelData::iD2);
    return sum;
}


// PERFORMING THE PROCESS.

// Performs the process on the given system.  The responsible rate term is given
// by index.  Returns 0 on success, otherwise negative.
int Condensation::Perform(const real t, Cell &sys, const unsigned int iterm) const
{
    // Select particle based on which term was called.
    int i = -1;
    switch(iterm) {
        case 1:
            i = sys.Particles().SelectParticle(ParticleData::iDcol);
            break;
        case 2:
            i = sys.Particles().SelectParticle(CoagModel_ID, CoagModelData::iD2);
            break;
        case 0:
        default:
            i = sys.Particles().SelectParticle();
            break;
    }

    // Check for a valid particle (i>=0).
    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);

        real majr = MajorantRate(t, sys, *sp);

        if (m_mech->AnyDeferred()) {
            // Update particle with deferred processes.
            m_mech->UpdateParticle(*sp, sys, t);
        }

        // Check that the particle is still valid.
        if (sp->IsValid()) {
            // Get the true process rate (after updates).
            real truer = Rate(t, sys, *sp);

            // Check that the event is not ficticious by comparing the
            // majorant rate with the true rate.
            if (!Ficticious(majr, truer)) {
                // Adjust particle.
                sp->Adjust(m_dcomp, m_dvals);
                sys.Particles().Update(i);

                // Apply changes to gas-phase chemistry.
                adjustGas(sys);
            }
        } else {
            // If not valid then remove the particle.
            sys.Particles().Remove(i);
        }
    } else {
        // Failed to select a particle.
        return -1;
    }
    return 0;
}

// Performs the process on a given particle in the system.  Particle
// is given by index.  The process is performed n times.
int Condensation::Perform(real t, Cell &sys, Particle &sp,
                          unsigned int n) const
{
    unsigned int m = sp.Adjust(m_dcomp, m_dvals, n);
    adjustGas(sys, m);
    return 0;
}
