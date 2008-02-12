#include "swp_surfacereaction.h"
#include "swp_mechanism.h"
#include <cmath>

using namespace Sweep;
using namespace std;

const real SurfaceReaction::m_majfactor = 2.0;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SurfaceReaction::SurfaceReaction(void)
: m_arr(0.0,0.0,0.0), m_pid(0), m_modelid(BasicModel_ID)
{
    m_defer = true;
}

// Copy constructor.
SurfaceReaction::SurfaceReaction(const SurfaceReaction &copy)
{
    *this = copy;
}

// Default destructor.
SurfaceReaction::~SurfaceReaction(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
SurfaceReaction &SurfaceReaction::operator =(const Sweep::SurfaceReaction &rhs)
{
    if (this != &rhs) {
        ParticleProcess::operator =(rhs);
        m_arr     = rhs.m_arr;
        m_pid     = rhs.m_pid;
        m_modelid = rhs.m_modelid;
    }
    return *this;
}


// ARRHENIUS COEFFICIENTS.

// Returns the Arrhenius parameter.
Sprog::Kinetics::ARRHENIUS &SurfaceReaction::Arrhenius() {return m_arr;}
const Sprog::Kinetics::ARRHENIUS &SurfaceReaction::Arrhenius() const {return m_arr;}

// Sets the Arrhenius parameters.
void SurfaceReaction::SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr) {m_arr = arr;}




// PARTICLE PROPERTY ID.

// Returns the ID number of the particle property to which
// the rate of this process is proportional.
unsigned int SurfaceReaction::PropertyID(void) const {return m_pid;}

// Returns the ID number of the particle number for which the
// PropertyID is valid.  The mechanism should check that this
// model is enabled.
ModelType SurfaceReaction::ModelID(void) const {return m_modelid;}

// Sets the ID number of the particle property to which
// the rate of this process is proportional.
void SurfaceReaction::SetPropertyID(unsigned int i, ModelType modelid)
{
    m_pid     = i;
    m_modelid = modelid;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

// Returns rate of the process for the given system.
real SurfaceReaction::Rate(real t, const Cell &sys) const
{
    return Rate(t, sys, sys);
}

// Calculates the process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real SurfaceReaction::Rate(real t, const Sprog::Thermo::IdealGas &gas, 
                           const Cell &sys) const
{
    // Rate constant.
    real rate = m_arr.A;

    // Chemical species concentration dependence.
    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        for (unsigned int j=0; j<i->second; ++j) {
            rate *= gas.MolarConc(i->first);
        }
    }

    // Tempearature dependance.
    real T = gas.Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

    // Paticle dependence.
    rate *= sys.Particles().GetSum(m_modelid, m_pid);

    if (m_mech->AnyDeferred()) {
        return rate * m_majfactor;
    } else {
        return rate;
    }
}


// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
real SurfaceReaction::Rate(real t, const Cell &sys, const Particle &sp) const
{
    return Rate(t,sys, sys, sp);
}

// Returns rate of the process for the given particle using the
// given chemical conditions rather than those conditions in the
// the given system.
real SurfaceReaction::Rate(real t, const Sprog::Thermo::IdealGas &gas, 
                           const Cell &sys, const Particle &sp) const
{
    // Rate constant.
    real rate = m_arr.A;

    // Chemical species concentration dependence.
    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        for (unsigned int j=0; j<i->second; ++j) {
            rate *= gas.MolarConc(i->first);
        }
    }

    // Tempearature dependance.
    real T = gas.Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

    // Paticle dependence.
    if (m_modelid == BasicModel_ID) {
        rate *= sp.Property(static_cast<ParticleData::PropertyID>(m_pid));
    } else {
        rate *= sp.ModelCache(m_modelid)->Property(m_pid);
    }

    return rate;
}

// Returns majorant rate of the process for the given system.
real SurfaceReaction::MajorantRate(real t, const Cell &sys, 
                                   const Particle &sp) const
{
    return Rate(t, sys, sys, sp) * m_majfactor;
}

// Calculates the majorant process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real SurfaceReaction::MajorantRate(real t, const Sprog::Thermo::IdealGas &gas, 
                                   const Cell &sys, const Particle &sp) const
{
    return Rate(t, gas, sys, sp) * m_majfactor;
}


// RATE TERM CALCULATIONS.
//   These routines return the individual rate terms for a 
//   process, which may have multiple terms (e.g. condensation).

// Returns the number of rate terms for this process.
unsigned int SurfaceReaction::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.
real SurfaceReaction::RateTerms(real t, const Cell &sys, 
                                fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, sys, sys);
}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  The given chemical conditions are used instead of those
// in the given system object.
real SurfaceReaction::RateTerms(real t, const Sprog::Thermo::IdealGas &gas, 
                                const Cell &sys, fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, gas, sys);
}


// PERFORMING THE PROCESS.

// Performs the process on the given system.  The responsible rate term is given
// by index.  Returns 0 on success, otherwise negative.
int SurfaceReaction::Perform(real t, Cell &sys, unsigned int iterm) const
{
    int i = sys.Particles().SelectParticle(m_modelid, m_pid);

    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);

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
int SurfaceReaction::Perform(real t, Cell &sys, Particle &sp, 
                             unsigned int n) const
{
    unsigned int m = sp.Adjust(m_dcomp, m_dvals, n);
    adjustGas(sys, m);
    return 0;
}
