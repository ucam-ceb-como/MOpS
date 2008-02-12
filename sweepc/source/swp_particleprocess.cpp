#include "swp_particleprocess.h"
#include "swp_mechanism.h"

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ParticleProcess::ParticleProcess()
: m_defer(false)
{
}

// Copy constructor.
ParticleProcess::ParticleProcess(const Sweep::ParticleProcess &copy)
{
    *this = copy;
}

// Default destructor.
ParticleProcess::~ParticleProcess()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
ParticleProcess &ParticleProcess::operator=(const Sweep::ParticleProcess &rhs)
{
    if (this != &rhs) {
        Process::operator=(rhs);
        m_defer = rhs.m_defer;
    }
    return *this;
}


// DEFERRED PROCESSES.

// Returns TRUE if process should be deferred, otherwise false.
bool ParticleProcess::IsDeferred(void) const
{
    return m_defer;
}

// Sets the process to be deferred or not.
void ParticleProcess::SetDeferred(bool defer)
{
    m_defer = defer; 
    if (m_mech!=NULL) m_mech->CheckDeferred();
}


// RATE CALCULATION.

// Calculates the rates of multiple particle processes.
real ParticleProcess::CalcRates(real t, const Cell &sys, 
                                const PartProcPtrVector &proc,
                                fvector &rates, unsigned int start)
{
    PartProcPtrVector::const_iterator p;
    fvector::iterator i = (rates.begin()+start);
    real sum = 0.0;
    for (p=proc.begin(); p!=proc.end(); ++p,++i) {
        *i = (*p)->Rate(t, sys);    
        sum += *i;
    }
    return sum;
}
