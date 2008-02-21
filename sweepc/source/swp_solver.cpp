#include "swp_solver.h"
#include "rng.h"
#include <stdlib.h>
#include <cmath>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Solver::Solver(void)
: m_maxdt(0.0), m_splitratio(1.0e9)
{
    // Seed the random number generator.
    srnd(123);
}

// Default destructor.
Solver::~Solver(void)
{
    // Nothing special to destruct.
}


// RUNNING THE SOLVER.

// Performs stochastic stepping algorithm up to specified stop time using
// the given mechanism to define the stochastic processes.  Updates given
// system accordingly.  On error returns <0, otherwise returns 0.
int Solver::Run(real &t, real tstop, Cell &sys, const Mechanism &mech)
{
    int err = 0;
    real tsplit, dtg, dt, jrate;
    static fvector rates(mech.TermCount(), 0.0);

    // Ensure the process counters contain sufficient 
    // entries to track all rate terms.
    m_processcounter.resize(mech.TermCount(), 0);
    m_ficticiouscounter.resize(mech.TermCount(), 0);

    // Global maximum time step.
    dtg     = tstop - t;
    m_maxdt = dtg / 3.0;
    m_tstop = tstop;

    // Loop over time until we reach the stop time.
    while (t < tstop)
    {
        if (mech.AnyDeferred() && (sys.ParticleCount() > 0)) {
            // Get the process jump rates (and the total rate).
            jrate = mech.CalcJumpRates(t, sys, rates);

            // Calculate split end time.
            tsplit = calcSplitTime(t, tstop,jrate, sys.ParticleCount(), dtg);
        } else {
            // There are no deferred processes, therefore there
            // is no need to perform LPDA splitting steps.
            tsplit = tstop;
        }

        // Perform stochastic jump processes.
        while (t < tsplit) {
            jrate = mech.CalcJumpRates(t, sys, rates);
            dt = timeStep(t, sys, mech, rates, jrate);
            if (dt >= 0.0) {
                t += dt;
            } else {
                return -1;
            }
        }

        // Perform Linear Process Deferment Algorithm to
        // update all deferred processes.
        if (mech.AnyDeferred()) {
            mech.LPDA(t, sys);
        }
    }

    return err;
}


// TIME STEPPING ROUTINES.

// Calculates the splitting end time after which all particles
// shallbe updated using LPDA.
real Solver::calcSplitTime(real t, real tstop, real jrate, 
                           unsigned int n, real maxdt) const
{
    // Calculate the splitting time step, ensuring that it is
    // not longer than the maximum allowable time.
    real tsplit = n * m_splitratio / (jrate + 1.0);
    tsplit = min(tsplit, maxdt);

    // Now put the split end time into tsplit, again
    // checking that it is not beyond the stop time.
    return min(tsplit+t, tstop);
}

// Performs a single stochastic event on the ensemble from the given
// mechanism.  Returns the length of the time step if successful,
// otherwise returns negative.
real Solver::timeStep(real t, Cell &sys, const Mechanism &mech, 
                      const fvector &rates, real jrate)
{
    // The purpose of this routine is to perform a single stochastic jump process.  This
    // involves summing the total rate of all processes, generating a waiting time,
    // selecting a process and performing that process.

    int i, j;
    real dt;

    // Calculate exponentially distributed time step size.
    if (jrate > 0.0) {
        dt = rnd();
        dt = - log(dt) / jrate;
    } else {
        // Avoid divide by zero.
        dt = 1.0e30;
    }

    // Truncate if step is too long and select a process
    // to perform.
    if (dt > m_maxdt) {
        dt = m_maxdt;
        i = -1;
    } else {
        if (t+dt <= m_tstop) {
            i = chooseProcess(rates);
        } else {
            i = -1;
        }
    }

    // Perform process.
    if (i >= 0) {
        mech.DoProcess(i, t, sys);
/*        if (j==0) {
            m_processcounter[i]++;
        } else if (j==1) {
            m_ficticiouscounter[i]++;
        } else {
            // Return invalid time step on error.
            dt = -dt;
        }*/
    }

    return dt;
}

// Selects a process using a DIV algorithm and the process rates
// as weights.
int Solver::chooseProcess(const fvector &rates)
{
    // This routine implements a DIV algorithm to select a process from a 
    // discrete list of processes when each process's rate is given.

    // Add together all process rates.
    fvector::const_iterator i;
    real sum = 0.0;
    for (i=rates.begin(); i!=rates.end(); ++i) {
        sum += *i;
    }

    // Generate a uniform random number.
    real r = rnd() * sum ;

    // Select a process using DIV algorithm.
    int j=-1;
    i = rates.begin();
    while((r > 0.0) && (i != rates.end())) {
        r -= *i;
        ++i; ++j;
    }
    return j;
};
