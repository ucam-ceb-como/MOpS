#include "swpsolver.h"
#include "rng.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace Sweep;

Solver::Solver(void)
{
    Initialise();
}

Solver::~Solver(void)
{
}

void Solver::Initialise(void)
{
    // Seed the random number generator.
    srnd(123);

    // Set default values for solver object member variables.
    m_maxdt = 0.0;
    m_splitratio = 1.0e9;
};

int Solver::Run(real *t, const real tstop, System &sys, Mechanism &mech)
{
    int err = 0;
    real tsplit, dtg, dt;
    vector<real> rates;

    // Ensure the process counters contain sufficient entries to track all rate terms.
    m_processcounter.resize(mech.TermCount(), 0);
    m_ficticiouscounter.resize(mech.TermCount(), 0);

    // Global max time step.
    dtg = tstop - *t;
    m_maxdt = dtg / 3.0;

    while (*t < tstop)
    {
        // Calculate LPDA splitting time step.
        if ((mech.AnyDeferred()) && (sys.ParticleCount() > 0)) {
            mech.GetRates(rates, *t, sys);
            tsplit = sys.ParticleCount() * m_splitratio / (mech.JumpRate(rates) + 1.0);
            if (tsplit > dtg) tsplit = dtg;
            m_maxdt = tsplit / 3.0;
            tsplit += *t;
            if (tsplit > tstop) tsplit = tstop;
        } else {
            tsplit = tstop;
        }

        // Perform stochastic jump processes.
        while (*t < tsplit) {
            err = mech.GetStochasticRates(rates, *t, sys);
            if (err < 0) return err;
            dt = TimeStep(*t, sys, mech, rates);
            if (dt >= 0.0) {
                *t += dt;
            } else {
                return -1;
            }
        }

        // Perform Linear Process Deferment Algorithm to
        // update all deferred processes.
        if (mech.AnyDeferred()) {
            err = mech.LPDA(*t, sys) < 0;
            if (err < 0) {
                return err;
            }
        }
    }

    // Release memory.
    rates.clear();

    return err;
}

real Solver::TimeStep(const real t, System &sys, Mechanism &mech, const vector<real> &rates)
{
    /* The purpose of this routine is to perform a single stochastic jump process.  This
       involves summing the total rate of all processes, generating a waiting time,
       selecting a process and performing that process. */

    int i, j;
    real dt;
    real jrate = mech.JumpRate(rates);

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
        i = ChooseProcess(rates, mech.TermCount());
    }

    // Perform process.
    if (i >= 0) {
        j = mech.DoProcess(i, t, sys);
        if (j==0) {
            m_processcounter[i]++;
        } else if (j==1) {
            m_ficticiouscounter[i]++;
        } else {
            // Return invalid time step on error.
            dt = -dt;
        }
    }

    return dt;
}

int Solver::ChooseProcess(const vector<real> &rates, const unsigned int n)
{
    /* This routine implements a DIV algorithm to select a process from a 
       descrete list of processes when each process's rate is given. */

    // Add together all process rates.
    vector<real>::const_iterator i;
    real sum = 0.0;
    for (i=rates.begin(); i!=rates.end(); i++) {
        sum += *i;
    }

    // Generate a uniform random number.
    real r = rnd() * sum ;

    // Select a process using DIV algorithm.
    int j=-1;
    i = rates.begin();
    while((r > 0.0) && (i != rates.end())) {
        r -= *i;
        i++; j++;
    }
    return j;
};

