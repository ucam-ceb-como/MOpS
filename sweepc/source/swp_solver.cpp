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
int Solver::Run(real &t, real tstop, Cell &sys, Mechanism &mech)
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


// Performs stochastic stepping algorithm up to specified stop time using
// the given mechanism to define the stochastic processes.  Updates given
// system accordingly.  On error returns <0, otherwise returns 0.  In this
// flavour the gas-phase chemistry is interpolated from a vector of
// IdealGas objects rather than being taken from the given system object.
// However, the particles in the system object are updated accordingly.
int Solver::Run(real &t, real tstop, const GasProfile &gasphase, 
                Cell &sys, Mechanism &mech)
{
    int err = 0;
    real tsplit, dtg, dt, jrate;
    fvector rates(mech.TermCount(), 0.0);
    Sprog::Thermo::IdealGas gas(*mech.Species());

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
        // Calculate LPDA splitting time step.
        if (mech.AnyDeferred() && (sys.ParticleCount() > 0)) {
            // Calculate the chemical conditions.
            linInterpGas(t, gasphase, gas);

            // Get the process jump rates (and the total rate).
            jrate = mech.CalcJumpRates(t, gas, sys, rates);

            // Calculate the splitting end time.
            tsplit = calcSplitTime(t, tstop, jrate, sys.ParticleCount(), dtg);
        } else {
            // There are no deferred processes, therefore there
            // is no need to perform LPDA splitting steps.
            tsplit = tstop;
        }

        // Perform stochastic jump processes.
        while (t < tsplit) {
            // Calculate the chemical conditions.
            linInterpGas(t, gasphase, gas);

            // Calculate jump rates.
            jrate = mech.CalcJumpRates(t, gas, sys, rates);

            // Perform time step.
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


// HELPER ROUTINES.

void Solver::linInterpGas(Sweep::real t, 
                          const GasProfile &gasphase, 
                          Sprog::Thermo::IdealGas &gas) const
{
    // Get the time point after the required time.
    GasProfile::const_iterator j = gasphase.upper_bound(t);
    
    if (j == gasphase.begin()) {
        // This time is before the beginning of the profile.  Return
        // the first time point.
        gas = j->second;
    } else {       
        // Get the time point before the required time.
        GasProfile::const_iterator i = j; --i;

        // Assign the conditions to this point.
        gas = i->second;
        
        // Calculate time interval between points i and j.
        real dt_pro = j->first - i->first;

        // Calculate time interval between point i and current time.
        real dt = t - i->first;

        // Calculate the intermediate gas-phase mole fractions by linear
        // interpolation of the molar concentrations.
        real dens = 0.0;
        for (unsigned int k=0; k<gas.Species()->size(); ++k) {
            real dc = (j->second.MolarConc(k) - i->second.MolarConc(k)) * dt / dt_pro;
            gas.RawData()[k] = gas.MolarConc(k) + dc;
            dens += gas.RawData()[k];
        }
        gas.Normalise();

        // Now use linear interpolation to calculate the temperature.
        real dT = (j->second.Temperature() - i->second.Temperature()) * dt / dt_pro;
        gas.SetTemperature(gas.Temperature()+dT);

        // Now set the gas density, calculated using the values above.
        gas.SetDensity(dens);
    }
}