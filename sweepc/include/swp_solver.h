/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The Solver class defines the stochastic stepping algorithm and
    holds information about the progress of the solution.
*/

#ifndef SWEEP_SOLVER_H
#define SWEEP_SOLVER_H

#include "swp_mechanism.h"
#include "swp_cell.h"
#include "sprog.h"
#include <vector>
#include <map>

namespace Sweep
{
class Solver
{
public:
    typedef std::map<real,Sprog::Thermo::IdealGas> GasProfile;

    // Constructors.
    Solver(void); // Default constructor.

    // Destructor.
    ~Solver(void);

    // Performs stochastic stepping algorithm up to specified stop time using
    // the given mechanism to define the stochastic processes.  Updates given
    // system accordingly.  On error returns <0, otherwise returns 0.
    int Run(
        real &t,        // Simulation start time.  Will return the stop time.
        real tstop,     // Stop time for simulation.
        Cell &sys,      // System to solve.
        Mechanism &mech // Mechanism to use to solve system.
        );

    // Performs stochastic stepping algorithm up to specified stop time using
    // the given mechanism to define the stochastic processes.  Updates given
    // system accordingly.  On error returns <0, otherwise returns 0.  In this
    // flavour the gas-phase chemistry is interpolated from a vector of
    // IdealGas objects rather than being taken from the given system object.
    // However, the particles in the system object are updated accordingly.
    int Run(
        real &t,        // Simulation start time.  Will return the stop time.
        real tstop,     // Stop time for simulation.
        const GasProfile &gasphase, // Gas-phase profile.
        Cell &sys,      // System to solve.
        Mechanism &mech // Mechanism to use to solve system.
        );

protected:
    // Numerical parameters.
    real m_tstop;      // Simulation stop time.
    real m_maxdt;      // Maximum allowed time step size.
    real m_splitratio; // Parameter defining number of LPDA updates per particle events.

    // Process counters.  These should be moved to
    // the Mechanism class.
    std::vector<unsigned int> m_processcounter, m_ficticiouscounter; 

    // TIME STEPPING ROUTINES.

    // Calculates the splitting end time after which all particles
    // shallbe updated using LPDA.
    real calcSplitTime(
        real t,         // Current time.
        real tstop,     // Stop time.
        real jrate,     // Sum of all jump process rates.
        unsigned int n, // Number of particles.
        real maxdt      // Maximum time step.
        ) const;
    
    // Performs a single stochastic event on the ensemble from the given
    // mechanism.  Returns the length of the time step if successful,
    // otherwise returns negative.
    real timeStep(
        real t,                // Current solution time.
        Cell &sys,             // System to update.
        const Mechanism &mech, // Mechanism to use.
        const fvector &rates,  // Current process rates as an array.
        real jrate             // The total jump rate (non-deferred processes).
        );

    // Selects a process using a DIV algorithm and the process rates
    // as weights.
    int chooseProcess(const fvector &rates);


    // HELPER FUNCTIONS.

    // Uses linear interpolation to return the chemical conditions
    // at a given time using a profile of Idealgas objects.
    void linInterpGas(
        real t,                      // Time.
        const GasProfile &gasphase,  // Gas-phase profile.
        Sprog::Thermo::IdealGas &gas // Output gas conditions.
        ) const;
};
};

#endif
