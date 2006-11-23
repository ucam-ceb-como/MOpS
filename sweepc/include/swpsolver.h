/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The solver class defines the stochastic stepping algorithm and
    holds information about the progress of the solution.
*/

#ifndef SWEEP_SOLVER_H
#define SWEEP_SOLVER_H

#include "swpmechanism.h"
#include "swpsystem.h"
#include <vector>

namespace Sweep
{
class Solver
{
public:
    vector<unsigned int> m_processcounter, m_ficticiouscounter; 
protected:
    real m_maxdt;      // Maximum allowed time step size.
    real m_splitratio; // Parameter defining number of LPDA updates per particle events.
public:
    Solver(void);
    ~Solver(void);
    /* Initialises the solver. */
    void Initialise(void);
public:
    /* Performs stochastic stepping algorithm up to specified stop time using
       the given mechanism to define the stochastic processes.  Updates given
       system accordingly.  On error returns <0, otherwise returns 0.*/
    int Run(real *t, const real tstop, System &sys, Mechanism &mech);
protected:
    /* Performs a single stochastic event on the ensemble from the given
       mechanism.  Returns the length of the time step if successful,
       otherwise returns negative. */
    real TimeStep(const real t,               // Current solution time.
                  System &sys,                // System to update.
                  Mechanism &mech,            // Mechanism to use.
                  const vector<real> &rates); // Current process rates as an array.
    /* Selects a process using a DIV algorithm and the process rates
       as weights. */
    int ChooseProcess(const vector<real> &rates, // Process rates.
                      const unsigned int n);     // Process count.
};
};

#endif