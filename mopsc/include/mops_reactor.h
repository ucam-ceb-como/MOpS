/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The Reactor class is the base class for all types of reactor
    which can be solved with mops.  The base reactor solves a
    constant pressure batch reactor.  CVODE is used to do the
    ODE calculations.
*/

#ifndef MOPS_REACTOR_H
#define MOPS_REACTOR_H

#include "mops_params.h"
#include "mops_mixture.h"
#include "mops_mechanism.h"

// CVODE includes.
#include "nvector\nvector_serial.h"

namespace Mops
{
class Reactor
{
public:
    // Constructors.
    Reactor(void); // Default constructor.

    // Destructor.
    virtual ~Reactor(void); // Default destructor.

    // Enumeration of energy models.
    enum EnergyModel {ConstT, Adiabatic};

    // REACTOR SOLUTION.
    
    //Returns the current reaction time.
    real Time() const;

    // Initialises the reactor at the given time.
    virtual void Initialise(real time);

    // Solves the reactor equations up to the given time, assuming
    // that it is in future to the current time.
    virtual void Solve(real time);


    // REACTOR CONTENTS.

    // Returns a pointer to the mixture currently occupying
    // the reactor.
    Mops::Mixture *const Mixture() const;

    // Fills the reactor with the given mixture.
    void Fill(
        Mops::Mixture *const mix, // The mixture with which to fill the reactor.
        bool clearfirst = false   // Set to true if the reactor should clear current
                                  // mixture from memory first.
        );


    // REACTOR MECHANISM.

    // Returns the current mechanism.
    const Mops::Mechanism *const Mechanism() const;

    // Returns the current mechanism.
    void SetMechanism(const Mops::Mechanism *const mech);


    // ENERGY MODEL.

    // Returns the current energy model.
    EnergyModel EnergyEquation() const;

    // Sets the energy model.
    void SetEnergyEquation(EnergyModel model);

protected:
    // Definition of RHS function for constant temperature energy model.
    virtual void RHS_ConstT(
        real t,              // Flow time.
        const real *const y, // Solution values.
        real *ydot           // Derivatives to return.
        );

    // Definition of RHS function for adiabatic energy model.
    virtual void RHS_Adiabatic(
        real t,              // Flow time.
        const real *const y, // Solution values.
        real *ydot           // Derivatives to return.
        );

private:
    // Reactor variables.
    real m_time;                   // The current reaction time.
    Mops::Mixture *m_mix;          // The mixture contained in the reactor.
    const Mops::Mechanism *m_mech; // The mechanism which defines 
                                   // what happens in the reactor.
    EnergyModel m_emodel;          // The energy model used to describe the reactor.

    // CVODE variables.
    void *m_odewk;        // CVODE workspace.
    real m_rtol, m_atol;  // Relative and absolute tolerances.
    N_Vector m_solvec;    // Internal solution array for CVODE interface.
    int m_neq;            // Number of equations solved.
    int m_iT;             // Index of temperature in solution vectors.
    int m_iDens;          // Index of density in solution vectors.
    real *m_deriv;        // Array to hold current solution derivatives.

    // VERY IMPORTANT!

    // The right-hand side evaluator.  This function calculates the RHS of
    // the reactor differential equations.  CVODE uses a void* pointer to
    // allow the calling code to pass whatever information it wants to
    // the RHS function.  In this case the void* pointer should be cast
    // into a Reactor object.
    static int rhsFn(
        double t,      // Current flow time.
        N_Vector y,    // The current solution variables.
        N_Vector ydot, // Derivatives to return.
        void* reactor  // A Reactor object (to be cast).
        );
};
};

#endif
