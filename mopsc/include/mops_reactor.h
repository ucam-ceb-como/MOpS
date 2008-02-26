/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The Reactor class is the base class for all types of reactor
    which can be solved with mops.  The base reactor solves a
    batch reactor.  The default reactor is constant pressure and
    constant temperature. CVODE is used to do the ODE calculations.
*/

#ifndef MOPS_REACTOR_H
#define MOPS_REACTOR_H

#include "mops_params.h"
#include "mops_mixture.h"
#include "mops_mechanism.h"
#include "mops_reactor_type.h"

// CVODE includes.
#include "nvector\nvector_serial.h"

#include <istream>

namespace Mops
{
class Reactor
{
public:
    // Constructors.
    Reactor(const Mechanism &mech); // Default constructor.
    Reactor(const Reactor &copy);   // Copy constructor.
    Reactor(                        // Stream-reading constructor.
        std::istream &in,           //   - Input stream.
        const Mechanism &mech       //   - Mechanism which defines the reactor.
        );

    // Destructor.
    virtual ~Reactor(void); // Default destructor.

    // Enumeration of energy models.
    enum EnergyModel {ConstT, Adiabatic};

    // REACTOR SOLUTION.
    
    // Returns the current reaction time.
    real Time() const;

    // Sets the current reaction time.
    void SetTime(real t);

    // Initialises the reactor at the given time.
    virtual void Initialise(real time);

    // Reset the solver.  Need to do this if the the reactor
    // contents has been changed between calls to Solve().
    virtual void ResetSolver(void);

    // Solves the reactor equations up to the given time, assuming
    // that it is in future to the current time.
    virtual void Solve(real time);


    // REACTOR CONTENTS.

    // Returns a pointer to the mixture currently occupying
    // the reactor.
    Mops::Mixture *const Mixture() const;

    // Fills the reactor with the given mixture.
    void Fill(
        Mops::Mixture &mix,     // The mixture with which to fill the reactor.
        bool clearfirst = false // Set to true if the reactor should clear current
                                // mixture from memory first.
        );


    // REACTOR MECHANISM.

    // Returns the current mechanism.
    const Mops::Mechanism *const Mechanism() const;

    // Returns the current mechanism.
    void SetMechanism(const Mops::Mechanism &mech);


    // ENERGY MODEL.

    // Returns the current energy model.
    EnergyModel EnergyEquation() const;

    // Sets the energy model.
    void SetEnergyEquation(EnergyModel model);


    // EQUATION-OF-STATE MODEL.

    // Sets the reactor to solve using a constant pressure assumption.
    void SetConstP(void);

    // Sets the reactor solve using a constant volume assumption.
    void SetConstV(void);


    // ERROR TOLERANCES.

    // Returns the absolute error tolerance used for ODE
    // calculations.
    real ATOL() const;

    // Sets the absolute error tolerance used for ODE
    // calculations.
    void SetATOL(real atol);

    // Returns the relative error tolerance used for ODE
    // calculations.
    real RTOL() const;

    // Sets the relative error tolerance used for ODE
    // calculations.
    void SetRTOL(real rtol);


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the reactor object.
    virtual Reactor* Clone() const;

    // Writes the reactor to a binary data stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the reactor data from a binary data stream.
    virtual void Deserialize(
        std::istream &in,           // Input stream.
        const Mops::Mechanism &mech // Mechanism which defines reactor.
        );

    // Identifies the reactor type for serialisation.
    virtual Serial_ReactorType SerialType() const;

protected:
    // Reactor variables.
    real m_time;                   // The current reaction time.
    Mops::Mixture *m_mix;          // The mixture contained in the reactor.
    const Mops::Mechanism *m_mech; // The mechanism which defines 
                                   // what happens in the reactor.
    EnergyModel m_emodel;          // The energy model used to describe the reactor.
    bool m_constv; // true=constant volume model, false=constant pressure model.

    // ODE solution variables.
    real m_rtol, m_atol;  // Relative and absolute tolerances.
    unsigned int m_neq;   // Number of equations solved.
    unsigned int m_nsp;   // Number of species in current mechanism.
    int m_iT;             // Index of temperature in solution vectors.
    int m_iDens;          // Index of density in solution vectors.
    real *m_deriv;        // Array to hold current solution derivatives.

    // Reactors should not be defined without knowledge of a Mechanism
    // object.  Therefore the default constructor is declared as protected.
    Reactor(void);


    // GOVERNING EQUATIONS.

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
    // CVODE variables.
    void *m_odewk;        // CVODE workspace.
    N_Vector m_solvec;    // Internal solution array for CVODE interface.

    // VERY IMPORTANT FOR CVODE INTEGRATION!

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


    // INITIALISATION AND DESTRUCTION.
    
    // Initialises the reactor to the default state.
    void init(void);

    // Releases all memory used by the reactor object.
    void releaseMemory(void);
};

typedef Reactor Batch;
};

#endif
