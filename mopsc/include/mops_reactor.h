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
#include "mops_src_terms.h"

// CVODE includes.
#include "nvector/nvector_serial.h"

//#include "fortran_interface.h"

#include <istream>

namespace Mops
{
class Reactor
{
public:
    // Function pointer typedef which defines the routine used
    // to calculate external source terms in the ODE RHSs.

    // Constructors.
    Reactor(const Mechanism &mech); // Default constructor.
    Reactor(const Reactor &copy);   // Copy constructor.
    Reactor(                        // Stream-reading constructor.
        std::istream &in,           //   - Input stream.
        const Mechanism &mech       //   - Mechanism which defines the reactor.
        );

    // Destructor.
    virtual ~Reactor(void); // Default destructor.

    // Operators.
    virtual Reactor &operator=(const Reactor &rhs);


    // Enumeration of energy models.
    enum EnergyModel {ConstT, Adiabatic};

    /*
    // Enumeration of ODE solvers.
    enum ODE_Solver {CVODE_Solver, RADAU5_Solver};
    */

    // REACTOR SOLUTION.
    
    // Returns the current reaction time.
    real Time() const;

    // Sets the current reaction time.
    void SetTime(real t);

    /*
    // Initialises the reactor at the given time.
    virtual void Initialise(real time);

    // Reset the solver.  Need to do this if the the reactor
    // contents has been changed between calls to Solve().
    virtual void ResetSolver(void);

    // Solves the reactor equations up to the given time, assuming
    // that it is in future to the current time.
    virtual void Solve(real time);
    */

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
    const Mops::Mechanism *const Mech() const;

    // Returns the current mechanism.
    void SetMech(const Mops::Mechanism &mech);


    // ENERGY MODEL.

    // Returns the current energy model.
    EnergyModel EnergyEquation() const;

    // Sets the energy model.
    void SetEnergyEquation(EnergyModel model);


    // EQUATION-OF-STATE MODEL.

    // Returns true if the reactor is at constant pressure.
    bool IsConstP(void) const;

    // Sets the reactor to solve using a constant pressure assumption.
    void SetConstP(void);

    // Returns true if the reactor is at constant volume.
    bool IsConstV(void) const;

    // Sets the reactor solve using a constant volume assumption.
    void SetConstV(void);

    /*
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


    // EXTERNAL SOURCE TERMS.

    // Returns the vector of external source terms.
    SrcProfile &ExtSrcTerms(void);

    // Returns the vector of external source terms (const version).
    const SrcProfile &ExtSrcTerms(void) const;

    // Sets the external source terms.
    void SetExtSrcTerms(const SrcProfile &src);

    // Returns the external source term function.
    SrcTermFnPtr ExtSrcTermFn(void) const;

    // Sets the source term function pointer.
    void SetExtSrcTermFn(SrcTermFnPtr fn);
    */

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


    // GOVERNING EQUATIONS.

    // Returns the number of governing equations which describe
    // the reactor.  Usually this will be SpeciesCount+2, for temperature
    // and density.
    virtual unsigned int ODE_Count() const;

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

protected:
    // Reactor variables.
    real m_time;                   // The current reaction time.
    Mops::Mixture *m_mix;          // The mixture contained in the reactor.
    const Mops::Mechanism *m_mech; // The mechanism which defines 
                                   // what happens in the reactor.
    EnergyModel m_emodel;          // The energy model used to describe the reactor.
    bool m_constv; // true=constant volume model, false=constant pressure model.

    // Derived reactor properties.
    unsigned int m_neq;     // Number of equations solved.
    unsigned int m_nsp;     // Number of species in current mechanism.
    int m_iT;               // Index of temperature in solution vectors.
    int m_iDens;            // Index of density in solution vectors.
    real *m_deriv;          // Array to hold current solution derivatives.

    /*
    // ODE solution variables.
    real m_rtol, m_atol;    // Relative and absolute tolerances.
    SrcProfile m_srcterms;  // Vector of externally defined source terms on  the RHS.
    SrcTermFnPtr _srcTerms; // Source term function pointer.
    */

    // Reactors should not be defined without knowledge of a Mechanism
    // object.  Therefore the default constructor is declared as protected.
    Reactor(void);

private:
    /*
    // CVODE variables.
    void *m_odewk;        // CVODE workspace.
    N_Vector m_solvec;    // Internal solution array for CVODE interface.

    // RADAU variables.
    fvector m_rwk;          // Real workspace.
    std::vector<int> m_iwk; // Integer workspace.
    */

    // VERY IMPORTANT FOR CVODE INTEGRATION!

    /*
    // The right-hand side evaluator.  This function calculates the RHS of
    // the reactor differential equations.  CVODE uses a void* pointer to
    // allow the calling code to pass whatever information it wants to
    // the RHS function.  In this case the void* pointer should be cast
    // into a Reactor object.
    static int rhsFn_CVODE(
        double t,      // Current flow time.
        N_Vector y,    // The current solution variables.
        N_Vector ydot, // Derivatives to return.
        void* reactor  // A Reactor object (to be cast).
        );

    static void rhsFn_RADAU5(
        int   *N,             // System dimension.
        Fortran::dreal *X,    // Independent variable.
        Fortran::dreal *Y,    // Current solution.
        Fortran::dreal *F,    // Returns vector of dy/dx = F(x,y).
        Fortran::dreal *RPAR, // Real parameters.
        int *IPAR             // Integer parameters.
        );
    */

    // INITIALISATION AND DESTRUCTION.
    
    // Initialises the reactor to the default state.
    void init(void);

    // Releases all memory used by the reactor object.
    void releaseMemory(void);
};

typedef Reactor Batch;
};

#endif
