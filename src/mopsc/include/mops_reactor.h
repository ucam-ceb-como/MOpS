/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Reactor class is the base class for all types of reactor
    which can be solved with mops.  The base reactor solves a
    batch reactor.  The default reactor is constant pressure and
    constant temperature. CVODE is used to do the ODE calculations.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef MOPS_REACTOR_H
#define MOPS_REACTOR_H

#include "mops_params.h"
#include "mops_mixture.h"
#include "mops_mechanism.h"
#include "mops_reactor_type.h"
#include "mops_src_terms.h"
#include "sweep.h"

namespace Mops
{
class Reactor
{
public:
    // Type definition of RHS functional form to return
    // single derivative.
    typedef real (*RHS1_FnPtr) (
        real t,                 // Flow time.
        const real *const y,    // Solution values.
        const real *const ydot, // Derivatives to return.
        const Reactor &r        // Calling reactor.
        );

    // Enumeration of energy models.  Note that an imposed
    // temperature gradient can still be defined by an 
    // external function for both of these cases.  This energy
    // model merely dictates how the energy from reaction
    // is treated (ignored or not).
    enum EnergyModel {
        ConstT,    // Constant temperature.
        Adiabatic  // No heat transfer (full heating/cooling).
    };

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


    // REACTOR TIME.
    
    // Returns the current reaction time.
    real Time() const;

    // Sets the current reaction time.
    void SetTime(real t);


    // REACTOR CONTENTS.

    // Returns a pointer to the mixture currently occupying
    // the reactor.
    const Mops::Mixture *const Mixture() const;
    Mops::Mixture *const Mixture();

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


    // IMPOSED dT/dt PROFILE.

    // Sets the function used to calculate the dT/dt in the
    // RHS functions.
    void SetTempFunc(RHS1_FnPtr fn);

    // Tells the Reactor object to use its default dT/dt
    // calculation function (from internal profile).
    void UseDefaultTempFunc(void);

    // Turns off the imposed dT/dt RHS function.
    void DisableTempFunc(void);

    // Adds a dT/dt functional to the internal profile.  This
    // automatically tells the Reactor object to use its
    // default dT/dt calculation function.
    void Add_dTdt(
        real t, // Time from which functional is valid.
        const Sweep::Maths::Functional &fun // Functional.
        );


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
        ) const;

    // Definition of RHS function for adiabatic energy model.
    virtual void RHS_Adiabatic(
        real t,              // Flow time.
        const real *const y, // Solution values.
        real *ydot           // Derivatives to return.
        ) const;

    // Definition of Jacobian evaluator function for constant
    // temperature model.
    virtual void Jacobian(
        real t,                 // Flow time.
        real *const y,          // Solution values.
        real **J,               // Jacobian Matrix J[j][i] = dFi/dYj.
        real uround             // Perturbation size parameter.
        ) const;    
    
    //! Calculates Jacobian domegai/dcj instead of d/dxj[dxi/dt], as is done above.
    void RateJacobian(
        real t,                 // Flow time.
        real *const y,          // Solution values.
        real **J,               // Jacobian Matrix J[j][i] = domegai/dcj.
        real uround             // Perturbation size parameter.
        ) const;

    //!Create the Jacobian memory space and initialise to Identity Matrix
    double** CreateJac(int n_species) const;

    //!Destroy the Jacobian memory space
    void DestroyJac(double** J, int n_species) const;

protected:
    // Reactor variables.
    real m_time;                   // The current reaction time.
    Mops::Mixture *m_mix;          // The mixture contained in the reactor.
    const Mops::Mechanism *m_mech; // The mechanism which defines 
                                   // what happens in the reactor.

    // Reactor model variables.
    EnergyModel m_emodel; // The energy model used to describe the reactor.
    bool m_constv;        // true=const. volume model, false=const. pressure model.

    // Imposed temperature profile.
    RHS1_FnPtr m_Tfunc; // External function for calculating dT/dt.
    
    // The default imposed temperature profile can be set in
    // the Reactor class by defining a map of functional against
    // time (from which functional is valid).  There is also a
    // default function (see below) for using this profile.
    std::map<real,Sweep::Maths::Functional*> m_dTdt_profile;
    
    // Derived reactor properties.
    unsigned int m_neq;     // Number of equations solved.
    unsigned int m_nsp;     // Number of species in current mechanism.
    int m_iT;               // Index of temperature in solution vectors.
    int m_iDens;            // Index of density in solution vectors.
    real *m_deriv;          // Array to hold current solution derivatives.

    // Reactors should not be defined without knowledge of a Mechanism
    // object.  Therefore the default constructor is declared as protected.
    Reactor(void);


    // IMPOSED dT/dt PROFILE.

    // Definition of RHS function for adiabatic energy model.
    static real _RHS_dTdt_profile(
        real t,                 // Flow time.
        const real *const y,    // Solution values.
        const real *const ydot, // Derivatives to return.
        const Reactor &r        // Calling reactor.
        );

private:
    // INITIALISATION AND DESTRUCTION.
    
    // Initialises the reactor to the default state.
    void init(void);

    // Releases all memory used by the reactor object.
    void releaseMemory(void);
};

typedef Reactor Batch;
} //namespace Mops

#endif
