/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Reactor class declared in the
    mops_reactor.h header file.

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

#include "mops_reactor.h"
#include "mops_mixture.h"
#include "swp_ensemble.h"

#include <vector>
#include <list>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <memory>

using namespace Mops;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Reactor::Reactor(void)
{
    init();
}

// Default constructor (public, requires mechanism).
Reactor::Reactor(const Mops::Mechanism &mech)
{
    init();
    SetMech(mech);
}

// Copy constructor.
Reactor::Reactor(const Mops::Reactor &copy)
{
    init();
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
Reactor::Reactor(std::istream &in, const Mops::Mechanism &mech)
{
    init();
    Deserialize(in, mech);
}

// Default destructor.
Reactor::~Reactor(void)
{
    releaseMemory();
}


// OPERATORS.

// Assignment operator.
Reactor &Reactor::operator=(const Mops::Reactor &rhs)
{
    if (this != &rhs) {
        releaseMemory();

        // Copy the reactor properties.
        m_time    = rhs.m_time;

        // Clone the mixture, if it exists in rhs
        if(rhs.m_mix)
            m_mix     = rhs.m_mix->Clone(); // Need to clone mixture!
        else
            m_mix = NULL;

        m_emodel  = rhs.m_emodel;
        m_nsp     = rhs.m_nsp;
        m_neq     = rhs.m_neq;
        m_iT      = rhs.m_iT;
        m_iDens   = rhs.m_iDens;

        // Initialise the reactor with the mechanism.
        SetMech(*rhs.m_mech);

        // Copy ODE workspace, incl. derivatives.
        memcpy(m_deriv, rhs.m_deriv, sizeof(real)*m_neq);
    }
    return *this;
}


// REACTOR SOLUTION.

// Returns the current reactor time.
real Reactor::Time() const
{
    return m_time;
}

// Sets the current reactor time.
void Reactor::SetTime(real t)
{
    m_time = t;
}

// REACTOR CONTENTS.

// Returns a pointer to the current reactor contents.
const Mops::Mixture *const Reactor::Mixture() const
{
    return m_mix;
}

Mops::Mixture *const Reactor::Mixture()
{
    return m_mix;
}

// Fills the reactor with the given mixture.  If the clearfirst flag
// is set then the current mixture is first deleted from memory.
void Reactor::Fill(Mops::Mixture &mix, bool clearfirst)
{
    if (m_mix!=NULL) {
        // First delete current mixture from memory.
        delete m_mix;
        m_mix = NULL;
    }

    // Set new mixture (also take ownership and responibility for its deletion).
    m_mix = &mix;

    // Ensure that the reactor and mixture are using the same
    // mechanism.
    m_mix->SetSpecies(m_mech->Species());
}


// REACTOR MECHANISM.

// Returns the current mechanism.
const Mops::Mechanism *const Reactor::Mech() const
{
    return m_mech;
}

// Sets the reactor mechanism.
void Reactor::SetMech(const Mops::Mechanism &mech)
{
    m_mech  = &mech;

    // Set up vector indices.
    m_nsp   = m_mech->SpeciesCount();
    m_neq   = m_nsp + 2; // Equations of state; these are not differential equations. See gpc_mech.h m_rxns
    m_iT    = m_nsp;
    m_iDens = m_iT + 1;

    // Allocate the derivative array.
    if (m_deriv != NULL) delete [] m_deriv;
    m_deriv = new real[m_neq];
}


// ENERGY MODEL.

// Returns the current energy model.
Reactor::EnergyModel Reactor::EnergyEquation() const
{
    return m_emodel;
}

// Sets the energy model.
void Reactor::SetEnergyEquation(Reactor::EnergyModel model)
{
    m_emodel = model;
}


// EQUATION-OF-STATE MODEL.

// Returns true if the reactor is at constant pressure.
bool Reactor::IsConstP(void) const {return !m_constv;};

// Sets the reactor to solve using a constant pressure assumption.
void Reactor::SetConstP(void)
{
    m_constv = false;
}

// Returns true if the reactor is at constant volume.
bool Reactor::IsConstV(void) const {return m_constv;};

// Sets the reactor solve using a constant volume assumption.
void Reactor::SetConstV(void)
{
    m_constv = true;
}


// IMPOSED dT/dt PROFILE.

// Sets the function used to calculate the dT/dt in the
// RHS functions.
void Reactor::SetTempFunc(RHS1_FnPtr fn) {m_Tfunc = fn;}

// Tells the Reactor object to use its default dT/dt
// calculation function (from internal profile).
void Reactor::UseDefaultTempFunc(void) {m_Tfunc = &_RHS_dTdt_profile;}

// Turns off the imposed dT/dt RHS function.
void Reactor::DisableTempFunc(void) {m_Tfunc = NULL;}

// Adds a dT/dt functional to the internal profile.  This
// automatically tells the Reactor object to use its
// default dT/dt calculation function.
void Reactor::Add_dTdt(real t, const Sweep::Maths::Functional &fun)
{
    UseDefaultTempFunc();
    m_dTdt_profile[t] = fun.Clone();
}

// Definition of RHS function for adiabatic energy model.
real Reactor::_RHS_dTdt_profile(real t, const real *const y,
                                const real *const ydot, 
                                const Reactor &r)
{
    // Locate the first time point after t, or end() if
    // the profile is not long enough.
    map<real,Sweep::Maths::Functional*>::const_reverse_iterator i;
    for (i=r.m_dTdt_profile.rbegin(); i!=r.m_dTdt_profile.rend(); ++i) {
        if (i->first <= t) {
            break;
        }
    }

    // Evaluate the function to calculate dT/dt.
    if (i != r.m_dTdt_profile.rend()) {
        return i->second->Evaluate(t);
    } else {
        // t is before the start of the profile.    
        return 0.0;
    }
}


// READ/WRITE/COPY FUNCTIONS.

// Creates a copy of the reactor object.
Reactor* Reactor::Clone() const
{
    return new Reactor(*this);
}

// Writes the reactor to a binary data stream.
void Reactor::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;
    
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output the time.
        double val = (double)m_time;
        out.write((char*)&val, sizeof(val));

        // Output the mixture.
        if (m_mix != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            m_mix->Serialize(out);
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Output the energy model.
        unsigned int n = (unsigned int)m_emodel;
        out.write((char*)&n, sizeof(n));

        // Output the EoS model.
        if (m_constv) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Output equation count + special indices.
        n = (unsigned int)m_nsp;
        out.write((char*)&n, sizeof(n));
        n = (unsigned int)m_neq;
        out.write((char*)&n, sizeof(n));
        n = (unsigned int)m_iT;
        out.write((char*)&n, sizeof(n));
        n = (unsigned int)m_iDens;
        out.write((char*)&n, sizeof(n));

        // Output derivatives array.
        if (m_deriv != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            for (unsigned int i=0; i<m_neq; i++) {
                val = (double)m_deriv[i];
                out.write((char*)&val, sizeof(val));
            }
        } else {
            out.write((char*)falseval, sizeof(falseval));
        }
    } else {
        throw invalid_argument("Output stream not ready (Mops, Reactor::Serialize).");
    }
}

// Reads the Reactor data from a binary data stream.
void Reactor::Deserialize(std::istream &in, const Mops::Mechanism &mech)
{
    // Clear current reactor data.
    if (m_mix != NULL) delete m_mix;
    if (m_deriv != NULL) delete [] m_deriv;
    init();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;
        unsigned int n = 0;

        switch (version) {
            case 0:
                // Read the time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_time = (real)val;

                // Read the mixture.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_mix = new Mops::Mixture(in, mech.ParticleMech());
                }

                // Read the energy model.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_emodel = (EnergyModel)n;

                // Read the EoS model.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_constv = (n==1);

                // Read equation count + special indices.
                in.read(reinterpret_cast<char*>(&m_nsp), sizeof(m_nsp));
                in.read(reinterpret_cast<char*>(&m_neq), sizeof(m_neq));
                in.read(reinterpret_cast<char*>(&m_iT), sizeof(m_iT));
                in.read(reinterpret_cast<char*>(&m_iDens), sizeof(m_iDens));

                // Store the mechanism.
                SetMech(mech);

                // Read derivatives array.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_deriv = new real[m_neq];
                    for (unsigned int i=0; i<m_neq; i++) {
                        in.read(reinterpret_cast<char*>(&val), sizeof(val));
                        m_deriv[i] = (real)val;
                    }
                }
                
                break;
            default:
                throw runtime_error("Reactor serialized version number "
                                    "is invalid (Mops, Reactor::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Mops, Reactor::Deserialize).");
    }
}

// Identifies the mixture type for serialisation.
Serial_ReactorType Reactor::SerialType() const
{
    return Serial_Reactor;
}


// RHS FUNCTION AND GOVERNING EQUATIONS.

// Returns the number of governing equations which describe
// the reactor.  Usually this will be SpeciesCount+2, for temperature
// and density.
unsigned int Reactor::ODE_Count() const
{
    return m_neq;
}

// Definition of RHS form for constant temperature energy equation.
void Reactor::RHS_ConstT(real t, const real *const y,  real *ydot) const
{
    static fvector wdot;
    real wtot = 0.0;

    // Calculate molar production rates.
    wtot = m_mech->Reactions().GetMolarProdRates(y[m_iT], y[m_iDens], y, 
                                                 m_nsp, *m_mix, wdot);

    // Calculate mole fraction derivatives.
    for (unsigned int i=0; i!=m_neq-2; ++i) {
        ydot[i] = ((wdot[i] - (y[i]*wtot)) / y[m_iDens]);
    }

    // Temperature derivative.
    if (m_Tfunc) {
        // Add imposed temperature gradient, if defined.
        ydot[m_iT] = m_Tfunc(t, y, ydot, *this);
    } else {
        // Constant temperature.
        ydot[m_iT] = 0.0;
    }

    // Density derivative.
    if (m_constv) {
        ydot[m_iDens] = wtot; // Constant volume.
    } else {
        ydot[m_iDens] = 0.0;  // Constant pressure.
    }
}

// Definition of RHS form for adiabatic energy equation.
void Reactor::RHS_Adiabatic(real t, const real *const y,  real *ydot) const
{
    static fvector wdot, Hs;
    real wtot = 0.0, Cp = 0.0;

    // Calculate mixture thermodynamic properties.
    m_mix->CalcHs_RT(y[m_iT], Hs);
    Cp = m_mix->ThermoInterface::CalcBulkCp_R(y[m_iT], y, m_nsp);
    
    // Calculate molar production rates.
    wtot = m_mech->Reactions().GetMolarProdRates(y[m_iT], y[m_iDens], y, 
                                                 m_nsp, *m_mix, wdot);


    // Calculate mole fraction and temperature derivatives.
    ydot[m_iT] = 0.0;
    for (unsigned int i=0; i!=m_nsp; ++i) {
        ydot[i] = ((wdot[i] - (y[i]*wtot)) / y[m_iDens]);
        ydot[m_iT] += wdot[i] * Hs[i];
    }

    // Complete temperature derivative.
    ydot[m_iT] *= - y[m_iT] / (Cp * y[m_iDens]);

    // Add imposed temperature gradient, if defined.
    if (m_Tfunc) ydot[m_iT] += m_Tfunc(t, y, ydot, *this);

    // Calculate density derivative.
    if (m_constv) {
        // Constant volume.
        ydot[m_iDens] = wtot;
    } else {
        // Constant pressure (Use EoS to calculate).
        ydot[m_iDens] = - (y[m_iDens] * ydot[m_iT] / y[m_iT]);
    }
}

// Definition of Jacobian evaluator function for constant
// temperature model.
void Reactor::Jacobian(real t, real *const y, 
                       real **J,
                       real uround) const
{
    m_mech->Reactions().CalcJacobian(y[m_iT], y[m_iDens], y, 
                                     m_nsp, *m_mix, uround, J, 
                                     m_constv, m_emodel==ConstT);
}

/*!
@param[in]          t       Time step
@param[in]          y       solution vector with mole fractions and density and temperature
@param[in, out]     J       Jacobian array
@param[in]          uround  The value of the perturbation factor for finite differencing.
*/
void Reactor::RateJacobian(real t, real *const y, 
                       real **J,
                       real uround) const
{
    m_mech->Reactions().RateJacobian(y[m_iT], y[m_iDens], y, 
                                     m_nsp, *m_mix, uround, J, 
                                     m_constv, m_emodel==ConstT);
}


/*!
@param[in]         n_species    number of species in the reaction
@return            J            Jacobian array, initialised to zero
*/
double** Reactor::CreateJac( int n_species) const{
    
    double** J;
    int Dim = n_species + 2; //Temperature and Pressure
    J = new double*[Dim];
    for (int i = 0; i < Dim; i++){
        J[i] = new double[Dim];
    }
    for (int i = 0; i < Dim; i++){
        for (int j = 0; j < Dim; j++){
            if (i == j)
                 J[i][j] = 1.0;
            else
                J[i][j] = 0.0;
        }
    }
    return J;
}

/*!
@param[in]         n_species   number of species in the reaction
*/

void Reactor::DestroyJac(double** J, int n_species) const{

    for (int i = 0; i < n_species; i++){
        delete [] J[i];
    }
    delete J;
}

// INITIALISATION.

// Initialises the reactor to a default state, this is called
// by the constructors.
void Reactor::init(void)
{
    // Reactor variables.
    m_time    = 0.0;
    m_mix     = NULL;
    m_mech    = NULL;
    // Reactor model variables.
    m_emodel  = ConstT;
    m_constv  = false;
    m_Tfunc   = NULL;
    // Derived reactor properties.
    m_neq     = 0;
    m_nsp     = 0;
    m_iT      = -1;
    m_iDens   = -1;
    m_deriv   = NULL;
}

// Releases all object memory.
void Reactor::releaseMemory(void)
{
    if (m_mix != NULL) delete m_mix;
    m_mix = NULL;
    if (m_deriv != NULL) delete [] m_deriv;
    m_deriv = NULL;
}
