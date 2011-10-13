/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the PSR class declared in the
    mops_psr.h header file.

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

#include "mops_psr.h"
#include "mops_mixture.h"

#include <vector>
#include <cmath>
#include <stdexcept>

using namespace Mops;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
PSR::PSR(void)
{
    init();
}

// Default constructor (public, requires mechanism).
PSR::PSR(const Mops::Mechanism &mech)
: Reactor(mech)
{
    init();
    m_infHs.resize(m_neq);
}

// Copy constructor.
PSR::PSR(const Mops::PSR &copy)
{
    init();
    *this = copy;
}

// Stream-reading constructor.
PSR::PSR(std::istream &in, const Mops::Mechanism &mech)
: Reactor(mech)
{
    init();
    Deserialize(in, mech);
}

// Default destructor.
PSR::~PSR(void)
{
    releaseMemory();
}

// OPERATORS.

Mops::PSR &Mops::PSR::operator=(const Mops::PSR &rhs)
{
    if (this != &rhs) {
        // Copy base class.
        Reactor::operator=(rhs);

        // Residence time.
        m_restime = rhs.m_restime;

        // Inflow stream.
        delete m_in;
        m_in = NULL;
        if (rhs.m_in) m_in = rhs.m_in->Clone();

        // Precalculated terms.
        m_infH = rhs.m_infH;
        m_infHs.assign(rhs.m_infHs.begin(), rhs.m_infHs.end());
    }
    return *this;
}


// RESIDENCE TIME.

// Returns the reactor residence time.
real PSR::ResidenceTime(void) const
{
    return m_restime;
}

// Sets the reactor residence time.
void PSR::SetResidenceTime(real t)
{
    if (t > 0.0) {
        m_restime = t;
        m_invrt = 1.0 / m_restime;
        m_mix->AddOutflow(m_invrt, m_mech->ParticleMech());
    } else {
        throw out_of_range("Residence time must be greater "
                           "than zero (Mops, PSR::SetResidenceTime).");
    }
}


// INFLOW CONDITIONS.

// Returns the mixture which describes the inflow conditions.
const Mops::FlowStream *const PSR::Inflow(void) const {return m_in;}

// Sets the mixture which describes the inflow conditions.
void PSR::SetInflow(Mops::FlowStream &inf)
{
    // Delete current inflow stream.
    delete m_in;

    // Currently clone stream, as PSR retains ownership
    // of the inflow.
    m_in = inf.Clone();

    // Calculate inflow enthalpies to speed up
    // later calculations.
    m_infH = m_in->Mixture()->GasPhase().BulkH() / (Sprog::R * m_in->Mixture()->GasPhase().Temperature());
    m_in->Mixture()->GasPhase().Hs_RT(m_infHs);
}


// READ/WRITE/COPY FUNCTIONS.

// Creates a copy of the PSR object.
PSR* PSR::Clone() const
{
    return new PSR(*this);
}

// Writes the PSR to a binary data stream.
void PSR::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;
    
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output base class Reactor.
        this->Reactor::Serialize(out);

        // Output the residence time.
        double val = (double)m_restime;
        out.write((char*)&val, sizeof(val));

        // Output the inflow stream.
        if (m_in != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            m_in->Serialize(out);
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

    } else {
        throw invalid_argument("Output stream not ready (Mops, PSR::Serialize).");
    }
}

// Reads the Reactor data from a binary data stream.
void PSR::Deserialize(std::istream &in, const Mops::Mechanism &mech)
{
    // Clear current reactor data.
    releaseMemory();

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
                // Read the base class Reactor.
                this->Reactor::Deserialize(in, mech);

                // Read the residence time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_restime = (real)val;

                // Calculate inverse RT.
                m_invrt = m_restime;

                // Read the inflow mixture.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_in = new Mops::FlowStream(in, mech);
                }

                // Recalculate enthalpy.
                m_infH = m_in->Mixture()->GasPhase().BulkH();
                m_in->Mixture()->GasPhase().Hs(m_infHs);
                
                break;
            default:
                throw runtime_error("Reactor serialized version number "
                                    "is invalid (Mops, PSR::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Mops, PSR::Deserialize).");
    }
}

// Identifies the mixture type for serialisation.
Serial_ReactorType PSR::SerialType() const
{
    return Serial_PSR;
}


// GOVERNING EQUATIONS.

// Definition of RHS form for constant temperature energy equation.
void PSR::RHS_ConstT(real t, const real *const y,  real *ydot) const
{
    static fvector wdot;
    real wtot = 0.0;

    // Calculate molar production rates.
    wtot = m_mech->Reactions().GetMolarProdRates(y[m_iT], y[m_iDens], y, 
                                                 m_nsp, m_mix->GasPhase(), wdot);

    // Calculate mole fraction derivatives.
    for (unsigned int i=0; i!=m_nsp; ++i) {
        ydot[i] = ((wdot[i] - (y[i]*wtot)) +
                  // Inflow/Outflow term:
                  (m_in->Mixture()->GasPhase().Density() * m_invrt *
                   (m_in->Mixture()->GasPhase().MoleFraction(i) - y[i]))) / y[m_iDens];
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
        // Constant volume.
        ydot[m_iDens] = wtot + (m_invrt * (m_in->Mixture()->GasPhase().Density() - y[m_iDens]));
    } else {
        // Constant pressure.
        ydot[m_iDens] = 0.0;
    }
}

// Definition of RHS form for adiabatic energy equation.
void PSR::RHS_Adiabatic(real t, const real *const y,  real *ydot) const
{
    static fvector wdot, Hs, Cps;
    real wtot = 0.0, Cp = 0.0, H = 0.0;

    // Calculate mixture thermodynamic properties.
    m_mix->GasPhase().CalcHs_RT(y[m_iT], Hs);
    H = m_mix->GasPhase().BulkH();
    Cp = m_mix->GasPhase().CalcBulkCp_R(y[m_iT], y, m_nsp);

    // Calculate molar production rates of species (mol/m3s).
    wtot = m_mech->Reactions().GetMolarProdRates(y[m_iT], y[m_iDens], y, 
                                                 m_nsp, m_mix->GasPhase(), wdot);

    // Calculate mole fraction and temperature derivatives.
    ydot[m_iT] = 0.0;
    for (unsigned int i=0; i!=m_nsp; ++i) {
        // Mole fraction derivative.
        ydot[i] = ((wdot[i] - (y[i]*wtot)) +
                  // Inflow/Outflow term:
                  (m_in->Mixture()->GasPhase().Density() * m_invrt *
                   (m_in->Mixture()->GasPhase().MoleFraction(i) - y[i]))) / y[m_iDens];

        // Temperature derivative.
        ydot[m_iT] += Hs[i] * wdot[i];
    }

    // Complete temperature derivative (including inflow/outflow term).
    ydot[m_iT] *= - y[m_iT] / (Cp * y[m_iDens]);
    ydot[m_iT] += (m_in->Mixture()->GasPhase().Density() / (y[m_iDens] * Cp * m_restime)) *
                  ((H*y[m_iT]) - (m_infH * m_in->Mixture()->GasPhase().Temperature()));

    // Add imposed temperature gradient, if defined.
    if (m_Tfunc) ydot[m_iT] += m_Tfunc(t, y, ydot, *this);

    // Calculate density derivative.
    if (m_constv) {
        // Constant volume.
        ydot[m_iDens] = wtot + (m_invrt * (m_in->Mixture()->GasPhase().Density() - y[m_iDens]));
    } else {
        // Constant pressure (use EoS to evaluate).
        ydot[m_iDens] = - y[m_iDens] * ydot[m_iT] / y[m_iT];
    }
}


// INITIALISATION.

// Initialises the PSR to a default state, this is called
// by the constructors.
void PSR::init(void)
{
    m_restime = 0.0;
    m_in      = NULL;
    m_out     = NULL;
    m_invrt   = 0.0;
    m_infH    = 0.0;
    m_infHs.clear();
}

// Releases all object memory.
void PSR::releaseMemory(void)
{
    if (m_in != NULL) delete m_in;
    init();
}
