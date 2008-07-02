/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The PSR class derives from the Reactor class and defines a 
    perfectly-stirred tank reactor with one inflow and one
    outflow with identical volumetric flow rates.

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

#ifndef MOPS_PSR_H
#define MOPS_PSR_H

#include "mops_params.h"
#include "mops_mechanism.h"
#include "mops_reactor.h"
#include "mops_reactor_type.h"

#include <istream>

namespace Mops
{
class PSR : public Reactor
{
public:
    // Constructors.
    PSR(const Mops::Mechanism &mech); // Default constructor.
    PSR(const PSR &copy);             // Copy constructor.
    PSR(                              // Stream-reading constructor.
        std::istream &in,             //   - Input stream.
        const Mops::Mechanism &mech   //   - Mechanism which defines the reactor.
        );

    // Destructor.
    ~PSR(void); // Default destructor.


    // RESIDENCE TIME.

    // Returns the reactor residence time.
    real ResidenceTime(void) const;

    // Sets the reactor residence time.
    void SetResidenceTime(real t);


    // INFLOW CONDITIONS.

    // Returns the mixture which describes the inflow conditions.
    const Mops::Mixture *const Inflow(void) const;

    // Sets the mixture which describes the inflow conditions.
    void SetInflow(Mops::Mixture &inf);


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the reactor object.
    virtual PSR* Clone() const;

    // Writes the PSR to a binary data stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the PSR data from a binary data stream.
    virtual void Deserialize(
        std::istream &in,           // Input stream.
        const Mops::Mechanism &mech // Mechanism which defines reactor.
        );

    // Identifies the reactor type for serialisation.
    virtual Serial_ReactorType SerialType() const;

protected:
    // Reactors should not be defined without knowledge of a Mechanism
    // object.  Therefore the default constructor is declared as protected.
    PSR(void);


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
    // PSR variables.
    real m_restime;          // Residence time.
    Mops::Mixture *m_inflow; // Inflow mixture.

    // Precalculated values.
    real m_invrt;    // Inverse residence time.
    real m_infH;     // Inflow enthalpy.
    fvector m_infHs; // Inflow species enthalpies.

    // INITIALISATION AND DESTRUCTION.
    
    // Initialises the reactor to the default state.
    void init(void);

    // Releases all memory used by the reactor object.
    void releaseMemory(void);
};
};

#endif
