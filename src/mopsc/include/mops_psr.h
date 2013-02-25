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
#include "mops_flow_stream.h"

#include <istream>

namespace Mops
{

// Forward-declare the flowstream class
class FlowStream;

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

    // Operators.
    PSR &operator=(const PSR &rhs);

    //! Overload of the << operator
    friend std::ostream& operator<<(
            std::ostream &os,
            const Mops::PSR &r);


    // RESIDENCE TIME.

    // Returns the reactor residence time.
    double ResidenceTime(void) const;

    // Sets the reactor residence time.
    void SetResidenceTime(double t);


    // FLOW CONDITIONS.

    //! Returns the flow-stream which describes the inflow conditions.
    Mops::FlowStream *const Inflow(void) const;

    //! Returns the flow-stream which describes the outflow conditions.
    Mops::FlowStream *const Outflow(void) const;

    //! Intialise the inflow's birth process
    void InitialiseInflow();

    //! Initialise the outflow's death process
    void InitialiseOutflow();

    //! Sets the mixture which describes the inflow conditions.
    void SetInflow(Mops::FlowStream &inf);

    //! Sets the outflow
    void SetOutflow(Mops::FlowStream &out);

    //! Does the reactor have an outflow set?
    bool HasOutflow() const {if (m_out!=NULL) return true; else return false;}

    //! Does the reactor have an inflow set?
    bool HasInflow() const {if (m_in!=NULL) return true; else return false;}


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
        double t,              // Flow time.
        const double *const y, // Solution values.
        double *ydot           // Derivatives to return.
        ) const;

    // Definition of RHS function for adiabatic energy model.
    virtual void RHS_Adiabatic(
        double t,              // Flow time.
        const double *const y, // Solution values.
        double *ydot           // Derivatives to return.
        ) const;

private:
    // PSR variables.
    double m_restime; // Residence time.

    // PSR inflow and outflow.
    Mops::FlowStream *m_in;  // Inflow stream.
    Mops::FlowStream *m_out; // Outflow stream.

    // Precalculated values.
    double m_invrt;    // Inverse residence time.
    double m_infH;     // Inflow enthalpy.
    fvector m_infHs; // Inflow species enthalpies.

    // INITIALISATION AND DESTRUCTION.
    
    // Initialises the reactor to the default state.
    void init(void);

    // Releases all memory used by the reactor object.
    void releaseMemory(void);

    //! Returns the bulk enthalpy of the inflow
    double InflowBulkEnthalpy() const;
};
};

#endif
