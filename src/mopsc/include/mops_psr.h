 /*!
  * @file   mops_psr.h
  * @author Matthew Celnik, William Menz
  * @brief  Declaration of the Pefectly Stirred Reactor (PSR) class
  *
  *   About:
  *      The PSR class derives from the Reactor class and defines a perfectly
  *      stirred reactor with multiple inflows and outflows. A constant total
  *      volumetic flowrate is assumed to pass through the reactor.
  *
  *      Streams connected to each reactor have a 'flow fraction' which
  *      represent the portion of the total flowrate that that stream makes-up.
  *      2 Apr 2013: PSRs no longer own the flowstreams. These should be
  *      administered outside the reactor class.
  *
  *   Licence:
  *      mops is free software; you can redistribute it and/or
  *      modify it under the terms of the GNU Lesser General Public License
  *      as published by the Free Software Foundation; either version 2
  *      of the License, or (at your option) any later version.
  *
  *      This program is distributed in the hope that it will be useful,
  *      but WITHOUT ANY WARRANTY; without even the implied warranty of
  *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *      GNU Lesser General Public License for more details.
  *
  *      You should have received a copy of the GNU Lesser General Public
  *      License along with this program; if not, write to the Free Software
  *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  *      02111-1307, USA.
  *
  *   Contact:
  *      Prof Markus Kraft
  *      Dept of Chemical Engineering
  *      University of Cambridge
  *      New Museums Site
  *      Pembroke Street
  *      Cambridge
  *      CB2 3RA, UK
  *
  *      Email:       mk306@cam.ac.uk
  *      Website:     http://como.cheng.cam.ac.uk
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
typedef std::vector<Mops::FlowStream*> FlowPtrVector;

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
    Mops::FlowStream *const Inflow(unsigned int i) const;

    //! Returns the inflow stream pointers
    // aab64: changed definition from 
    // Mops::FlowPtrVector Inflows() const; 
    const Mops::FlowPtrVector &Inflows() const;

    //! Returns the flow-stream which describes the outflow conditions.
    Mops::FlowStream *const Outflow(unsigned int i) const;

    //! Returns the Outflow stream pointers
    // aab64: changed definition from 
    // Mops::FlowPtrVector Outflows() const;
    const Mops::FlowPtrVector &Outflows() const;

    //! Intialise the inflow's birth processes
    void InitialiseInflows();

    //! Initialise the outflow's death processes
    void InitialiseOutflows();

    //! Sets the mixture which describes the inflow conditions.
    void SetInflow(Mops::FlowStream &inf);

    //! Sets the outflow
    void SetOutflow(Mops::FlowStream &out);

    //! Does the reactor have an outflow set?
    bool HasOutflow() const;

    //! Does the reactor have an inflow set?
    bool HasInflow() const;

    //! Set the inflow process type
    void SetInflowType(Sweep::Processes::BirthProcess::BirthType btype);

    //! Set the inflow process type
    void SetOutflowType(Sweep::Processes::DeathProcess::DeathType dtype);

    //! Normalise the particle birth/death process rates.
    void NormaliseIOProcessRates();

    //! Clear the memory associated with any flow streams
    void ClearStreamMemory();

    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the reactor object.
    PSR* Clone() const;

    // Writes the PSR to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the PSR data from a binary data stream.
    void Deserialize(
        std::istream &in,           // Input stream.
        const Mops::Mechanism &mech // Mechanism which defines reactor.
        );

    // Identifies the reactor type for serialisation.
    Serial_ReactorType SerialType() const;

protected:
    // Reactors should not be defined without knowledge of a Mechanism
    // object.  Therefore the default constructor is declared as protected.
    PSR(void);


    // GOVERNING EQUATIONS.

    //! Definition of RHS function for constant temperature energy model
    void RHS_ConstT(
        double t,              // Flow time.
        const double *const y, // Solution values.
        double *ydot           // Derivatives to return.
        ) const;

    //! Definition of RHS function for adiabatic energy model
    void RHS_Adiabatic(
        double t,              // Flow time.
        const double *const y, // Solution values.
        double *ydot           // Derivatives to return.
        ) const;

private:

    //! Initialise the birth process of the specific inflow
    void InitialiseInflow(Mops::FlowStream& inf);

    //! Initialise the birth process of the specific outf
    void InitialiseOutflow(Mops::FlowStream& outf);

    //! Helper function for the ODE solver
    void RHS_Complete(double t, const double *const y, double *ydot) const;

    //! Residence time
    double m_restime;

    //! Inflow streams
    Mops::FlowPtrVector m_inflow_ptrs;

    //! Outflow streams
    Mops::FlowPtrVector m_outflow_ptrs;

    //! Inverse residence time
    double m_invrt;

    //! Scaling factor for inflow streams' flow fractions
    double m_iscaling;

    // (Not used for the gas-phase solver)
    //! Scaling factor for the outflow streams' flow fractions
    double m_oscaling;

    // Cached inflow and outflow types
    Sweep::Processes::BirthProcess::BirthType m_default_birth;
    Sweep::Processes::DeathProcess::DeathType m_default_death;
};
};

#endif
