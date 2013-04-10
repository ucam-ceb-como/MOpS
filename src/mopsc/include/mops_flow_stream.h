/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The FlowStream class defines an object which connects
    to the input or output of a reactor.  If the stream is
    connected to an output, then its properties are those of the
    reactor to which it is connected, otherwise the properties
    are defined in the FlowStream object.

    No reactions occur in a flow stream.

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

#ifndef MOPS_FLOW_STREAM_H
#define MOPS_FLOW_STREAM_H

#include "mops_psr.h"
#include "mops_mixture.h"
#include <iostream>

namespace Mops
{
class PSR; // Forward declare.

class FlowStream
{
public:
    // Constructors.
    FlowStream(const Mechanism &mech);  // Default constructor.
    FlowStream(const FlowStream &copy); // Copy constructor.
    FlowStream(               // Stream-reading constructor.
        std::istream &in,     //   - Input stream.
        const Mechanism &mech //   - Mechanism which defines the stream.
        );

    //! Destructor.
    ~FlowStream(void); // Default destructor.

    // Operators.
    FlowStream &operator=(const FlowStream &rhs);


    // FLOW STREAM PROPERTIES.

    //! Returns a pointer to the Mixture in the flowstream
    Mops::Mixture *const Mixture() const {return m_mix;}

    //! Sets the Mixture in the flowstream
    void SetConditions(const Mops::Mixture &mix);

    // FLOW STREAM CONNECTIONS.

    //! Returns the reactor which is the inflow to this stream
    PSR *const Inflow() const {return m_in;}

    //! Returns the reactor which receives the outflow of this stream
    PSR *const Outflow() const {return m_out;}

    //! Connects the inflow of the stream to a reactor outflow
    void ConnectInflow(PSR &r);

    //! Connects the ouflow of the stream to a reactor inflow
    void ConnectOutflow(PSR &r);

    //! Has an inflow? Keep different from PSR one for clarity
    bool HasReacInflow() const {if (m_in!=NULL) return true; else return false;}

    //! Has an outflow? Keep different from PSR one for clarity
    bool HasReacOutflow() const {if (m_out!=NULL) return true; else return false;}

    //! Set the flow rate
    void SetFlowFraction(double ff) {m_flow_frac = ff;}

    //! Get the flow rate
    double GetFlowFraction() const {return m_flow_frac;}

    //! Returns the current mechanism.
    const Mops::Mechanism *const Mech() const {return m_mech;}


    // READ/WRITE/COPY FUNCTIONS.

    //! Creates a copy of the flow-stream object
    FlowStream* Clone() const;

    //! Writes the stream to a binary data stream
    void Serialize(std::ostream &out) const;

    //! Reads the stream data from a binary data stream
    void Deserialize(
        std::istream &in,           // Input stream.
        const Mops::Mechanism &mech // Mechanism which defines stream.
        );


private:
    //! Private default constructor
    FlowStream();
    // Flow stream inflow and outflow.
    Mops::PSR *m_in, *m_out;

    // Stream mixture conditions.  The flow-stream object
    // may control this mixture (if no inflow is specified)
    // or it may point to reactor conditions (specified inflow).
    Mops::Mixture *m_mix;

    // Defining mechanism.
    const Mops::Mechanism *m_mech;

    //! The fraction of the total reactor flowrate that this stream represents
    double m_flow_frac;
};
};

#endif
