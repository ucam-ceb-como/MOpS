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

#include "mops_params.h"
#include "mops_mixture.h"
#include "mops_mechanism.h"
#include <iostream>

namespace Mops
{
class Reactor; // Forward declare.

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

    // Destructor.
    ~FlowStream(void); // Default destructor.

    // Operators.
    FlowStream &operator=(const FlowStream &rhs);


    // FLOW STREAM PROPERTIES.

    // Returns a pointer to the mixture currently occupying
    // the reactor.
    Mops::Mixture *const Mixture() const;

    // Sets the flow stream mixture conditions, if they are
    // not dictated by the inflow.
    void SetConditions(const Mops::Mixture &mix);


    // FLOW STREAM CONNECTIONS.

    // Returns the reactor which is the inflow to this stream.
    Reactor *const Inflow(void) const;

    // Returns the reactor which receives the outflow 
    // of this stream.
    Reactor *const Outflow(void) const;

    // Connects the flow stream as an output of the given
    // reactor (i.e. the stream conditions become those
    // of the reactor).
    void ConnectInflow(Reactor &r);

    // Connects the flow stream as an input to the given 
    // reactor (i.e. the reactor receives the output from
    // the stream).
    void ConnectOutflow(Reactor &r);

    // Disconnects the flow-stream inflow reactor.
    void DisconnectInflow(void);

    // Disconnects the flow-stream outflow reactor.
    void DisconnectOutflow(void);


    // FLOW STREAM MECHANISM.

    // Returns the current mechanism.
    const Mops::Mechanism *const Mech() const;

    // Returns the current mechanism.
    void SetMech(const Mops::Mechanism &mech);


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the flow-stream object.
    FlowStream* Clone() const;

    // Writes the stream to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the stream data from a binary data stream.
    void Deserialize(
        std::istream &in,           // Input stream.
        const Mops::Mechanism &mech // Mechanism which defines stream.
        );

protected:
    // Flow stream inflow and outflow.
    Reactor *m_in, *m_out;

    // Stream mixture conditions.  The flow-stream object
    // may control this mixture (if no inflow is specified)
    // or it may point to reactor conditions (specified inflow).
    Mops::Mixture *m_mix;

    // Defining mechanism.
    const Mops::Mechanism *m_mech;

    // FlowStreams should not be defined without knowledge of a Mechanism
    // object.  Therefore the default constructor is declared as protected.
    FlowStream(void);

private:
    // INITIALISATION AND DESTRUCTION.
    
    // Initialises the flow-stream to the default state.
    void init(void);

    // Releases all memory used by the flow-stream object.
    void releaseMemory(void);
};
};

#endif
