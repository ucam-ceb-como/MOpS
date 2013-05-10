/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the FlowStream class declared in the
    mops_flow_stream.h header file.

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

#include "mops_flow_stream.h"
#include "mops_reactor.h"
#include <stdexcept>

using namespace Mops;

// CONSTRUCTORS AND DESTRUCTORS.

// Initialising constructor.
Mops::FlowStream::FlowStream(const Mops::Mechanism &mech)
: m_in(NULL),
  m_out(NULL),
  m_mix(NULL),
  m_mech(&mech),
  m_flow_frac(1.0)
{
    // Initialise the mixture.
    m_mix = new Mops::Mixture(mech.ParticleMech());
}

// Copy constructor.
Mops::FlowStream::FlowStream(const Mops::FlowStream &copy)
: m_in(NULL),
  m_out(NULL),
  m_mix(NULL),
  m_mech(copy.m_mech),
  m_flow_frac(copy.m_flow_frac)
{*this = copy;}

// Stream-reading constructor.
Mops::FlowStream::FlowStream(
        std::istream &in,
        const Mops::Mechanism &mech)
: m_in(NULL),
  m_out(NULL),
  m_mix(NULL),
  m_mech(&mech),
  m_flow_frac(1.0)
{
    Deserialize(in, mech);
}

// Default destructor.
Mops::FlowStream::~FlowStream()
{
    if (!m_in) delete m_mix;
}

// OPERATORS.

// Assignment operator.
Mops::FlowStream &Mops::FlowStream::operator=(const Mops::FlowStream &rhs)
{
    if (this != &rhs) {
        m_in = rhs.m_in;
        m_out = rhs.m_out;
        if (m_in) {
            // If stream inflow is defined then this flow-stream
            // does not own the mixture.
            m_mix = rhs.m_mix;
        } else {
            // The flow stream does own the mixture, so create a copy.
            m_mix = rhs.m_mix->Clone();
        }
        m_mech = rhs.m_mech;
    }
    return *this;
}

// FLOW STREAM PROPERTIES.

// Sets the flow stream mixture conditions, if they are
// not dictated by the inflow.
void Mops::FlowStream::SetConditions(const Mops::Mixture &mix)
{
    if (!m_in) {
        // Only set conditions if the stream inflow is
        // undefined, as this means that the flow-stream
        // owns the mixture.
        delete m_mix;
        m_mix = mix.Clone();
    }
}

// FLOW STREAM CONNECTIONS.

void Mops::FlowStream::ConnectInflow(Mops::PSR &r)
{
    // If the inflow is currently undefined then we 
    // need to delete the mixture memory.
    if (!m_in) {
        delete m_mix;
        m_mix = NULL;
    }
    // Set the inflow and store pointer to inflow mixture.
    m_in  = &r;
    m_mix = r.Mixture();
}


void Mops::FlowStream::ConnectOutflow(Mops::PSR &r)
{
    m_out = &r;
}

// READ/WRITE/COPY FUNCTIONS.

// Creates a copy of the flow-stream object.
Mops::FlowStream* Mops::FlowStream::Clone() const
{
    return new FlowStream(*this);
}

// Writes the stream to a binary data stream.
void Mops::FlowStream::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;
    
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output whether or not the flow-stream owns
        // its mixture.  Serialize the mixture if it
        // does.
        if (!m_in) {
             out.write((char*)&trueval, sizeof(trueval));
             m_mix->Serialize(out);
       } else {
            out.write((char*)&falseval, sizeof(falseval));
       }
    } else {
        throw std::invalid_argument("Output stream not ready (Mops, FlowStream::Serialize).");
    }
}

// Reads the stream data from a binary data stream.
void Mops::FlowStream::Deserialize(std::istream &in, const Mops::Mechanism &mech)
{

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        const unsigned int trueval  = 1;

        switch (version) {
            case 0:
                // Read if the mixture was serialized.          
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n==trueval) {
                    // Read the mixture.
                    m_mix = new Mops::Mixture(in, mech.ParticleMech());
                } else {
                    // Define a default mixture.
                    m_mix = new Mops::Mixture(mech.ParticleMech());
                }
                break;
            default:
                throw std::runtime_error("Reactor serialized version number "
                                    "is invalid (Mops, FlowStream::Deserialize).");
        }
    } else {
        throw std::invalid_argument("Input stream not ready "
                               "(Mops, FlowStream::Deserialize).");
    }
}

