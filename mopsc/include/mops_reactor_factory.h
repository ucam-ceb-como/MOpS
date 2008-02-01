/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc.

  File purpose:
    The ReactorFactory class is a factory class for mops reactor
    objects.  It contains routines for creating specific reactor types.
    As there are only a limited number of possible Reactor-derived classes, 
    there is no requirement for an extensible factory class.  Instead the 
    type IDs are maintained using an enumeration Serial_ReactorType which 
    is defined in its own header.
*/

#ifndef MOPS_REACTOR_FACTORY_H
#define MOPS_REACTOR_FACTORY_H

#include "mops_params.h"
#include "mops_reactor.h"
#include "mops_reactor_type.h"

#include <string>
#include <fstream>

namespace Mops
{
class ReactorFactory
{
public:
    // REACTOR CREATION.
    // Use these routines if a generic Reactor pointer is required.

    // Creates a new reactor object of the given type.
    static Reactor *const Create(
        Serial_ReactorType type,    // Type of reactor to create.
        const Mops::Mechanism &mech // Mechanism which defines the reactor.
        );

    // Reads a Reactor object from a binary stream.  The first thing
    // read from the stream is a ReactorType to define what type
    // of reactor to read from the stream.
    static Reactor *const Read(
        std::istream &in,           // Input stream from which to read mixture.
        const Mops::Mechanism &mech // Mechanism which defines the reactor.
        );


    // REACTOR STREAM OUTPUT.

    static void Write(const Reactor &r, std::ostream &out);
};
};

#endif
