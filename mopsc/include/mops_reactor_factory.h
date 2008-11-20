/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc.
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ReactorFactory class is a factory class for mops reactor
    objects.  It contains routines for creating specific reactor types.
    As there are only a limited number of possible Reactor-derived classes, 
    there is no requirement for an extensible factory class.  Instead the 
    type IDs are maintained using an enumeration Serial_ReactorType which 
    is defined in its own header.

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
