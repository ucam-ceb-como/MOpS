/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ProcessFactory is a factory class for sweep processes.  It
    provides routines for creating, reading and writing process
    data objects and specialised routines for inception, particle-process
    and coagulation processes.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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

#ifndef SWEEP_PROCESS_FACTORY_H
#define SWEEP_PROCESS_FACTORY_H

#include "swp_params.h"
#include "swp_process_type.h"
#include "swp_process.h"
#include "swp_inception.h"
#include "swp_particle_process.h"
#include "swp_coagulation.h"
#include "swp_fragmentation.h"
#include "swp_death_process.h"
#include <iostream>

namespace Sweep
{
namespace Processes
{
class ProcessFactory
{
public:
    // Reads an inception from a binary stream.  The first item read
    // is the inception ID which tells the ModelFactory what type
    // of inception to read.
    static Inception *const ReadInception(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

    // Reads a particle process from a binary stream.  The first item read
    // is the process ID which tells the ModelFactory what type
    // of process to read.
    static ParticleProcess *const ReadPartProcess(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

    // Reads an coagulation from a binary stream.  The first item read
    // is the coagulation ID which tells the ModelFactory what type
    // of coagulation to read.
    static Coagulation *const ReadCoag(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

    // Reads an coagulation from a binary stream.  The first item read
    // is the coagulation ID which tells the ModelFactory what type
    // of coagulation to read.
    static Fragmentation *const ReadFrag(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

    // STREAM OUTPUT.

    // Writes a process, along with its ID to an output stream.
    static void Write(
        const Process &proc, // Process to write.
        std::ostream &out    // Output stream.
        );
};
};
};

#endif
