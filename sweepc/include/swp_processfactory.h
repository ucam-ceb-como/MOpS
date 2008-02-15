/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ProcessFactory is a factory class for sweep processes.  It
    provides routines for creating, reading and writing process
    data objects and specialised routines for inception, particle-process
    and coagulation processes.
*/

#ifndef SWEEP_PROCESS_FACTORY_H
#define SWEEP_PROCESS_FACTORY_H

#include "swp_params.h"
#include "swp_processtype.h"
#include "swp_process.h"
#include "swp_inception.h"
#include "swp_particleprocess.h"
#include "swp_coagulation.h"
#include <iostream>

namespace Sweep
{
class ProcessFactory
{
public:
    // PROCESS CREATION.

    // Creates a new process data object of the given type.
    static Process *const Create(ProcessType id);


    // STREAM INPUT.

    // Reads a process from a binary stream.  The first item read
    // is the process ID which tells the ModelFactory what type
    // of process to read.
    static Process *const Read(std::istream &in);

    // Reads an inception from a binary stream.  The first item read
    // is the inception ID which tells the ModelFactory what type
    // of inception to read.
    static Inception *const ReadInception(std::istream &in);

    // Reads a particle process from a binary stream.  The first item read
    // is the process ID which tells the ModelFactory what type
    // of process to read.
    static ParticleProcess *const ReadPartProcess(std::istream &in);

    // Reads an coagulation from a binary stream.  The first item read
    // is the coagulation ID which tells the ModelFactory what type
    // of coagulation to read.
    static Coagulation *const ReadCoag(std::istream &in);


    // STREAM OUTPUT.

    // Writes a process, along with its ID to an output stream.
    static void Write(
        const Process &proc, // Process to write.
        std::ostream &out    // Output stream.
        );
};
};

#endif
