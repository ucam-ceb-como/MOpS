/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc.

  File purpose:
    The Reactor_IO class is used to read/write Reactor objects from binary
    data files.  It also includes routines to convert these data files into
    comma-separated files, thereby provided some post-processing facilities
    for mops.
*/

#ifndef MOPS_REACTOR_IO_H
#define MOPS_REACTOR_IO_H

#include "mops_params.h"
#include "mops_reactor.h"
#include "mops_mechanism.h"

#include <string>
#include <fstream>

namespace Mops
{
class Reactor_IO
{
public:
    // Constructors.
    Reactor_IO(void); // Default constructor.

    // Destructors.
    virtual ~Reactor_IO(void); // Default destructor.

    // Opens a file.
    void Open(std::string name);

    // Closes the current file.
    void Close();

    // Writes a reactor object to the file.
    void Write(const Reactor &r);

    // Post-processes a binary data file to produce a series
    // of CSV files which describe the reactor time sequence therein.
    void ConvertToCSV(const std::string &filename, const Mechanism &mech);
private:
    std::string m_filename; // File name of current output file.
    fstream m_file; // Current file stream.
};
};

#endif