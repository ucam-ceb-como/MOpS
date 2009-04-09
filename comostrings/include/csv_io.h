/*
  Author(s):      Matthew Celnik (msc37)
  Project:        comostrings
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    A class which facillitates reading and writing data to CSV files.  The
    class maintains a filestream for input and output.  Functions are provided
    for opening and closing the file, and for reading/writing data to the file.

    Currently the class to only write complete lines of a single type.  The
    supported types are:
        std::string
        float
        double
        long double.

  Licence:
    This file is part of the "comostrings" library.

    comostrings is free software; you can redistribute it and/or
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

#ifndef COMOSTRINGS_CSV_IO_H
#define COMOSTRINGS_CSV_IO_H

#include <string>
#include <vector>
#include <fstream>

class CSV_IO
{
public:
    // Constructors.
    CSV_IO(void);  // Default constructors.
    CSV_IO(const std::string &name, bool del = false); // Initialising constructor.

    // Destructors.
    ~CSV_IO(void); // Default destructors.

    // FILE OPEN/CLOSE.

    // Opens a csv file.  Returns -1 on failure.  Returns 0 if 
    // successful.
    int Open(
        const std::string &name, // Name of file to open.
        bool del = false         // Set del=true to delete current file contents.
        );

    // Closes the current csv file.
    void Close();

    // FILE READING.

    // Read a line of comma-separated values from the file into a vector
    // of strings.  Throws exception on error.
    void Read(std::vector<std::string> &values);

    // Reads a line of comma-separated floats from the file.
    void Read(std::vector<float> &values);

    // Reads a line of comma-separated doubles from the file.
    void Read(std::vector<double> &values);

    // Reads a line of comma-separated long doubles from the file.
    void Read(std::vector<long double> &values);

    // FILE WRITING.

    // Writes a line of comma-separated values to the file from a
    // vector of strings.  Throws exception on error.
    void Write(const std::vector<std::string> &values);

    // Writes a line of comma-separated floats to the file.
    void Write(const std::vector<float> &values);

    // Writes a line of comma-separated doubles to the file.
    void Write(const std::vector<double> &values);

    // Writes a line of comma-separated long doubles to the file.
    void Write(const std::vector<long double> &values);



private:
    std::string m_name;  // The file name.
    std::fstream m_file; // A file stream for reading/writing files

    // Reads a line from the file into a string.  Returns blank string
    // if the file is not open.
    void readLine(std::string &line);

    // Reads a line of values of type T from the file stream.  Returns
    // a blank vector if no values are read.
    template<class T>
    void readLine(std::vector<T> &values);

    // Writes a line of values of type T to the file stream.
    template<class T>
    void writeLine(const std::vector<T> &values);


};

#endif
