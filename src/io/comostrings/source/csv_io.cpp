/*

  Author(s):      Matthew Celnik (msc37)
  Project:        comostrings (mopssuite)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the CSV_IO class which is declared in
    the csv_io.h header file.

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

#include <iomanip>


#include "csv_io.h"
#include "string_functions.h"
#include <string>
#include <vector>
#include <fstream>

using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
CSV_IO::CSV_IO(void)
{
    m_name = "";
}

// Initialising constructor.
CSV_IO::CSV_IO(const std::string &name, bool del)
{
    // Open the file.
    Open(name, del);
}

// Default destructor.
CSV_IO::~CSV_IO(void)
{
    // Close the file, if it is not already closed.
    if (m_file.is_open()) m_file.close();
}


// FILE OPENING AND CLOSING.

// Opens a CSV file.
int CSV_IO::Open(const std::string &name, bool del)
{
    // Attempt to open the file for input and output.
    if (del) {
        // Use ios_base::trunc to delete current file contents.
        m_file.open(name.c_str(), ios_base::in | ios_base::out | ios_base::trunc);
    } else {
        // Use ios_base::app to append to current file contents.
        m_file.open(name.c_str(), ios_base::in | ios_base::out | ios_base::app);
    }

    if (m_file.good()) {
        // Save the file name.
        m_name = name;
        return 0;
    } else {
        // Failed to open file.
        return -1;
    }
}

// Close the current file.
void CSV_IO::Close()
{
    // Close the file, if it is not already closed.
    if (m_file.is_open()) m_file.close();
    // Clear file name.
    m_name = "";
}


// FILE READING.

// Reads a line into a vector of strings.
void CSV_IO::Read(std::vector<std::string> &values)
{
    string line;

    // Check that the file is open before reading.
    if (m_file.is_open() && m_file.good()) {
        // Read the next line.
        readLine(line);

        // Now split the line using commas as delimiters.  This routine
        // will first clear the values vector, so we don't need to
        // do that here.
        split(line, values, ",");
    } else {
        // The file is not open, so just return an empty vector.
        values.clear();
    }
}

// Reads a line of floats into a vector.
void CSV_IO::Read(std::vector<int> &values)
{
    readLine<int>(values);
}

// Reads a line of floats into a vector.
void CSV_IO::Read(std::vector<float> &values)
{
    readLine<float>(values);
}

// Reads a line of doubles into a vector.
void CSV_IO::Read(std::vector<double> &values)
{
    readLine<double>(values);
}

// Reads a line of long doubles into a vector.
void CSV_IO::Read(std::vector<long double> &values)
{
    readLine<long double>(values);
}


// FILE WRITING.

// Writes a line from a vector of strings.
void CSV_IO::Write(const std::vector<std::string> &values)
{
    writeLine<string>(values);
}

// Writes a line from a vector of ints.
void CSV_IO::Write(const std::vector<int> &values)
{
    writeLine<int>(values);
}

// Writes a line from a vector of floats.
void CSV_IO::Write(const std::vector<float> &values)
{
    writeLine<float>(values);
}

// Writes a line from a vector of doubles.
void CSV_IO::Write(const std::vector<double> &values)
{
    writeLine<double>(values);
}

// Writes a line from a vector of long doubles.
void CSV_IO::Write(const std::vector<long double> &values)
{
    writeLine<long double>(values);
}

//wtite tab seperated strings


// INTERNAL READ/WRITE FUNCTIONS.

// Flush CSV_IO file.
void CSV_IO::flush_file()
{
    m_file.flush();
}

// Reads a string line from the file.
void CSV_IO::readLine(std::string &line)
{
    char c;

    // Clear string.
    line.clear();

    // Check that the file is open before reading.
    if (m_file.is_open() && m_file.good()) {
        // Read the next line (terminated by return or EOF).
        m_file.get(c);
        while ((c!='\n') && (c!='\r') && !m_file.eof()) {
            line.append(&c, 1); // Append character to line.
            m_file.get(c);
        }
    }
}

// Reads a line of values into a vector.
template<class T>
void CSV_IO::readLine(std::vector<T> &values)
{
    char c, cnext;
    T val;

    // Clear the current vector.
    values.clear();

    // Check that the file is open before reading.
    if (m_file.is_open() && m_file.good() && !m_file.eof()) {
        // We need to keep checking the file stream to ensure we
        // only read one line.  The terminating characters are '\n' and
        // '\r', which are the new line and carriage return characters
        // respectively.

        cnext = m_file.peek();
        if ((cnext!='\n') && (cnext!='\r')) {
            // Read the first value from the file and add it
            // to the vector.
            m_file >> val;
            values.push_back(val);
        }

        while (!m_file.eof()) {
            // Strip out the next character, which should be the comma,
            // but we check it to see if it is end-of-line.
            m_file.get(c);
            if ((c!='\n') && (c!='\r') && !m_file.eof()) {
                // Read the next value and add it to the vector.
                m_file >> val;
                values.push_back(val);
            } else {
                break;
            }
        }
    } else {
        // The file is not open, so just return an empty vector.
        values.clear();
    }
}

// Writes values to the file.
template<class T>
void CSV_IO::writeLine(const std::vector<T> &values)
{
    // If this is not the first line, then output a newline
    // character.
    ios::streamoff i = m_file.tellp(); 
    if (i != 0) {
        m_file << '\n';
    }

    // Check that the file is open before writing.  Also check that
    // the values vector has some values to write.
    if (m_file.is_open() && m_file.good() && (values.size()>0)) {
        // Declare vector iterator and set it to beginning of values.
        typename vector<T>::const_iterator it = values.begin();
		
		//Added for better comparison with Chemkin.
		m_file.precision(15);

        // Output the comma-separated values.
        m_file << *(it++);
        while (it!=values.end()) {
            m_file << ',' << *(it++);
        }
    }
}


