/*

  Author(s):      Matthew Celnik (msc37)
  Project:        comostrings (mopssuite)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Console_io class, declared in the
    console_io.h header file.

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

#include "console_io.h"
#include <vector>
#include <iostream>
#include <string>
#include <stdio.h>

using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Console_IO::Console_IO(void)
{
    m_rowlength = 80; // Standard console is 80 characters wide.
    m_ncolumns = 6;   // Default is 6 columns
    m_formats.resize(6, Scientific); // Default scientific output.
    m_ifmt = " %10.10d |";
    m_ffmt = " %10.2f |";
    m_scifmt = " %10.3e |";
    m_strfmt = " %10.10s |";
    m_divider = "|------------|------------|------------|------------|------------|------------|\n";
    m_divint = 20;
    m_divnum = 20;
    m_autodivide = false;

#ifdef _WIN32
    // Ensure exponents only have two digits.
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
}

// Copy constructor.
Console_IO::Console_IO(const Console_IO &cio)
{
    if (&cio != this) {
        m_rowlength = cio.m_rowlength;
        m_ncolumns = cio.m_ncolumns;
        m_formats.assign(cio.m_formats.begin(), cio.m_formats.end());
        m_ifmt = cio.m_ifmt;
        m_ffmt = cio.m_ffmt;
        m_scifmt = cio.m_scifmt;
        m_strfmt = cio.m_strfmt;
        m_divider = cio.m_divider;
        m_divint = cio.m_divint;
        m_divnum = cio.m_divnum;
        m_autodivide = cio.m_autodivide;
    }
}

// Default destructor.
Console_IO::~Console_IO(void)
{
}


// TABULAR DATA.

// Returns the current column count.
unsigned int Console_IO::ColumnCount() const
{
    return m_ncolumns;
}
/*
// Sets the column count for tabular data.
void Console_IO::SetColumnCount(unsigned int ncol)
{
    m_ncolumns = ncol;
    m_formats.resize(ncol);
}
*/
// Returns the format of the ith column.
Console_IO::ColFormat Console_IO::ColumnFormat(unsigned int i) const
{
    if (i < m_ncolumns) {
        return m_formats[i];
    } else {
        return None;
    }
}

// Sets the format of the ith column.
void Console_IO::SetColumnFormat(unsigned int i, Console_IO::ColFormat format)
{
    if (i < m_ncolumns) {
        m_formats[i] = format;
    }
}

// Sets the formats of all columns.
void Console_IO::SetColumnFormats(const std::vector<ColFormat> &formats)
{
    for (unsigned int i=0; i<min(m_ncolumns, (unsigned int)formats.size()); ++i) {
        m_formats[i] = formats[i];
    }

    // If the formats vector does not conain sufficient values, then
    // set the remaining values to Scientific.
    if (m_ncolumns > formats.size()) {
        for (unsigned int i=formats.size(); i<m_ncolumns; ++i) {
            m_formats[i] = Scientific;
        }
    }
}

// Prints a data row given a vector of floats.
void Console_IO::PrintRow(const std::vector<float> &data, 
                          const std::vector<unsigned int> &mask) const
{
    // Auto divider lines.
    printAutoDivider();

    // Determine the maximum number of values that can be output.
    unsigned int i, n = min(m_ncolumns, (unsigned int)data.size());
    if (mask.size() > 0) {
        n = min(n, (unsigned int)mask.size());
    }

    // Write the values.
    printf("|");

    if (mask.size() > 0) {
        // Print values using the mask vector to get the indices.
        for (i=0; i<n; i++) {
            switch (m_formats[i]) {
                case Integer:
                    printf(m_ifmt.c_str(), (int)data[mask[i]]);
                    break;
                case Float:
                    printf(m_ffmt.c_str(), data[mask[i]]);
                    break;
                case Scientific: // Note: Scientific is the default.
                default:
                    printf(m_scifmt.c_str(), data[mask[i]]);
                    break;
            }
        }
    } else {
        // Print values starting from the beginning of 
        // the data vector.
        for (i=0; i<n; i++) {
            switch (m_formats[i]) {
                case Integer:
                    printf(m_ifmt.c_str(), (int)data[i]);
                    break;
                case Float:
                    printf(m_ffmt.c_str(), data[i]);
                    break;
                case Scientific: // Note: Scientific is the default.
                default:
                    printf(m_scifmt.c_str(), data[i]);
                    break;
            }
        }
    }

    // If we haven't printed all the columns because the data vector
    // is too small then fill the remaining columns with zeros.
    if (n < m_ncolumns) {
        for (i=n; i<m_ncolumns; i++) {
            switch (m_formats[i]) {
                case Integer:
                    printf(m_ifmt.c_str(), 0);
                    break;
                case Float:
                    printf(m_ffmt.c_str(), 0.0);
                    break;
                case Scientific: // Note: Scientific is the default.
                default:
                    printf(m_scifmt.c_str(), 0.0);
                    break;
            }
        }
    }

    // Start new line.
    printf("\n");
}

// Prints a data row given a vector of doubles.
void Console_IO::PrintRow(const std::vector<double> &data, 
                          const std::vector<unsigned int> &mask) const
{
    // Auto divider lines.
    printAutoDivider();

    // Determine the maximum number of values that can be output.
    unsigned int i, n = min(m_ncolumns, (unsigned int)data.size());
    if (mask.size() > 0) {
        n = min(n, (unsigned int)mask.size());
    }

    // Write the values.
    printf("|");

    if (mask.size() > 0) {
        // Print values using the mask vector to get the indices.
        for (i=0; i<n; i++) {
            switch (m_formats[i]) {
                case Integer:
                    printf(m_ifmt.c_str(), (int)data[mask[i]]);
                    break;
                case Float:
                    printf(m_ffmt.c_str(), data[mask[i]]);
                    break;
                case Scientific: // Note: Scientific is the default.
                default:
                    printf(m_scifmt.c_str(), data[mask[i]]);
                    break;
            }
        }
    } else {
        // Print values starting from the beginning of 
        // the data vector.
        for (i=0; i<n; i++) {
            switch (m_formats[i]) {
                case Integer:
                    printf(m_ifmt.c_str(), (int)data[i]);
                    break;
                case Float:
                    printf(m_ffmt.c_str(), data[i]);
                    break;
                case Scientific: // Note: Scientific is the default.
                default:
                    printf(m_scifmt.c_str(), data[i]);
                    break;
            }
        }
    }

    // If we haven't printed all the columns because the data vector
    // is too small then fill the remaining columns with zeros.
    if (n < m_ncolumns) {
        for (i=n; i<m_ncolumns; i++) {
            switch (m_formats[i]) {
                case Integer:
                    printf(m_ifmt.c_str(), 0);
                    break;
                case Float:
                    printf(m_ffmt.c_str(), 0.0);
                    break;
                case Scientific: // Note: Scientific is the default.
                default:
                    printf(m_scifmt.c_str(), 0.0);
                    break;
            }
        }
    }

    // Start new line.
    printf("\n");
}

// Prints a data row given a vector of long doubles.
void Console_IO::PrintRow(const std::vector<long double> &data, 
                          const std::vector<unsigned int> &mask) const
{
    // Auto divider lines.
    printAutoDivider();

    // Determine the maximum number of values that can be output.
    unsigned int i, n = min(m_ncolumns, (unsigned int)data.size());
    if (mask.size() > 0) {
        n = min(n, (unsigned int)mask.size());
    }

    // Write the values.
    printf("|");

    if (mask.size() > 0) {
        // Print values using the mask vector to get the indices.
        for (i=0; i<n; i++) {
            switch (m_formats[i]) {
                case Integer:
                    printf(m_ifmt.c_str(), (int)data[mask[i]]);
                    break;
                case Float:
                    printf(m_ffmt.c_str(), data[mask[i]]);
                    break;
                case Scientific: // Note: Scientific is the default.
                default:
                    printf(m_scifmt.c_str(), data[mask[i]]);
                    break;
            }
        }
    } else {
        // Print values starting from the beginning of 
        // the data vector.
        for (i=0; i<n; i++) {
            switch (m_formats[i]) {
                case Integer:
                    printf(m_ifmt.c_str(), (int)data[i]);
                    break;
                case Float:
                    printf(m_ffmt.c_str(), data[i]);
                    break;
                case Scientific: // Note: Scientific is the default.
                default:
                    printf(m_scifmt.c_str(), data[i]);
                    break;
            }
        }
    }

    // If we haven't printed all the columns because the data vector
    // is too small then fill the remaining columns with zeros.
    if (n < m_ncolumns) {
        for (i=n; i<m_ncolumns; i++) {
            switch (m_formats[i]) {
                case Integer:
                    printf(m_ifmt.c_str(), 0);
                    break;
                case Float:
                    printf(m_ffmt.c_str(), 0.0);
                    break;
                case Scientific: // Note: Scientific is the default.
                default:
                    printf(m_scifmt.c_str(), 0.0);
                    break;
            }
        }
    }

    // Start new line.
    printf("\n");
}

// Prints a data row given a vector of strings.
void Console_IO::PrintRow(const std::vector<std::string> &data, 
                          const std::vector<unsigned int> &mask) const
{
    // Determine the maximum number of values that can be output.
    unsigned int i, n = min(m_ncolumns, (unsigned int)data.size());
    if (mask.size() > 0) {
        n = min(n, (unsigned int)mask.size());
    }

    // Write the values.
    printf("|");

    if (mask.size() > 0) {
        // Print values using the mask vector to get the indices.
        for (i=0; i<n; i++) {
            printf(m_strfmt.c_str(), data[mask[i]].c_str());
        }
    } else {
        // Print values starting from the beginning of 
        // the data vector.
        for (i=0; i<n; i++) {
            printf(m_strfmt.c_str(), data[i].c_str());
        }
    }

    // If we haven't printed all the columns because the data vector
    // is too small then fill the remaining columns with blanks.
    if (n < m_ncolumns) {
        for (i=n; i<m_ncolumns; i++) {
            printf(m_strfmt.c_str(), "");
        }
    }

    // Start new line.
    printf("\n");
}

// Prints a dividing row, "---------------- etc", to the console.
void Console_IO::PrintDivider() const
{
    printf("%s", m_divider.c_str());
}


// AUTOMATIC DIVIDERS.

// Sets the automatic divider interval.  This is the number of
// rows printed between dividing lines.
void Console_IO::SetDividerInterval(unsigned int i)
{
    m_divint = max(i,(unsigned int)1);
}

// Turns on automatic dividers.
void Console_IO::EnableAutoDividers()
{
    m_autodivide = true;
    m_divnum = m_divint;
}

// Turns off automatic dividers.
void Console_IO::DisableAutoDividers()
{
    m_autodivide = false;
}

// Sets the auto header row to print with auto dividers.
void Console_IO::SetAutoHeader(const std::vector<std::string> &header)
{
    m_autoheader.assign(header.begin(), header.end());
}

// Prints the auto divider.
void Console_IO::printAutoDivider() const
{
    if (m_autodivide) {
        if (m_divnum == 0) {
            // Print dividing line.
            PrintDivider();
            PrintRow(m_autoheader);
            PrintDivider();
            m_divnum = m_divint-1;
        } else {
            --m_divnum;
        }
    }
}
