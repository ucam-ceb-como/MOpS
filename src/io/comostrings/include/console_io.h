/*

  Author(s):      Matthew Celnik (msc37)
  Project:        comostrings (mopssuite)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    A class which facillitates reading and writing formatted data to the
    console.  This class provides routines for sending tabular data to the
    console as well as standard text strings.

    Currently this class supports output only of vector of a single type.  The
    supported types are:
      - std::string
      - float
      - double
      - long double

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

#ifndef COMOSTRINGS_CONSOLE_IO_H
#define COMOSTRINGS_CONSOLE_IO_H    

#include <vector>
#include <string>

class Console_IO
{
public:
    // Constructors.
    Console_IO(void);                  // Default constructor.
    Console_IO(const Console_IO &cio); // Copy constructor.

    // Destructors.
    ~Console_IO(void); // Default destructor.

    // Enumeration of allowable column formats.
    enum ColFormat {None=-1, Integer, Float, Scientific};


    // TABULAR DATA.

    // Returns the current column count.
    unsigned int ColumnCount() const;
    
    /*
    // Sets the number of columns to be outputted using tabular format.
    void SetColumnCount(unsigned int ncol);
    */

    // Returns the current format of the ith column.
    ColFormat ColumnFormat(unsigned int i) const;

    // Sets the column format of the ith column.
    void SetColumnFormat(unsigned int i, ColFormat format);

    // Sets the format of all columns using a std::vector.  If the vector
    // is not long enough, then Scientific is used for the unspecified
    // columns.
    void SetColumnFormats(const std::vector<ColFormat> &formats);

    // Prints a vector of numbers to the console using the currently specified
    // formatting.  The number of columns is specified by the SetColumnCount()
    // function, and only this many values shall be printed, regardless
    // of the length of the vector.  Optionally a "mask" vector with indices of
    // the output values in the output vector can be provided, though again the
    // number of output values shall match the column count.
    void PrintRow(
        const std::vector<float> &data,       // Data vector to output.
        const std::vector<unsigned int> &mask // Mask of output indices.
          = std::vector<unsigned int>()       //  - Default empty vector for mask.
        ) const;
    void PrintRow(
        const std::vector<double> &data,      // Data vector to output.
        const std::vector<unsigned int> &mask // Mask of output indices.
          = std::vector<unsigned int>()       //  - Default empty vector for mask.
        ) const;
    void PrintRow(
        const std::vector<long double> &data, // Data vector to output.
        const std::vector<unsigned int> &mask // Mask of output indices.
          = std::vector<unsigned int>()       //  - Default empty vector for mask.
        ) const;

    // Prints a data row in tabular format given a vector of strings.  This
    // flavour of the PrintRow function ignores the current column formatting.
    void PrintRow(
        const std::vector<std::string> &data, // List of strings to output.
        const std::vector<unsigned int> &mask // Mask of output indices.
          = std::vector<unsigned int>()       //  - Default empty vector for mask.
        ) const;

    // Prints a dividing row, "---------------- etc", to the console.
    void PrintDivider() const;


    // AUTOMATIC DIVIDERS.

    // Sets the automatic divider interval.  This is the number of
    // rows printed between dividing lines.
    void SetDividerInterval(unsigned int i);

    // Sets the header row to print with the automatic dividers.
    void SetAutoHeader(const std::vector<std::string> &header);

    // Turns on automatic dividers.
    void EnableAutoDividers();

    // Turns off automatic dividers.
    void DisableAutoDividers();

private:
    // The number of characters per row of output.
    unsigned int m_rowlength;

    // The number of columns to output to the console.
    unsigned int m_ncolumns;

    // Vector of column formats.
    std::vector<ColFormat> m_formats;

    // printf format specifications.
    std::string m_ifmt, m_ffmt, m_scifmt, m_strfmt;

    // Divider string.
    std::string m_divider;

    // Auto-divide interval.  This is the number of rows printed between
    // automatically-printed dividers.
    unsigned int m_divint;

    // The number of lines remaining until next auto divider.
    mutable unsigned int m_divnum;

    // String data row to print with the automatic dividers.
    std::vector<std::string> m_autoheader;

    // Set to true if divider lines are to be printed automatically.
    bool m_autodivide;

    // Prints the auto divider.
    void printAutoDivider() const;
};

#endif
