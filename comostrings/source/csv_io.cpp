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
    // Select the file stream open flags depending on whether the
    // current file contents should be deleted or not.
    int flags;
    if (del) {
        // Use ios_base::trunc to delete current file contents.
        flags = ios_base::in | ios_base::out | ios_base::trunc;
    } else {
        // Use ios_base::app to append to current file contents.
        flags = ios_base::in | ios_base::out | ios_base::app;
    }

    // Attempt to open the file for input and output.
    m_file.open(name.c_str(), flags);

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
    char c;
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


// INTERNAL READ/WRITE FUNCTIONS.

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
    streamoff i = m_file.tellp();
    if (i != 0) {
        m_file << '\n';
    }

    // Check that the file is open before writing.  Also check that
    // the values vector has some values to write.
    if (m_file.is_open() && m_file.good() && (values.size()>0)) {
        // Declare vector iterator and set it to beginning of values.
        vector<T>::const_iterator i = values.begin();

        // Output the comma-separated values.
        m_file << *i++;
        while (i!=values.end()) {
            m_file << ',' << *i++;
        }
    }
}