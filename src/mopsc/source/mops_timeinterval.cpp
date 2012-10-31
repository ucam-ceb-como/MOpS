/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the TimeInterval class declared in the
    mops_timeinterval.h header file.

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

#include "mops_timeinterval.h"
#include <vector>
#include <iostream>
#include <stdexcept>

using namespace Mops;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
TimeInterval::TimeInterval(void)
{
    m_start  = 0.0;
    m_end    = 0.0;
    m_steps  = 1;
    m_splits = 1;
    m_subsplits = 1;
}

// Copy constructor.
TimeInterval::TimeInterval(const Mops::TimeInterval &ti)
{
    // Check for self assignment.
    if (this != &ti) {
        m_start  = ti.m_start;
        m_end    = ti.m_end;
        m_steps  = ti.m_steps;
        m_splits = ti.m_splits;
        m_subsplits = ti.m_subsplits;
    }
}

// Stream-reading constructor.
TimeInterval::TimeInterval(std::istream &in)
{
    Deserialize(in);
}

// Default destructor.
TimeInterval::~TimeInterval(void)
{
    // Nothing special to destruct.
}


// PROPERTIES.

// Returns the start time.
double TimeInterval::StartTime() const
{
    return m_start;
}

// Sets the start time.
void TimeInterval::SetStartTime(double t)
{
    m_start = t;
}

// Returns the end time.
double TimeInterval::EndTime() const
{
    return m_end;
}

// Sets the end time.
void TimeInterval::SetEndTime(double t)
{
    m_end = t;
}

// Returns the number of output steps.
unsigned int TimeInterval::StepCount() const
{
    return m_steps;
}

// Sets the number of output steps.
void TimeInterval::SetStepCount(unsigned int n)
{
    m_steps = n;
}

// Returns the number of splitting steps per output step.
unsigned int TimeInterval::SplittingStepCount() const
{
    return m_splits;
}

// Sets the number of splitting steps per output step.
void TimeInterval::SetSplittingStepCount(unsigned int n)
{
    m_splits = n;
}

/*!
 *@return   Number of sub steps to take within one splitting step
 */
unsigned int TimeInterval::SubSplittingStepCount() const
{
    return m_subsplits;
}

/*!
 * @param[in]   n   Number of sub steps to take within one splitting step
 */
void TimeInterval::SetSubSplittingStepCount(unsigned int n)
{
    m_subsplits = n;
}


// Calculates the output step size.
double TimeInterval::StepSize() const
{
    return (m_end - m_start) / (double)m_steps;
}

// Calculates the splitting step size.
double TimeInterval::SplitStepSize() const
{
    return StepSize() / (double)m_splits;
}


// READ/WRITE/COPY.

// Writes the time interval to a binary stream.
void TimeInterval::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialisation version number to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write the start time.
        double t = (double)m_start;
        out.write((char*)&t, sizeof(t));

        // Write the end time.
        t = (double)m_end;
        out.write((char*)&t, sizeof(t));

        // Write the number of steps.
        out.write((char*)&m_steps, sizeof(m_steps));

        // Write the number of splitting steps.
        out.write((char*)&m_splits, sizeof(m_splits));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Mops, TimeInterval::Serialize).");
    }
}

// Reads the time interval from a binary stream.
void TimeInterval::Deserialize(std::istream &in)
{
    // Clear the time interval of its current data.
    m_start = m_end = 0.0;
    m_steps = m_splits = 0;

    if (in.good()) {
        // Read the serialized version number.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;

        switch (version) {
            case 0:

                // Read the start time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_start = (double)val;

                // Read the end time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_end = (double)val;

                // Read the number of steps.
                in.read(reinterpret_cast<char*>(&m_steps), sizeof(m_steps));

                // Read the number of splitting steps.
                in.read(reinterpret_cast<char*>(&m_splits), sizeof(m_splits));

                break;
            default:
                throw runtime_error("TimeInterval serialized version "
                                    "number is unsupported (Mops, TimeInterval::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Mops, TimeInterval::Deserialize).");
    }
}
