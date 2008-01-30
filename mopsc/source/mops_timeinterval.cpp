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
real TimeInterval::StartTime() const
{
    return m_start;
}

// Sets the start time.
void TimeInterval::SetStartTime(real t)
{
    m_start = t;
}

// Returns the end time.
real TimeInterval::EndTime() const
{
    return m_end;
}

// Sets the end time.
void TimeInterval::SetEndTime(real t)
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

// Calculates the output step size.
real TimeInterval::StepSize() const
{
    return (m_end - m_start) / (real)m_steps;
}

// Calculates the splitting step size.
real TimeInterval::SplitStepSize() const
{
    return StepSize() / (real)m_splits;
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
                m_start = (real)val;
                
                // Read the end time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_end = (real)val;

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
