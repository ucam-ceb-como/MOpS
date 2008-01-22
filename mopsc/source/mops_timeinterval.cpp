#include "mops_timeinterval.h"

using namespace Mops;

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

// Default destructor.
TimeInterval::~TimeInterval(void)
{
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
