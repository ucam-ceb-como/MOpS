/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The Settings class holds simulation settings for mops.
*/

#ifndef MOPS_TIMEINTERVAL_H
#define MOPS_TIMEINTERVAL_H

#include "mops_params.h"

namespace Mops
{
class TimeInterval
{
public:
    // Constructors.
    TimeInterval(void);                   // Default constructor.
    TimeInterval(const TimeInterval &ti); // Copy constructor.

    // Destructors.
    virtual ~TimeInterval(void); // Default destructors.

    // Returns the start time.
    real StartTime() const;
    // Sets the start time.
    void SetStartTime(real t);

    // Returns the end time.
    real EndTime() const;
    // Sets the end time.
    void SetEndTime(real t);

    // Returns the number of output steps.
    unsigned int StepCount() const;
    // Sets the number of output steps.
    void SetStepCount(unsigned int n);

    // Returns the number of splitting steps per output step.
    unsigned int SplittingStepCount() const;
    // Sets the number of splitting steps per output step.
    void SetSplittingStepCount(unsigned int n);

    // Calculates the output step size.
    real StepSize() const;

    // Calculates the splitting step size.
    real SplitStepSize() const;

private:
    // Start and end times (s).
    real m_start, m_end;

    // Number of output steps to take before this time.
    unsigned int m_steps;

    // Number of splitting steps to perform per output step.
    unsigned int m_splits;
};
};

#endif
