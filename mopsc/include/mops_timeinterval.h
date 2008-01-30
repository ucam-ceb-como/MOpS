/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The TimeInterval class defines a time interval used by the mops
    solver to control output.
*/

#ifndef MOPS_TIMEINTERVAL_H
#define MOPS_TIMEINTERVAL_H

#include "mops_params.h"
#include <vector>
#include <iostream>

namespace Mops
{
class TimeInterval
{
public:
    // Constructors.
    TimeInterval(void);                   // Default constructor.
    TimeInterval(const TimeInterval &ti); // Copy constructor.
    TimeInterval(std::istream &in);       // Stream-reading constructor.

    // Destructors.
    ~TimeInterval(void); // Default destructors.


    // START TIME.

    // Returns the start time.
    real StartTime() const;
    
    // Sets the start time.
    void SetStartTime(real t);


    // END TIME.

    // Returns the end time.
    real EndTime() const;
    // Sets the end time.
    void SetEndTime(real t);


    // STEP COUNT.

    // Returns the number of output steps.
    unsigned int StepCount() const;
    // Sets the number of output steps.
    void SetStepCount(unsigned int n);


    // SPLITTING STEP COUNT.

    // Returns the number of splitting steps per output step.
    unsigned int SplittingStepCount() const;
    // Sets the number of splitting steps per output step.
    void SetSplittingStepCount(unsigned int n);


    // STEP SIZE CALCULATION.

    // Calculates the output step size.
    real StepSize() const;

    // Calculates the splitting step size.
    real SplitStepSize() const;


    // READ/WRITE/COPY.

    // Writes the time interval to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the time interval from a binary stream.
    void Deserialize(std::istream &in);

private:
    // Start and end times (s).
    real m_start, m_end;

    // Number of output steps to take before this time.
    unsigned int m_steps;

    // Number of splitting steps to perform per output step.
    unsigned int m_splits;
};

// A typedef for a vector of TimeIntervals.  This is defined as
// TimeInterval objects will generally be processed in groups.
typedef std::vector<TimeInterval> timevector;
};

#endif
