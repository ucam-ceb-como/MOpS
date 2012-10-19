/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The TimeInterval class defines a time interval used by the mops
    solver to control output.

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

#ifndef MOPS_TIMEINTERVAL_H
#define MOPS_TIMEINTERVAL_H

#include "mops_params.h"
#include <vector>
#include <iostream>

namespace Mops
{
/*!
 * Specification of the time intervals over which calculations are to be
 * performed.
 */
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
    double StartTime() const;

    // Sets the start time.
    void SetStartTime(double t);


    // END TIME.

    // Returns the end time.
    double EndTime() const;
    // Sets the end time.
    void SetEndTime(double t);


    // STEP COUNT.

    // Returns the number of output steps.
    unsigned int StepCount() const;
    // Sets the number of output steps.
    void SetStepCount(unsigned int n);


    // SPLITTING STEP COUNT.

    //! Returns the number of splitting steps per output step.
    unsigned int SplittingStepCount() const;
    //! Sets the number of splitting steps per output step.
    void SetSplittingStepCount(unsigned int n);

    //! Returns the number of sub steps per splitting step.
    unsigned int SubSplittingStepCount() const;
    //! Sets the number of sub steps steps per splitting step.
    void SetSubSplittingStepCount(unsigned int n);

    // STEP SIZE CALCULATION.

    // Calculates the output step size.
    double StepSize() const;

    // Calculates the splitting step size.
    double SplitStepSize() const;


    // READ/WRITE/COPY.

    // Writes the time interval to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the time interval from a binary stream.
    void Deserialize(std::istream &in);

private:
    // Start and end times (s).
    double m_start, m_end;

    //! Number of output steps to take before this time.
    unsigned int m_steps;

    //! Number of splitting steps to perform per output step.
    unsigned int m_splits;

    //! Number of pieces into which to break each splitting step.
    unsigned int m_subsplits;
};

// A typedef for a vector of TimeIntervals.  This is defined as
// TimeInterval objects will generally be processed in groups.
typedef std::vector<TimeInterval> timevector;
};

#endif
