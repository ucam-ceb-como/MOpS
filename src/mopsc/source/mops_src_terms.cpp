/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the SrcPoint class declared in the
    mops_src_terms.h header file.

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

#include "mops_src_terms.h"
#include <algorithm>

using namespace Mops;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SrcPoint::SrcPoint()
: Time(0.0)
{
}

// Initialising constructor.
SrcPoint::SrcPoint(unsigned int nterms)
: Time(0.0), Terms(nterms, 0.0)
{
}

// Copy constructor.
SrcPoint::SrcPoint(const SrcPoint &copy)
: Time(copy.Time), Terms(copy.Terms)
{
}

// Default destructor.
SrcPoint::~SrcPoint(void)
{
}


// OPERATORS.

// Assignment operator.
SrcPoint &SrcPoint::operator=(const SrcPoint &rhs)
{
    if (this != &rhs) {
        Time = rhs.Time;
        Terms.assign(rhs.Terms.begin(), rhs.Terms.end());
    }
    return *this;
}


// POINT COMPARISONS.

// Returns true if the current point is before the given point
// in time.
bool SrcPoint::IsBefore(const SrcPoint &rhs)
{
    return Time < rhs.Time;
}

// Returns true if the LHS point is before the RHS point
// in time.
bool SrcPoint::IsBeforePoint(const SrcPoint &lhs, const SrcPoint &rhs)
{
    return lhs.Time < rhs.Time;
}

// Returns true if the current point is before the given time.
bool SrcPoint::IsBeforeTime(const SrcPoint &lhs, double t)
{
    return lhs.Time < t;
}

// Returns true if the current point is after the given point
// in time.
bool SrcPoint::IsAfter(const SrcPoint &rhs)
{
    return Time > rhs.Time;
}

// Returns true if the LHS point is after the RHS point
// in time.
bool SrcPoint::IsAfterPoint(const SrcPoint &lhs, const SrcPoint &rhs)
{
    return lhs.Time > rhs.Time;
}

// Returns true if the current point is after the given time.
bool SrcPoint::IsAfterTime(const SrcPoint &lhs, double t)
{
    return lhs.Time > t;
}

// GLOBAL FUNCTIONS IN HEADER FILE.

// Sort a gas-profile in order of ascending time.
void Mops::SortSrcProfile(Mops::SrcProfile &prof)
{
    std::sort(prof.begin(), prof.end(), SrcPoint::IsAfterPoint);
}

// Returns the first SrcPoint defined after the given time.  If the time
// is out-of-range then returns the end() of the vector.
SrcProfile::const_iterator Mops::LocateSrcPoint(const SrcProfile &prof, double t)
{
    for (SrcProfile::const_iterator i=prof.begin(); i!=prof.end(); ++i) {
        if (SrcPoint::IsAfterTime(*i, t)) return i;
    }
    return prof.end();
}
