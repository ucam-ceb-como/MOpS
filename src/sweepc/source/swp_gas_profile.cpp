/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the GasPoint class declared in the
    swp_gas_profile.h header file.

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

#include "swp_gas_profile.h"

#include <algorithm>
#include <stdexcept>

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
GasPoint::GasPoint(const Sprog::SpeciesPtrVector &sp)
: Time(0.0), Gas(sp)
{
}

// Copy constructor.
GasPoint::GasPoint(const GasPoint &copy)
: Time(copy.Time), Gas(copy.Gas)
{
}

// Default destructor.
GasPoint::~GasPoint(void)
{
}


// OPERATORS.

// Assignment operator.
GasPoint &GasPoint::operator=(const GasPoint &rhs)
{
    if (this != &rhs) {
        Time = rhs.Time;
        Gas  = rhs.Gas;
    }
    return *this;
}


// POINT COMPARISONS.

// Returns true if the current point is before the given point
// in time.
bool GasPoint::IsBefore(const GasPoint &rhs)
{
    return Time < rhs.Time;
}

// Returns true if the LHS point is before the RHS point
// in time.
bool GasPoint::IsBeforePoint(const GasPoint &lhs, const GasPoint &rhs)
{
    return lhs.Time < rhs.Time;
}

//Checks that two entries have the same time point
bool GasPoint::IsEqualTime(const GasPoint &lhs, const GasPoint &rhs)
{
    return lhs.Time == rhs.Time;
}

// Returns true if the current point is before the given time.
bool GasPoint::IsBeforeTime(const GasPoint &lhs, double t)
{
    return lhs.Time < t;
}

// Returns true if the current point is after the given point
// in time.
bool GasPoint::IsAfter(const GasPoint &rhs)
{
    return Time > rhs.Time;
}

// Returns true if the LHS point is after the RHS point
// in time.
bool GasPoint::IsAfterPoint(const GasPoint &lhs, const GasPoint &rhs)
{
    return lhs.Time > rhs.Time;
}

// Returns true if the current point is after the given time.
bool GasPoint::IsAfterTime(const GasPoint &lhs, double t)
{
    return lhs.Time > t;
}

// GLOBAL FUNCTIONS IN HEADER FILE.

/*!
 * Sort a gas-profile in order of ascending time.
 * Duplicate time entries are removed
 *
 *@param[in,out]    prof    Profile to sort
 *
 *@exception        std::logic_error    Failed to produce a sequence with strictly increasing times
 */

void Sweep::SortGasProfile(Sweep::GasProfile &prof)
{   
    if(!prof.empty()) {
        //copy the original profile to temp
        Sweep::GasProfile temp(prof);
        //delete the original profile
        prof.clear();
        //sort the temp profile according to the time stamps
        std::stable_sort(temp.begin(), temp.end(), GasPoint::IsBeforePoint);
        //move duplicate entries at the end, new_end is the last entry before the duplicate entries
        GasProfile::iterator new_end=std::unique(temp.begin(), temp.end(), GasPoint::IsEqualTime);
        //copy the temp to prof without the duplicate entries
        prof.insert(prof.begin(),temp.begin(),new_end);

//        GasProfile::const_iterator it = prof.begin();
//        const GasProfile::const_iterator itEnd = prof.end();
//        while(++it != itEnd) {
//            if(!((it-1)->Time < it->Time)) {
//                std::cerr << (it-1)->Time << ", " << it->Time << '\n';
//                throw std::logic_error("Failed to put gas profile points in strictly increasing time order (Sweep::SortGasProfile)");
//            }
//        }
    }
}

// Returns the first GasPoint defined after the given time.  If the time
// is out-of-range then returns the end() of the vector.
GasProfile::const_iterator Sweep::LocateGasPoint(const GasProfile &prof, double t)
{
    for (GasProfile::const_iterator i=prof.begin(); i!=prof.end(); ++i) {
        if (GasPoint::IsAfterTime(*i, t)) return i;
    }
    // return the last element
    return prof.end()-1;
}
