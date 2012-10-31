/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The SrcPoint and SrcProfile objects define the external source terms
    of a reactor object.  The idea is that a source term profile is held
    by the reactor object and an external function is used to interpolate 
    the points, depending on how many points are defined and the structure
    of the interpolation scheme.  As the source terms are external to the 
    reactor, the reactor class has no knowledge of how to use them, merely
    that they exist.

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

#ifndef MOPS_SRC_TERMS_H
#define MOPS_SRC_TERMS_H

#include "mops_params.h"
#include <vector>

namespace Mops
{
class SrcPoint
{
public:
    // Constructors.
    SrcPoint();                     // Default constructor.
    SrcPoint(const SrcPoint &copy); // Copy constructor.
    SrcPoint(unsigned int nterms);  // Initialising constructor.

    // Destructors.
    ~SrcPoint(void);

    // Operators.
    SrcPoint &operator=(const SrcPoint &rhs);

    // COMPARISON.
    
    // Returns true if the current point is before the given point
    // in time.
    bool IsBefore(const SrcPoint &rhs);

    // Returns true if the LHS point is before the RHS point
    // in time.
    static bool IsBeforePoint(const SrcPoint &lhs, const SrcPoint &rhs);

    // Returns true if the current point is before the given time.
    static bool IsBeforeTime(const SrcPoint &lhs, double t);

    // Returns true if the current point is after the given point
    // in time.
    bool IsAfter(const SrcPoint &rhs);

    // Returns true if the LHS point is after the RHS point
    // in time.
    static bool IsAfterPoint(const SrcPoint &lhs, const SrcPoint &rhs);

    // Returns true if the current point is after the given time.
    static bool IsAfterTime(const SrcPoint &lhs, double t);

    // DATA.

    // Time at which source terms are valid.
    double Time;

    // Source terms.
    fvector Terms;
};

// GasProfile typedef defines a vector of gas points.
typedef std::vector<SrcPoint> SrcProfile;

// Sort a gas-profile in order of ascending time.
void SortSrcProfile(SrcProfile &prof);

// Returns the last GasPoint defined before the given time.  If the time
// is out-of-range then returns the end() of the vector.
SrcProfile::const_iterator LocateSrcPoint(const SrcProfile &prof, double t);

// Typedef of a function pointer which defines how source terms are used
// to calculate the reactor ODE RHSs.
typedef void (*SrcTermFnPtr)(
    double *rhs,             // ODE right-hand sides (to be updated)
    unsigned int n,        // Number of values in rhs array.
    double t,                // Current time.
    const SrcProfile &prof // Source term profile.
    );
};

#endif
