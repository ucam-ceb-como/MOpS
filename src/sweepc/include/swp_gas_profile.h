/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.


  File purpose:
    The GasPoint and GasProfile objects allow a gas system to be
    described over a series of time points.  They are used principally
    by the Sweep::FlameSolver to store the result of a premixed flame
    calculation of the gas-phase.  Additionally they are used by the
    predictor-corrector solver Mops::PredCorSolver to store the gas
    profile over a time step.

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

#ifndef SWEEP_GAS_PROFILE_H
#define SWEEP_GAS_PROFILE_H

#include "swp_params.h"

#include "gpc_species.h"
#include "gpc_idealgas.h"

#include <vector>

namespace Sweep
{
class GasPoint
{
public:
    // Constructors.
    GasPoint(const Sprog::SpeciesPtrVector &species); // Default constructor.
    GasPoint(const GasPoint &copy); // Copy constructor.

    // Destructors.
    ~GasPoint(void);

    // Operators.
    GasPoint &operator=(const GasPoint &rhs);
/*
    // Comparison operators consider only the time coordinate.
    bool operator==(const GasPoint &rhs) const;
    bool operator!=(const GasPoint &rhs) const;
    bool operator>(const GasPoint &rhs) const;
    bool operator<(const GasPoint &rhs) const;
    bool operator>=(const GasPoint &rhs) const;
    bool operator<=(const GasPoint &rhs) const;
*/

    // COMPARISON.
    
    // Returns true if the current point is before the given point
    // in time.
    bool IsBefore(const GasPoint &rhs);

    // Returns true if the LHS point is before the RHS point
    // in time.
    static bool IsBeforePoint(const GasPoint &lhs, const GasPoint &rhs);

    // Returns true if the current point is before the given time.
    static bool IsBeforeTime(const GasPoint &lhs, double t);

    // Returns true if the current point is after the given point
    // in time.
    bool IsAfter(const GasPoint &rhs);

    // Returns true if the LHS point is after the RHS point
    // in time.
    static bool IsAfterPoint(const GasPoint &lhs, const GasPoint &rhs);

    // Returns true if the current point is after the given time.
    static bool IsAfterTime(const GasPoint &lhs, double t);

    //Checks that two entries have the same time point
    static bool IsEqualTime(const GasPoint &lhs, const GasPoint &rhs);

    // DATA.

    // Time at which gas conditions are valid.
    double Time;

    // Gas conditions.
    Sprog::Thermo::IdealGas Gas;
};

// GasProfile typedef defines a vector of gas points.
typedef std::vector<GasPoint> GasProfile;

// Sort a gas-profile in order of ascending time.
void SortGasProfile(GasProfile &prof);

// Returns the last GasPoint defined before the given time.  If the time
// is out-of-range then returns the end() of the vector.
GasProfile::const_iterator LocateGasPoint(const GasProfile &prof, double t);
};

#endif
