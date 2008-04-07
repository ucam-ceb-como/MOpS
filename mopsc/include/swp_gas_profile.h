/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The GasPoint and GasProfile objects allow a gas system to be
    described over a series of time points.  They are used principally
    by the Sweep::FlameSolver to store the result of a premixed flame
    calculation of the gas-phase.  Additionally they are used by the
    predictor-corrector solver Mops::PredCorSolver to store the gas
    profile over a time step.
*/

#ifndef SWEEP_GAS_PROFILE_H
#define SWEEP_GAS_PROFILE_H

#include "swp_params.h"
#include "sprog.h"
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
    static bool IsBeforeTime(const GasPoint &lhs, real t);

    // Returns true if the current point is after the given point
    // in time.
    bool IsAfter(const GasPoint &rhs);

    // Returns true if the LHS point is after the RHS point
    // in time.
    static bool IsAfterPoint(const GasPoint &lhs, const GasPoint &rhs);

    // Returns true if the current point is after the given time.
    static bool IsAfterTime(const GasPoint &lhs, real t);

    // DATA.

    // Time at which gas conditions are valid.
    real Time;

    // Gas conditions.
    Sprog::Thermo::IdealGas Gas;
};

// GasProfile typedef defines a vector of gas points.
typedef std::vector<GasPoint> GasProfile;

// Sort a gas-profile in order of ascending time.
void SortGasProfile(GasProfile &prof);

// Returns the last GasPoint defined before the given time.  If the time
// is out-of-range then returns the end() of the vector.
GasProfile::const_iterator LocateGasPoint(const GasProfile &prof, real t);
};

#endif
