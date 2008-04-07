#include "swp_gas_profile.h"
#include <algorithm>

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

// Returns true if the current point is before the given time.
bool GasPoint::IsBeforeTime(const GasPoint &lhs, real t)
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
bool GasPoint::IsAfterTime(const GasPoint &lhs, real t)
{
    return lhs.Time > t;
}

// GLOBAL FUNCTIONS IN HEADER FILE.

// Sort a gas-profile in order of ascending time.
void Sweep::SortGasProfile(Sweep::GasProfile &prof)
{
    std::sort(prof.begin(), prof.end(), GasPoint::IsAfterPoint);
}

// Returns the first GasPoint defined after the given time.  If the time
// is out-of-range then returns the end() of the vector.
GasProfile::const_iterator Sweep::LocateGasPoint(const GasProfile &prof, real t)
{
    for (GasProfile::const_iterator i=prof.begin(); i!=prof.end(); ++i) {
        if (GasPoint::IsAfterTime(*i, t)) return i;
    }
    return prof.end();
}
