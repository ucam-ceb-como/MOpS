#include "swp_tracker.h"

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Tracker::Tracker()
: m_name("")
{
}

// Copy constructor.
Tracker::Tracker(const Sweep::Tracker &copy)
{
    *this = copy;
}

// Initialising constructor.
Tracker::Tracker(const std::string &name)
{
    m_name = name;
}

// Default destructor.
Tracker::~Tracker()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
Tracker &Tracker::operator =(const Sweep::Tracker &rhs)
{
    if (this != &rhs) {
        m_name = rhs.m_name;
    }
    return *this;
}