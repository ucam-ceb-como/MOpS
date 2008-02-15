/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    
*/

#ifndef SWEEP_TRACKER_H
#define SWEEP_TRACKER_H

#include "swp_params.h"
#include <vector>
#include <string>
#include <iostream>

namespace Sweep
{
class Tracker
{
public:
    // Constructors.
    Tracker(void);                    // Default constructor.
    Tracker(const Tracker &copy);     // Copy constructor.
    Tracker(const std::string &name); // Initialising constructor.
    Tracker(std::istream &in);        // Stream-reading constructor.

    // Destructor.
    ~Tracker(void);

    // Operators.
    Tracker &operator=(const Tracker &rhs);


    // NAME.

    // Returns the tracker name.
    const std::string &Name(void) const;

    // Sets the tracker name.
    void SetName(const std::string &name);


    // READ/WRITE/COPY.

    // Returns a copy of the tracker variable.
    Tracker *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

private:
    std::string m_name;
};

// Typedef of a vector of pointers to Component objects.
typedef std::vector<Tracker*> TrackPtrVector;

// Include inline function definitions.
#include "swp_tracker_inl.h"
};

#endif
