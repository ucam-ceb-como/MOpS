/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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
