/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Coagulation class declared in the
    swp_coagulation.h header file.

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

#include "swp_coagulation.h"
#include "swp_mechanism.h"
#include <stdexcept>
#include <iostream>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// Base coagulation class.

/**
  Default rate scaling to 1 for backwards compatibilty
 */
 Coagulation::Coagulation(const Sweep::Mechanism &mech)
 : Process(mech)
 , m_a(1.0)
 {}

// Writes the object to a binary stream.
void Coagulation::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output rate scaling factor
        out.write(reinterpret_cast<const char*>(&m_a), sizeof(m_a));

        // Serialize base class.
        Process::Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Coagulation::Serialize).");
    }
}

// Reads the object from a binary stream.
void Coagulation::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                real a;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                SetA(a);

                // Deserialize base class.
                Process::Deserialize(in, mech);

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Coagulation::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Coagulation::Deserialize).");
    }
}


//==============================================================
