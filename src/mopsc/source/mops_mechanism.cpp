/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Mechanism class declared in the
    mops_mechanism.h header file.

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

#include "mops_mechanism.h"
#include <stdexcept>
#include <fstream>

using namespace Mops;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Mechanism::Mechanism(void)
{
}

// Default destructor.
Mechanism::~Mechanism(void)
{
}


// PARTICLE MECHANISM.

// READ/WRITE/COPY FUNCTIONS.

// Writes the mechanism to a binary data stream.
void Mechanism::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialize version to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));
        
        // Write the gas mechanism
        m_gmech.Serialize(out);
        // and the particle mechanism.
        m_pmech.Serialize(out);
		/*
		// added by mm864 (09 June 2012): the surface mechanism. 
        m_smech.Serialize(out);
		*/
	} else {
        throw std::invalid_argument("Output stream not ready "
                               "(Mops, Mechanism::Serialize).");
    }
}

// Reads the mechanism data from a binary data stream.
void Mechanism::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the serialized mechanism version.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                m_gmech.Deserialize(in);
				
				/*
				// Added by mm864 - read the surface mechanism
				m_smech.Deserialize(in);
				m_smech.SetSpecies(SurfMech().Species()); 
				*/
				
                // Read the particle mechanism.
                m_pmech.Deserialize(in);
                m_pmech.SetSpecies(GasMech().Species()); // takes all surface 

                break;
            default:
                throw std::runtime_error("Mechanism serialized version number "
                                    "is unsupported (Mops, Mechanism::Deserialize).");
        }
    } else {
        throw std::invalid_argument("Input stream not ready (Mops, Mechanism::Deserialize).");
    }
}
