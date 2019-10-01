/*
	Author(s):      Casper Lindberg (csl37)
	Project:        sweep (population balance solver)
	Sourceforge:    http://sourceforge.net/projects/mopssuite

	Copyright (C) 2019 Casper Lindberg

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
		Prof. Markus Kraft
		Dept of Chemical Engineering
		University of Cambridge
		Philippa Fawcett Drive
		Cambridge
		CB3 0AS
		UK

	Email:       mk306@cam.ac.uk
	Website:     http://como.cheng.cam.ac.uk
*/

#include "swp_phase.h"
#include <stdexcept>

using namespace Sweep;

// CONSRUCTORS AND DESTRUCTORS.

// Default constructor.
Phase::Phase()
	: m_name(""), m_liquid(false)
{
}

// Initialising constructor.
Phase::Phase(const std::string &name, bool liquid)
{
    // Initialise the component properties.
    m_name    = name;
	m_liquid = liquid;
}

// Copy constructor.
Phase::Phase(const Phase &copy)
{
    *this = copy;
}

// Stream-reading constructor.
Phase::Phase(std::istream &in)
{
    Deserialize(in);
}

// Default destructor.
Phase::~Phase()
{
}

// OPERATOR OVERLOADS.
Phase &Phase::operator=(const Phase &rhs)
{
    if (this != &rhs) {
        m_name    = rhs.m_name;
		m_liquid  = rhs.m_liquid;
		// Copy component indices.
		for (std::vector<unsigned int>::const_iterator it = rhs.m_compIndex.begin(); it != rhs.m_compIndex.end(); ++it) {
			m_compIndex.push_back((*it));
		}
	}
    return *this;
}

// Add index
void Phase::AddComponent(unsigned int id){
	m_compIndex.push_back(id);
}

// Get components
std::vector<unsigned int> Phase::GetComponents() const{
	return m_compIndex;
}

// Set liquid phase
void Phase::SetLiquid(){
	m_liquid = true;
}

// Return if phase is liquid
bool Phase::GetLiquid() const{
	return m_liquid;
}

// Set name 
void Phase::SetName(const std::string &name) { m_name = name; }

// Return name 
std::string Phase::Name() { return m_name; }

// READ/WRITE/COPY.

// Creates a copy of the component.
Phase *const Phase::Clone(void) const { return new Phase(*this); }

// Writes the object to a binary stream.
void Phase::Serialize(std::ostream &out) const
{

	const unsigned int trueval = 1;
	const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

		// Write the length of the component name to the stream.
        unsigned int n = m_name.length();
        out.write((char*)&n, sizeof(n));

        // Write the component name to the stream.
        out.write(m_name.c_str(), n);

		// Write if liquid
		if (m_liquid) {
			out.write((char*)&trueval, sizeof(trueval));
		}
		else {
			out.write((char*)&falseval, sizeof(falseval));
		}

    } else {
        throw std::invalid_argument("Output stream not ready "
                                    "(Sweep, Phase::Serialize).");
    }
}

// Reads the object from a binary stream.
void Phase::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;
        unsigned int n = 0;
        char *name = NULL;

        switch (version) {
            case 0:
                // Read the length of the species name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                
                // Read the species name.
                name = new char[n];
                in.read(name, n);
                m_name.assign(name, n);
                delete [] name;

				// read if liquid
				in.read(reinterpret_cast<char*>(&n), sizeof(n));
				m_liquid = (n == 1);

                break;
            default:
                throw std::runtime_error("Serialized version number is invalid "
                                         "(Sweep, Phase::Deserialize).");
        }
    } else {
        throw std::invalid_argument("Input stream not ready "
                                    "(Sweep, Phase::Deserialize).");
    }
}
