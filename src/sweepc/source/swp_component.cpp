/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Component class declared in the
    swp_component.h header file.

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

#include "swp_component.h"
#include <stdexcept>

using namespace Sweep;

// CONSRUCTORS AND DESTRUCTORS.

// Default constructor.
Component::Component()
: m_density(0.0), m_molwt(0.0), m_minValid(0.0), m_name(""),m_coalesc_thresh(1.0), m_growthfact(1.0), m_minPAH(0), 
m_phase(""), m_element("")
{
}

// Initialising constructor.
Component::Component(double molwt, 
                     double dens,
                     double min,
                     const std::string &name)
{
    // Initialise the component properties.
    m_density = dens;
    m_molwt   = molwt;
    m_minValid = min;
    m_name    = name;
    m_coalesc_thresh = 1.0;     //added by ms785, do not coalesce particles by default
    m_growthfact=1.0;           //added by ms785
    m_minPAH=0;                 //added by ms785
	m_phase = "";
	m_element = "";
}

// Copy constructor.
Component::Component(const Component &copy) 
{
    *this = copy;
}

// Stream-reading constructor.
Component::Component(std::istream &in) 
{
    Deserialize(in);
}

// Default destructor.
Component::~Component()
{
}


// OPERATOR OVERLOADS.
Component &Component::operator=(const Component &rhs)
{
    if (this != &rhs) {
        m_density = rhs.m_density;
        m_molwt   = rhs.m_molwt;
        m_minValid = rhs.m_minValid;
        m_name    = rhs.m_name;
        m_coalesc_thresh= rhs.m_coalesc_thresh;
        m_growthfact= rhs.m_growthfact;
        m_minPAH= rhs.m_minPAH;
		m_phase = rhs.m_phase;
		m_element = rhs.m_element;
    }
    return *this;
}

/*! 
 * Method intended for use in Particle::IsValid()
 *
 *@param[in]    r   Amount of this component in a particle
 *
 *@return       True if it is possible for a valid particle to have the 
 *              specifed amount of this component
 */
bool Component::IsValidValue(const double r) const {
    return r >= m_minValid;
}

// READ/WRITE/COPY.

// Creates a copy of the component.
Component *const Component::Clone(void) const {return new Component(*this);}

// Writes the object to a binary stream.
void Component::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write molecular weight.
        double v = (double)m_molwt;
        out.write((char*)&v, sizeof(v));

        // Write density
        v = (double)m_density;
        out.write((char*)&v, sizeof(v));

        // Write minimum value for a valid particle
        v = static_cast<double>(m_minValid);
        out.write((char*)&v, sizeof(v));

        // Write the length of the component name to the stream.
        unsigned int n = m_name.length();
        out.write((char*)&n, sizeof(n));

        // Write the component name to the stream.
        out.write(m_name.c_str(), n);

		// Write phase to the stream.
		n = m_phase.length();
		out.write((char*)&n, sizeof(n));
		out.write(m_phase.c_str(), n);

		// Write element to the stream.
		n = m_element.length();
		out.write((char*)&n, sizeof(n));
		out.write(m_element.c_str(), n);

    } else {
        throw std::invalid_argument("Output stream not ready "
                                    "(Sweep, Component::Serialize).");
    }
}

// Reads the object from a binary stream.
void Component::Deserialize(std::istream &in)
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
                // Read molecular weight.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_molwt = (double)val;

                // Read density
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_density = (double)val;

                // Read minimum value for a valid particle
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_minValid = static_cast<double>(val);

                // Read the length of the species name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                
                // Read the species name.
                name = new char[n];
                in.read(name, n);
                m_name.assign(name, n);
                delete [] name;

				// Read the phase.
				in.read(reinterpret_cast<char*>(&n), sizeof(n));
				name = new char[n];
				in.read(name, n);
				m_phase.assign(name, n);
				delete[] name;

				// Read the element.
				in.read(reinterpret_cast<char*>(&n), sizeof(n));
				name = new char[n];
				in.read(name, n);
				m_element.assign(name, n);
				delete[] name;

                break;
            default:
                throw std::runtime_error("Serialized version number is invalid "
                                         "(Sweep, Component::Deserialize).");
        }
    } else {
        throw std::invalid_argument("Input stream not ready "
                                    "(Sweep, Component::Deserialize).");
    }
}
