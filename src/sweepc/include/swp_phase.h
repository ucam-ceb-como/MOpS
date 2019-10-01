/*
  Author(s):      Casper Lindberg (csl37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2019 Casper Lindberg

  File purpose:
	Definition of a phase. The class identifies components that belong to a specific phase. 

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

#ifndef SWEEP_PHASE_H
#define SWEEP_PHASE_H

#include "swp_params.h"
#include <string>
#include <vector>
#include <iostream>

namespace Sweep
{
class Phase
{
public:
    // Constructors.
	Phase(void); // Default constructor.
	Phase(                  // Initialising constructor.
        const std::string &name, //   - name or symbol.
		bool liquid
        );
	Phase(const Phase &copy); // Copy constructor.
	Phase(std::istream &in); // Stream-reading constructor.

    // Destructor.
	~Phase(void);

    // Operators.
	Phase &operator=(const Phase &rhs);

    // Returns phase symbol or name.
    const std::string &Name() const;

    // Sets the symbol or name.
    void SetName(const std::string &name);  

	// Get the name
	std::string Name();

	// Add component index to phase
	void AddComponent(unsigned int id);

	// Return component indices
	std::vector<unsigned int> GetComponents() const;

	// Set if phase is liquid (for the melting model)
	void SetLiquid();

	// Return true if the phase is liquid (for the melting model)
	bool GetLiquid() const;

    // READ/WRITE/COPY.

    // Creates a copy of the phase
	Phase *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

private:
  
    std::string m_name; // Phase name.

	std::vector<unsigned int> m_compIndex; // Vector of component indices that belong to phase

	bool m_liquid; // Liquid phase

};

// Typedef of a vector of pointers to phase objects.
typedef std::vector<Phase*> PhasePtrVector;

};

#endif
