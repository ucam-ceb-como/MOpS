/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of a particle component.  A component is a chemical element
    or molecule which is the smallest repeating unit in a particle.  A 
    particle can consist of more than one type of component, though the assumption
    is generally that each component exists in discrete phases rather than being
    well mixed.  This makes calculating particle volume and mass much easier.

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

#ifndef SWEEP_COMPONENT_H
#define SWEEP_COMPONENT_H

#include "swp_params.h"
#include <string>
#include <vector>
#include <iostream>

namespace Sweep
{
class Component
{
public:
    // Constructors.
    Component(void); // Default constructor.
    Component(                  // Initialising constructor.
        real molwt,             //   - Component molecular weight.
        real dens,              //   - Component density.
        real min,               //   - Min value for a valid particle
        const std::string &name //   - Component name or symbol.
        );
    Component(const Component &copy); // Copy constructor.
    Component(std::istream &in); // Stream-reading constructor.

    // Destructor.
    ~Component(void);

    // Operators.
    Component &operator=(const Component &rhs);


    // MOLECULAR WEIGHT
    
    // Returns component molecular weight (g/mol).
    real MolWt() const;

    // Sets the molecular weight (g/mol).
    void SetMolWt(real molwt);


    // DENSITY.

    // Returns component density (g/cm3).
    real Density() const;

    // Sets the density (g/cm3).
    void SetDensity(real dens);


    // COMPONENT NAME.

    // Returns component symbol or name.
    const std::string &Name() const;

    // Sets the symbol or name.
    void SetName(const std::string &name);

    // MINIMUM VALID VALUE
    //! Get the minimum amount that a valid particle may have of this component
    real MinValid() const {return m_minValid;}

    //! Set the minimum amount that a valid particle must have of this component
    void SetMinValid(const real min) {m_minValid = min;}

    //! Check is a valid particle can have specified amount of this component
    bool IsValidValue(const real r) const;

    // READ/WRITE/COPY.

    // Creates a copy of the component.
    Component *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

private:
    real m_density;     // Component density (g/cm3).
    real m_molwt;       // Component molecular weight (g/mol).

    //! Minimum value a valid particle can have for this component
    real m_minValid;

    std::string m_name; // Component symbol or name.
};

// Typedef of a vector of pointers to Component objects.
typedef std::vector<Component*> CompPtrVector;
};

// Include inline function definitions.
#include "swp_component_inl.h"
#endif
