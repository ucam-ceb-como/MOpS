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
        double molwt,             //   - Component molecular weight.
        double dens,              //   - Component density.
        double min,               //   - Min value for a valid particle
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
    double MolWt() const;

    // Sets the molecular weight (g/mol).
    void SetMolWt(double molwt);


    // DENSITY.

    // Returns component density (g/cm3).
    double Density() const;

    double CoalescThresh() const;

    double GrowthFact() const;

    double MinPAH() const;

	// PHASE

	// Return phase
	const std::string &Phase() const;

	// Set phase
	void SetPhase(const std::string &phase);

	// Return element
	const std::string &Element() const;

	// Set element
	void SetElement(const std::string &element);

    //! If a jump process reduces the total number of 6-member rings (excludes 5-member rings) in a PAH below a certain threshold it is removed. Return this threshold in terms of the total number of 6-member rings.
    double ThresholdOxidation() const;

    //! Allow PAHs in soot particles to point to the same memory location after a doubling event.
    double SharedPointers() const;

	//! Allow particles composed of only single PAHs to be respresented with weighted particles
	double WeightedPAHs() const;

    // Sets the density (g/cm3).
    void SetDensity(double dens);

    void SetCoalescThresh(double ct);

    void SetGrowthFact(double gf);

    void SetMinPAH(int mp);

    //! If a jump process reduces the total number of 6-member rings (excludes 5-member rings) in a PAH below a certain threshold it is removed. Sets this threshold in terms of the total number of 6-member rings.
    void SetThresholdOxidation(int to);

    //! Allow PAHs in soot particles to point to the same memory location after a doubling event.
    void SetSharedPointers(int sp);

	//! Allow particles composed of only single PAHs to be respresented with weighted particles
	void SetWeightedPAHs(int wpah);

    // Returns component symbol or name.
    const std::string &Name() const;

    // Sets the symbol or name.
    void SetName(const std::string &name);

    // MINIMUM VALID VALUE
    //! Get the minimum amount that a valid particle may have of this component
    double MinValid() const {return m_minValid;}

    //! Set the minimum amount that a valid particle must have of this component
    void SetMinValid(const double min) {m_minValid = min;}

    //! Check is a valid particle can have specified amount of this component
    bool IsValidValue(const double r) const;

    // READ/WRITE/COPY.

    // Creates a copy of the component.
    Component *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

private:
    double m_density;     // Component density (g/cm3).
    double m_molwt;       // Component molecular weight (g/mol).

    //! Minimum value a valid particle can have for this component
    double m_minValid;

    std::string m_name; // Component symbol or name.

    //! threshold after which two individual primaries are replaced by one
    double m_coalesc_thresh;

    //! factor with which the PAH growth is multiplied
    double m_growthfact;

    //! minimum number of PAHs withinin the primary that are needed to apply m_growthfact
    int m_minPAH;

    //! If a jump process reduces the total number of 6-member rings (excludes 5-member rings) in the PAH (in a particle) below this threshold it is removed.
    double m_thresholdOxidation;
    
    //! Allow PAHs in soot particles to point to the same memory location after a doubling event.
    double m_sharedPointers;

	//! Allow particles composed of only single PAHs to be respresented with weighted particles
	double m_weightedPAHs;

	//! Component phase
	std::string m_phase;

	//! Component element
	std::string m_element;
};

// Typedef of a vector of pointers to Component objects.
typedef std::vector<Component*> CompPtrVector;
};

// Include inline function definitions.
#include "swp_component_inl.h"
#endif
