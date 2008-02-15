/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a particle component.  A component is a chemical element
    or molecule which is the smallest repeating unit in a particle.  A 
    particle can consist of more than one type of component, though the assumption
    is generally that each component exists in discrete phases rather than being
    well mixed.  This makes calculating particle volume and mass much easier.
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


    // READ/WRITE/COPY.

    // Creates a copy of the component.
    Component *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

private:
    real m_molwt;       // Component molecular weight (g/mol).
    real m_density;     // Component density (g/cm3).
    std::string m_name; // Component symbol or name.
};

// Typedef of a vector of pointers to Component objects.
typedef std::vector<Component*> CompPtrVector;
};

// Include inline function definitions.
#include "swp_component_inl.h"
#endif
