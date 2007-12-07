/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This file contains the definition of a structure for a chemical species.  File
    also contains typdefs and structure definitions related to Species objects.
*/

#ifndef GPC_SPECIES_H
#define GPC_SPECIES_H

#include "gpc_el_comp.h"
#include "gpc_thermo_params.h"
#include <vector>
#include <string>

namespace Sprog
{
class Mechanism; // Forward declaration of mechanism.

class Species
{
public:
    // Constructors.
    Species(void);  // Default constructor.
    Species(const Species &sp); // Copy constructor.

    // Default destructor.
    virtual ~Species(void);

    // Operator overloads.
    Species &operator=(const Species &sp);
    bool operator==(const Species &sp) const;
    bool operator==(const std::string &name) const;
    bool operator!=(const Species &sp) const;
    bool operator!=(const std::string &name) const;

    // Species name.
    const std::string &Name(void) const;  // Returns the species name.   
    void SetName(const std::string &name); // Sets the species name.

    // Species elemental composition.
    const ElCompVector &Composition(void) const; // Returns the elemental composition of the species.
    
    void AddElement(const ElComp &elcomp);           // Adds an element to the species composition.
    void AddElement(unsigned int i, unsigned int n); // Adds an element given the index and count.
    void AddElement(const std::string &name, unsigned int n); // Adds an element given the name.
    
    bool ContainsElement(unsigned int i) const;          // Returns true if species contains the 
                                                         // element (given by index).
    bool ContainsElement(const std::string &name) const; // Returns true if species contains the 
                                                         // element (given by name).
    bool ContainsElement(const Element &el) const;       // Returns true if species contains the 
                                                         // element (given by object).

    // Species molecular weight.
    const real MolWt(void) const; // Returns the species molecular weight.
    real CalcMolWt(void);         // Recalculates the molecular weight of the species using the elements.

    // Elements.
    const ElementPtrVector *const Elements(void); // Returns the vector of elements used to define species.
    void SetElements(const ElementPtrVector *const els); // Sets the vector of elements.

    // Thermodynamic fitting parameters.
    unsigned int ThermoRangeCount(void) const;    // Returns the number of parameter ranges.
    void SetThermoStartTemperature(const real T); // Sets the start temperature for the range.  
    const Sprog::Thermo::THERMO_PARAMS &ThermoParams(const real T) const; // Returns the set of parameters valid for 
                                                           // the given temperature.
    void AddThermoParams(const real T, const Sprog::Thermo::THERMO_PARAMS &params); // Adds a set of parameters with the 
                                                                     // given end point temperature.
    void RemoveThermoParams(const real T); // Removes the parameters from the list valid for the given temperature.

    // Parent mechanism.
    Sprog::Mechanism *const Mechanism(void) const;   // Returns pointer to parent mechanism.
    void SetMechanism(Sprog::Mechanism *const mech); // Sets the parent mechanism.

    // Cloning.
    virtual Species *const Clone(void) const; // Returns a pointer to a copy of the Species object.

protected:
    std::string m_name;        // Name/symbol.
    ElCompVector m_elcomp;     // Elemental composition.
    real m_molwt;              // Molecular weight (kg/mol).

    const ElementPtrVector *m_elements; // Elements used to define species.
    Sprog::Mechanism *m_mech;           // Parent mechanism.

    // Thermo parameters for different temperature ranges.  The map key is the
    // end point temperature up to which the parameters are valid.
    Sprog::Thermo::ThermoMap m_thermoparams;
    real m_T1; // Start temperature for range.
};

// Inline function definitions.
#include "gpc_species_inl.h"

// A typedef for a STL vector of species.
typedef std::vector<Species> SpeciesVector;
typedef std::vector<Species*> SpeciesPtrVector;
};

#endif