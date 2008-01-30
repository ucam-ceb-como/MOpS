/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This file contains the definition of a structure for a chemical species.  File
    also contains typdefs and structure definitions related to Species objects.
*/

#ifndef GPC_SPECIES_H
#define GPC_SPECIES_H

#include "gpc_element.h"
#include "gpc_el_comp.h"
#include "gpc_thermo_params.h"
#include <vector>
#include <string>
#include <iostream>

namespace Sprog
{
class Mechanism; // Forward declaration of mechanism.

class Species
{
public:
    // Constructors.
    Species(void);              // Default constructor.
    Species(const Species &sp); // Copy constructor.
    Species(std::istream &in);  // Stream-reading constructor.

    // Default destructor.
    ~Species(void);

    // Operator overloads.
    Species &operator=(const Species &sp);
    bool operator==(const Species &sp) const;
    bool operator==(const std::string &name) const;
    bool operator!=(const Species &sp) const;
    bool operator!=(const std::string &name) const;


    // SPECIES NAME.

    // Returns the species name.
    const std::string &Name(void) const;

    // Sets the species name.
    void SetName(const std::string &name);


    // ELEMENTAL COMPOSITION.
    
    // Returns the elemental composition of the species.
    const ElCompVector &Composition(void) const;
    
    // Returns the number of elements required to define the species.
    unsigned int ComponentCount(void) const;

    // Returns the total number of atoms in the species.
    unsigned int AtomCount(void) const;

    // Returns the number of the given element in the species.
    unsigned int AtomCount(unsigned int iel) const;

    // Adds an element to the species composition using an ElComp object.
    void AddElement(const ElComp &elcomp);

    // Adds an element given the index and count.
    void AddElement(unsigned int i, unsigned int n);

    // Adds an element given the name.  Element found using parent mechanism.
    void AddElement(const std::string &name, unsigned int n);
    
    // Returns true if species contains the element (given by index).
    bool ContainsElement(unsigned int i) const;

    // Returns true if species contains the element (given by name).
    bool ContainsElement(const std::string &name) const;
    
    // Returns true if species contains the element (given by object).
    bool ContainsElement(const Element &el) const;


    // MOLECULAR WEIGHT.

    // Returns the species molecular weight.
    real MolWt(void) const;

    // Recalculates the molecular weight of the 
    // species using the elements in the parent mechanism.
    real CalcMolWt(void);


    // Elements.
    /* Removed as the elements can be defined using the parent mechanism.
    const ElementPtrVector *const Elements(void); // Returns the vector of elements used to define species.
    void SetElements(const ElementPtrVector *const els); // Sets the vector of elements.
    */


    // THERMODYNAMIC FITTING PARAMETERS.

    // Returns the number of thermo parameter ranges.
    unsigned int ThermoRangeCount(void) const;

    // Sets the start temperature for the thermo parameter range.  
    void SetThermoStartTemperature(const real T);

    // Returns the set of thermo parameters valid for the given temperature.
    const Sprog::Thermo::THERMO_PARAMS &ThermoParams(const real T) const;

    // Adds a set of thermo parameters with the given end point temperature.
    void AddThermoParams(
        const real T, // Maximum temperature for which parameters are valid.
        const Sprog::Thermo::THERMO_PARAMS &params // Thermo params to add.
        );

    // Removes the parameters from the list valid for the given temperature.
    void RemoveThermoParams(const real T);


    // PARENT MECHANISM.

    // Returns pointer to parent mechanism.
    const Sprog::Mechanism *const Mechanism(void) const;

    // Sets the parent mechanism.
    void SetMechanism(Sprog::Mechanism &mech);


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the species object.
    Species *const Clone(void) const;

    // Writes the species to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the species data from a binary data stream.
    void Deserialize(std::istream &in);

protected:
    // Species data.
    std::string m_name;       // Name/symbol.
    ElCompVector m_elcomp;    // Elemental composition.
    real m_molwt;             // Molecular weight (kg/mol).
    Sprog::Mechanism *m_mech; // Parent mechanism.

    //const ElementPtrVector *m_elements; // Elements used to define species.

    // Thermo parameters for different temperature ranges.  The map key is the
    // end point temperature up to which the parameters are valid.
    Sprog::Thermo::ThermoMap m_thermoparams;
    real m_T1; // Start temperature for range.
};

// Inline function definitions.
#include "gpc_species_inl.h"

// A typedef for a STL vector of species.
typedef std::vector<Species> SpeciesVector;

// A typedef for a STL vector of pointers to species.
typedef std::vector<Species*> SpeciesPtrVector;
};

#endif
