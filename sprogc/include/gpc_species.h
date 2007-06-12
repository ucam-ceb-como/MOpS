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
    void AddElement(const ElComp &elcomp);       // Adds an element to the species composition.

    // Species molecular weight.
    const real MolWt(void) const;    // Returns the species molecular weight.
    void SetMolWt(const real molwt); // Sets the molecular weight.

    // Thermodynamic fitting parameters.
    unsigned int ThermoRangeCount(void) const;    // Returns the number of parameter ranges.
    void SetThermoStartTemperature(const real T); // Sets the start temperature for the range.  
    const THERMO_PARAMS &ThermoParams(const real T) const; // Returns the set of parameters valid for 
                                                           // the given temperature.
    void AddParams(const real T, const THERMO_PARAMS &params); // Adds a set of parameters with the given end 
                                                               // point temperature.
    void RemoveParams(const real T); // Removes the parameters from the list valid for the given temperature.

protected:
    std::string m_name;    // Name/symbol.
    ElCompVector m_elcomp; // Elemental composition.
    real m_molwt;          // Molecular weight (kg/mol).

    // Thermo parameters for different temperature ranges.  The map key is the
    // end point temperature up to which the parameters are valid.
    ThermoMap m_thermoparams;
    real m_T1; // Start temperature for range.
};

// Inline function definitions.
#include "gpc_species_inl.h"

// A typedef for a STL vector of species.
typedef std::vector<Species> SpeciesVector;
typedef std::vector<Species*> SpeciesPtrVector;

// A data structure used for internal data referencing and calculations.  A pointer
// to a Species object is maintained with an index for this species, which will 
// normally refer to its location in a vector or array.
struct IndexedSpecies
{
public:
    // Data.
    Species *Pointer;
    unsigned int Index;

    // Constructors.
    IndexedSpecies(void); // Default constructor.
    IndexedSpecies(const IndexedSpecies &sp); // Copy constructor.

    // Destructor.
    ~IndexedSpecies(void);

    // Operator overloads.
    IndexedSpecies &operator=(const IndexedSpecies &sp);
    bool operator==(const IndexedSpecies &sp) const;
    bool operator!=(const IndexedSpecies &sp) const;
};
};

#endif