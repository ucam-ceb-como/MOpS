/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a chemical mechanism.
*/

#ifndef GPC_MECH_H
#define GPC_MECH_H

#include "gpc_params.h"
#include "gpc_element.h"
#include "gpc_species.h"
#include "gpc_reaction_set.h"
#include "gpc_unit_systems.h"
#include <vector>
#include <string>
#include <map>

namespace Sprog
{
class Mechanism
{
public:
    // Structure for species/reaction stoichiometry cross-referencing.
    typedef std::map<unsigned int, real> RxnStoichMap;
    typedef std::pair<unsigned int, real> RxnStoichPair;
    struct StoichXRef
    {
        unsigned int Species;
        RxnStoichMap RxnStoich;
    };
    typedef std::vector<StoichXRef> StoichXRefVector;

    // Constructors.
    Mechanism(void); // Default constructor.
    Mechanism(const Mechanism &mech); // Copy constructor.

    // Destructors.
    virtual ~Mechanism(void);

    // Operator overloads.
    Mechanism &operator=(const Mechanism &mech);

    // Units.
    UnitSystem Units(void); // Returns the current unit system.
    void SetUnits(UnitSystem u); // Converts the mechanism to a new units system.

    // Chemical elements.
    const ElementPtrVector &Elements(void) const; // Returns the vector of elements.
    Element *const AddElement(void);              // Adds a default element to the mechanism.
    Element *const AddElement(const Element &el); // Copies the given element into the mechanism.
    int FindElement(const std::string &name);  // Returns index of element.  Returns -1 if not found.
    int FindElement(const Element &el);        // Returns index of element.  Returns -1 if not found.

    // Element updates.
    void CheckElementChanges(const Element &el); // Updates mechanism with changes applied to an element.

    // Species.
    const SpeciesPtrVector &Species(void) const; // Returns the vector of species.
    Sprog::Species *const AddSpecies(void);      // Adds an empty species to the mechanism.
    Sprog::Species *const AddSpecies(const Sprog::Species &sp); // Copies given species into the mechanism.
    int FindSpecies(const Sprog::Species &sp) const; // Returns index of species.  Returns -1 if not found.
    int FindSpecies(const std::string &name) const;  // Returns index of species.  Returns -1 if not found.
    Sprog::Species *const GetSpecies(const unsigned int i) const; // Returns pointer to species at index i.  NULL if not found.
    Sprog::Species *const GetSpecies(const std::string &name) const;  // Returns pointer to species with given name.  NULL if not found.

    // Reactions.
    const Kinetics::ReactionSet &Reactions(void) const; // Returns the reaction set.
    Kinetics::Reaction *const AddReaction(void);        // Adds an empty reaction to the mechanism.
    Kinetics::Reaction *const AddReaction(const Kinetics::Reaction *const rxn); // Copies a reaction into the mechanism.

protected:
    // Mechanism data.
    UnitSystem m_units;             // The system of units used by this mechanism.
    ElementPtrVector m_elements;    // Vector of chemical elements defined by mechanism.
    SpeciesPtrVector m_species;     // Vector of chemical species defined by mechanism.
    Kinetics::ReactionSet m_rxns;   // Set of reactions defined by mechanism.
    StoichXRefVector m_stoich_xref; // Reaction stoichiometry cross-referenced for each species.

    // Copying routines.
    void copyInElements(const ElementPtrVector &els); // Copies elements from given array into this mechanism.
    void copyInSpecies(const SpeciesPtrVector &els);  // Copies species from given array into this mechanism.

    // Builds the species-reaction stoichiometry cross-reference table.
    void buildStoichCrossRef(void);

    // Memory management.
    void releaseMemory(void); // Clears memory used by the mechanism object.
};
};

#endif