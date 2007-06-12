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

using namespace std;

namespace Sprog
{


class Mechanism
{
public:
    // Structure for species/reaction stoichiometry cross-referencing.
    struct RxnStoichf
    {
        IndexedRxn Reaction;
        real Stoich;
    };

    // Constructors.
    Mechanism(void); // Default constructor.
    Mechanism(const Mechanism &mech); // Copy constructor.

    // Destructors.
    ~Mechanism(void);

    // Operator overloads.
    Mechanism &operator=(const Mechanism &mech);
    Mechanism &operator+=(const Mechanism &mech);
    const Mechanism operator+(const Mechanism &mech) const;

    // Chemical elements.
    const ElementVector &Elements(void) const; // Returns the vector of elements.
   
    // Species.
    const SpeciesVector &Species(void) const; // Returns the vector of species.
   
    // Reactions.
    const ReactionSet &Reactions(void) const; // Returns the reaction set.

protected:
    // Mechanism data.
    UnitSystem m_units;       // The system of units used by this mechanism.
    ElementVector m_elements; // Vector of chemical elements defined by mechanism.
    SpeciesVector m_species;  // Vector of chemical species defined by mechanism.
    ReactionSet m_rxns;       // Set of reactions defined by mechanism.
    map<IndexedSpecies, RxnStoichf> m_stoich_xref; // Reaction stoichiometry cross-referenced for each species.

    // Builds the species-reaction stoichiometry cross-reference table.
    void buildStoichCrossRef(void);
};
};

#endif