/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a chemical mechanism.
*/

#ifndef GPC_MECH_H
#define GPC_MECH_H

#include <vector>
#include <string>
#include "gpc_params.h"
#include "gpc_element.h"
#include "gpc_species.h"
#include "gpc_reaction_set.h"
#include "gpc_unit_systems.h"

using namespace std;

namespace Sprog
{

struct FRXN_STOICH
{
    IndexedRxn Reaction;
    real Stoich;
};

class Mechanism
{
protected:
    UnitSystem m_units;       // The system of units used by this mechanism.
    ElementVector m_elements; // Vector of chemical elements defined by mechanism.
    SpeciesVector m_species;  // Vector of chemical species defined by mechanism.
    ReactionSet m_rxns;       // Set of reactions defined by mechanism.
    /* Reaction stoichiometry cross-referenced for each species. */
    map<IndexedSpecies, FRXN_STOICH> m_stoich_xref;
public:
    Mechanism(void);  // Default constructor.
    ~Mechanism(void); // Default destructor.
    Mechanism(const Mechanism &mech); // Copy constructor.
public:
    Mechanism &operator=(const Mechanism &mech);
    Mechanism &operator+=(const Mechanism &mech);
    const Mechanism operator+(const Mechanism &mech) const;
public:
    /* Returns the vector of elements. */
    const ElementVector &Elements(void) const;
    /*Returns the vector of species. */
    const SpeciesVector &Species(void) const;
    /* Returns the reaction set. */
    const ReactionSet &Reactions(void) const;
protected:
    /* Builds the species-reaction stoichiometry cross-reference table. */
    void BuildStoichCrossRef(void);
};
};

#endif