/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a structure for a set of chemical reactions, 
    including code to provide speed enhancements when working with many reactions.
*/

#ifndef GPC_REACTION_SET_H
#define GPC_REACTION_SET_H

#include "gpc_params.h"
#include "gpc_reaction.h"
#include "gpc_gasphase.h"
#include <vector>
#include <map>
#include <iostream>

namespace Sprog
{
class Mechanism; // Forward declare Mechanism class.

namespace Kinetics
{
class ReactionSet
{
public:
    // Constructors.
    ReactionSet(void);                   // Default constructor.
    ReactionSet(const ReactionSet &rxn); // Copy constructor.
    ReactionSet(std::istream &in);       // Stream-reading constructor.

    // Destructor.
    ~ReactionSet(void);

    // Operator overloads.
    ReactionSet &operator=(const ReactionSet &rxns);
    ReactionSet &operator+=(const ReactionSet &rxns);
    const ReactionSet operator+(const ReactionSet &rxns) const;
    Reaction *const operator[](unsigned int i);
    const Reaction *const operator[](unsigned int i) const;
   

    // REACTIONS.

    // Returns the number of reactions in the set.
    unsigned int Count(void) const;

    // Returns the list of reactions.
    const RxnPtrVector &Reactions(void) const;

    // Returns a pointer to the ith reaction.  Returns NULL if i is invalid.
    const Reaction *const Reactions(unsigned int i) const;

    // Adds a reaction to the set.
    Reaction *const AddReaction(const Reaction &rxn);


    // TIDYING UP.

    // Clears all reactions from the set.
    void Clear(void);


    // SPECIES MOLAR PRODUCTION RATES.

    // Calculates the molar production rates of all species.
    void GetMolarProdRates(
        const fvector &rop, // Rate of Progress of each reaction.
        fvector &wdot       // Return vector for molar prod. rates.
        ) const;


    // REACTION RATES OF PROGRESS.

    // Calculates the rate of progress of each reaction.
    void GetRatesOfProgress(
        const Sprog::Thermo::GasPhase &mix, // The mixture for which to calculate the rates.
        const fvector &kforward,  // Forward rate constants of all reactions.
        const fvector &kreverse,  // Reverse rate constants of all reactions.
        fvector & rop             // Return vector for rates of progress.
        ) const;

    // Calculates the rate of progress of each reaction.
    void GetRatesOfProgress(
        real density,            // Mixture molar density.
        const real *const x,     // Species mole fractions.
        unsigned int n,          // Number of values in x array.
        const fvector &kforward, // Forward rate constants of all reactions.
        const fvector &kreverse, // Reverse rate constants of all reactions.
        fvector &rop             // Return vector for rates of progress.
        ) const;


    // RATE CONSTANTS.

    // Calculates the forward and reverse rate constants 
    // of all reactions given a mixture object.
    void GetRateConstants(
        const Sprog::Thermo::GasPhase &mix, // The mixture for which to calculate the rate constants.
        const fvector &Gs, // Dimensionless Gibbs free energy of each species (1/mol).
        fvector &kforward, // Return vector for forward rate constants.
        fvector &kreverse  // Return vector for reverse rate constants.
        ) const;

    // Calculates the forward and reverse rate constants 
    // of all reactions given mixture temperature, density
    // and species mole fractions.
    void GetRateConstants( 
        real T,              // The mixture temperature.
        real density,        // Mixture molar density.
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x array.
        const fvector &Gs,   // Dimensionless Gibbs free energy of each species (1/mol).
        fvector &kforward,   // Return vector for forward rate constants.
        fvector &kreverse    // Return vector for reverse rate constants.
        ) const;


    // PARENT MECHANISM.

    // Returns a pointer to the parent mechanism.
    const Sprog::Mechanism *const Mechanism() const;

    // Sets the parent mechanism.
    void SetMechanism(Sprog::Mechanism &mech);


    // READ/WRITE/COPY FUNCTIONS.

    // Writes the reaction set to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the reaction set data from a binary data stream.
    void Deserialize(std::istream &in);

protected:
    // In the following map the key is the index in the vector of all reactions.
    typedef std::map<unsigned int,const Reaction*> RxnMap;

    // Reaction set data.
    RxnPtrVector m_rxns; // Vector of all reactions in the set.
    RxnMap m_rev_rxns;   // Map of reactions which have explicit reverse Arrhenius parameters.
    RxnMap m_tb_rxns;    // Map of third body reactions.
    RxnMap m_fo_rxns;    // Map of fall-off reactions.
    RxnMap m_lt_rxns;    // Map of reactions with Landau Teller parameters.
    RxnMap m_revlt_rxns; // Map of reactions with reverse Landau Teller parameters.


    // MEMORY MANAGEMENT.

    // Clears all memory used by the set.
    void releaseMemory(void);

private:
    // Pointer to mechanism to which this ReactionSet belongs.
    Sprog::Mechanism *m_mech;
};
};
};

#endif
