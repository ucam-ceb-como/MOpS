/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a structure for chemical reactions of different
    types.
*/

#ifndef GPC_REACTION_H
#define GPC_REACTION_H

#include "gpc_params.h"
#include "gpc_species.h"
#include <vector>
#include <string>

namespace Sprog
{
// Structure to hold integer species stoichiometry for a reaction.
struct STOICH
{
    IndexedSpecies Species;
    int Stoich;
};

// Structre to hold real species stoichiometry for a reaction.
struct FSTOICH
{
    IndexedSpecies Species;
    real Stoich;
};

// Reaction Arrhenius parameters.
struct ARRHENIUS
{
    real A; // Pre-exponential factor.
    real n; // Temperature exponent.
    real E; // Activation energy.
};

// Landau Teller reaction parameters.
struct LTCOEFFS
{
    real B, C;
};


class Reaction
{
public:
    // Constructors.
    Reaction(void);  // Default constructor.
    Reaction(const Reaction &rxn); // Copy constructor.

    // Destructor.
    ~Reaction(void);

    // Operator overloads.
    Reaction &operator=(const Reaction &rxn);

    // Reaction name.
    const std::string &Name(void) const; // Returns the reaction name.
    void SetName(const std::string &name); // Sets the reaction name.

    // Reaction reversibility.
    bool IsReversible(void) const;        // Returns true if the reaction is reversible.
    void SetReversible(const bool isrev); // Sets the reaction to be reversible or not.

    // Reactants.
    const std::vector<STOICH> &Reactants(void) const;   // Returns the vector of integer stoichiometric
                                                        // reactant coefficients.
    const std::vector<FSTOICH> &FReactants(void) const; // Returns the vector of real stoichiometric 
                                                        // reactant coefficients.
    void AddReactant(const STOICH &reac);  // Adds an integer reactant to the reaction.  
    void AddReactant(const FSTOICH &reac); // Adds a real reactant to the reaction. 
    void RemoveReactant(const std::string &name); // Removes a reactant, given by name, from 
                                                  // the reaction (integer or real).

    // Products.
    const std::vector<STOICH> &Products(void) const; // Returns the vector of integer stoichiometric 
                                                     // product coefficients.
    const std::vector<FSTOICH> &FProducts(void) const; // Returns the vector of real stoichiometric 
                                                       // product coefficients.
    void AddProduct(const STOICH &prod);  // Adds an integer product to the reaction.   
    void AddProduct(const FSTOICH &prod); // Adds a real product to the reaction.
    void RemoveProduct(const std::string &name); // Removes a product, given by name, from the 
                                                 // reaction (integer or real).

    // Forward Arrhenius coefficients.
    const ARRHENIUS &Arrhenius(void) const;  // Returns the forward Arrhenius rate parameters.
    void SetArrhenius(const ARRHENIUS &arr); // Sets the forward Arrhenius rate parameters.

    // Reverse Arrhenius coefficients.
    const ARRHENIUS *const RevArrhenius(void) const; // Returns a pointer to the reverse Arrhenius parameters.
    void SetRevArrhenius(const ARRHENIUS &arr);      // Sets the reverse Arrhenius parameters.

    // Forward Landau Teller parameters.
    const LTCOEFFS *const LTCoeffs(void) const; // Returns a pointer to the forward Landau Teller coefficients.   
    void SetLTCoeffs(const LTCOEFFS &lt);       // Sets the forward Landau Teller coefficients.

    // Reverse Landau Teller parameters.
    const LTCOEFFS *const RefLTCoeffs(void) const; // Returns a pointer to the reverse Landau Teller coefficients.   
    void SetRevLTCoeffs(const LTCOEFFS &lt);       // Sets the reverse Landau Teller coefficients.

protected:
    // Reaction data.
    std::string m_name;                    // Reaction description.
    bool m_reversible;                     // Is the reaction reversible or not?
    std::vector<STOICH> m_reac, m_prod;    // Integer reactant & product stoichiometry.
    std::vector<FSTOICH> m_freac, m_fprod; // Real reactant & product stoichiometry.
    real m_dstoich, m_dreac, m_dprod;      // Total stoichiometry changes.
    ARRHENIUS m_arrf, *m_arrr;             // Forward and reverse Arrhenius parameters.
    LTCOEFFS *m_lt, *m_revlt;              // Landau-Teller forward and reverse coefficients.
   
    // Memory management.
    void releaseMemory(void); // Releases all memory used by the reaction object.
};

// A typedef for a STL vector of reactions.
typedef std::vector<Reaction> RxnVector;
typedef std::vector<Reaction*> RxnPtrVector;

// A data structure used for internal data referencing and calculations.  A pointer
// to a Reaction object is maintained with an index for this reaction, which will 
// normally refer to its location in a vector or array.
struct IndexedRxn
{
public:
    // Data.
    Reaction *Pointer;
    unsigned int Index;

    // Constructors.
    IndexedRxn(void); // Default constructor.
    IndexedRxn(const IndexedRxn &rxn); // Copy constructor.

    // Destructor.
    ~IndexedRxn(void);

    // Operator overloads.
    IndexedRxn &operator=(const IndexedRxn &rxn);
    bool operator==(const IndexedRxn &rxn) const;
    bool operator!=(const IndexedRxn &rxn) const;
};
};

#endif