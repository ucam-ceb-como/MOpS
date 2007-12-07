/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains the definition of a structure for chemical reactions.  Also
    contains typedefs and other data structures related to chemical reactions.
*/

#ifndef GPC_REACTION_H
#define GPC_REACTION_H

#include "gpc_params.h"
#include "gpc_species.h"
#include "gpc_stoich.h"
#include "gpc_rate_params.h"
#include <vector>
#include <string>

namespace Sprog
{
class Mechanism; // Forward declaration of mechanism.

namespace Kinetics
{
class Reaction
{
public:
    // Constructors.
    Reaction(void);  // Default constructor.
    Reaction(const Reaction &rxn); // Copy constructor.

    // Destructor.
    virtual ~Reaction(void);

    // Operator overloads.
    Reaction &operator=(const Reaction &rxn);

    // Reaction name.
    const std::string &Name(void) const; // Returns the reaction name.
    void SetName(const std::string &name); // Sets the reaction name.

    // Reaction reversibility.
    bool IsReversible(void) const;        // Returns true if the reaction is reversible.
    void SetReversible(const bool isrev); // Sets the reaction to be reversible or not.

    // Reactants.
    const std::vector<Stoich> &Reactants(void) const;   // Returns the vector of integer stoichiometric
                                                        // reactant coefficients.
    const std::vector<Stoichf> &FReactants(void) const; // Returns the vector of real stoichiometric 
                                                        // reactant coefficients.
    void AddReactant(const Stoich &reac);  // Adds an integer reactant to the reaction.  
    void AddReactant(const Stoichf &reac); // Adds a real reactant to the reaction.
    void AddReactant(const std::string &name, unsigned int stoich); // Adds an integer reactant given
                                                                    // the species name.
    void AddReactant(const std::string &name, real stoich); // Adds a real reactant given the species name.
    void RemoveReactant(const std::string &name); // Removes a reactant, given by name, from 
                                                  // the reaction (integer or real).

    // Products.
    const std::vector<Stoich> &Products(void) const; // Returns the vector of integer stoichiometric 
                                                     // product coefficients.
    const std::vector<Stoichf> &FProducts(void) const; // Returns the vector of real stoichiometric 
                                                       // product coefficients.
    void AddProduct(const Stoich &prod);  // Adds an integer product to the reaction.   
    void AddProduct(const Stoichf &prod); // Adds a real product to the reaction.
    void AddProduct(const std::string &name, unsigned int stoich); // Adds an integer product given
                                                                   // the species name.
    void AddProduct(const std::string &name, real stoich); // Adds a real product given the species name.
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
    const LTCOEFFS *const RevLTCoeffs(void) const; // Returns a pointer to the reverse Landau Teller coefficients.   
    void SetRevLTCoeffs(const LTCOEFFS &lt);       // Sets the reverse Landau Teller coefficients.

    // Species vector.
    const SpeciesPtrVector *const Species(void) const; // Returns pointer to the vector of species used to 
                                                       // define reaction.
    void SetSpecies(const SpeciesPtrVector *const sp); // Sets the species vector.

    // Parent mechanism.
    Sprog::Mechanism *const Mechanism(void) const;   // Returns a pointer to the parent mechanism.
    void SetMechanism(Sprog::Mechanism *const mech); // Sets the parent mechanism.

    // Cloning.
    virtual Reaction *Clone(void) const; // Returns a pointer to a clone of the reaction.

protected:
    // Reaction data.
    std::string m_name;                    // Reaction description.
    bool m_reversible;                     // Is the reaction reversible or not?
    std::vector<Stoich> m_reac, m_prod;    // Integer reactant & product stoichiometry.
    std::vector<Stoichf> m_freac, m_fprod; // Real reactant & product stoichiometry.
    real m_dstoich, m_dreac, m_dprod;      // Total stoichiometry changes.
    ARRHENIUS m_arrf, *m_arrr;             // Forward and reverse Arrhenius parameters.
    LTCOEFFS *m_lt, *m_revlt;              // Landau-Teller forward and reverse coefficients.
   
    const SpeciesPtrVector *m_species; // Vector of species used to define the reaction.
    Sprog::Mechanism *m_mech;          // Parent mechanism.

    // Memory management.
    virtual void releaseMemory(void); // Releases all memory used by the reaction object.
};

// Inline function definitions.
#include "gpc_reaction_inl.h"

// A typedef for a STL vector of reactions.
typedef std::vector<Reaction> RxnVector;
typedef std::vector<Reaction*> RxnPtrVector;
};
};

#endif