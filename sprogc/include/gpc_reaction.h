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
    Reaction(void); // Default constructor.
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
    const Stoich Reactant(unsigned int k) const;   // Returns the stoichiometry of the kth integer reactant.
    const Stoichf FReactant(unsigned int k) const; // Returns the stoichiometry of the kth real reactant.
    int ReactantCount() const;  // Returns the number of integer reactants.
    int FReactantCount() const; // Returns the number of real reactants.

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
    const Stoich Product(unsigned int k) const;    // Returns the stoichiometry of the kth integer product.
    const Stoichf FProduct(unsigned int k) const;  // Returns the stoichiometry of the kth real product.
    int ProductCount() const;   // Returns the number of integer products.
    int FProductCount() const;  // Returns the number of real product.

    real TotalStoich() const; // Returns the total stoichiometry of the reaction.
    real ReactantStoich() const; // Returns the reactant stoichiometry of the reaction.
    real ProductStoich() const;  // Returns the product stoichiometry of the reaction.

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

    // Third bodies.
    bool UseThirdBody() const;        // Returns true if this reaction uses third bodies.
    void SetUseThirdBody(bool usetb); // Sets whether or not this reaction uses third bodies.
    const std::vector<Stoichf> &ThirdBodies(void) const;    // Returns the vector of third-body coefficients.
    Stoichf ThirdBody(unsigned int i) const;                // Returns the coefficient for the ith third body.
    int ThirdBodyCount() const; // Returns the number of third body coefficients defined for this reaction.
    void AddThirdBody(const Stoichf &tb);                   // Adds a third body to the reaction.
    void AddThirdBody(unsigned int sp, real coeff);         // Adds a third body to the reaction.
    void AddThirdBody(const std::string &name, real coeff); // Adds a third body given the species name.
    void RemoveThirdBody(const std::string &name);          // Removes a third body, given by name, from 
                                                            // the reaction.

    // Low pressure limit.
    const ARRHENIUS &LowPressureLimit(void) const;   // Returns the low pressure limit Arrhenius coefficients.   
    void SetLowPressureLimit(const ARRHENIUS &lowp); // Sets the low pressure limit Arrhenius coefficients.

    // Fall-off third body.
    const Sprog::Species *const FallOffThirdBody(void) const; //   Returns a pointer to the species used as a 
                                                              // third body for fall-off calculations.
    void SetFallOffThirdBody(int sp);                   //   Sets the species used as a third body for 
                                                        // fall-off calculations.
    void SetFallOffThirdBody(const std::string &name);  //   Sets the species used as a third body for 
                                                        // fall-off calculations given the species name.

    // Fall-off parameters.
    FALLOFF_FORM FallOffType() const;                // Returns the fall-off type.
    const FALLOFF_PARAMS &FallOffParams(void) const; // Returns the fall-off parameter type.
    void SetFallOffParams(const FALLOFF_FORM form,   // Sets the fall-off type and parameters.
                          const real params[FALLOFF_PARAMS::MAX_FALLOFF_PARAMS]);

    // Fall-off functions.
    real FTROE3(real T, real logpr) const; // 3-parameter Troe form.
    real FTROE4(real T, real logpr) const; // 4-parameter Troe form.
    real FSRI(real T, real logpr) const;   // SRI form.
    FallOffFnPtr FallOffFn() const;        // Custom function.

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

    // Third bodies.
    bool m_usetb; // Set to true if this reaction requires third bodies.
    std::vector<Stoichf> m_thirdbodies; // Reaction third bodies and their coefficients.

    // Fall-off data.
    FALLOFF_FORM m_fotype;     // The type of the fall-off reaction.
    FALLOFF_PARAMS m_foparams; // Fall-off parameters.
    FallOffFnPtr m_fofn;       // Custom fall-off function, if required.

    // Useful data to put reaction in context.
    const SpeciesPtrVector *m_species; // Vector of species used to define the reaction.
    Sprog::Mechanism *m_mech;          // Parent mechanism.

    // Memory management.
    virtual void releaseMemory(void); // Releases all memory used by the reaction object.
};

// Inline function definitions.
#include "gpc_reaction_inl.h"

// A typedef for a STL vector of reactions.
typedef std::vector<Reaction> RxnVector;
// A typedef fo a STL vector of pointers to reactions.
typedef std::vector<Reaction*> RxnPtrVector;
};
};

#endif
