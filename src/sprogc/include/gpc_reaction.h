/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This files contains the definition of a structure for chemical reactions.  Also
    contains typedefs and other data structures related to chemical reactions.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef GPC_REACTION_H
#define GPC_REACTION_H

#include "gpc_params.h"
#include "gpc_species.h"
#include "gpc_stoich.h"
#include "gpc_delta_stoich.h"
#include "gpc_rate_params.h"
#include <vector>
#include <string>
#include <iostream>

namespace Sprog
{
class Mechanism; // Forward declaration of mechanism.

namespace Kinetics
{
class Reaction
{
public:
    // Constructors.
    Reaction(void);                // Default constructor.
    Reaction(const Reaction &rxn); // Copy constructor.
    Reaction(std::istream &in);

    // Destructor.
    ~Reaction(void);

    // Operator overloads.
    Reaction &operator=(const Reaction &rxn);

    // REACTION NAME.

    // Returns the reaction name.
    const std::string &Name(void) const;

    // Sets the reaction name.
    void SetName(const std::string &name);


    // REVERSIBILITY.

    // Returns true if the reaction is reversible.
    bool IsReversible(void) const;

    // Sets the reaction to be reversible or not.
    void SetReversible(const bool isrev);


    // REACTANTS.

    // Returns the vector of integer stoichiometric reactant coefficients.
    const std::vector<Stoich> &Reactants(void) const;

    // Returns the vector of integer stoichiometric reactant coefficients of a given phaseName (type).
    // const std::vector<Stoich> &Reactants(const Mechanism::Phase &ph) const;

    // Adds a reactant to the reaction.
    void AddReactant(const Stoich &reac, const std::string &spName);

    // Adds a reactant given the species name.
    void AddReactant(const std::string &name, double stoich);

    // Removes a reactant, given by name, from the reaction.
    void RemoveReactant(const std::string &name);

    // Returns the stoichiometry of the kth reactant.
    const Stoich Reactant(unsigned int k) const;

    // Returns the number of integer reactants.
    int ReactantCount() const;


    // PRODUCTS.

    // Returns the vector of integer stoichiometric product coefficients.
    const std::vector<Stoich> &Products(void) const;

    // Adds a product to the reaction.
    void AddProduct(const Stoich &prod, const std::string &spName);

    // Adds a product given the species name.
    void AddProduct(const std::string &name, double stoich);

    // Removes a product, given by name, from the reaction.
    void RemoveProduct(const std::string &name);

    // Returns the stoichiometry of the kth integer product.
    const Stoich Product(unsigned int k) const;

    // Returns the number of integer products.
    int ProductCount() const;


    // STOICHIOMETRY.

    
    
    // Returns the total stoichiometry of the reaction.
    double TotalStoich() const;

    // Returns the reactant stoichiometry of the reaction.
    double ReactantStoich() const;

    // Returns the product stoichiometry of the reaction.
    double ProductStoich() const;
    

	// STOICH CHANGES in each phase
	/*
	* Added by mm864
	*
	*/

// Returns the number of delta stoich.
    unsigned int PhaseCount(void) const;

 // Returns phase at index i.  NULL if not found.  This
    // function returns a delta stoich object.
    const std::string GetPhaseName(const unsigned int i) const;

// Adds an PHASE vector to the reaction.
 void AddPhaseName(const std::string &phName);

// Find an PHASE vector to the reaction.
 int FindPhaseName(const std::string &phName) const;


	
    // Returns the number of phase in the reaction.
    unsigned int DeltaStoichCount(void) const;

	 // Returns the vector of delta stoich.
    const DeltaStoichPtrVector &DeltaStoich(void) const;

    // Returns a pointer to the ith delta stoich.  Returns NULL if i is invalid.
    const Sprog::Kinetics::DeltaStoich *const DeltaStoich(unsigned int i) const;

    // Returns pointer to delta stoich with given name.  NULL if not found.
    Sprog::Kinetics::DeltaStoich *const DeltaStoich(const std::string &name) const;


    // Returns index of delta stoich.  Returns -1 if not found.
    int FindDeltaStoich(const std::string &name) const;

    // Adds an empty  delta stoich to the reaction.
    Sprog::Kinetics::DeltaStoich *const AddDeltaStoich(const std::string &spName);

    // Copies given delta stoich into the reaction.
    Sprog::Kinetics::DeltaStoich *const AddDeltaStoich(const Sprog::Kinetics::DeltaStoich &delta_st, const std::string &spName);

    // Returns pointer to delta stoich at index i.  NULL if not found.  This
    // function returns a modifiable (non-const) delta stoich object.
    Sprog::Kinetics::DeltaStoich *const GetDeltaStoich(const unsigned int i) const;

    // Returns pointer to delta stoich with given name.  NULL if not found.  Thi
    // function returns a modifiable (non-const) delta stoich object.
    Sprog::Kinetics::DeltaStoich *const GetDeltaStoich(const std::string &name) const;
    
	
    // FORWARD ARRHENIUS COEFFICIENTS.

    // Returns the forward Arrhenius rate parameters.
    const ARRHENIUS &Arrhenius(void) const;

    // Sets the forward Arrhenius rate parameters.
    void SetArrhenius(const ARRHENIUS &arr);


    // REVERSE ARRHENIUS COEFFICIENTS.

    // Returns a pointer to the reverse Arrhenius parameters.
    const ARRHENIUS *const RevArrhenius(void) const;

    // Sets the reverse Arrhenius parameters.
    void SetRevArrhenius(const ARRHENIUS &arr);

	
	// FORD
	// Returns true if this reaction uses ford.
	bool IsFORD() const;
	
	// Sets whether or not this reaction uses ford.
	void SetUseFORD(const bool isFord); 
	
	// Returns a pointer to the Ford parameters.
	const std::vector<FORD> &FordReac(void) const;

	// Sets the Ford parameters.
	void SetFord(const double c, const std::string &name);

	// Returns the coefficient for the ith ford of the reaction.
	FORD FORDElement(unsigned int i) const;

	
	
	// Count the number of ford params 
	int FORDCount() const;
	
	// SURFACE 
	// Returns true if one of the reactants or products is surface type
	bool IsSURF() const;

	// Sets whether or not this reaction is surface types
	void SetIsSurface(const std::string &phName);  

	// STICK
	// Returns true if this reaction uses stick.
	bool IsSTICK() const;
	
	// Sets whether or not this reaction uses stick.
	void SetUseSTICK(const bool isSticking);
	
	// MWON
	// Returns true if this reaction uses Mott-Wise.
	bool IsMottWise() const;
	// Sets whether or not this reaction uses Mott-Wise.
	void SetUseMottWise(const bool isMott);
	
	
	// COVERAGE 
	// Returns true if this reaction uses coverage.
	bool IsCOVERAGE() const; // (inline)
	
	// Sets whether or not this reaction uses coverage.
	void SetUseCOV(const bool isCov); // (inline)
	
	// Returns the vector of the Coverage parameters.
	const std::vector<COVERAGE> &CoverageReac(void) const;

	// Returns the coefficient for the ith coverage.
	COVERAGE CoverageElement(unsigned int i) const;

	// Sets the Coverage parameters.
	void SetCoverage(const double e, const double m, const double eps, const std::string &name);

	// Count the number of coverage params 
	int COVERAGECount() const;
	
    // FORWARD LANDAU-TELLER PARAMETERS.

    // Returns a pointer to the forward Landau Teller coefficients.
    const LTCOEFFS *const LTCoeffs(void) const;

    // Sets the forward Landau Teller coefficients.
    void SetLTCoeffs(const LTCOEFFS &lt);


    // REVERSE LANDAU-TELLER PARAMETERS.

    // Returns a pointer to the reverse Landau Teller coefficients.
    const LTCOEFFS *const RevLTCoeffs(void) const;

    // Sets the reverse Landau Teller coefficients.
    void SetRevLTCoeffs(const LTCOEFFS &lt);


    // THIRD BODIES.

    // Returns true if this reaction uses third bodies.
    bool UseThirdBody() const; // (inline)

    // Sets whether or not this reaction uses third bodies.
    void SetUseThirdBody(bool usetb); // (inline)

    // Returns the vector of third-body coefficients.
    const std::vector<Stoich> &ThirdBodies(void) const;

    // Returns the coefficient for the ith third body.
    Stoich ThirdBody(unsigned int i) const;

    // Returns the number of third body coefficients
    // defined for this reaction.
    int ThirdBodyCount() const;

    // Adds a third body to the reaction.
    void AddThirdBody(const Stoich &tb);

    // Adds a third body to the reaction.
    void AddThirdBody(unsigned int sp, double coeff);

    // Adds a third body given the species name.
    void AddThirdBody(const std::string &name, double coeff);

    // Removes a third body, given by name, from the reaction.
    void RemoveThirdBody(const std::string &name);


    // LOW-PRESSURE LIMIT.

    // Returns the low pressure limit Arrhenius coefficients.
    const ARRHENIUS &LowPressureLimit(void) const;

    // Sets the low pressure limit Arrhenius coefficients.
    void SetLowPressureLimit(const ARRHENIUS &lowp);


    // FALL-OFF THIRD BODY.

    // Returns a pointer to the species used as a third
    // body for fall-off calculations.
    const Sprog::Species *const FallOffThirdBody(void) const;

    // Sets the species used as a third body for fall-off calculations.
    void SetFallOffThirdBody(int sp);

    // Sets the species used as a third body for fall-off calculations
    // given the species name.
    void SetFallOffThirdBody(const std::string &name);


    // FALL-OFF PARAMETERS.

    // Returns the fall-off type.
    FALLOFF_FORM FallOffType() const;

    // Returns the fall-off parameter type.
    const FALLOFF_PARAMS &FallOffParams(void) const;

    // Sets the fall-off type and parameters.
    void SetFallOffParams(
        const FALLOFF_FORM form,
        const double params[FALLOFF_PARAMS::MAX_FALLOFF_PARAMS]
        );


    // FALL-OFF FUNCTIONS.

    double FTROE3(double T, double logpr) const; // 3-parameter Troe fall-off form.
    double FTROE4(double T, double logpr) const; // 4-parameter Troe fall-off form.
    double FSRI(double T, double logpr) const;   // SRI fall-off form.
    //FallOffFnPtr FallOffFn() const;        // Custom fall-off function.


    // PARENT MECHANISM.

    // Returns a pointer to the parent mechanism.  The mechanism
    // is stored in order for the reaction to be aware of the
    // chemical species available.
    const Sprog::Mechanism *const Mechanism(void) const;

    // Sets the parent mechanism.  The mechanism
    // is stored in order for the reaction to be aware of the
    // chemical species available.
    void SetMechanism(Sprog::Mechanism &mech);


    // RATE CALCULATION.

    // Calculates the rate of progress of this reaction.
    double RateOfProgress(
        double density,        // Mixture molar density.
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x array.
        double kforward,       // Forward rate constant.
        double kreverse        // Reverse rate constant.
        ) const;


    // READ/WRITE/COPY FUNCTIONS.

	 //! Indicates whether the pre-exponential should be converted to cgs
    bool ConvertPreexponential(void) const;
	
    // Creates a copy of the species object.
    Reaction *Clone(void) const;

    // Writes the element to a binary data stream.
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
        ar.template register_type<ARRHENIUS>();
	//ar.template register_type<COVERAGE>();
	//ar.template register_type<FORD>();
        ar.template register_type<Sprog::Stoichiometry>();
        ar.register_type(static_cast<LTCOEFFS *>(NULL));

        ar & m_name & m_reversible
           & m_reac & m_prod & m_dstoich & m_dreac & m_dprod
           & m_arrf & m_arrr & m_lt & m_revlt & m_fo & m_covr & m_usetb
	   & m_thirdbodies & m_fotype & m_foparams & m_isSurface & m_sticking & m_mottwise & m_isCoverage & m_coverage & m_isFord & m_ford & m_mech & m_deltaStoich & m_phaseVector;
    }

    // Writes the element to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the element data from a binary data stream.
    void Deserialize(std::istream &in);

    // Prints a diagnostic output file containing all the
    // reaction data.  This is used to debug.
    void WriteDiagnostics(std::ostream &out) const;

   

    //! Writes the reduced mechanism reactions to the output file
    void WriteReducedMechReacs(std::ostream &out, std::vector<std::string> RejectSpecies) const;

protected:
    // Reaction data.
    std::string m_name;                    // Reaction description.
    bool m_reversible;                     // Is the reaction reversible or not?
    bool m_sticking;                       // Is the reaction involves sticking reaction 
    bool m_mottwise;                       // Is the MottWise evaluation used to solve the sticking reaction  
    bool  m_isCoverage;
    bool m_isFord;
    bool m_isSurface;                      // Is the reaction surface reaction of any types.  
    std::vector<Stoich> m_reac, m_prod;    // reactant & product stoichiometry.
    std::vector<std::string> m_phaseVector;       // Vector of phase name for delta stoich
    DeltaStoichPtrVector m_deltaStoich;  // Vector of stoichiometric changes defined by mechanism. 
    double m_dstoich, m_dreac, m_dprod;      // Total stoichiometry changes.
    ARRHENIUS m_arrf, *m_arrr;             // Forward and reverse Arrhenius parameters.
    LTCOEFFS *m_lt, *m_revlt;              // Landau-Teller forward and reverse coefficients.
    COVERAGE m_covr; // for Deserialisation and Serialisation 

    FORD m_fo; // for Deserialisation and Serialisation     
    // Third bodies.
    bool m_usetb; // Set to true if this reaction requires third bodies.
    std::vector<Stoich> m_thirdbodies; // Reaction third bodies and their coefficients.

    

    // Ford reaction 

    std::vector<FORD> m_ford; // Vector of FORD struct, storing species name and its ford coefficient 
    
    // Sticking reaction 

    std::vector<COVERAGE> m_coverage; // Vector of COVERAGE struct, storing species name and their parameters 
    
    // Fall-off data.
    FALLOFF_FORM m_fotype;     // The type of the fall-off reaction.
    FALLOFF_PARAMS m_foparams; // Fall-off parameters.
    //FallOffFnPtr m_fofn;       // Custom fall-off function, if required.

    // Useful data to put reaction in context.
    Sprog::Mechanism *m_mech; // Parent mechanism.

    // MEMORY MANAGEMENT.

    // Releases all memory used by the reaction object.
    void releaseMemory(void);
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
