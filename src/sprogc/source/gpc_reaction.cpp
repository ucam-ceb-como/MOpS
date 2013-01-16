/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Reaction class declared in the
    gpc_reaction.h header file.

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

#include "gpc_params.h"
#include "gpc_stoich.h"
#include "gpc_rate_params.h"
#include "gpc_reaction.h"
#include "gpc_mech.h"
#include <exception>
#include <stdexcept>
#include <string>
#include <math.h>
#include "string_functions.h"

using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Reaction::Reaction(void)
{
    // Reaction data.
    m_name       = "";
    m_reversible = false;
    m_dstoich    = m_dreac = m_dprod = 0.0;
    m_arrf.A     = m_arrf.n = m_arrf.E = 0.0;
    m_arrr       = NULL;
    m_lt         = NULL;
    m_revlt      = NULL;
    m_fo.F_k     = 0.0;
    m_fo.spName  = "";
    m_covr.Eta   = m_covr.Miu = m_covr.Epsilon =0.0;  
    m_covr.spName = "";
    // Third-bodies.
    m_usetb = false;
    m_thirdbodies.clear();

    // Fall-off data.
    m_fotype   = None;
    m_foparams = FALLOFF_PARAMS();
    //m_fofn     = NULL;
    // Surface
    m_isSurface = false; 
    // Stick
    m_sticking = false;
    // Mottwise
    m_mottwise = false; 

    // Coverage data
    m_isCoverage = false;
    m_coverage.clear();

    // Ford data 
    m_isFord = false; 
    m_ford.clear();

    // Reaction context.
    m_mech = NULL;
    m_deltaStoich.clear();
    m_phaseVector.clear();
}

// Stream-reading constructor.
Reaction::Reaction(std::istream &in)
{
    // Reaction data.
    m_name       = "";
    m_reversible = false;
    m_dstoich    = m_dreac = m_dprod = 0.0;
    m_arrf.A     = m_arrf.n = m_arrf.E = 0.0;
    m_arrr       = NULL;
    m_lt         = NULL;
    m_revlt      = NULL;
    m_fo.F_k     = 0.0;
    m_fo.spName  = "";
    m_covr.Eta   = m_covr.Miu = m_covr.Epsilon =0.0;  
    m_covr.spName = "";

    // Third-bodies.
    m_usetb = false;
    m_thirdbodies.clear();

    // Fall-off data.
    m_fotype   = None;
    m_foparams = FALLOFF_PARAMS();
    //m_fofn     = NULL;

    // Surface
    m_isSurface = false; 
    // Stick
    m_sticking = false;
    // Mottwise
    m_mottwise = false; 

    // Coverage data
    m_isCoverage = false;
    m_coverage.clear();

    // Ford data 
    m_isFord = false; 
    m_ford.clear();

    // Reaction context.
    m_mech = NULL;
    m_deltaStoich.clear();
    m_phaseVector.clear();
    Deserialize(in);
}

// Copy constructor.
Reaction::Reaction(const Sprog::Kinetics::Reaction &rxn)
{
    // Nullify pointer members.
    m_arrr    = NULL;
    m_lt      = NULL;
    m_revlt   = NULL;
    m_mech    = NULL;

    // Use assignment operator to copy.
    *this = rxn;
}

// Destructor.
Reaction::~Reaction(void)
{
    // Clear memory associated with the reaction object.
    releaseMemory();
}

// OPERATOR OVERLOADING.

// Assignment operator.
Reaction &Reaction::operator=(const Sprog::Kinetics::Reaction &rxn)
{
    // Check for self assignment!
    if (this != &rxn) {
        // Clear current memory used.
        releaseMemory();

        // Copy name and reversibility.
        m_name = rxn.m_name;
        m_reversible = rxn.m_reversible;

        // Copy products and reactants.
        m_reac.assign(rxn.m_reac.begin(), rxn.m_reac.end());
        m_prod.assign(rxn.m_prod.begin(), rxn.m_prod.end());

        // Copy total stoichiometries.
        m_dstoich = rxn.m_dstoich;
        m_dreac = rxn.m_dreac;
        m_dprod = rxn.m_dprod;

        // Copy Arrhenius coefficients (forward and reverse).
        m_arrf = rxn.m_arrf;
        if (rxn.m_arrr != NULL) m_arrr = new ARRHENIUS(*rxn.m_arrr);

        // Copy Landau Teller parameters (forward and reverse).
        if (rxn.m_lt != NULL) m_lt = new LTCOEFFS(*rxn.m_lt);
        if (rxn.m_revlt != NULL) m_revlt = new LTCOEFFS(*rxn.m_revlt);

	// Copy FORD and COVERAGE 
        m_covr = rxn.m_covr;
	m_fo = rxn.m_fo;

        // Copy third bodies.
        m_usetb = rxn.m_usetb;
        m_thirdbodies.assign(rxn.m_thirdbodies.begin(), rxn.m_thirdbodies.end());

        // Copy fall-off data.
        m_fotype = rxn.m_fotype;
        m_foparams = rxn.m_foparams;
        //m_fofn = rxn.m_fofn;
	// All surface reactions
	m_isSurface = rxn.m_isSurface;
	 // Stick
	m_sticking = rxn.m_sticking;
	// Mottwise
	m_mottwise = rxn.m_mottwise; 

	// Coverage data
	m_isCoverage = rxn.m_isCoverage;
	m_coverage.assign(rxn.m_coverage.begin(), rxn.m_coverage.end());

	// Ford data 
	m_isFord = rxn.m_isFord; 
	m_ford.assign(rxn.m_ford.begin(), rxn.m_ford.end());
	
        // Copy pointer to mechanism.
        m_mech = rxn.m_mech;
	m_deltaStoich.assign(rxn.m_deltaStoich.begin(), rxn.m_deltaStoich.end());
	m_phaseVector.assign(rxn.m_phaseVector.begin(), rxn.m_phaseVector.end());
    }

    return *this;
}


// REACTANTS.

// Adds an integer stoichiometry reactant to the reaction.
void Reaction::AddReactant(const Sprog::Stoich &reac, const std::string &spName)
{
    // We must check if the reactant is already defined in this reaction.  Need to
    // check both integer and double reactants.

    // Add the new contribution to the total stoichiometry sums.
    m_dstoich -= reac.Mu();
    m_dreac += reac.Mu();
    AddDeltaStoich(spName);
    const double min_val = -(reac.Mu());
    DeltaStoich(spName)->IncrementTotalStoich(min_val);
    const double val = reac.Mu();
    DeltaStoich(spName)->IncrementReacStoich(val);

    // reactants.
    vector<Stoich>::iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); i++) {
        if ((*i).Index() == reac.Index()) {
            // We have found this species in the reactants already.  Append
            // the new contribution.
            (*i).IncMu(reac.Mu());
            return;
        }
    }

    // If we have got here then the reactant is not defined for this reaction, so
    // we must add it.
    m_reac.push_back(reac);
}


// Adds a reactant to the reaction given the species name.
void Reaction::AddReactant(const std::string &name, double stoich)
{
    if (m_mech != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_mech->SpeciesCount(); i++) {
            if (name.compare(m_mech->Species(i)->Name()) == 0) {
                // Found species in the list, now add it as a reactant.
                AddReactant(Stoich(i, stoich), m_mech->Species(i)->Name());
                AddPhaseName(m_mech->Species(i)->PhaseName());
		SetIsSurface(m_mech->Species(i)->PhaseName()); 
		return;
            }
        }

        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in species "
                           "list (Sprog, Reaction::AddReactant)."));
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Trying to add reactant before assigning "
                          "parent mechanism (Sprog, Reaction::AddReactant).");
    }
}


// Removes the reactant species with the given name from the reaction.  If the
// species is not a reactant then does nothing.
void Reaction::RemoveReactant(const std::string &name)
{
    if (m_mech != NULL) {
        // Search through the list of species to find that with
        // the given name.
        int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species:  Loop though integer stoichiometry
            // reactant to find that with this index.

            vector<Stoich>::iterator j;
            for (j=m_reac.begin(); j!=m_reac.end(); j++) {
                if ((*j).Index() == i) {
                    // Remove contribution from the total stoichiometry sums.
                    m_dstoich += (*j).Mu();
                    m_dreac -= (*j).Mu();

		    const double val = (*j).Mu();
		    DeltaStoich(m_mech->Species(i)->Name())->IncrementTotalStoich(val);
		    const double min_val = -((*j).Mu());
		    DeltaStoich(m_mech->Species(i)->Name())->IncrementReacStoich(min_val);

                    // We have found the species in the list.
                    m_reac.erase(j);
                }
            }

        }
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Can't search for species before assigning parent "
                          "mechanism (Sprog, Reaction::RemoveReactant).");
    }
}

// Returns the stoichiometry of the kth integer reactant.
const Stoich Reaction::Reactant(unsigned int k) const
{
    if (k < m_reac.size()) {
        return m_reac[k];
    } else {
        return Stoich(0,0);
    }
}

// Returns the number of integer reactants.
int Reaction::ReactantCount() const
{
    return m_reac.size();
}

// PRODUCTS.

// Adds an integer stoichiometry product to the reaction.
void Reaction::AddProduct(const Sprog::Stoich &prod, const std::string &spName)
{
    // We must check if the product is already defined in this reaction.  Need to
    // check both integer and double products.

    // Add the new contribution to the total stoichiometry sums.
    m_dstoich += prod.Mu();
    m_dprod += prod.Mu();
    AddDeltaStoich(spName);
    const double val = (prod.Mu());
    DeltaStoich(spName)->IncrementTotalStoich(val);
    DeltaStoich(spName)->IncrementProdStoich(val);

    // Integer products.
    vector<Stoich>::iterator i;
    for (i=m_prod.begin(); i!=m_prod.end(); i++) {
        if ((*i).Index() == prod.Index()) {
            // We have found this species in the products already.  Append
            // the new contribution.
            (*i).IncMu(prod.Mu());
            return;
        }
    }

    // If we have got here then the product is not defined for this reaction, so
    // we must add it.
    m_prod.push_back(prod);
}



// Adds an integer product to the reaction given the species name.
void Reaction::AddProduct(const std::string &name, double stoich)
{
    if (m_mech != NULL) {
        // Find the species in the list by name.
        unsigned int i;
        for (i=0; i<m_mech->SpeciesCount(); i++) {
            if (name.compare(m_mech->Species(i)->Name()) == 0) {
                // Found species in the list, now add it as a product.
                // AddProduct(Stoich(i, stoich)); (changed to below by mm864)
		AddProduct(Stoich(i, stoich), m_mech->Species(i)->Name());
		AddPhaseName(m_mech->Species(i)->PhaseName()); 
		SetIsSurface(m_mech->Species(i)->PhaseName()); 
                return;
            }
        }

        // We have got here because the species name was not found
        // in the species list.
        throw out_of_range(string(name).append(" was not found in "
                           "species list (Sprog, Reaction::AddProduct)."));
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Trying to add reactant before assigning "
                          "parent mechanism (Sprog, Reaction::AddProduct).");
    }
}

// Removes the product species with the given name from the reaction.  If the
// species is not a product then does nothing.
void Reaction::RemoveProduct(const std::string &name)
{
    if (m_mech != NULL) {
        // Search through the list of species to find that with
        // the given name.
        int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species:  Loop though integer stoichiometry
            // product to find that with this index.
            vector<Stoich>::iterator j;
            for (j=m_prod.begin(); j!=m_prod.end(); j++) {
                if ((*j).Index() == i) {
                    // Remove contribution from the total stoichiometry sums.
                    m_dstoich -= (*j).Mu();
                    m_dprod -= (*j).Mu();
		    const double min_val = -((*j).Mu());
		    DeltaStoich(m_mech->Species(i)->Name())->IncrementTotalStoich(min_val);
		    DeltaStoich(m_mech->Species(i)->Name())->IncrementProdStoich(min_val);

                    // We have found the species in the list.
                    m_prod.erase(j);
                }
            }
        }
    } else {
        // Parent mechanism has not been assigned.
        throw logic_error("Can't search for species before assigning "
                          "parent mechanism (Sprog, Reaction::RemoveProduct).");
   }
}

// Returns the stoichiometry of the kth integer product.
const Stoich Reaction::Product(unsigned int k) const
{
    if (k < m_prod.size()) {
        return m_prod[k];
    } else {
        return Stoich(0,0);
    }
}

// Returns the number of integer products.
int Reaction::ProductCount() const
{
    return m_prod.size();
}

	// STOICH CHANGES in each phase
	/*
	* Added by mm864
	*
	*/

// Returns the number of delta stoich.
unsigned int Reaction::PhaseCount() const
{
    return m_phaseVector.size();
}

 // Returns phase at index i.  NULL if not found.  This
    // function returns a modifiable (non-const) delta stoich object.
const std::string Reaction::GetPhaseName(const unsigned int i) const {
  if (i < m_phaseVector.size()) {
        return m_phaseVector[i];
  } else {
        return "";
  }
}

// Adds an PHASE vector to the reaction.
void Reaction::AddPhaseName(const std::string &phName)
{
    // We must check if the phase is already defined in this reaction. 
  int i = FindPhaseName(phName);
  
     if (i < 0) {
        // A phase with this name has been recorded. 

        m_phaseVector.push_back(phName);
    }  
   
}

// Adds an PHASE vector to the reaction.
int Reaction::FindPhaseName(const std::string &phName) const
{
  // Loop over phases to find index.
    unsigned int i;
    for (i=0; i<m_phaseVector.size(); i++) {
      if (phName.compare(m_phaseVector[i])== 0) {
            // Found phase!
            return i;
        }
    }

    // We are here because the phase wasn't found.
    return -1;

}


// Returns the number of delta stoich.
unsigned int Reaction::DeltaStoichCount() const
{
    return m_deltaStoich.size();
}

// Returns the vector of delta stoich.
const DeltaStoichPtrVector &Reaction::DeltaStoich(void) const
{
return m_deltaStoich;
}
    // Returns a pointer to the ith delta stoich.  Returns NULL if i is invalid.
const Sprog::Kinetics::DeltaStoich *const Reaction::DeltaStoich(unsigned int i) const
{
  if (i < m_deltaStoich.size()) {
        return m_deltaStoich[i];
    } else {
        return NULL;
    }
}
    // Returns pointer to delta stoich with given name.  NULL if not found.
Sprog::Kinetics::DeltaStoich *const Reaction::DeltaStoich(const std::string &name) const
{
  int i = FindDeltaStoich(name);
    if (i >= 0) {
        return m_deltaStoich[i];
    } else {
        return NULL;
    }
}

// Returns index of delta stoich.  Returns -1 if not found.
int Reaction::FindDeltaStoich(const std::string &name) const
{
    // Loop over delta stoich to find index.
    unsigned int i;
    for (i=0; i<m_deltaStoich.size(); i++) {
        if (*m_deltaStoich[i] == name) {
            // Found delta stoich!
            return i;
        }
    }
	
	// We are here because the phase wasn't found.
    return -1;
	
	
}

    // Adds an empty delta stoich to the reaction.
Sprog::Kinetics::DeltaStoich *const Reaction::AddDeltaStoich(const string &spName)
    {
    // Adds delta stoich to vector.
      Sprog::Kinetics::DeltaStoich delta_st;
    return AddDeltaStoich(delta_st, spName);
    }

    // Copies given delta stoich into the reaction.
Sprog::Kinetics::DeltaStoich *const Reaction::AddDeltaStoich(const Sprog::Kinetics::DeltaStoich &delta_st, const string &spName)
    {
      // First check  delta stoich vector to see if a  delta stoich with the
      // same name is already defined.  We can only have one delta stoich
      // per name.
    int i = FindDeltaStoich(delta_st.Name());
    if (i >= 0) {
        // A species with this name has been found.
        return m_deltaStoich.at(i);
    }

    // Adds delta stoich to vector.
    Sprog::Kinetics::DeltaStoich *delta_stoich_new = delta_st.Clone();
    delta_stoich_new->SetSpeciesName(spName);
    m_deltaStoich.push_back(delta_stoich_new);
    
    // Set up delta stoich.
    delta_stoich_new->SetReaction(*this); // You need a set reaction  function 
    

    // Return delta stoich reference.
    return delta_stoich_new;

    }
    
     

	// Returns pointer to delta stoich at index i.  NULL if not found.  This
    // function returns a modifiable (non-const) delta stoich object.
Sprog::Kinetics::DeltaStoich *const Reaction::GetDeltaStoich(const unsigned int i) const {
    if (i < m_deltaStoich.size()) {
        return m_deltaStoich[i];
    } else {
        return NULL;
    }
    }

    // Returns pointer to delta stoich with given name.  NULL if not found.  This
    // function returns a modifiable (non-const) delta stoich object.
Sprog::Kinetics::DeltaStoich *const Reaction::GetDeltaStoich(const std::string &name) const
    {
     int i = FindDeltaStoich(name);
    if (i >= 0) {
        return m_deltaStoich[i];
    } else {
        return NULL;
    }
    }


// ARRHENIUS COEFFICIENTS.

// Sets the forward Arrhenius parameters.
void Reaction::SetArrhenius(const ARRHENIUS &arr)
{
    m_arrf = arr;
}


// Sets the explicit reverse Arrhenius parameters.
void Reaction::SetRevArrhenius(const Sprog::Kinetics::ARRHENIUS &arr)
{
    if (m_arrr == NULL) {
        m_arrr = new ARRHENIUS(arr);
    } else {
        *m_arrr = arr;
    }
}



// LANDAU TELLER COEFFICIENTS.

// Sets the forward Landau Teller rate parameters.
void Reaction::SetLTCoeffs(const Sprog::Kinetics::LTCOEFFS &lt)
{
    if (m_lt == NULL) {
        m_lt = new LTCOEFFS(lt);
    } else {
        *m_lt = lt;
    }
}

// Sets the reverse Landau Teller rate parameters.
void Reaction::SetRevLTCoeffs(const Sprog::Kinetics::LTCOEFFS &lt)
{
    if (m_revlt == NULL) {
        m_revlt = new LTCOEFFS(lt);
    } else {
        *m_revlt = lt;
    }
}

// SURFACE
// Sets whether or not this reaction is surface types
void Reaction::SetIsSurface(const std::string &phName)
{
  if (phName.compare("gas") != 0)
    {
      m_isSurface = true; 
    }
 
}   

// COVERAGE

// Returns a constant reference to the vector of coverage
const vector<COVERAGE> &Reaction::CoverageReac() const
{
    return m_coverage;
}

// Returns the number of coverage params defined for this reaction.
int Reaction::COVERAGECount() const
{
  // return m_coverage.size(); 

    // Added new capability since coverage can be split from COV and STICK reaction then if this reaction is stick only count == 0

    if (m_isCoverage == false){

      return 0;
    }
    else{ 
      return m_coverage.size();
    }
}

// Returns the coefficient for the ith coverage.
COVERAGE Reaction::CoverageElement(unsigned int i) const
{
    if (i < m_coverage.size()) {
        return m_coverage[i];
    } else {
        // Invalid index.
      return COVERAGE(); // Default constructor
    }
}

// Sets the Coverage (always forward) parameters.
void Reaction::SetCoverage(const double e, const double m, const double eps, const std::string &name)
{
    m_isCoverage = true;

    COVERAGE m_cover; 
    
    m_cover.Eta = e;
    m_cover.Miu = m;
    m_cover.Epsilon = eps;
    m_cover.spName = name;

    m_coverage.push_back(m_cover);
}

// FORD
// Returns a constant reference to the vector of Ford
const vector<FORD> &Reaction::FordReac() const
{
    return m_ford;
}

// Returns the number of ford coefficients defined for this reaction.
int Reaction::FORDCount() const
{
    return m_ford.size();
}

// Returns the coefficient for the ith ford of the reaction.
FORD Reaction::FORDElement(unsigned int i) const
{
    if (i < m_ford.size()) {
        return m_ford[i];
    } else {
        // Invalid index.
      return FORD(); // Default constructor
    }
}

// Sets the Ford (always forward) parameters.
void Reaction::SetFord(const double c, const std::string &name)
{
    m_isFord = true;

    FORD m_forw; 
    
    m_forw.F_k = c;
    m_forw.spName = name;
    cout << "Set Ford: " << m_forw.spName << endl;
    m_ford.push_back(m_forw);
}


// THIRD-BODIES.

// Returns a constant reference to the vector of third bodies with their
// coefficients.
const vector<Stoich> &Reaction::ThirdBodies() const
{
    return m_thirdbodies;
}

// Returns the coefficient for the ith third body.
Stoich Reaction::ThirdBody(unsigned int i) const
{
    if (i < m_thirdbodies.size()) {
        return m_thirdbodies[i];
    } else {
        // Invalid index.
        return Stoich(0, 1.0);
    }
}

// Returns the number of third body coefficients defined for this reaction.
int Reaction::ThirdBodyCount() const
{
    return m_thirdbodies.size();
}

// Adds a third body to the reaction given the stoichiometric structure.
void Reaction::AddThirdBody(const Sprog::Stoich &tb)
{
    m_usetb = true;
    m_thirdbodies.push_back(tb);
}

// Adds a third body to the reaction given the species and the
// coefficient.
void Reaction::AddThirdBody(const unsigned int sp, double coeff)
{
    // Add a new Stoichf to the array of third bodies.
    m_usetb = true;
    m_thirdbodies.push_back(Stoich(sp, coeff));
}

// Adds a third body to the reaction given the species name.
void Reaction::AddThirdBody(const std::string &name, double coeff)
{
    if (m_mech != NULL) {
        // Find the species in the mechanism with the given name.
        unsigned int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species!
            m_usetb = true;
            m_thirdbodies.push_back(Stoich(i, coeff));
        }
    } else {
        throw logic_error("Can't add third body before assigning "
                          "parent mechanism (Sprog, Reaction::AddThirdBody).");
    }
}

// Removes a third body, given by name, from the reaction.  If the third
// body is not defined for this reaction then does nothing.
void Reaction::RemoveThirdBody(const std::string &name)
{
    if (m_mech != NULL) {
        // Search through the list of species to find that with
        // the given name.
        int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species:  Loop though third bodies and find
            // that with this index.
            vector<Stoich>::iterator j;
            for (j=m_thirdbodies.begin(); j!=m_thirdbodies.end(); j++) {
                if ((*j).Index() == i) {
                    // We have found the species in the list.
                    m_thirdbodies.erase(j);
                }
            }
        }
    } else {
        throw logic_error("Can't find third body name before assigning "
                          "parent mechanism (Sprog, Reaction::RemoveThirdBody).");
    }
}


// LOW PRESSURE LIMIT.

// Returns the Arrhenius parameters for the low pressure limit.
const ARRHENIUS &Reaction::LowPressureLimit() const
{
    return m_foparams.LowP_Limit;
}

// Sets the low pressure limit Arrhenius parameters.
void Reaction::SetLowPressureLimit(const Sprog::Kinetics::ARRHENIUS &lowp)
{
    if (m_fotype == None) m_fotype = Lindemann;
    m_foparams.LowP_Limit = lowp;
}


// FALL-OFF THIRD BODY.

// Returns a pointer to the species used as a third body in the fall-off reaction.
const Species *const Reaction::FallOffThirdBody() const
{
    if ((m_mech != NULL) && (m_foparams.ThirdBody >= 0)) {
        return m_mech->Species(m_foparams.ThirdBody);
    }
    return NULL;
}

// Sets the species to use as a third body for fall-off calculations.
void Reaction::SetFallOffThirdBody(int sp)
{
    m_foparams.ThirdBody = sp;
}

// Sets the species to use as a third body for fall-off calculations
// given the species name.
void Reaction::SetFallOffThirdBody(const std::string &name)
{
    if (m_mech != NULL) {
        // Locate the species by name in the vector.
        unsigned int i = m_mech->FindSpecies(name);
        if (i >= 0) {
            // Found the species!
            m_foparams.ThirdBody = i;
        } else {
            // If the name is not a valid species then
            // use all species as third body.
            m_foparams.ThirdBody = -1;
        }
    } else {
        throw logic_error("Can't find set fall-off third body before assigning "
                          "parent mechanism (Sprog, Reaction::SetFallOffThirdBody).");
    }
}


// FALL-OFF PARAMETERS.

// Returns the fall-off type.
Kinetics::FALLOFF_FORM Reaction::FallOffType() const
{
    return m_fotype;
}

// Returns the fall-off parameters.
const Kinetics::FALLOFF_PARAMS &Reaction::FallOffParams() const
{
    return m_foparams;
}

// Sets the fall-off type and parameters.
void Reaction::SetFallOffParams(const FALLOFF_FORM form,
                                const double params[FALLOFF_PARAMS::MAX_FALLOFF_PARAMS])
{
    m_fotype = form;

    for (int i=0; i<FALLOFF_PARAMS::MAX_FALLOFF_PARAMS; i++) {
        m_foparams.Params[i] = params[i];
    }
}


// FALL-OFF FUNCTIONAL FORMS.

// 3-parameter Troe fall-off form.
double Reaction::FTROE3(double T, double logpr) const
{
    double fcent, c, n, F;
    const double d = 0.14;

    fcent = log10(((1.0 - m_foparams.Params[0]) * exp(-T / m_foparams.Params[1])) +
                  (m_foparams.Params[0] * exp(-T / m_foparams.Params[2])));
    c = logpr - 0.4 - (0.67 * fcent);
    n = 0.75 - (1.27 * fcent) - (d * c);
    c = c / n;
    F = pow(10.0, fcent / (1 + (c*c)));
    return F;
}

// 4-parameter Troe fall-off form.
double Reaction::FTROE4(double T, double logpr) const
{
    double fcent, c, n, F;
    const double d = 0.14;

    fcent = log10(((1.0 - m_foparams.Params[0]) * exp(-T / m_foparams.Params[1])) +
                  (m_foparams.Params[0] * exp(-T / m_foparams.Params[2])) +
                  exp(-m_foparams.Params[3] / T));
    c = logpr - 0.4 - (0.67 * fcent);
    n = 0.75 - (1.27 * fcent) - (d * c);
    c = c / n;
    F = pow(10.0, fcent / (1 + (c*c)));
    return F;
}

// SRI fall-off form.
double Reaction::FSRI(double T, double logpr) const
{
    double x = 1.0 / (1.0 + (logpr * logpr));
    double F = m_foparams.Params[3] * pow(T, m_foparams.Params[4]) *
             pow((m_foparams.Params[0]*exp(-m_foparams.Params[1]/T)) +
                 exp(-T/m_foparams.Params[2]), x);
    return F;
}

// Custom functional form for fall-off.
/*FallOffFnPtr Reaction::FallOffFn() const
{
    return m_fofn;
}*/

// PARENT MECHANISM.

// Returns a pointer to the parent mechanism.
const Mechanism *const Reaction::Mechanism() const
{
    return m_mech;
}

// Sets the parent mechanism.
void Reaction::SetMechanism(Sprog::Mechanism &mech)
{
    m_mech = &mech;
}


// RATE CALCULATION.

// Calculates the rate of progress of this reaction.
double Reaction::RateOfProgress(double density, const double *const x,
                              unsigned int n, double kforward,
                              double kreverse) const
{
    int j=0;
    unsigned int k=0;
    double rop=0.0, rev=0.0;

    if (n >= m_mech->Species().size()) {
        // Use rop to store forward rates of production,
        // and rev to store reverse rates.
        rop = kforward;
        rev = kreverse;

        // Integer reactants.
        for (k=0; k!=m_reac.size(); ++k) {
            // As the stoichiometry is integer, it is more computationally efficient
            // to multiply the values together than to use the pow() function.
            for (j=0; j!=m_reac[k].Mu(); ++j) {
                rop *= density * x[m_reac[k].Index()];
            }
        }

        // Integer products.
        for (k=0; k!=m_prod.size(); ++k) {
            // As the stoichiometry is integer, it is more computationally efficient
            // to multiply the values together than to use the pow() function.
            for (j=0; j!=m_prod[k].Mu(); ++j) {
                rev *= density * x[m_prod[k].Index()];
            }
        }

        // Calculate the net rate of production.
        rop -= rev;
    }

    return rop;
}


// READ/WRITE/COPY FUNCTIONS.

// Creates a copy of the reaction object.
Reaction *Reaction::Clone(void) const
{
    return new Reaction(*this);
}

// Prints a diagnostic output file containing all the
// reaction data.  This is used to debug.
void Reaction::WriteDiagnostics(std::ostream &out) const
{
 
    string data = "";
    double val = 0.0;
    int ival = 0;

    if (out.good()) {
        // Name.
        out.write(string(m_name+" ").c_str(), m_name.length());

        // Reaction reversibility.
        if (m_reversible) {
            data = "R ";
        } else {
            data = "- ";
        }
        out.write(data.c_str(), data.length());

        //  reactant stoichiometry.
        for (unsigned int i=0; i!=m_reac.size(); ++i) {
            // Species name.
            data = m_mech->Species(m_reac[i].Index())->Name() + " ";
            out.write(data.c_str(), data.length());
            // Stoichiometry.
            val = m_reac[i].Mu();
            data = cstr(val) + " ";
            out.write(data.c_str(), data.length());
        }

        // Integer product stoichiometry.
        for (unsigned int i=0; i!=m_prod.size(); ++i) {
            // Species name.
            data = m_mech->Species(m_prod[i].Index())->Name() + " ";
            out.write(data.c_str(), data.length());
            // Stoichiometry.
            val = m_prod[i].Mu();
            data = cstr(val) + " ";
            out.write(data.c_str(), data.length());
        }

	

        // Stoichiometry changes.
        data = cstr(m_dstoich) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_dreac) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_dprod) + " ";
        out.write(data.c_str(), data.length());

        // Forward Arrhenius coefficients.
        data = cstr(m_arrf.A) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_arrf.n) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_arrf.E) + " ";
        out.write(data.c_str(), data.length());

        // Reverse Arrhenius coefficients.
        if (m_arrr != NULL) {
            data = cstr(m_arrr->A) + " ";
            out.write(data.c_str(), data.length());
            data = cstr(m_arrr->n) + " ";
            out.write(data.c_str(), data.length());
            data = cstr(m_arrr->E) + " ";
            out.write(data.c_str(), data.length());
        }

        // Forward LT parameters.
        if (m_lt != NULL) {
            data = cstr(m_lt->B) + " ";
            out.write(data.c_str(), data.length());
            data = cstr(m_lt->C) + " ";
            out.write(data.c_str(), data.length());
        }

        // Reverse LT parameters.
        if (m_revlt != NULL) {
            data = cstr(m_revlt->B) + " ";
            out.write(data.c_str(), data.length());
            data = cstr(m_revlt->C) + " ";
            out.write(data.c_str(), data.length());
        }

        // Third body flag.
        if (m_usetb) {
            data = "TB ";
        } else {
            data = "-- ";
        }
        out.write(data.c_str(), data.length());

        // Write third bodies.
        for (unsigned int i=0; i!=m_thirdbodies.size(); ++i) {
            // Write species name.
            data = m_mech->Species(m_thirdbodies[i].Index())->Name() + " ";
            out.write(data.c_str(), data.length());
            // Write mu.
            data = cstr(m_thirdbodies[i].Mu()) + " ";
            out.write(data.c_str(), data.length());
        }


	// Surface flag.
        if (m_isSurface) {
            data = "SURF ";
        } else {
            data = "-- ";
        }
        out.write(data.c_str(), data.length());

	// STICK flag.
        if (m_sticking) {
            data = "ST ";
        } else {
            data = "-- ";
        }
        out.write(data.c_str(), data.length());

	

	// Mott-Wise flag.
        if (m_mottwise) {
            data = "MW ";
        } else {
            data = "-- ";
        }
        out.write(data.c_str(), data.length());


	// Ford flag.
        if (m_isFord) {
            data = "FD ";
        } else {
            data = "-- ";
        }
        out.write(data.c_str(), data.length());


	// Write ford.
        for (unsigned int i=0; i!=m_ford.size(); ++i) {
	data = cstr(m_ford[i].F_k) + " ";
        out.write(data.c_str(), data.length());
	out.write(string(m_ford[i].spName+" ").c_str(), m_ford[i].spName.length());
	}
	

	// Coverage flag.
        if (m_isCoverage) {
            data = "CV ";
        } else {
            data = "-- ";
        }
        out.write(data.c_str(), data.length());


	// Write coverage.
        for (unsigned int i=0; i!=m_coverage.size(); ++i) {
	data = cstr(m_coverage[i].Eta) + " ";
        out.write(data.c_str(), data.length());
	data = cstr(m_coverage[i].Miu) + " ";
        out.write(data.c_str(), data.length());
	data = cstr(m_coverage[i].Epsilon) + " ";
        out.write(data.c_str(), data.length());
	out.write(string(m_coverage[i].spName+" ").c_str(), m_coverage[i].spName.length());
	}

        // Write fall-off type.
        ival = (int)m_fotype;
        data = cstr(ival) + " ";
        out.write(data.c_str(), data.length());

        // Write fall-off low pressure limit.
        data = cstr(m_foparams.LowP_Limit.A) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_foparams.LowP_Limit.n) + " ";
        out.write(data.c_str(), data.length());
        data = cstr(m_foparams.LowP_Limit.E) + " ";
        out.write(data.c_str(), data.length());

        // Write fall-off third body.
        if (m_foparams.ThirdBody >= 0) {
            data = m_mech->Species(m_foparams.ThirdBody)->Name() + " ";
            out.write(data.c_str(), data.length());
        }

        // Write fall-off parameters.
        for (unsigned int i=0; i!=(unsigned)FALLOFF_PARAMS::MAX_FALLOFF_PARAMS; ++i) {
            data = cstr(m_foparams.Params[i]) + " ";
            out.write(data.c_str(), data.length());
        }

				
	// Phase Involved.
	data = "Phase involved:\n";
	out.write(data.c_str(), data.length());
        for (unsigned int i=0; i!=m_phaseVector.size(); ++i) {
            // Phase name in this reaction.
	  out.write(string(m_phaseVector[i]+"\t").c_str(), m_phaseVector[i].length());
          
        }

	// Write the DeltaStoich to the file.
	data = " DeltaStoich:\n";
	out.write(data.c_str(), data.length());
	for (unsigned int i=0; i!=m_deltaStoich.size(); ++i) {
        m_deltaStoich[i]->WriteDiagnostics(out);
	}

        // New line.
        data = "\n";
        out.write(data.c_str(), data.length());
    }
}




/*!
 * @brief       Should the pre-exponential factor be converted to cgs?
 *
 * Only three-body reactions include the [M] term explicitly in their rate
 * expression, that is, falloff reactions DO NOT require conversion of the
 * high-pressure pre-exponential's units $A_\infty$.
 *
 * @return
 */
bool Reaction::ConvertPreexponential(void) const {
    bool ans(false);
    if ((m_usetb) && (m_fotype == Sprog::Kinetics::None)) {
        ans = true;
    }
    return ans;
}


/*!
@param[in]      out             Output stream for file creation.
@param[in]      RejectSpecies   String vector containing the names of LOI rejected species.
*/
void Reaction::WriteReducedMechReacs(std::ostream &out, std::vector<std::string> RejectSpecies) const
{

    if (out.good()) {

        //Check to see if the reaction contains any of the rejected species.
        //If so, return without writing anything.
        for (unsigned int n = 0; n < RejectSpecies.size(); n++){
            for (unsigned int j = 0; j < m_reac.size(); j++){
                if (RejectSpecies[n] == m_mech->Species(m_reac[j].Index())->Name())
                    return;
            }
            for (unsigned int k = 0; k < m_prod.size(); k++){
                if (RejectSpecies[n] == m_mech->Species(m_prod[k].Index())->Name())
                    return;
            }
        }

		std::string tbstr;
        if (m_fotype != Sprog::Kinetics::None){
            tbstr = "(+M)";
        } else if (m_usetb) {
            tbstr = "+M";
        };

        // Integer reactant stoichiometry.
        if (m_reac.size() > 0 && m_prod.size() > 0){
            for (unsigned int i = 0; i != m_reac.size(); ++i) {
                // Species name.
                for (int j = 0; j < m_reac[i].Mu(); j++) {
                    out << m_mech->Species(m_reac[i].Index())->Name();
                    if (i < m_reac.size() - 1 || j < m_reac[i].Mu() - 1)
                        out << " + ";
                }
            }


            if (m_usetb || m_fotype != Sprog::Kinetics::None){
                out << tbstr;
            }

                // Reaction reversibility.
                if (m_reversible) {
                    out << " =  ";
                } else {
                    out << " => ";
                }

            // Integer product stoichiometry.
            for (unsigned int i = 0; i != m_prod.size(); ++i) {
                // Species name.
                for (int j = 0; j < m_prod[i].Mu(); j++) {
                    out << m_mech->Species(m_prod[i].Index())->Name();
                    if (i < m_prod.size() - 1 || j < m_prod[i].Mu() - 1)
                        out << " + ";
                }
            }
            if (m_usetb || m_fotype != Sprog::Kinetics::None){
                out << tbstr;
            }

            // Forward Arrhenius coefficients.
            out << " " << m_arrf.A / pow(1.0e-6, ReactantStoich() - 1.0
                    + (ConvertPreexponential()?1.0:0.0)) << " " << m_arrf.n  << " " << m_arrf.E  / 4.184E7 / 1.0e-7 << "\n";

            // Reverse Arrhenius coefficients.
            if (m_arrr != NULL)
                 out << "Rev / " << m_arrr->A / pow(1.0e-6, ProductStoich() - 1.0
                    + (ConvertPreexponential()?1.0:0.0)) << " " << m_arrr->n << " " << m_arrr->E / 4.184E7 / 1.0e-7 << " /\n";
        }

        // Forward LT parameters.
        if (m_lt != NULL) {
            out << m_lt->B << " " << m_lt->C << " ";
        }

        // Reverse LT parameters.
        if (m_revlt != NULL) {
            out << m_revlt->B << " " << m_revlt->C << " ";
        }

        // Write fall-off low pressure limit.
        if (m_fotype != Sprog::Kinetics::None){
            out << "LOW /" << cstr((m_foparams.LowP_Limit.A)/(pow(1.0e-6, ReactantStoich()))) << " " << m_foparams.LowP_Limit.n << " " << (m_foparams.LowP_Limit.E/ 4.184E7 / 1.0e-7) << " ";
            out << " / " << "\n";

             //Write fall-off third body.
            if (m_foparams.ThirdBody >= 0) {
                out << m_mech->Species(m_foparams.ThirdBody)->Name() << " ";
            }

            // Write fall-off parameters.
            if (m_fotype == Sprog::Kinetics::Troe3 || m_fotype == Sprog::Kinetics::Troe4) {
                out << "TROE /";
                for (unsigned int i = 0; i != (unsigned) (FALLOFF_PARAMS::MAX_FALLOFF_PARAMS - 1); ++i) {
                    out << m_foparams.Params[i] << " ";
                }
                out << "/ " << "\n";
            }

        }

        // Write third-body efficiencies.
        for (int i=0; i<ThirdBodyCount(); i++) {
            out << m_mech->Species(ThirdBody(i).Index())->Name();
            out << "/" << ThirdBody(i).Mu() << "/";
            if (i>=ThirdBodyCount()-1) {out << "\n";} else {out << " ";};
        };

        // New line.
        out << "\n";

    }
    out.flush();
}


// MEMORY MANAGEMENT.

void Reaction::releaseMemory(void)
{
    m_name.clear();
    m_reac.clear();
    m_prod.clear();
    m_dstoich = m_dreac = m_dprod = 0.0;
    m_arrf.A = m_arrf.n = m_arrf.E = 0.0;
    if (m_arrr != NULL) delete m_arrr;
    m_arrr = NULL;
    if (m_lt != NULL) delete m_lt;
    m_lt = NULL;
    if (m_revlt != NULL) delete m_revlt;
    m_revlt = NULL;
    m_fo.F_k     = 0.0;
    m_fo.spName  = "";
    m_covr.Eta   = m_covr.Miu = m_covr.Epsilon =0.0;  
    m_covr.spName = "";
    m_usetb = false;
    m_thirdbodies.clear();
    m_fotype = None;
    m_foparams.LowP_Limit.A = m_foparams.LowP_Limit.n = m_foparams.LowP_Limit.E = 0.0;
    m_foparams.ThirdBody = -1;
    //m_fofn = NULL;
    m_isSurface = false; 
    // Stick
    m_sticking = false;
    // Mottwise
    m_mottwise = false; 
    // Coverage data
    m_isCoverage = false;
    m_coverage.clear();
    // Ford data 
    m_isFord = false; 
    m_ford.clear();

    m_mech = NULL;
    m_deltaStoich.clear();
    m_phaseVector.clear();
}

// Writes the reaction to a binary data stream.
void Reaction::Serialize(std::ostream &out) const
{
    if (out.good()) {
        const unsigned int trueval  = 1;
        const unsigned int falseval = 0;
	unsigned int u;
        // Write the serialisation version number to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write the length of the reaction name to the stream.
        unsigned int n = m_name.length();
        out.write((char*)&n, sizeof(n));

        // Write the reaction name to the stream.
        if (n > 0) {
            out.write(m_name.c_str(), n);
        }

        // Write reaction reversibility.
        if (m_reversible) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write integer reactant stoichiometry count.
        n = m_reac.size();
        out.write((char*)&n, sizeof(n));

        // Write integer reactant stoichiometry.
        for (unsigned int i=0; i<n; i++) {
            // Write species index.
            int ix = m_reac[i].Index();
            out.write((char*)&ix, sizeof(ix));

            // Write mu.
            int mu = m_reac[i].Mu();
            out.write((char*)&mu, sizeof(mu));
        }

        // Write integer product stoichiometry count.
        n = m_prod.size();
        out.write((char*)&n, sizeof(n));

        // Write integer product stoichiometry.
        for (unsigned int i=0; i<n; i++) {
            // Write species index.
            int ix = m_prod[i].Index();
            out.write((char*)&ix, sizeof(ix));

            // Write mu.
            int mu = m_prod[i].Mu();
            out.write((char*)&mu, sizeof(mu));
        }

	 // Write the PhaseName involved.
        u = m_phaseVector.size();
        out.write((char*)&u, sizeof(u));
	
	for (unsigned int i=0; i<u; i++) {
            // Write phase name.
	  n = m_phaseVector[i].length();
	  out.write((char*)&n, sizeof(n));
	  out.write(m_phaseVector[i].c_str(), m_phaseVector[i].length());

        }

	// Write the number of delta stoich to the stream.
        u = m_deltaStoich.size();
        out.write((char*)&u, sizeof(u));

        // Write the delta stoich to the stream.
        for (DeltaStoichPtrVector::const_iterator ids=m_deltaStoich.begin(); ids!=m_deltaStoich.end(); ids++) {
            (*ids)->Serialize(out);
        }

        // Write stoichiometry changes.
        double mu = (double)m_dstoich;
        out.write((char*)&mu, sizeof(mu));
        mu = (double)m_dreac;
        out.write((char*)&mu, sizeof(mu));
        mu = (double)m_dprod;
        out.write((char*)&mu, sizeof(mu));

        // Write forward Arrhenius coefficients.
        double A  = (double)m_arrf.A;
        double nn = (double)m_arrf.n;
        double E  = (double)m_arrf.E;
        out.write((char*)&A, sizeof(A));
        out.write((char*)&nn, sizeof(nn));
        out.write((char*)&E, sizeof(E));

        // Write reverse Arrhenius coefficients.
        if (m_arrr != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            A  = (double)m_arrr->A;
            nn = (double)m_arrr->n;
            E  = (double)m_arrr->E;
            out.write((char*)&A, sizeof(A));
            out.write((char*)&nn, sizeof(nn));
            out.write((char*)&E, sizeof(E));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write forward LT parameters.
        if (m_lt != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            A = (double)m_lt->B;
            E = (double)m_lt->C;
            out.write((char*)&A, sizeof(A));
            out.write((char*)&E, sizeof(E));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write reverse LT parameters.
        if (m_revlt != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            A = (double)m_revlt->B;
            E = (double)m_revlt->C;
            out.write((char*)&A, sizeof(A));
            out.write((char*)&E, sizeof(E));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write third body flag.
        if (m_usetb) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write third body count.
        n = m_thirdbodies.size();
        out.write((char*)&n, sizeof(n));

        // Write third bodies.
        for (unsigned int i=0; i<n; i++) {
            // Write species index.
            int ix = m_thirdbodies[i].Index();
            out.write((char*)&ix, sizeof(ix));

            // Write mu.
            double mu = m_thirdbodies[i].Mu();
            out.write((char*)&mu, sizeof(mu));
        }

	// Write SURFACE flag.
        if (m_isSurface) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }


	// Write STICK flag.
        if (m_sticking) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

	 // Write Mott-Wise flag.
        if (m_mottwise) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

	 // Write ford flag.
        if (m_isFord) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

	// Write ford count.
        n = m_ford.size();
        out.write((char*)&n, sizeof(n));


	// Write ford.
        for (unsigned int i=0; i<n; i++) {
	double coef  = (double)m_ford[i].F_k;
        out.write((char*)&coef, sizeof(coef));
	// Write the length of the species name to the stream.
        unsigned int sp = m_ford[i].spName.length();
        out.write((char*)&sp, sizeof(sp));
	cout << m_ford[i].spName << endl; // for debugging
        // Write the species name to the stream.
        if (sp > 0) {
            out.write(m_ford[i].spName.c_str(), sp);
        }
	
        }


	 // Write coverage flag.
        if (m_isCoverage) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

	// Write coverage count.
        n = m_coverage.size();
        out.write((char*)&n, sizeof(n)); 

	// Write coverage.
        for (unsigned int i=0; i<n; i++) {
	double e  = (double)m_coverage[i].Eta;
        double m = (double)m_coverage[i].Miu;
        double epsil  = (double)m_coverage[i].Epsilon;
        out.write((char*)&e, sizeof(e));
        out.write((char*)&m, sizeof(m));
        out.write((char*)&epsil, sizeof(epsil));
	// Write the length of the species name to the stream.
        unsigned int sp = m_coverage[i].spName.length();
        out.write((char*)&sp, sizeof(sp));

        // Write the speciesn name to the stream.
        if (sp > 0) {
	  out.write(m_coverage[i].spName.c_str(), sp);
        }
	
        }

        // Write fall-off type.
        n = (unsigned int)m_fotype;
        out.write((char*)&n, sizeof(n));

        // Write fall-off low pressure limit.
        A  = (double)m_foparams.LowP_Limit.A;
        nn = (double)m_foparams.LowP_Limit.n;
        E  = (double)m_foparams.LowP_Limit.E;
        out.write((char*)&A, sizeof(A));
        out.write((char*)&nn, sizeof(nn));
        out.write((char*)&E, sizeof(E));

        // Write fall-off third body.
        int ix = m_foparams.ThirdBody;
        out.write((char*)&ix, sizeof(ix));

        // Write fall-off parameter count.
        n = FALLOFF_PARAMS::MAX_FALLOFF_PARAMS;
        out.write((char*)&n, sizeof(n));

        // Write fall-off parameters.
        for (unsigned int i=0; i<n; i++) {
            A = (double)m_foparams.Params[i];
            out.write((char*)&A, sizeof(A));
        }

        // There is no way of outputting custom fall-off functions at the moment.
        // This needs to corrected in the future.

    } else {
        throw invalid_argument("Output stream not ready (Sprog, Reaction::Serialize).");
    }
}

// Reads the reaction data from a binary data stream.
void Reaction::Deserialize(std::istream &in)
{
    // Clear the reaction of all current data.
    releaseMemory();

    if (in.good()) {
        // Read the serialized species version number.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0, u =0, spN = 0; // Need for reading name length.
        char *name = NULL;
        double A = 0.0, nn = 0.0, E = 0.0;
        int itb = 0;

        switch (version) {
            case 0:
                // Read the length of the reaction name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the reaction name.
                if (n > 0) {
                    name = new char[n];
                    in.read(name, n);
                    m_name.assign(name, n);
                    delete [] name;
                } else {
                    m_name = "";
                }

                // Read reaction reversibility.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n==1) {
                    m_reversible = true;
                } else {
                    m_reversible = false;
                }

                // Read integer reactant stoichiometry count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_reac.reserve(n);

                // Read integer reactant stoichiometry.
                for (unsigned int i=0; i<n; i++) {
                    int ix = 0;
                    int mu = 0;

                    // Read species index.
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Read mu.
                    in.read(reinterpret_cast<char*>(&mu), sizeof(mu));

                    // Push a new Stoich object into the vector.
                    m_reac.push_back(Stoich(ix, mu));
                }

                // Read integer product stoichiometry count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_prod.reserve(n);

                // Read integer product stoichiometry.
                for (unsigned int i=0; i<n; i++) {
                    int ix = 0;
                    int mu = 0;

                    // Read species index.
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Read mu.
                    in.read(reinterpret_cast<char*>(&mu), sizeof(mu));

                    // Push a new Stoich object into the vector.
                    m_prod.push_back(Stoich(ix, mu));
                }

		// Read PhaseName involved. 
		in.read(reinterpret_cast<char*>(&n), sizeof(n));
		m_phaseVector.reserve(n);

		  for (unsigned int i=0; i<n; i++) {
                    
		    std::string NAME;
		     // Read the length of the phase name.
		    in.read(reinterpret_cast<char*>(&spN), sizeof(spN));

                // Read the reaction name.
		    if (spN > 0) {
		      name = new char[spN];
		      in.read(name, spN);
		      NAME.assign(name, spN);
		      delete [] name;
		    } else {
		      NAME = "";
		    }
                  

                    // Push a new string into the vector.
                    m_phaseVector.push_back(NAME);
		  }

		  // Read the number of delta stoich and reserve memory.
                in.read(reinterpret_cast<char*>(&u), sizeof(u));
                m_deltaStoich.reserve(u);

                // Read the delta stoich.
                try {
                    for (unsigned int i=0; i<u; i++) {
                        // Read the delta stoich from the stream using the
                        // appropriate constructor.
		      Sprog::Kinetics::DeltaStoich *ds = new Sprog::Kinetics::DeltaStoich(in);
                        ds->SetReaction(*this);

                        // Add the delta stoich to the vector.
                        m_deltaStoich.push_back(ds);
                    }
                } catch (exception &e) {
                    // Ensure the mechanism is cleared before throwing
                    // the exception to the next layer up.
                    releaseMemory();
                    throw;
                }
			
		
                // Read stoichiometry changes.
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                m_dstoich = (double)A;
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                m_dreac = (double)A;
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                m_dprod = (double)A;

                // Read forward Arrhenius coefficients.
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
                in.read(reinterpret_cast<char*>(&E), sizeof(E));
                m_arrf.A = (double)A;
                m_arrf.n = (double)nn;
                m_arrf.E = (double)E;

                // Read reverse Arrhenius coefficients.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    in.read(reinterpret_cast<char*>(&A), sizeof(A));
                    in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
                    in.read(reinterpret_cast<char*>(&E), sizeof(E));
                    m_arrr = new Kinetics::ARRHENIUS();
                    m_arrr->A = (double)A;
                    m_arrr->n = (double)nn;
                    m_arrr->E = (double)E;
                }

                // Read forward LT parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    in.read(reinterpret_cast<char*>(&A), sizeof(A));
                    in.read(reinterpret_cast<char*>(&E), sizeof(E));
                    m_lt = new LTCOEFFS();
                    m_lt->B = (double)A;
                    m_lt->C = (double)E;
                }

                // Read reverse LT parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    in.read(reinterpret_cast<char*>(&A), sizeof(A));
                    in.read(reinterpret_cast<char*>(&E), sizeof(E));
                    m_revlt = new LTCOEFFS();
                    m_revlt->B = (double)A;
                    m_revlt->C = (double)E;
                }

                // Read third body flag.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_usetb = true;
                } else {
                    m_usetb = false;
                }

                // Read third body count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_thirdbodies.reserve(n);

                // Read third bodies.
                for (unsigned int i=0; i<n; i++) {
                    int ix = 0;
                    double mu = 0.0;

                    // Read species index.
                   in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Read mu.
                   in.read(reinterpret_cast<char*>(&mu), sizeof(mu));

                   // Add third body to vector.
                   m_thirdbodies.push_back(Stoich(ix,mu));
                }

		// Read SURFACE flag (include all surface reactions) 
		in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_isSurface = true;
                } else {
                    m_isSurface = false;
                }

		// Read STICK flag 
		in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_sticking = true;
                } else {
                    m_sticking = false;
                }

		// Read Mott-Wise flag
		in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_mottwise = true;
                } else {
                    m_mottwise = false;
                }
		
		// Read FORD flag
		in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_isFord = true;
                } else {
                    m_isFord = false;
                }

		// Read FORD count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_ford.reserve(n);


		 // Read FORD.
                for (unsigned int i=0; i<n; i++) {
                    
		  double coef = 0.0;
		  std::string speciesName; 
                    in.read(reinterpret_cast<char*>(&coef), sizeof(coef));
		    in.read(reinterpret_cast<char*>(&spN), sizeof(spN));

		    // Read the species name.
                    name = new char[spN];
                    in.read(name, spN);
                    speciesName.assign(name, spN);
                    delete [] name;
		    
		    
                    m_fo.F_k = (double)coef;
		    cout << (double)coef << endl; 
		    m_fo.spName = speciesName;
		    
		    // Add covr to vector.
		    m_ford.push_back(m_fo);
                }


		// Read COVERAGE flag 
		in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_isCoverage = true;
                } else {
                    m_isCoverage = false;
                }

		// Read COVERAGE count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_coverage.reserve(n);

		 // Read COVERAGE.
                for (unsigned int i=0; i<n; i++) {
                    
		  double e, m, epsil = 0.0;
		  std::string speciesName; 
                    in.read(reinterpret_cast<char*>(&e), sizeof(e));
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));
                    in.read(reinterpret_cast<char*>(&epsil), sizeof(epsil));
		    
		    in.read(reinterpret_cast<char*>(&spN), sizeof(spN));

		    // Read the species name.
                    name = new char[spN];
                    in.read(name, spN);
                    speciesName.assign(name, spN);
                    delete [] name;
		   

                    m_covr.Eta = (double)e;
                    m_covr.Miu = (double)m;
                    m_covr.Epsilon = (double)epsil;
		    m_covr.spName = speciesName;
		    
		    // Add covr to vector.
		    m_coverage.push_back(m_covr);
                }


                // Read fall-off type.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_fotype = (FALLOFF_FORM)n;

                // Read fall-off low pressure limit.
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
                in.read(reinterpret_cast<char*>(&E), sizeof(E));
                m_foparams.LowP_Limit.A = (double)A;
                m_foparams.LowP_Limit.n = (double)nn;
                m_foparams.LowP_Limit.E = (double)E;

                // Read fall-off third body.
                in.read(reinterpret_cast<char*>(&itb), sizeof(itb));
                m_foparams.ThirdBody = itb;

                // Read fall-off parameter count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read fall-off parameters.
                for (unsigned int i=0; i<n; i++) {
                    in.read(reinterpret_cast<char*>(&A), sizeof(A));
                    m_foparams.Params[i] = (double)A;
                }

                break;
            default:
                throw runtime_error("Reaction serialized version number is "
                                    "unsupported (Sprog, Reaction::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, Reaction::Deserialize).");
    }
}
