/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Mechanism class declared in the
    gpc_mech.h header file.

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

#include "gpc_mech.h"
#include "gpc_unit_systems.h"
#include <string>
#include <math.h>
#include <stdexcept>
#include <fstream>
#include "string_functions.h"

using namespace Sprog;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Mechanism::Mechanism()
{
    m_units = SI;
    m_stoich_xref_valid = false;
    m_rxns.SetMechanism(*this);
}

// Copy constructor.
Mechanism::Mechanism(const Sprog::Mechanism &mech)
{
    *this = mech;
}

// Destructor.
Mechanism::~Mechanism()
{
    releaseMemory();
}


// OPERATOR OVERLOADING.

// Assignment operator
Mechanism &Mechanism::operator=(const Sprog::Mechanism &mech)
{
    // Check for self assignment!
    if (this != &mech) {
        // Clear current memory.
        releaseMemory();

        // Copy mechanism data.
        m_units = mech.m_units;

        // Copy new elements and species into mechanism.
        copyInPhase(mech.m_phase); // Added by mm864
        copyInElements(mech.m_elements);
        copyInSpecies(mech.m_species);

        // Copy reaction set and stoich cross-referencing.
        m_rxns = mech.m_rxns;
        m_stoich_xref.assign(mech.m_stoich_xref.begin(), mech.m_stoich_xref.end());
        m_stoich_xref_valid = mech.m_stoich_xref_valid;

        // Inform species of new elements vector and mechanism.
        SpeciesPtrVector::iterator sp;
        for (sp=m_species.begin(); sp!=m_species.end(); sp++) {
            (*sp)->SetMechanism(*this);
        }

       // Inform phase of new species vector and mechanism. 
	PhasePtrVector::iterator ph;
        for (ph=m_phase.begin(); ph!=m_phase.end(); ph++) {
            (*ph)->SetMechanism(*this);
        }
		
        // Inform reactions of new species vector and mechanism.
        unsigned int i;
        for (i=0; i<m_rxns.Count(); i++) {
            m_rxns[i]->SetMechanism(*this);
        }
    }
    return *this;
}


// Clears the mechanism of all elements, species and reactions.
void Mechanism::Clear()
{
    releaseMemory();
}


// UNITS.

// Returns the current system of units of this mechanism.
UnitSystem Mechanism::Units() const
{
    return m_units;
}

// Sets the system of units used by this mechanism by converting all
// element, species and reaction properties.
void Mechanism::SetUnits(Sprog::UnitSystem u)
{
    Kinetics::ARRHENIUS arr;

    // Is new system SI?
    if (u == SI) {
        // Is current system CGS?
        if (m_units == CGS) {
            // Convert elements' mol. weights.
            ElementPtrVector::iterator iel;
            for (iel=m_elements.begin(); iel!=m_elements.end(); iel++) {
                // Convert from g/mol to kg/mol.
                (*iel)->SetMolWt((*iel)->MolWt()*1.0e-3);
            }

            // Scale species' mol. weights.
            SpeciesPtrVector::iterator isp;
            for (isp=m_species.begin(); isp!=m_species.end(); isp++) {
                // Recalculate mol. weight.
                (*isp)->CalcMolWt();
            }

	     
	     // Scale phase site density to mol/m2. (Added by mm864)
            PhasePtrVector::iterator iph;
            for (iph=m_phase.begin(); iph!=m_phase.end(); iph++) {
		                 
		// Recalculate site density.
                (*iph)->SetSiteDensity((*iph)->SiteDen()*1.0e4, (*iph)->Name());
            }	
	     			
            // Scale reaction coefficients.
            for (unsigned int irxn=0; irxn!=m_rxns.Count(); ++irxn) {
                // Convert volumes from cm3 to m3, and convert energies from ergs/mol to
                // J/mol.
				double gasReactantStoich = 0.0, surfReactantStoich = 0.0, gasProductStoich = 0.0, surfProductStoich = 0.0;
				
				for (unsigned int dS = 0; dS!= m_rxns[irxn]->DeltaStoichCount(); ++dS){
					
					string spName = m_rxns[irxn]->DeltaStoich(dS)->Name();
					string phName = m_rxns[irxn]->Mechanism()->GetSpecies(spName)->PhaseName(); 
					
					if ((phName.compare("gas") == 0) && (m_rxns[irxn]->DeltaStoich(dS)->ReacStoich() > 0)){
					gasReactantStoich += m_rxns[irxn]->DeltaStoich(dS)->ReacStoich();
					}
					else if ((phName.compare("gas") == 0) && (m_rxns[irxn]->DeltaStoich(dS)->ProdStoich() > 0)){
					gasProductStoich += m_rxns[irxn]->DeltaStoich(dS)->ProdStoich();
					}
					else if ((phName.compare("gas") != 0) && (m_rxns[irxn]->DeltaStoich(dS)->ReacStoich() > 0)){
					surfReactantStoich += m_rxns[irxn]->DeltaStoich(dS)->ReacStoich();
					}
					else{
					surfProductStoich += m_rxns[irxn]->DeltaStoich(dS)->ProdStoich();
					}
				}
				
                // Forward rate coefficients.
                arr = m_rxns[irxn]->Arrhenius();
                if (m_rxns[irxn]->IsSURF()){
					
				
					if (m_rxns[irxn]->IsFORD()) {
						
						double total_surf_Reactant_stoich_to_Replace = 0.0;
						double total_gas_Reactant_stoich_to_Replace = 0.0;
						double total_ford_surf = 0.0;
						double total_ford_gas = 0.0;
						
						for (int count_Ford = 0; count_Ford!= m_rxns[irxn]->FORDCount();count_Ford++){
						
						string Fordspecies = m_rxns[irxn]->FORDElement(count_Ford).spName;
						string FordSpeciesPhase =  FindPhaseName(Fordspecies);
						
							// calculate the total stoich of the reactants replaced by Ford	(this is if the COVERAGE with FORD or ONLY FORD )					
							double single_sp_stoich_to_Replace = m_rxns[irxn]->DeltaStoich(Fordspecies)->ReacStoich();
							
								if (FordSpeciesPhase.compare("gas") == 0){
								total_gas_Reactant_stoich_to_Replace += single_sp_stoich_to_Replace; 
								total_ford_gas += m_rxns[irxn]->FORDElement(count_Ford).F_k;
								}
								else {
								total_surf_Reactant_stoich_to_Replace += single_sp_stoich_to_Replace; 
								total_ford_surf += m_rxns[irxn]->FORDElement(count_Ford).F_k;
								}
							
						}
						
						if((m_rxns[irxn]->IsSTICK()) || (m_rxns[irxn]->IsMottWise())) {
							
							for (int count_Ford = 0; count_Ford!= m_rxns[irxn]->FORDCount();count_Ford++){
							string Fordspecies = m_rxns[irxn]->FORDElement(count_Ford).spName;
							string FordSpeciesPhase =  FindPhaseName(Fordspecies);
							
								if (FordSpeciesPhase.compare("gas") == 0){
								throw runtime_error("FORD gas species modification cannot exist with STICK/MOTTWISE ");
								}
								else {
								arr.A = arr.A ; // No Units
								}
								
							}
						}
						
						
						else if (m_rxns[irxn]->IsCOVERAGE()) { // FORD and COVERAGE
						// This assume that if FORD species should be all gas or all surface if they are more than one (just use the first ford species as a way to indicate the phase)
						string FordSp = m_rxns[irxn]->FORDElement(0).spName; 
						string FordPh = FindPhaseName(FordSp);
						
							if (FordPh.compare("gas") != 0){
							arr.A *= pow(1.0e-6, gasReactantStoich) * pow (1.0e-4, (surfReactantStoich - total_surf_Reactant_stoich_to_Replace -1 +total_ford_surf ));
							}
							else{
							arr.A *= pow(1.0e-4, surfReactantStoich - 1) * pow (1.0e-6, (gasReactantStoich - total_gas_Reactant_stoich_to_Replace +total_ford_gas ));
							}
						}
						
						else { // If it is FORD only reaction
						
						// This assume that if FORD species should be all gas or all surface if they are more than one (just use the first ford as a way to indicate the phase)
						string FordSp = m_rxns[irxn]->FORDElement(0).spName; 
						string FordPh = FindPhaseName(FordSp);
						
							if (FordPh.compare("gas") != 0){
							arr.A *= pow(1.0e-6, gasReactantStoich) * pow (1.0e-4, (surfReactantStoich - total_surf_Reactant_stoich_to_Replace -1 +total_ford_surf ));
							}
							else{
							arr.A *= pow(1.0e-4, surfReactantStoich - 1) * pow (1.0e-6, (gasReactantStoich - total_gas_Reactant_stoich_to_Replace +total_ford_gas ));
							}
							  /*

							  if ( m_rxns[irxn]->IsSURF() == false){ // Gas phase reaction
							  arr.A *= pow (1.0e-6, (gasReactantStoich - total_gas_Reactant_stoich_to_Replace +total_ford_gas ));  

							  }
							  */

						
						}
					
					}
					else if ((m_rxns[irxn]->IsSTICK()) || (m_rxns[irxn]->IsMottWise())){
					arr.A = arr.A ; // No Units
					}
					else { // COV ONLY or NORMAL SURFACE REACTIONS 
					arr.A *= pow(1.0e-6, gasReactantStoich) * pow(1.0e-4, surfReactantStoich-1.0); 
					}
				}
				else{
				arr.A *= pow(1.0e-6, m_rxns[irxn]->ReactantStoich()-1.0);
                }
				arr.E *= 1.0e-7; // Ergs to Joules
                m_rxns[irxn]->SetArrhenius(arr);

                // Reverse rate coefficients.
                if (m_rxns[irxn]->RevArrhenius() != NULL) {
                    arr = *(m_rxns[irxn]->RevArrhenius());
					
                    if (m_rxns[irxn]->IsSURF()){
					
						if ((m_rxns[irxn]->IsSTICK()) || (m_rxns[irxn]->IsMottWise())){
						arr.A = arr.A ; // No Units
						}
						else { // FORD cannot be applied to Reversible (Unless Rord also specified)
						arr.A *= pow(1.0e-6, gasProductStoich) * pow(1.0e-4, surfProductStoich-1.0); 
						}
					}
					else{
					arr.A *= pow(1.0e-6, m_rxns[irxn]->ProductStoich()-1.0);
					}
					
                    arr.E *= 1.0e-7;
                    m_rxns[irxn]->SetRevArrhenius(arr);
                }

                // Fall-off parameters.
                if (m_rxns[irxn]->FallOffType() != Kinetics::None) {
                    // Low-pressure limit.
                    arr = m_rxns[irxn]->LowPressureLimit();
                    // Note there is no -1 term here because there is also
                    // a third-body concentration.
                    arr.A *= pow(1.0e-6, m_rxns[irxn]->ReactantStoich());
                    arr.E *= 1.0e-7;
                    m_rxns[irxn]->SetLowPressureLimit(arr);
                } else {
                    // Third-body concentrations also need scaling.
                    if (m_rxns[irxn]->UseThirdBody()) {
                        // Forward rate coefficients.
                        arr = m_rxns[irxn]->Arrhenius();
                        arr.A *= 1.0e-6;
                        m_rxns[irxn]->SetArrhenius(arr);

                        // Reverse rate coefficients.
                        if (m_rxns[irxn]->RevArrhenius() != NULL) {
                            arr = *(m_rxns[irxn]->RevArrhenius());
                            arr.A *= 1.0e-6;
                            m_rxns[irxn]->SetRevArrhenius(arr);
                        }
                    }
                }
            }

            m_units = SI;
        }
    } else if (u == CGS) {
//        throw invalid_argument("Cannot currently convert mechanism to "
//                               "CGS units.  Consult your programmer.");
        m_units = CGS;
    }
}


// CHEMICAL ELEMENTS.

// Returns the number of chemical elements.
unsigned int Mechanism::ElementCount() const
{
    return m_elements.size();
}

// Returns the vector of chemical elements.
const ElementPtrVector &Mechanism::Elements() const
{
    return m_elements;
}

// Returns a pointer to the ith element.  NULL if i invalid.
const Sprog::Element *const Mechanism::Elements(unsigned int i) const
{
    if (i < m_elements.size()) {
        return m_elements[i];
    } else {
        return NULL;
    }
}

// Returns iterator to first element.
Mechanism::el_iterator Mechanism::ElBegin() {return m_elements.begin();}

// Returns const iterator to first element.
Mechanism::const_el_iterator Mechanism::ElBegin() const {return m_elements.begin();}

    // Returns iterator to position after last element.
Mechanism::el_iterator Mechanism::ElEnd() {return m_elements.end();}

// Returns const iterator to position after last element.
Mechanism::const_el_iterator Mechanism::ElEnd() const {return m_elements.end();}

// Adds a default element to the mechanism and returns it.
Element *const Mechanism::AddElement()
{
    // Add a new element.
    Element el;
    return AddElement(el);
}

// Copies the given element into the mechanism and returns a reference to
// the copy.
Element *const Mechanism::AddElement(const Sprog::Element &el)
{
    // First check current list for this element's name.  We
    // cannot have two identically named elements.
    int i = FindElement(el.Name());
    if (i >= 0) {
        // Element with this name already defined!
        return m_elements.at(i);
    }

    // Add the element.
    Element *elnew = el.Clone();
    m_elements.push_back(elnew);

    // Set up the element.
    elnew->SetMechanism(*this);

    // Return the element.
    return elnew;
}

// Returns index of element in list.  Returns -1 if not found.
int Mechanism::FindElement(const Sprog::Element &el) const
{
    // Loop over elements to find index.
    unsigned int i;
    for (i=0; i<m_elements.size(); i++) {
        if (el == *m_elements[i]) {
            // Found element!
            return i;
        }
    }

    // We are here because the element wasn't found.
    return -1;
}

// Returns index of element in list.  Returns -1 if not found.
int Mechanism::FindElement(const std::string &name) const
{
    // Loop over elements to find index.
    unsigned int i;
    for (i=0; i<m_elements.size(); i++) {
        if (*m_elements[i] == name) {
            // Found element!
            return i;
        }
    }

    // We are here because the element wasn't found.
    return -1;
}


// ELEMENT UPDATES.

// Updates the mechanism with changes to the given element.
void Mechanism::CheckElementChanges(const Element &el)
{
    // Find the index of the element in the list.
    int i = FindElement(el);
    if (i >= 0) {
        // Element found!

        // Loop over all species.  If a species contains the element, then
        // we must recalculate its mol. wt.
        SpeciesPtrVector::iterator sp;
        for (sp=m_species.begin(); sp!=m_species.end(); sp++) {
            if ((*sp)->ContainsElement(i)) {
                // Species contains element, recalculate mol. wt.
                (*sp)->CalcMolWt();
            }
        }
    }
}


// SPECIES.

// Returns the total number of species in the mechanism.
unsigned int Mechanism::SpeciesCount(void) const
{
    return m_species.size();
}

// Returns the total number of gas species in the mechanism.
unsigned int Mechanism::GasSpeciesCount(void) const
{

  unsigned int gas_sp = 0;
  for (unsigned int j = 0; j < m_species.size(); ++j ){
    std::string nm = FindPhaseName(m_species[j]->Name()); 
    
    if (FindID(nm).compare("g")==0){
      gas_sp++; 			     
    }
 }
  
  return gas_sp;

}

// Returns the vector of chemical species.
const SpeciesPtrVector &Mechanism::Species() const
{
    return m_species;
}

// Returns a pointer to the ith species.  Returns NULL if i is invalid.
const Sprog::Species *const Mechanism::Species(unsigned int i) const
{
    if (i < m_species.size()) {
        return m_species[i];
    } else {
        return NULL;
    }
}

// Returns pointer to species with given name.  NULL if not found.
const Sprog::Species *const Mechanism::Species(const std::string &name) const
{
    int i = FindSpecies(name);
    if (i >= 0) {
        return m_species[i];
    } else {
        return NULL;
    }
}

// Returns iterator to first element.
Mechanism::sp_iterator Mechanism::SpBegin() {return m_species.begin();}

// Returns const iterator to first element.
Mechanism::const_sp_iterator Mechanism::SpBegin() const {return m_species.begin();}

    // Returns iterator to position after last element.
Mechanism::sp_iterator Mechanism::SpEnd() {return m_species.end();}

// Returns const iterator to position after last element.
Mechanism::const_sp_iterator Mechanism::SpEnd() const {return m_species.end();}

// Adds an empty species to the mechansism and returns a reference to it.
Species *const Mechanism::AddSpecies()
{
    // Adds species to vector.
    Sprog::Species sp;
    return AddSpecies(sp);
}

// Copies species into mechanism and returns a reference to
// the copy.
Species *const Mechanism::AddSpecies(const Sprog::Species &sp)
{
    // First check species vector to see if a species with the
    // same name is already defined.  We can only have one species
    // per name.
    int i = FindSpecies(sp.Name());
    if (i >= 0) {
        // A species with this name has been found.
        return m_species.at(i);
    }

    // Adds species to vector.
    Sprog::Species *spnew = sp.Clone();
    m_species.push_back(spnew);

    // Set up species.
    spnew->SetMechanism(*this);

    // Return species reference.
    return spnew;
}

// Returns index of species.  Returns -1 if not found.
int Mechanism::FindSpecies(const Sprog::Species &sp) const
{
    // Loop over species to find index.
    unsigned int i;
    for (i=0; i<m_species.size(); i++) {
        if (sp == *m_species[i]) {
            // Found species!
            return i;
        }
    }

    // We are here because the species wasn't found.
    return -1;
}

// Returns index of species.  Returns -1 if not found.
int Mechanism::FindSpecies(const std::string &name) const
{
    // Loop over species to find index.
    unsigned int i;
    for (i=0; i<m_species.size(); i++) {
        if (*m_species[i] == name) {
            // Found species!
            return i;
        }
    }

    // We are here because the species wasn't found.
    return -1;
}

// Returns a pointer to the species at the given index.  Returns NULL if not found.
Sprog::Species *const Mechanism::GetSpecies(const unsigned int i) const
{
    if (i < m_species.size()) {
        return m_species[i];
    } else {
        return NULL;
    }
}

// Returns a pointer to the species with the given name.  Returns NULL if not found.
Sprog::Species *const Mechanism::GetSpecies(const std::string &name) const
{
    int i = FindSpecies(name);
    if (i >= 0) {
        return m_species[i];
    } else {
        return NULL;
    }
}

// Returns site occupancy of species.  Returns 0 if not found.
int Mechanism::FindSiteOccup(const std::string &name) const

{
   int i = FindSpecies(name);
    if (i >= 0) {
      return m_species[i]->SiteOccupancy();
    } else {
        return 0;
    }

}

// Returns phase Name of species.  Returns "" if not found.
std::string Mechanism::FindPhaseName(const std::string &name) const

{
   int i = FindSpecies(name);
    if (i >= 0) {
      return m_species[i]->PhaseName();
    } else {
        return "";
    }

}



// CHEMICAL PHASES.

// Returns the number of phases.
unsigned int Mechanism::PhaseCount(void) const
{
    return m_phase.size();
}

// Returns the vector of phases.
const PhasePtrVector &Mechanism::Phase() const
{
    return m_phase;
}

// Returns a pointer to the ith phase.  NULL if i invalid.
const Sprog::Phase *const Mechanism::Phase(unsigned int i) const
{
    if (i < m_phase.size()) {
        return m_phase[i];
    } else {
        return NULL;
    }
}

// Returns pointer to phase with given name.  NULL if not found.
const Sprog::Phase *const Mechanism::Phase(const std::string &name) const
{
    int i = FindPhase(name);
    if (i >= 0) {
        return m_phase[i];
    } else {
        return NULL;
    }
}

// Returns iterator to first phase.
Mechanism::phase_iterator Mechanism::PhaseBegin() 
{return m_phase.begin();}

// Returns const iterator to first phase.
Mechanism::const_phase_iterator Mechanism::PhaseBegin() const 
{return m_phase.begin();}

// Returns iterator to position after last phase.
Mechanism::phase_iterator Mechanism::PhaseEnd() 
{return m_phase.end();}

// Returns const iterator to position after last phase.
Mechanism::const_phase_iterator Mechanism::PhaseEnd() const 
{return m_phase.end();}

// Adds an empty phase to the mechansism and returns a reference to it.
Phase *const Mechanism::AddPhase()
{
    // Adds species to vector.
    Sprog::Phase ph;
    return AddPhase(ph);
}

// Copies phase into mechanism and returns a reference to
// the copy.
Phase *const Mechanism::AddPhase(const Sprog::Phase &phase)
{
    // First check phase vector to see if a phase with the
    // same name is already defined.  We can only have one phase
    // per name.
    int i = FindPhase(phase.Name());
    if (i >= 0) {
        // A phase with this name has been found.
        return m_phase.at(i);
    }

    // Adds phase to vector.
    Sprog::Phase *phasenew = phase.Clone();
    m_phase.push_back(phasenew);

    // Set up phase.
    phasenew->SetMechanism(*this);

    // Return phase reference.
    return phasenew;
}

// Returns index of phase.  Returns -1 if not found.
int Mechanism::FindPhase(const Sprog::Phase &phase) const
{
    // Loop over phases to find index.
    unsigned int i;
    for (i=0; i<m_phase.size(); i++) {
        if (phase == *m_phase[i]) {
            // Found phase!
            return i;
        }
    }

    // We are here because the phase wasn't found.
    return -1;
}



// Returns index of phase.  Returns -1 if not found.
int Mechanism::FindPhase(const std::string &name) const
{
    // Loop over phases to find index.
    unsigned int i;
    for (i=0; i<m_phase.size(); i++) {
        if (*m_phase[i] == name) {
            // Found phase!
            return i;
        }
    }

    // We are here because the phase wasn't found.
    return -1;
}


// Returns the phase id given the phase name.  Returns NULL if not found.
string Mechanism::FindID(const std::string &name) const
{
    int i = FindPhase(name);
    if (i >= 0) {
      return m_phase[i]->ID();
    } else {
        return "";
    }
}

// Returns the phase site density given the phase name.  Returns NULL if not found.
double Mechanism::FindSiteDensity(const std::string &name) const
{
    int i = FindPhase(name);
    if (i >= 0) {
      return m_phase[i]->SiteDen();
    } else {
        return 0.0;
    }
}


// Returns a pointer to the phase at the given index.  Returns NULL if not found.
Sprog::Phase *const Mechanism::GetPhase(const unsigned int i) const
{
    if (i < m_phase.size()) {
        return m_phase[i];
    } else {
        return NULL;
    }
}

// Returns a pointer to the phase with the given name.  Returns NULL if not found.
Sprog::Phase *const Mechanism::GetPhase(const std::string &name) const
{
    int i = FindPhase(name);
    if (i >= 0) {
        return m_phase[i];
    } else {
        return NULL;
    }
}




// REACTIONS.

// Returns the number of reactions in the mechanism.
unsigned int Mechanism::ReactionCount(void) const
{
    return m_rxns.Count();
}

// Returns the set of chemical reactions.
const Sprog::Kinetics::ReactionSet &Mechanism::Reactions() const
{
    return m_rxns;
}

// Returns a pointer to the ith reaction. Returns NULL if i is invalid.
const Kinetics::Reaction *const Mechanism::Reactions(unsigned int i) const
{
    return m_rxns[i];
}
// Returns a pointer to the ith reaction. Returns NULL if i is invalid.
Kinetics::Reaction * Mechanism::GetReactions(unsigned int i)
{
    return m_rxns[i];
}

// Adds an empty reaction to the mechanism.
Sprog::Kinetics::Reaction *const Mechanism::AddReaction()
{
    return AddReaction(new Sprog::Kinetics::Reaction());
}

// Copies the reaction to the mechanism.
Sprog::Kinetics::Reaction *const Mechanism::AddReaction(const Sprog::Kinetics::Reaction *const rxn)
{
    Sprog::Kinetics::Reaction *prxn = m_rxns.AddReaction(*rxn);
    prxn->SetMechanism(*this);
    return prxn;
}


// STOICHIOMETRY CROSS REFERENCE.

// Builds the species-reaction stoichiometry cross-reference table.
void Mechanism::BuildStoichXRef()
{
    unsigned int i, j;
    int k;
    RxnStoichMap::iterator ij;
    double mu;

    // Clear current table.
    m_stoich_xref.clear();

    // Set up empty table.
    for (i=0; i!=m_species.size(); ++i) {
        m_stoich_xref.push_back(StoichXRef());
        m_stoich_xref[i].Species = i;
    }

    // Loop over all reactions.
    for (j=0; j!=m_rxns.Count(); ++j) {
        // Sum up integer reactant stoich.
        for (k=0; k!=m_rxns[j]->ReactantCount(); ++k) {
            // Get the species index and the stoichiometry.
            i  = m_rxns[j]->Reactant(k).Index();
            mu = m_rxns[j]->Reactant(k).Mu();

            // Add up the contribution of this reaction to this species.
            ij = m_stoich_xref[i].RxnStoich.find(j);
            if (ij != m_stoich_xref[i].RxnStoich.end()) {
                (*ij).second -= mu;
            } else {
                m_stoich_xref[i].RxnStoich.insert(RxnStoichPair(j,-mu));
            }
        }

        // Sum up integer product stoich.
        for (k=0; k!=m_rxns[j]->ProductCount(); ++k) {
            // Get the species index and the stoichiometry.
            i  = m_rxns[j]->Product(k).Index();
            mu = m_rxns[j]->Product(k).Mu();

            // Add up the contribution of this reaction to this species.
            ij = m_stoich_xref[i].RxnStoich.find(j);
            if (ij != m_stoich_xref[i].RxnStoich.end()) {
                (*ij).second += mu;
            } else {
                m_stoich_xref[i].RxnStoich.insert(RxnStoichPair(j,mu));
            }
        }

    }

    m_stoich_xref_valid = true;
}

// Returns true if the stoichiometry xref map is valid, otherwise false.
bool Mechanism::IsStoichXRefValid()
{
    return m_stoich_xref_valid;
}

// Returns the stoichiometry for all reactions which
// involve the species with the given index.  Throws error if index is
// invalid.
const Sprog::RxnStoichMap &Mechanism::GetStoichXRef(unsigned int isp) const
{
    if (isp < m_species.size()) {
        return m_stoich_xref[isp].RxnStoich;
    }
    // Species index is invalid.
    throw invalid_argument("Invalid species index given when "
                           "finded RxnStoichMap (Sprog, Mechanism::GetStoichXRef.");
}

// COPYING ROUTINES.

// Copies elements from given array into this mechanism.
void Mechanism::copyInElements(const Sprog::ElementPtrVector &els)
{
    ElementPtrVector::const_iterator el;
    for (el=els.begin(); el!=els.end(); el++) {
        // Use Clone() function to create a copy of the element object.
        m_elements.push_back((*el)->Clone());
    }
}

// Copies species from given array into this mechanism.
void Mechanism::copyInSpecies(const Sprog::SpeciesPtrVector &sps)
{
    SpeciesPtrVector::const_iterator sp;
    for (sp=sps.begin(); sp!=sps.end(); sp++) {
        // Use Clone() function to create a copy of the species object.
        m_species.push_back((*sp)->Clone());
    }
}

// Copies phase from given array into this mechanism.
void Mechanism::copyInPhase(const Sprog::PhasePtrVector &phs)
{
    PhasePtrVector::const_iterator ph;
    for (ph=phs.begin(); ph!=phs.end(); ph++) {
        // Use Clone() function to create a copy of the phase object.
        m_phase.push_back((*ph)->Clone());
    }
}


// MEMORY MANAGEMENT.

// Clears memory used by the mechanism object.
void Mechanism::releaseMemory()
{
    // Clear elements.
    ElementPtrVector::iterator el;
    for (el=m_elements.begin(); el!=m_elements.end(); el++) {
        delete *el;
    }
    m_elements.clear();

    // Clear species.
    SpeciesPtrVector::iterator sp;
    for (sp=m_species.begin(); sp!=m_species.end(); sp++) {
        delete *sp;
    }
    m_species.clear();

	// Clear phases.
    PhasePtrVector::iterator ph;
    for (ph=m_phase.begin(); ph!=m_phase.end(); ph++) {
        delete *ph;
    }
    m_phase.clear();
	
    // Clear reactions.
    m_rxns.Clear();

    // Clear other variables.
    m_stoich_xref.clear();
}


// OUTPUT FUNCTIONS.

// Prints a diagnostic output file containing all the
// mechanism data.  This is used to debug.
void Mechanism::WriteDiagnostics(const std::string &filename) const
{
    string data = "";

    // Open the output file.
    ofstream fout;
    fout.open(filename.c_str());

    // Write the Elements to the file.
    data = "Elements:\n";
    fout.write(data.c_str(), data.length());
    for (unsigned int i=0; i!=m_elements.size(); ++i) {
        m_elements[i]->WriteDiagnostics(fout);
    }
    data = "End of elements.\n";
    fout.write(data.c_str(), data.length());

    // Write the Species to the file.
    data = "Species:\n";
    fout.write(data.c_str(), data.length());
    for (unsigned int i=0; i!=m_species.size(); ++i) {
        m_species[i]->WriteDiagnostics(fout);
    }
    data = "End of species.\n";
    fout.write(data.c_str(), data.length());

    // Write the Phase to the file.
    data = "Phase:\n";
    fout.write(data.c_str(), data.length());
    for (unsigned int i=0; i!=m_phase.size(); ++i) {
        m_phase[i]->WriteDiagnostics(fout);
    }
    data = "End of phase.\n";
    fout.write(data.c_str(), data.length());
	
    // Write the reactions to the file.
    data = "Reactions:\n";
    fout.write(data.c_str(), data.length());
    for (unsigned int i=0; i!=m_rxns.Count(); ++i) {
        m_rxns.Reactions(i)->WriteDiagnostics(fout);
    }
    data = "End of reactions.\n";
    fout.write(data.c_str(), data.length());

    // Close the output file.
    fout.close();
}

/*!
@param[in]      filename            The name of the reduced mechanism output file.
@param[in]      RejectSpecies       The rejected species as specified by LOI reduction.
@param[in]      KeptSpecies         The kept species as specified by LOI reduction.
*/
void Mechanism::WriteReducedMech(const std::string &filename, std::vector<std::string> RejectSpecies) const {

    // Open the output file.
    ofstream fout;
    fout.open(filename.c_str());

    // Write the Elements to the file.
    fout << "ELEMENTS\n";
    for (unsigned int i = 0; i != m_elements.size(); ++i) {
        m_elements[i]->WriteElements(fout);
    }
    fout << "END\n\n";

    // Write the Species to the file.
    fout << "SPECIES\n";
    for (unsigned int i = 0; i != m_species.size(); ++i) {
        if (m_species[i]->Name() != RejectSpecies[i])
            m_species[i]->WriteSpecies(fout);
    }
    fout << "END\n\n";

	/*
	// Write the Phase to the file.
    fout << "PHASE\n";
    for (unsigned int i = 0; i != m_phase.size(); ++i) {
        if (m_phase[i]-> !(ContainsSpecies(RejectSpecies[i]))) // If the phase doesn't contain rejected species
            m_phase[i]->WritePhase(fout);
    }
    fout << "END\n\n";
	*/
	
    // Write the reactions to the file.
    fout << "REAC\n";
    for (unsigned int i = 0; i != m_rxns.Count(); ++i) {
        m_rxns.Reactions(i)->WriteReducedMechReacs(fout, RejectSpecies);
    }
    fout << "END\n\n";

    // Close the output file.
    fout.close();
}

void Mechanism::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialize version to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write the units system to the stream.
        unsigned int u = (unsigned int)m_units;
        out.write((char*)&u, sizeof(u));

        // Write the number of elements to the stream.
        u = m_elements.size();
        out.write((char*)&u, sizeof(u));

        // Write the elements to the stream.
        for (ElementPtrVector::const_iterator iel=m_elements.begin(); iel!=m_elements.end(); iel++) {
            (*iel)->Serialize(out);
        }

        // Write the number of species to the stream.
        u = m_species.size();
        out.write((char*)&u, sizeof(u));

        // Write the species to the stream.
        for (SpeciesPtrVector::const_iterator isp=m_species.begin(); isp!=m_species.end(); isp++) {
            (*isp)->Serialize(out);
        }

	// Write the number of phase to the stream.
        u = m_phase.size();
        out.write((char*)&u, sizeof(u));

        // Write the phase to the stream.
        for (PhasePtrVector::const_iterator iph=m_phase.begin(); iph!=m_phase.end(); iph++) {
            (*iph)->Serialize(out);
        }
		
        // Write the reaction set to the stream.
        m_rxns.Serialize(out);

        // We don't need to write the stoich xref vector to the stream, as
        // this can be rebuilt when the mechanism is deserialised.

    } else {
        throw invalid_argument("Output stream not ready (Sprog, Mechanism::Serialize).");
    }
}

// Reads the mechanism data from a binary data stream.
void Mechanism::Deserialize(std::istream &in)
{
    // Clear the current mechanism.  We do this before checking
    // the stream condition to avoid confusion in the calling code.
    // Even if the possible exception is handled incorrectly, the
    // mechanism will still be empty.
    releaseMemory();

    if (in.good()) {
        // Read the serialized mechanism version.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int u = 0; // We'll need to to read unsigned ints.

        switch (version) {
            case 0:
                // Read the units.
                in.read(reinterpret_cast<char*>(&u), sizeof(u));
                m_units = (UnitSystem)u;

                // Read the number of elements and reserve memory.
                in.read(reinterpret_cast<char*>(&u), sizeof(u));
                m_elements.reserve(u);

                // Read the elements.
                try {
                    for (unsigned int i=0; i<u; i++) {
                        // Read the element from the stream using the
                        // appropriate constructor.
                        Element *el = new Element(in);
                        el->SetMechanism(*this);

                        // Add the element to the vector.
                        m_elements.push_back(el);
                    }
                } catch (exception &e) {
                    // Ensure the mechanism is cleared before throwing
                    // the exception to the next layer up.
                    releaseMemory();
                    throw;
                }

                // Read the number of species and reserve memory.
                in.read(reinterpret_cast<char*>(&u), sizeof(u));
                m_species.reserve(u);
                // Read the species.
                try {
                    for (unsigned int i=0; i<u; i++) {
                        // Read the species from the stream using the
                        // appropriate constructor.
                        Sprog::Species *sp = new Sprog::Species(in);
                        sp->SetMechanism(*this);
                        // Add the species to the vector.
                        m_species.push_back(sp);
                    }
                } catch (exception &e) {
                    // Ensure the mechanism is cleared before throwing
                    // the exception to the next layer up.
                    releaseMemory();
                    throw;
                }
		
	       	// Read the number of phase and reserve memory.
                in.read(reinterpret_cast<char*>(&u), sizeof(u));
                m_phase.reserve(u);
	
                // Read the phase.
                try {
                    for (unsigned int i=0; i<u; i++) {
                        // Read the phase from the stream using the
                        // appropriate constructor.
                        Sprog::Phase *ph = new Sprog::Phase(in);
                        ph->SetMechanism(*this);

                        // Add the species to the vector.
                        m_phase.push_back(ph);
                    }
                } catch (exception &e) {
                    // Ensure the mechanism is cleared before throwing
                    // the exception to the next layer up.
                    releaseMemory();
                    throw;
                }
			
                // Read the reaction set.
                try {
                    m_rxns.SetMechanism(*this);
                    m_rxns.Deserialize(in);
                } catch (exception &e) {
                    // Ensure the mechanism is cleared before throwing
                    // the exception to the next layer up.
                    releaseMemory();
                    throw;
                }
		  
                // Rebuild the stoich xref.
                BuildStoichXRef();

                break;
            default:
                throw runtime_error("Mechanism serialized version number "
                                    "is unsupported (Sprog, Mechanism::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, Mechanism::Deserialize).");
    }
}
