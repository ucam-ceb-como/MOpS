/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This files contains the definition of a chemical mechanism.  A chemical mechanism consists
    of 3 parts:  definition of chemical elements, definition of chemical species and a set of
    reactions which define the systems kinetics.  By combining these items into a single
    class the code can be made more efficient.

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

#ifndef GPC_MECH_H
#define GPC_MECH_H

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include "gpc_params.h"
#include "gpc_element.h"
#include "gpc_phase.h"
#include "gpc_species.h"
#include "gpc_reaction_set.h"
#include "gpc_unit_systems.h"
#include "gpc_stoich.h"
#include <vector>
#include <string>
#include <map>



namespace Sprog
{
class Mechanism
{
public:
    // Element iterators.
    typedef ElementPtrVector::iterator el_iterator;
    typedef ElementPtrVector::const_iterator const_el_iterator;

    // Species iterators.
    typedef SpeciesPtrVector::iterator sp_iterator;
    typedef SpeciesPtrVector::const_iterator const_sp_iterator;

    // Phase iterators.
    typedef PhasePtrVector::iterator phase_iterator;
    typedef PhasePtrVector::const_iterator const_phase_iterator;

    // Constructors.
    Mechanism();                  // Default constructor.
    Mechanism(const Mechanism &mech); // Copy constructor.
    Mechanism(std::istream &in);      // Stream-reading constructor.

    // Destructors.
    ~Mechanism();

    // Operator overloads.
    Mechanism &operator=(const Mechanism &mech);

    // Empties the mechanism of all elements, species and reactions.
    void Clear();


    // UNITS.

    // Returns the current unit system.
    UnitSystem Units(void) const;

    // Converts the mechanism to a new units system.
    void SetUnits(UnitSystem u);


    // CHEMICAL ELEMENTS.

    // Returns the number of elements in the mechanism.
    unsigned int ElementCount(void) const;

    // Returns the vector of elements.
    const ElementPtrVector &Elements(void) const;

    // Returns a pointer to the ith element.  NULL if i invalid.
    const Sprog::Element *const Elements(unsigned int i) const;

    // Returns iterator to first element.
    el_iterator ElBegin();

    // Returns const iterator to first element.
    const_el_iterator ElBegin() const;

        // Returns iterator to position after last element.
    el_iterator ElEnd();

    // Returns const iterator to position after last element.
    const_el_iterator ElEnd() const;

    // Adds a default element to the mechanism.
    Sprog::Element *const AddElement(void);

    // Copies the given element into the mechanism.
    Sprog::Element *const AddElement(const Sprog::Element &el);

    // Returns index of element.  Returns -1 if not found.
    int FindElement(const std::string &name) const;

    // Returns index of element.  Returns -1 if not found.
    int FindElement(const Sprog::Element &el) const;

    // Updates mechanism with changes applied to an element.
    void CheckElementChanges(const Sprog::Element &el);


    // 

    // SPECIES.

    // Returns the number of species in the mechanism.
    unsigned int SpeciesCount(void) const;

    // Returns the number of gas phase species in the mechanism. 
    unsigned int GasSpeciesCount(void) const;

    // Returns the vector of species.
    const SpeciesPtrVector &Species(void) const;

    // Returns a pointer to the ith species.  Returns NULL if i is invalid.
    const Sprog::Species *const Species(unsigned int i) const;

    // Returns pointer to species with given name.  NULL if not found.
    const Sprog::Species *const Species(const std::string &name) const;

    // Returns iterator to first species.
    sp_iterator SpBegin();

    // Returns const iterator to first species.
    const_sp_iterator SpBegin() const;

    // Returns iterator to position after last species.
    sp_iterator SpEnd();

    // Returns const iterator to position after last species.
    const_sp_iterator SpEnd() const;

    // Adds an empty species to the mechanism.
    Sprog::Species *const AddSpecies(void);

    // Copies given species into the mechanism.
    Sprog::Species *const AddSpecies(const Sprog::Species &sp);

    // Returns index of species.  Returns -1 if not found.
    int FindSpecies(const Sprog::Species &sp) const;

    // Returns index of species.  Returns -1 if not found.
    int FindSpecies(const std::string &name) const;

    // Returns pointer to species at index i.  NULL if not found.  This
    // function returns a modifiable (non-const) species object.
    Sprog::Species *const GetSpecies(const unsigned int i) const;

    // Returns pointer to species with given name.  NULL if not found.  This
    // function returns a modifiable (non-const) species object.
    Sprog::Species *const GetSpecies(const std::string &name) const;


    // Returns site occupancy of species.  Returns 0 if not found.
    int FindSiteOccup(const std::string &name) const;

    // Returns phase Name of species.  Returns "" if not found.
    std::string FindPhaseName(const std::string &name) const;

 // PHASE.

    // Returns the number of phase in the mechanism.
    unsigned int PhaseCount(void) const;

    // Returns the vector of phase.
    const PhasePtrVector &Phase(void) const;

    // Returns a pointer to the ith phase.  Returns NULL if i is invalid.
    const Sprog::Phase *const Phase(unsigned int i) const;

    // Returns pointer to phase with given name.  NULL if not found.
    const Sprog::Phase *const Phase(const std::string &name) const;

    // Returns iterator to first phase.
    phase_iterator PhaseBegin();

    // Returns const iterator to first phase.
    const_phase_iterator PhaseBegin() const;

    // Returns iterator to position after last phase.
    phase_iterator PhaseEnd();

    // Returns const iterator to position after last phase.
    const_phase_iterator PhaseEnd() const;

    // Adds an empty phase to the mechanism.
    Sprog::Phase *const AddPhase(void);

    // Copies given phase into the mechanism.
    Sprog::Phase *const AddPhase(const Sprog::Phase &phase);

    // Returns index of phase.  Returns -1 if not found.
    int FindPhase(const Sprog::Phase &phase) const;

    // Returns index of phase.  Returns -1 if not found.
    int FindPhase(const std::string &name) const;

    // Returns phaseid given the name.  Returns NULL if not found.
    std::string FindID(const std::string &name) const;

	// Returns phase site density given the name.  Returns NULL if not found.
    double FindSiteDensity(const std::string &name) const;

	
    // Returns pointer to phase at index i.  NULL if not found.  This
    // function returns a modifiable (non-const) phase object.
    Sprog::Phase *const GetPhase(const unsigned int i) const;

    // Returns pointer to phase with given name.  NULL if not found.  This
    // function returns a modifiable (non-const) phase object.
    Sprog::Phase *const GetPhase(const std::string &name) const;


    // REACTIONS.

    // Returns the number of reactions in the mechanism.
    unsigned int ReactionCount(void) const;

    // Returns the reaction set.
    const Kinetics::ReactionSet &Reactions(void) const;

    // Returns a pointer to the ith reaction. Returns NULL if i is invalid.
    const Kinetics::Reaction *const Reactions(unsigned int i) const;
    // Returns a pointer to the ith reaction. Returns NULL if i is invalid.
    Kinetics::Reaction * GetReactions(unsigned int i);

    // Adds an empty reaction to the mechanism.
    Kinetics::Reaction *const AddReaction(void);

    // Copies a reaction into the mechanism.
    Kinetics::Reaction *const AddReaction(const Kinetics::Reaction *const rxn);


    // SPECIES-REACTIONS STOICHIOMETRY CROSS-REFERENCE.

    // Builds the species-reaction stoichiometry cross-reference table.
    void BuildStoichXRef();

    // Returns true if the stoich xref map is valid, otherwise false.
    bool IsStoichXRefValid();

    // Returns the stoichiometry for all reactions which
    // involve the species with the given index.  Throws error
    // if index is invalid.
    const RxnStoichMap &GetStoichXRef(unsigned int isp) const;


    // READ/WRITE/COPY FUNCTIONS.

    // Writes the element to a binary data stream.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
          ar.register_type(static_cast<Sprog::Element *>(NULL));
          ar.register_type(static_cast<Sprog::Species *>(NULL));
	  ar.register_type(static_cast<Sprog::Phase *>(NULL));
          ar.template register_type<Sprog::StoichXRef>();

          ar & m_units & m_elements & m_species & m_phase & m_rxns & m_stoich_xref
             & m_stoich_xref_valid & m_nec_spec;
    }

    friend class boost::serialization::access;

    // Writes the mixture to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the mixture data from a binary data stream.
    void Deserialize(std::istream &in);

    // OUTPUT FUNCTIONS.

    // Prints a diagnostic output file containing all the
    // mechanism data.  This is used to debug.
    void WriteDiagnostics(const std::string &filename) const;

    //! Creates and writes out a Chemkin file containing the reduced mechanism as specified by LOI analysis
    void WriteReducedMech(const std::string &filename, std::vector<std::string> RejectSpecies) const;

protected:
    // Mechanism data.
    UnitSystem m_units;             // The system of units used by this mechanism.
    ElementPtrVector m_elements;    // Vector of chemical elements defined by mechanism.
    SpeciesPtrVector m_species;     // Vector of chemical species defined by mechanism.
    PhasePtrVector m_phase;         // Vector of chemical phase defined by mechanism.
	Kinetics::ReactionSet m_rxns;   // Set of reactions defined by mechanism.
    StoichXRefVector m_stoich_xref; // Reaction stoichiometry cross-referenced for each species.
    bool m_stoich_xref_valid;       // Flag which tells whether or not the stoich xref map is valid.
    std::vector<std::string> m_nec_spec;    //String vector containing the species that must be present in the reduced mechanism.

 
    // COPYING ROUTINES.

    // Copies elements from given array into this mechanism.
    void copyInElements(const ElementPtrVector &els);

    // Copies species from given array into this mechanism.
    void copyInSpecies(const SpeciesPtrVector &sps);

    // Copies phase from given array into this mechanism.
    void copyInPhase(const PhasePtrVector &phs);
	
	
    // MEMORY MANAGEMENT.

    // Clears memory used by the mechanism object.
    void releaseMemory(void);
};
};

#endif
