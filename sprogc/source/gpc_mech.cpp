#include "gpc_mech.h"
#include "gpc_unit_systems.h"
#include <string>

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Mechanism::Mechanism()
{
    m_units = SI;
}

// Copy constructor.
Mechanism::Mechanism(const Sprog::Mechanism &mech)
{
    m_units = SI;
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
        copyInElements(mech.m_elements);
        copyInSpecies(mech.m_species);

        // Copy reaction set and stoich cross-referencing.
        m_rxns = mech.m_rxns;
        m_stoich_xref.assign(mech.m_stoich_xref.begin(), mech.m_stoich_xref.end());

        // Inform species of new elements vector and mechanism.
        SpeciesPtrVector::iterator sp;
        for (sp=m_species.begin(); sp!=m_species.end(); sp++) {
            (*sp)->SetElements(&m_elements);
            (*sp)->SetMechanism(this);
        }

        // Inform reactions of new species vector and mechanism.
        unsigned int i;
        for (i=0; i<m_rxns.Count(); i++) {
            m_rxns[i]->SetSpecies(&m_species);
            m_rxns[i]->SetMechanism(this);
        }
    }
    return *this;
}


// CHEMICAL ELEMENTS.

// Returns the vector of chemical elements.
const ElementPtrVector &Mechanism::Elements() const
{
    return m_elements;
}

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
    elnew->SetMechanism(this);

    // Return the element.
    return elnew;
}

// Returns index of element in list.  Returns -1 if not found.
int Mechanism::FindElement(const Sprog::Element &el)
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
int Mechanism::FindElement(const std::string &name)
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

// Returns the vector of chemical species.
const SpeciesPtrVector &Mechanism::Species() const
{
    return m_species;
}

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
    spnew->SetElements(&m_elements);
    spnew->SetMechanism(this);

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


// REACTIONS.

// Returns the set of chemical reactions.
const Sprog::Kinetics::ReactionSet &Mechanism::Reactions() const
{
    return m_rxns;
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
    prxn->SetSpecies(&m_species);
    prxn->SetMechanism(this);
    return prxn;
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

    // Clear reactions.
    m_rxns.Clear();

    // Clear other variables.
    m_stoich_xref.clear();
}