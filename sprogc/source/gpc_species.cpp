#include "gpc_species.h"
#include "gpc_el_comp.h"
#include "gpc_mech.h"
#include <exception>
#include <stdexcept>
#include <string>

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Species::Species(void)
{
    m_name = "";
    m_molwt = 0.0;
    m_T1 = 0.0;
    m_elements = NULL;
    m_mech = NULL;
}

// Copy constructor.
Species::Species(const Sprog::Species &sp)
{
    m_elements = NULL;
    m_mech = NULL;
    *this = sp;
}

// Destructor.
Species::~Species(void)
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
Species &Species::operator=(const Sprog::Species &sp)
{
    // Remember to check for self-assignment!
    if (this!=&sp) {
        // Copy simple member data.
        m_name = sp.m_name;
        m_elcomp.assign(sp.m_elcomp.begin(), sp.m_elcomp.end());
        m_molwt = sp.m_molwt;
        m_T1 = sp.m_T1;
        m_elements = sp.m_elements;
        m_mech = sp.m_mech;

        // Copy thermo parameters.
        map<real,Thermo::THERMO_PARAMS>::const_iterator i;
        for (i=sp.m_thermoparams.begin(); i!=sp.m_thermoparams.end(); i++) {
            m_thermoparams[i->first] = i->second;
        }
    }
    return *this;
}

// Comparison operator:  Compares two Species objects.  Returns true
// if they have the same name.
bool Species::operator==(const Sprog::Species &sp) const
{
    return m_name.compare(sp.m_name)==0;
}

// Comparison operator:  Compares a species object to a string.  Returns
// true if the species name is the same as the string.
bool Species::operator==(const std::string &name) const
{
    return m_name.compare(name)==0;
}

// Inequality operator:  Returns false of the two species objects have the
// same name.
bool Species::operator!=(const Sprog::Species &sp) const
{
    return !(*this==sp);
}

// Inequality operator:  Returns false if the species name matches
// the string.
bool Species::operator!=(const std::string &name) const
{
    return !(*this==name);
}


// SPECIES NAME.

// Sets the species name.
void Species::SetName(const std::string &name)
{
    // If this species is part of a mechanism then we must check
    // if a species with this name is already defined.
    if (m_mech != NULL) {
        if (m_mech->FindSpecies(name) >= 0) {
            // Oh dear:  Species with this name already defined.  We'd
            // better throw an error.
            throw invalid_argument("Cannot have two species with the same name (Species::SetName).");
        }
    }

    // It is ok to set the name now.
    m_name = name;
}


// SPECIES COMPOSITION.

// Adds an element to the species composition vector.
void Species::AddElement(const Sprog::ElComp &elcomp)
{
    if (m_elements != NULL) {
        if ((elcomp.Index() >= 0) && (elcomp.Index() < m_elements->size())) {
            // Must check if this element is already defined for this
            // species.
            ElCompVector::iterator i;
            for (i=m_elcomp.begin(); i!=m_elcomp.end(); i++) {
                if ((*i) == elcomp) {
                    // Found element:  Append value.
                    (*i) += elcomp.Count();
                    break;
                }
            }

            // Element not currently defined for species, so add it.
            m_elcomp.push_back(elcomp);

            // Update species molecular weight:  Sum of element weights.
            m_molwt += ((*m_elements)[elcomp.Index()]->MolWt() * (real)elcomp.Count());
        } else {
            // Element index in elcomp is out-of-range.
            throw out_of_range("Element index is out of range (Sprog::Species::AddElement).");
        }
    } else {
        // Element vector has not been assigned yet.
        throw exception("Attempt to add element before elements vector is assigned.");
    }
}

// Adds an element to the species given the index and count.
void Species::AddElement(unsigned int i, unsigned int n)
{
    if (m_elements != NULL) {
        if ((i >= 0) && (i < m_elements->size())) {
            // Must check if this element is already defined for this
            // species.
            ElCompVector::iterator j;
            for (j=m_elcomp.begin(); j!=m_elcomp.end(); j++) {
                if ((*j).Index() == i) {
                    // Found element:  Append value.
                    (*j) += n;
                    break;
                }
            }

            // Element not currently defined for species, so add it.
            m_elcomp.push_back(ElComp(i, n));

            // Update species molecular weight:  Sum of element weights.
            m_molwt += ((*m_elements)[i]->MolWt() * (real)n);
        } else {
            // Element index is out-of-range.
            throw out_of_range("Element index is out of range (Sprog::Species::AddElement).");
        }
    } else {
        // Element vector has not been assigned yet.
        throw exception("Attempt to add element before elements vector is assigned.");
    }
}

// Adds an element to the species given the element name and count.
void Species::AddElement(const std::string &name, unsigned int n)
{
    if (m_elements != NULL) {
        // Must find the index of the element in the list.
        unsigned int i;
        for (i=0; i<m_elements->size(); i++) {
            if (name.compare(m_elements->at(i)->Name()) == 0) {
                // Found element in the list!

                // Must check if this element is already defined for this
                // species.
                ElCompVector::iterator j;
                for (j=m_elcomp.begin(); j!=m_elcomp.end(); j++) {
                    if ((*j).Index() == i) {
                        // Found element:  Append value.
                        (*j) += n;
                        break;
                    }
                }

                // Element not currently defined for species, so add it.
                m_elcomp.push_back(ElComp(i, n));

                // Update species molecular weight:  Sum of element weights.
                m_molwt += ((*m_elements)[i]->MolWt() * (real)n);

                return;
            }
        }
        
        // We have got here because the element wasn't found in the list.
        throw invalid_argument(string(name).append(" not found in element list (Sprog::Species::AddElement)."));
    } else {
        // Element vector has not been assigned yet.
        throw exception("Attempt to add element before elements vector is assigned.");
    }
}

// Returns true if the species contains the element (given by index).
bool Species::ContainsElement(unsigned int i) const 
{
    // Loop over composition to find element.
    ElCompVector::const_iterator el;
    for (el=m_elcomp.begin(); el!=m_elcomp.end(); el++) {
        if (el->Index() == i) {
            return true;
        }
    }

    // We have arrived here because the element wasn't found.
    return false;
}

// Returns true if the species contains the element (given by reference).
bool Species::ContainsElement(const Sprog::Element &el) const
{
    if (m_elements != NULL) {
        // Loop over the elements to find the index.
        unsigned int i;
        for (i=0; i<m_elements->size(); i++) {
            if (*m_elements->at(i) == el) {
                return ContainsElement(i);
            }
        }

        // We have got here because the element isn't in the list, therefore
        // the species cannot contain it.
        return false;
    } else {
        // The elements vector is unassigned.
        throw exception("Elements vector is unassigned (Species::ContainsElement).");
    }
}

// Returns true if the species contains the element (given by name).
bool Species::ContainsElement(const std::string &name) const
{
    if (m_elements != NULL) {
        // Loop over the elements to find the index.
        unsigned int i;
        for (i=0; i<m_elements->size(); i++) {
            if (*m_elements->at(i) == name) {
                return ContainsElement(i);
            }
        }

        // We have got here because the element isn't in the list, therefore
        // the species cannot contain it.
        return false;
    } else {
        // The elements vector is unassigned.
        throw exception("Elements vector is unassigned (Species::ContainsElement).");
    }
}


// MOLECULAR WEIGHT.

// Recalculates the species molecular weight using the elements.
real Species::CalcMolWt()
{
    m_molwt = 0.0; // Reset.

    if (m_elements != NULL) {
        // Loop over composition vector, summing up the molecular weight.
        ElCompVector::const_iterator el;
        for (el=m_elcomp.begin(); el!=m_elcomp.end(); el++) {
            m_molwt += (*m_elements)[el->Index()]->MolWt() * (real)el->Count();
        }
        return m_molwt;
    } else {
        // Elements vector not assigned.
        throw exception("Cannot calculate molecular weight before assigning elements vector (Species::CalcMolWt).");
    }
}

// ELEMENTS VECTOR.

// Returns a pointer to the vector of elements used to define the species.
const ElementPtrVector *const Species::Elements()
{
    return m_elements;
}

// Sets the pointer to the vector of elements used to define the species.
void Species::SetElements(const Sprog::ElementPtrVector *const els)
{
    m_elements = els;
}


// THERMODYNAMIC FITTING PARAMETERS.

// Returns the thermo parameters which are valid for the given temperature.
const Thermo::THERMO_PARAMS &Species::ThermoParams(const Sprog::real T) const
{
    if (T >= m_T1) {
        // The thermo params are indexed by temperature.  The lower_bound() function
        // returns the first point at which the key (Temperature) is greater than
        // the supplied value.
        Thermo::ThermoMap::const_iterator i = m_thermoparams.lower_bound(T);

        if (i != m_thermoparams.end()) {
            // Success.
            return i->second;
        } else {
            // T is to large!
            throw out_of_range("Temperature above acceptable range. (Species::ThermoParams)");
        }
    } else {
        // The given temperature is lower than the thermo range start temperature.
        throw out_of_range("Temperature under acceptable range. (Species::ThermoParams)");
    }
}

// Adds a set of thermo parameters valid up to the given temperature to the
// species object.
void Species::AddThermoParams(const Sprog::real T, const Sprog::Thermo::THERMO_PARAMS &params)
{
    m_thermoparams[T] = params;
}

// Removes the thermo parameters associated with the given temperature, assuming
// it is in the list.
void Species::RemoveThermoParams(const Sprog::real T)
{
    m_thermoparams.erase(T);
}


// PARENT MECHANISM.

// Returns the pointer to the parent mechanism.
Sprog::Mechanism *const Species::Mechanism(void) const
{
    return m_mech;
}

// Sets the parent mechanism.
void Species::SetMechanism(Sprog::Mechanism *const mech)
{
    m_mech = mech;
}


// CLONING.

// Returns a pointer to a copy of the Species object.
Species *const Species::Clone(void) const
{
    return new Species(*this);
}