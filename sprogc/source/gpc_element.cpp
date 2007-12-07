#include "gpc_element.h"
#include "gpc_mech.h"
#include <stdexcept>

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Element::Element(void)
{
    m_name = "";
    m_molwt = 0.0;
    m_mech = NULL;
}

// Copy constructor.
Element::Element(const Element &e)
{
    m_name = e.m_name;
    m_molwt = e.m_molwt;
    m_mech = e.m_mech;
}

// Initialising constructor.
Element::Element(const std::string &name, const Sprog::real molwt)
{
    m_name = name;
    m_molwt = molwt;
    m_mech = NULL;
}

// Destructor.
Element::~Element(void)
{
    m_name.clear();
}


// OPERATOR OVERLOADING.

// Assignment operator.
Element &Element::operator=(const Sprog::Element &el)
{
    if (this != &el) {
        m_name = el.m_name;
        m_molwt = el.m_molwt;
        m_mech = el.m_mech;
    }
    return *this;
}

// Comparison operator:  Compares two elements.  Returns true
// if the element names are the same.
bool Element::operator==(const Sprog::Element &el) const
{
    return m_name.compare(el.m_name)==0;
}

// Comparison operator: Compares an element to a string (name).  Returns
// true if the element name matches the string.
bool Element::operator==(const std::string &name) const
{
    return m_name.compare(name)==0;
}

// Inequality operator:  Compares two elements.  Returns false if the
// element names are equal.
bool Element::operator!=(const Sprog::Element &el) const
{
    return !(*this==el);
}

// Inequality operator:  Compares an element to a string (name).  Returns
// false of the element name matches the string.
bool Element::operator !=(const std::string &name) const
{
    return !(*this==name);
}


// NAME.

// Sets the element name.
void Element::SetName(const std::string &name) 
{
    // If the element is part of a mechanism then we must check
    // that the mechanism does not already contain an element with
    // that name.
    if (m_mech != NULL) {
        if (m_mech->FindElement(name) >= 0) {
            // Oh dear:  mechanism already contains an element with
            // this name.
            throw invalid_argument("Cannot have two elements with the same name (Element::SetName).");
        }
    }

    // Element is not in a mech or name is not defined in mech.  Can set name here.
    m_name = name;
}


// MOLECULAR WEIGHT.

// Sets the element molecular weight.
void Element::SetMolWt(const real molwt)
{
    if (molwt > 0.0) {
        m_molwt = molwt;

        // Tell mechanism that element mol. wt. has changed:  It probably
        // wants to update its species.
        if (m_mech != NULL) {
            m_mech->CheckElementChanges(*this);
        }
    } else {
        // Attempting to set zero or negative mol. wt.  This 
        // is unphysical.
        throw out_of_range("Molecular weight must be positive and non-zero! (Element::SetMolWt)");
    }
}

// Sets the element molecular weight by search the list of known elements.
bool Element::SetMolWtFromLibrary()
{
    int i;
    for (i=0; i<m_nlib; i++) {
        if (m_name.compare(m_lib[i].Name()) == 0) {
            // We have found a matching element in the library.
            m_molwt = m_lib[i].MolWt();
            return true;
        }
    }

    // We have not found an element in the library.
    return false;
}


// PARENT MECHANISM.

// Returns a pointer to the parent mechanism.
Sprog::Mechanism *const Element::Mechanism()
{
    return m_mech;
}

// Sets the parent mechanism.
void Element::SetMechanism(Sprog::Mechanism *const mech)
{
    m_mech = mech;
}


// CLONING.

// Returns a pointer to a copy of the Element object.
Element *const Element::Clone(void) const
{
    return new Element(*this);
}


