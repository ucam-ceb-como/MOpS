#include "gpc_element.h"

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Element::Element(void)
{
    m_name = "";
    m_molwt = 0.0;
}

// Copy constructor.
Element::Element(const Element &e)
{
    m_name = e.m_name;
    m_molwt = e.m_molwt;
}

// Initialising constructor.
Element::Element(const std::string &name, const Sprog::real molwt)
{
    m_name = name;
    m_molwt = molwt;
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