#include "gpc_el_comp.h"
#include "gpc_element.h"

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ElComp::ElComp(void) 
{
    m_element = NULL;
    m_count = 0;
}

// Copy constructor.
ElComp::ElComp(const ElComp &el)
{
    *this = el;
}

// Destructor.
ElComp::~ElComp(void)
{
    m_element = NULL;
}


// OPERATOR OVERLOADING.

// Assignment operator.
ElComp &ElComp::operator=(const Sprog::ElComp &el)
{
    if (this!=&el) {
        m_element = el.m_element;
        m_count = el.m_count;
    }
    return *this;
}

// Increment operator:  Adds n to the element count.
ElComp &ElComp::operator+=(const unsigned int n)
{
    m_count += n;
    return *this;
}

// Decrement operator:  Removes n from the element count to a zero minimum.
ElComp &ElComp::operator-=(const unsigned int n)
{
    if (n > m_count) {
        m_count = 0;
    } else {
        m_count -= n;
    }
    return *this;
}

// Addition operator:  Returns a new ElComp object whose element count
// is this object's count + n.
const ElComp ElComp::operator+(const unsigned int n)
{
    ElComp el(*this);
    return el += n;
}

// Subtraction operator:  Returns a new ElComp object whose element count
// is this object's count - n to a minimum of zero.
const ElComp ElComp::operator-(const unsigned int n)
{
    ElComp el(*this);
    return el -= n;
}

// Comparison operator:  Returns true if both ElComp objects point to
// the same element.
bool ElComp::operator==(const ElComp &el) const
{
    return (m_element == el.m_element);
}

// Inequality operator:  Returns false if both ElComp objects point to
// the same element.
bool ElComp::operator!=(const ElComp &el) const
{
    return !(*this==el);
}


// ELEMENT POINTER.

// Returns a pointer to the element referred to by this ElComp object.
const Sprog::Element *const ElComp::Element(void) const
{
    return m_element;
}

// Sets the pointer to the element referred to by this ElComp object.
void Sprog::ElComp::SetElement(const Sprog::Element *const el)
{
    m_element = el;
}


// ELEMENT COUNT.

// Returns the element count for this composition.
unsigned int Sprog::ElComp::Count() const
{
    return m_count;
}

// Sets the element count for this composition.
void Sprog::ElComp::SetCount(unsigned int n)
{
    m_count = n;
}