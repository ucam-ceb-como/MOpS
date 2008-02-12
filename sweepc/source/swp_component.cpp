#include "swp_component.h"

using namespace Sweep;

// CONSRUCTORS AND DESTRUCTORS.

// Default constructor.
Component::Component()
: m_density(0.0), m_molwt(0.0), m_name("")
{
}

// Initialising constructor.
Component::Component(Sweep::real molwt, 
                     Sweep::real dens, 
                     const std::string &name)
{
    // Initialise the component properties.
    m_density = dens;
    m_molwt   = molwt;
    m_name    = name;
}

// Copy constructor.
Component::Component(const Component &copy) 
{
    *this = copy;
}

// Default destructor.
Component::~Component()
{
}


// OPERATOR OVERLOADS.
Component &Component::operator=(const Component &rhs)
{
    if (this != &rhs) {
        m_density = rhs.m_density;
        m_molwt   = rhs.m_molwt;
        m_name    = rhs.m_name;
    }
    return *this;
}