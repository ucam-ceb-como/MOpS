#include "swpcomponent.h"

using namespace Sweep;

Component::Component()
{
    // Give default values to component properties.
    m_density = 0.0;
    m_molwt = 0.0;
    m_name = "";
}

Component::Component(const Sweep::real molwt, const Sweep::real dens, const std::string &name)
{
    // Call the initialisation routine.
    Initialise(molwt, dens, name);
}

Component::~Component()
{
}

void Component::Initialise(const Sweep::real molwt, const Sweep::real dens, const std::string &name)
{
    // Initialise the component properties.
    m_density = dens;
    m_molwt = molwt;
    m_name = name;
}