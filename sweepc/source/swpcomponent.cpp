#include "swpcomponent.h"

using namespace Sweep;

Component::Component()
{
    m_density = 0.0;
    m_molwt = 0.0;
}

Component::Component(const Sweep::real molwt, const Sweep::real dens, const std::string &name)
{
    Initialise(molwt, dens, name);
}

Component::~Component()
{
}

void Component::Initialise(const Sweep::real molwt, const Sweep::real dens, const std::string &name)
{
    m_density = dens;
    m_molwt = molwt;
    m_name = name;
}