#include "swpparticlechanger.h"

using namespace Sweep;

ParticleChanger::ParticleChanger(void)
{
    m_components = NULL;
    m_comp.clear();
    m_values.clear();
}

ParticleChanger::~ParticleChanger(void)
{
    m_comp.clear();
    m_values.clear();
}

void ParticleChanger::SetComponents(std::vector<Component*> &components)
{
    m_components = &components;
    m_comp.resize(m_components->size());
}

void ParticleChanger::SetCompChange(const unsigned int i, const real dcomp)
{
    m_comp.resize(i+1);
    m_comp[i] = dcomp;
}

void ParticleChanger::SetValueChange(const unsigned int i, const real value)
{
    if ((int)m_values.size() == 0) {
        m_values.resize(i+1);
    } else {
        m_values.reserve(i+1);
    }
    m_values[i] = value;
}
