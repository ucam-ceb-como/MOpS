#include "mops_mixture.h"

using namespace Mops;

Mixture::Mixture(void)
{
}

Mixture::~Mixture(void)
{
}

void Mixture::Normalise()
{
    real xtot = 0.0;

    for (int i=0; i<m_species->size(); i++) {
        if (m_data[i] < 0.0) m_data[i] = 0.0;
        xtot += m_data[i];
    }

    if (xtot != 1.0) {
        for (int i=0; i<m_species->size(); i++) {
            m_data[i] /= xtot;
        }
    }
}