#include "swphomogas.h"

using namespace Sweep;

HomoGas::HomoGas(void)
{
    m_chem.clear();
    m_t = 0.0;
    m_p = 0.0;
    m_initchem.clear();
    m_initT = 0.0;
    m_initP = 0.0;
    m_inittime = 0.0;
}

HomoGas::~HomoGas(void)
{
    m_chem.clear();
    m_initchem.clear();
}

void HomoGas::SetInitConditions(const Sweep::real time, const std::vector<Sweep::real> &chem, 
                                const Sweep::real T, const Sweep::real P)
{
    m_inittime = time;
    m_initchem.assign(chem.begin(), chem.end());
    m_initT = T;
    m_initP = P;
    m_chem.assign(chem.begin(), chem.end());
    m_t = T;
    m_p = P;
}
