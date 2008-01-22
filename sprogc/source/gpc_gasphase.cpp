#include "gpc_gasphase.h"
#include<vector>

using namespace Sprog;
using namespace Sprog::Thermo;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

GasPhase::GasPhase(void)
{
    /*
    m_precalc = false;
    m_Us_valid = false;
    m_Hs_valid = false;
    m_Ss_valid = false;
    m_Gs_valid = false;
    m_Cps_valid = false;
    m_Cvs_valid = false;
    m_Us_RT_valid = false;
    m_Hs_RT_valid = false;
    m_Ss_R_valid = false;
    m_Gs_RT_valid = false;
    m_Cps_R_valid = false;
    m_Cvs_R_valid = false;
    */
}

GasPhase::~GasPhase(void)
{
}


// THERMODYNAMIC PROPERTIES.

void GasPhase::Us(Sprog::fvector &U) const
{
    CalcUs(*m_pT, U);
}

void GasPhase::Hs(Sprog::fvector &H) const
{
    CalcHs(*m_pT, H);
}

void GasPhase::Ss(Sprog::fvector &S) const
{
    CalcSs(*m_pT, S);
}

void GasPhase::Gs(Sprog::fvector &G) const
{
    CalcGs(*m_pT, G);
}

void GasPhase::Cps(Sprog::fvector &Cp) const
{
    CalcCps(*m_pT, Cp);
}

void GasPhase::Cvs(Sprog::fvector &Cv) const
{
    CalcCvs(*m_pT, Cv);
}


// THERMODYNAMIC PROPERTIES (DIMENSIONLESS).

void GasPhase::Us_RT(Sprog::fvector &U) const
{
    CalcUs_RT(*m_pT, U);
}

void GasPhase::Hs_RT(Sprog::fvector &H) const
{
    CalcHs_RT(*m_pT, H);
}

void GasPhase::Ss_R(Sprog::fvector &S) const
{
    CalcSs_R(*m_pT, S);
}

void GasPhase::Gs_RT(Sprog::fvector &G) const
{
    CalcGs_RT(*m_pT, G);
}

void GasPhase::Cps_R(Sprog::fvector &Cp) const
{
    CalcCps_R(*m_pT, Cp);
}

void GasPhase::Cvs_R(Sprog::fvector &Cv) const
{
    CalcCvs_R(*m_pT, Cv);
}

// BULK MIXTURE PROPERTIES.

// Calculates the bulk internal energies in current units.
real GasPhase::BulkU() const
{
    return CalcBulkU(*m_pT, &m_data[0], m_species->size());
}

// Calculates the bulk enthalpy in current units.
real GasPhase::BulkH() const
{
    return CalcBulkH(*m_pT, &m_data[0], m_species->size());
}

// Calculates the bulk entropy in current units.
real GasPhase::BulkS() const
{
    return CalcBulkS(*m_pT, &m_data[0], m_species->size());
}

// Calculates the bulk Gibbs free energies in current units.
real GasPhase::BulkG() const
{
    return CalcBulkG(*m_pT, &m_data[0], m_species->size());
}

// Calculates the mean molar heat capacity at const. P.
real GasPhase::BulkCp() const
{
    return CalcBulkCp(*m_pT, &m_data[0], m_species->size());
}

// Calculates the mean molar heat capacity at const. V.
real GasPhase::BulkCv() const
{
    return CalcBulkCv(*m_pT, &m_data[0], m_species->size());
}

