#include "gpc_idealgas.h"
#include "gpc_params.h"
#include <math.h>

using namespace Sprog;
using namespace Sprog::Thermo;
using namespace std;


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
IdealGas::IdealGas(void)
{
}

// Default destructor.
IdealGas::~IdealGas(void)
{
}


// EQUATION OF STATE FUNCTIONS.

// Calculates the mixture pressure.
real IdealGas::Pressure() const
{
    // P = rho R T.
    return m_dens * R * m_T;
}


// THERMODYNAMIC PROPERTY CALCULATION.

// Calculates enthalpies of all species.
void IdealGas::CalcHs(std::vector<real> &Hs) const
{
    int i, n;
    real t[Thermo::H_PARAM_COUNT];

    // Ensure output array has sufficient length.
    Hs.resize(m_species->size());

    // Precalculate temperature terms.
    t[0] = R * m_T;
    n = Thermo::H_PARAM_COUNT;
    for (i=1; i<n-1; i++) {
        t[i] = (real)i * t[i-1] * m_T / (real)(i+1);
    }
    t[n-1] = 1.0 / m_T;

    // Sum terms in polynomial for Hs.
    sumTerms(t, n, Hs);
}

// Calculates entropies of all species.
void IdealGas::CalcSs(std::vector<real> &Ss) const
{
    int i, n;
    real t[Thermo::S_PARAM_COUNT];

    // Ensure output array has sufficient length.
    Ss.resize(m_species->size());

    // Precalculate temperature terms.
    n = Thermo::S_PARAM_COUNT;
    t[0] = R * log(m_T);
    t[1] = R * m_T;
    for (i=2; i<n-2; i++) {
        t[i] = (real)(i-1) * t[i-1] * m_T / (real)i;
    }
    t[n-2] = 0.0;
    t[n-1] = R;

    // Sum terms in polynomial for Ss.
    sumTerms(t, n, Ss);
}

// Calculates molar heat capacity at const. P of all species.
void IdealGas::CalcCps(std::vector<real> &Cps) const
{
    int i, n;
    real t[Thermo::CP_PARAM_COUNT];

    // Ensure output array has sufficient length.
    Cps.resize(m_species->size());

    // Precalculate temperature terms.
    n = Thermo::CP_PARAM_COUNT;
    t[0] = R;
    for (i=1; i<n; i++) {
        t[i] = t[i-1] * m_T;
    }

    // Sum terms in polynomial for Cps.
    sumTerms(t, n, Cps);
}

// Calculates molar heat capacity at const. V of all species.
void IdealGas::CalcCvs(std::vector<real> &Cvs) const
{
    vector<real>::iterator i;

    // Calculate heat capacties at const. P.
    CalcCps(Cvs);

    // Convert to const. V.
    for(i=Cvs.begin(); i!=Cvs.end(); i++) {
        (*i) -= R;
    }
}

// Calculates molar internal energies of each species.
void IdealGas::CalcUs(std::vector<real> &Us) const
{
    // Calculate the enthalpies.
    CalcHs_RT(Us);

    // Convert to internal energies.
    vector<real>::iterator i;
    for (i=Us.begin(); i!=Us.end(); i++) {
        (*i) -= R * m_T;
    }
}

// Calculates molar Gibbs free energies of each species.
void IdealGas::CalcGs(std::vector<real> &Gs) const
{
    int i, n;
    real t[Thermo::H_PARAM_COUNT];

    // Ensure output array has sufficient length.
    Gs.resize(m_species->size());

    // Precalculate temperature terms.
    n = Thermo::H_PARAM_COUNT;
    t[1] = - R * m_T;
    t[0] = t[1] * (log(m_T) - 1.0);
    t[S_PARAM_COUNT-1] = t[1];
    t[1] *= 0.5 * m_T;
    for (i=2; i<Thermo::CP_PARAM_COUNT; i++) {
        t[i] = (real)(i-1) * t[i-1] * m_T / (real)(i+1);
    }
    t[H_PARAM_COUNT-1] = m_T;

    // Sum terms in polynomial for Cps.
    sumTerms(t, n, Gs);
}

// Calculates enthalpies of all species.
void IdealGas::CalcHs_RT(std::vector<real> &Hs) const
{
    int i, n;
    real t[Thermo::H_PARAM_COUNT];

    // Ensure output array has sufficient length.
    Hs.resize(m_species->size());

    // Precalculate temperature terms.
    t[0] = 1.0;
    n = Thermo::H_PARAM_COUNT;
    for (i=1; i<n-1; i++) {
        t[i] = (real)i * t[i-1] * m_T / (real)(i-1);
    }
    t[n-1] = 1.0 / m_T;

    // Sum terms in polynomial for Hs.
    sumTerms(t, n, Hs);
}


// Calculates entropies of all species.
void IdealGas::CalcSs_R(std::vector<real> &Ss) const
{
    int i, n;
    real t[Thermo::S_PARAM_COUNT];

    // Ensure output array has sufficient length.
    Ss.resize(m_species->size());

    // Precalculate temperature terms.
    n = Thermo::S_PARAM_COUNT;
    t[0] = log(m_T);
    t[1] = m_T;
    for (i=2; i<n-2; i++) {
        t[i] = (real)(i-1) * t[i-1] * m_T / (real)i;
    }
    t[n-2] = 0.0;
    t[n-1] = 1.0;

    // Sum terms in polynomial for Ss.
    sumTerms(t, n, Ss);
}

// Calculates molar heat capacity at const. P of all species.
void IdealGas::CalcCps_R(std::vector<real> &Cps) const
{
    int i, n;
    real t[Thermo::CP_PARAM_COUNT];

    // Ensure output array has sufficient length.
    Cps.resize(m_species->size());

    // Precalculate temperature terms.
    n = Thermo::CP_PARAM_COUNT;
    t[0] = 1.0;
    for (i=1; i<n; i++) {
        t[i] = t[i-1] * m_T;
    }

    // Sum terms in polynomial for Cps.
    sumTerms(t, n, Cps);
}

// Calculates molar heat capacity at const. V of all species.
void IdealGas::CalcCvs_R(std::vector<real> &Cvs) const
{
    vector<real>::iterator i;

    // Calculate heat capacties at const. P.
    CalcCps_R(Cvs);

    // Convert to const. V.
    for(i=Cvs.begin(); i!=Cvs.end(); i++) {
        (*i) -= 1.0;
    }
}

// Calculates molar internal energies of each species.
void IdealGas::CalcUs_RT(std::vector<real> &Us) const
{
    // Calculate the enthalpies.
    CalcHs_RT(Us);

    // Convert to internal energies.
    vector<real>::iterator i;
    for (i=Us.begin(); i!=Us.end(); i++) {
        (*i) -= 1.0;
    }
}

// Calculates molar Gibbs free energies of each species.
void IdealGas::CalcGs_RT(std::vector<real> &Gs) const
{
    int i, n;
    real t[Thermo::H_PARAM_COUNT];

    // Ensure output array has sufficient length.
    Gs.resize(m_species->size());

    // Precalculate temperature terms.
    n = Thermo::H_PARAM_COUNT;
    t[0] = 1.0 - log(m_T);
    t[1] = - 0.5 * m_T;
    for (i=2; i<Thermo::CP_PARAM_COUNT; i++) {
        t[i] = (real)(i-1) * t[i-1] * m_T / (real)(i+1);
    }
    t[H_PARAM_COUNT-1] = 1.0 / m_T;
    t[S_PARAM_COUNT-1] = - 1.0;

    // Sum terms in polynomial for Cps.
    sumTerms(t, n, Gs);
}


// BULK THERMODYNAMIC PROPERTIES.

// Calculates the bulk enthalpy in current units.
real IdealGas::BulkH(void) const
{
    // Get individual species enthalpies.
    vector<real> Hs;
    CalcHs(Hs);
    
    // Sum to get bulk enthalpy.
    real H = 0.0;
    for(int i=0; i<m_species->size(); i++) {
        H += m_frac[i] * Hs[i];
    }

    return H;
}

// Calculates the bulk entropy in current units.
real IdealGas::BulkS(void) const
{
    // Get individual species entropies.
    vector<real> Ss;
    CalcSs(Ss);
    
    // Sum to get bulk entropy.
    real S = 0.0;
    for(int i=0; i<m_species->size(); i++) {
        S += m_frac[i] * Ss[i];
    }

    return S;
}

// Calculates the mean molar heat capacity at const. P.
real IdealGas::BulkCp(void) const
{
    // Get individual species heat capacities.
    vector<real> Cps;
    CalcCps(Cps);
    
    // Sum to get bulk heat capacity.
    real Cp = 0.0;
    for(int i=0; i<m_species->size(); i++) {
        Cp += m_frac[i] * Cps[i];
    }

    return Cp;
}

// Calculates the mean molar heat capacity at const. V.
real IdealGas::BulkCv(void) const
{
    // Get individual species heat capacities.
    vector<real> Cvs;
    CalcCvs(Cvs);
    
    // Sum to get bulk heat capacity.
    real Cv = 0.0;
    for(int i=0; i<m_species->size(); i++) {
        Cv += m_frac[i] * Cvs[i];
    }

    return Cv;
}


// PRIVATE FUNCTIONS.

// Calculates a polynomial fit of any thermo property given the
// temperature terms.  The polynomial coefficients are found per
// species.
void IdealGas::sumTerms(real *t, int n, std::vector<real> &Xs) const
{
    int i, k;
    const THERMO_PARAMS *a;

    for (i=0; i<min(m_species->size(),Xs.size()); i++) {
        // Get the thermo fitting parameters for this species
        // at this temperature.
        a = &(*m_species)[i]->ThermoParams(m_T);

        // Calculate the thermo property for this species.
        for (k=0; k<n; k++) {
            Xs[i] += a->Params[k] * t[k];
        }
    }
}

