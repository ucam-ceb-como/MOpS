#include "gpc_params.h"
#include "gpc_mixture.h"
#include <vector>

using namespace Sprog;
using namespace Sprog::Thermo;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Mixture::Mixture(void)
{
    m_frac.clear();
    m_T = 0.0;
    m_dens = 0.0;
    m_species = NULL;
}

// Default desstructor.
Mixture::~Mixture(void)
{
    m_frac.clear();
}


// OPERATOR OVERLOADING.

// Assignment operator.
Mixture &Mixture::operator=(const Mixture &mix)
{
    // Check for self assignment.
    if (this != &mech) {
        m_frac.assign(mix.m_frac.begin(), mix.m_frac.end());
        m_T = mix.m_T;
        m_dens = mix.m_dens;
        m_species = mix.m_species;
    }

    return *this;
}

// CLONING.

Mixture* Mixture::Clone() const
{
    return NULL;
}


// TEMPERATURE.

// Returns the mixture temperature.
real Mixture::T() const
{
    return m_T;
}


// Set the mixture temperature.
void Mixture::SetT(Sprog::real T)
{
    m_T = T;
}


// CONCENTRATIONS/FRACTIONS.

// Returns the vector of mixture species mole fractions.
const vector<real> &Mixture::MoleFractions() const
{
    return m_frac;
}

// Returns a vector of species concentrations.
void Mixture::GetConcs(std::vector<real> &concs) const
{
    // Clear output vector.
    concs.clear();

    // Loop over all mole fractions and convert to concentrations.
    vector<real>::const_iterator i;
    for (i=m_frac.begin(); i!=m_frac.end(); i++) {
        concs.push_back((*i) * m_dens);
    }
}

// Returns a vector of species mass fractions.
void Mixture::GetMassFractions(std::vector<real> &fracs) const
{
    // Clear output vector.
    fracs.clear();

    // Loop over all mole fractions and convert to mass fractions:
    //   y = x * wt / sum(x*wt)
    int i;
    real val, tot = 0.0;
    for (i=0; i<m_frac.size(); i++) {
        val = m_frac[i] * (*m_species)[i]->MolWt();
        fracs.push_back(val);
        tot += val;
    }
    tot = 1.0 / tot;
    for (i=0; i<m_frac.size(); i++) {
        fracs[i] *= tot;
    }
}

// Returns the mole fraction of species i.
real Mixture::MoleFraction(unsigned int i) const
{
    if (i < m_species->size()) {
        return m_frac[i];
    } else {
        return 0.0;
    }
}

// Returns the molar concentration of species i.
real Mixture::MolarConc(unsigned int i) const
{
    if (i < m_species->size()) {
        return m_frac[i] * m_dens;
    } else {
        return 0.0;
    }
}

// Returns the mass fraction of species i.
real Mixture::MassFraction(unsigned int i) const
{
    if (i < m_species->size()) {
        // Get total x * W.
        real tot = 0.0;
        for (int j=0; j<m_frac.size(); j++) {
            tot += m_frac[j] * (*m_species)[j]->MolWt();
        }

        // Return mass fraction.
        return m_frac[i] * (*m_species)[i]->MolWt() / tot;
    } else {
        return 0.0;
    }
}

// Sets the vector of species mole fractions.
void Mixture::SetFracs(const std::vector<real> &fracs)
{
    m_frac.assign(fracs.begin(), fracs.end());
}

// Sets the species mole fractions using the supplied molar concentrations.
void Mixture::SetConcs(const std::vector<real> &concs)
{
    // Check that the concentration vector is of sufficient length.
    if (concs.size() >= m_species->size()) {
        int i;

        // Sum up the total concentration.
        m_dens = 0.0;
        for (i=0; i<m_species->size(); i++) {
            m_frac[i] = concs[i];
            m_dens += concs[i];
        }

        // Convert values to mole fractions.
        real invdens = 1.0 / m_dens;
        for (i=0; i<m_species->size(); i++) {
            m_frac[i] *= invdens;
        }
    }
}

// Sets the species mole fractions using the supplied mass fractions.
void Mixture::SetMassFracs(const std::vector<real> &fracs)
{
    // Check that the mass fraction vector is of sufficient length.
    if (fracs.size() >= m_species->size()) {
        int i;
        real val = 0.0, tot = 0.0;

        // Convert to mole fractions:
        //   x = y / (wt * sum(y/wt))
        for (i=0; i<m_species->size(); i++) {
            val = fracs[i] / (*m_species)[i]->MolWt();
            m_frac[i] = val;
            tot += val;
        }
        tot = 1.0 / tot;
        for (i=0; i<m_species->size(); i++) {
            m_frac[i] *= tot;
        }
    }
}


// MIXTURE DENSITY.

// Returns the mixture molar density.
real Mixture::Density() const
{
    return m_dens;
}

// Returns the mixture mass density.
real Mixture::MassDensity() const
{
    int i;
    real rho = 0.0;

    // Calcualate mass density:
    //   rho_mass = rho_mole * sum(x * wt)
    for (i=0; i<m_frac.size(); i++) {
        rho += m_frac[i] * (*m_species)[i]->MolWt();
    }
    rho *= m_dens;
    return rho;
}

// Sets the mixture molar density.
void Mixture::SetDensity(Sprog::real dens)
{
    m_dens = dens;
}

// Sets the molar density using the supplied mass density.
void Mixture::SetMassDensity(Sprog::real dens)
{
    int i;
    m_dens = 0.0;
    
    // Calcualate molar density:
    //   rho_mass = rho_mole * sum(x * wt)
    for (i=0; i<m_frac.size(); i++) {
        m_dens += m_frac[i] * (*m_species)[i]->MolWt();
    }
    m_dens = dens / m_dens;
}


// MIXTURE CONTEXT.

// Returns a pointer to the vector of species used to define the mixture.
const SpeciesPtrVector *const Mixture::Species() const
{
    return m_species;
}

// Sets the vector of species used to define the mixture.
void Mixture::SetSpecies(const Sprog::SpeciesPtrVector *const sp)
{
    m_species = sp;
}