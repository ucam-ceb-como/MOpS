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
    m_data.clear();
    m_pT = NULL;
    m_pdens = NULL;
    m_species = NULL;
}

// Default desstructor.
Mixture::~Mixture(void)
{
    m_data.clear();
}


// OPERATOR OVERLOADING.

// Assignment operator.
Mixture &Mixture::operator=(const Mixture &mix)
{
    // Check for self assignment.
    if (this != &mix) {
        m_data.assign(mix.m_data.begin(), mix.m_data.end());
        m_species = mix.m_species;
        m_pT = &m_data.at(m_species->size()-2);
        m_pdens = &m_data.at(m_species->size()-1);
    }

    return *this;
}

// CLONING.

Mixture* Mixture::Clone() const
{
    return new Mixture(*this);
}


// TEMPERATURE.

// Returns the mixture temperature.
real Mixture::Temperature() const
{
    return *m_pT;
}


// Set the mixture temperature.
void Mixture::SetTemperature(Sprog::real T)
{
    *m_pT = T;
}


// CONCENTRATIONS/FRACTIONS.

// Returns the vector of mixture species mole fractions.
const vector<real> &Mixture::MoleFractions() const
{
    return m_data;
}

// Returns a vector of species concentrations.
void Mixture::GetConcs(std::vector<real> &concs) const
{
    // Resize output vector.
    concs.reserve(m_species->size());

    // Loop over all mole fractions and convert to concentrations.
    vector<real>::const_iterator i;
    for (i=m_data.begin(); i!=m_data.end()-2; i++) {
        concs.push_back((*i) * (*m_pdens));
    }
}

// Returns a vector of species mass fractions.
void Mixture::GetMassFractions(std::vector<real> &fracs) const
{
    // Clear output vector.
    fracs.reserve(m_species->size());

    // Loop over all mole fractions and convert to mass fractions:
    //   y = x * wt / sum(x*wt)
    int i;
    real val, tot = 0.0;
    for (i=0; i<m_species->size(); i++) {
        val = m_data[i] * (*m_species)[i]->MolWt();
        fracs.push_back(val);
        tot += val;
    }
    tot = 1.0 / tot;
    for (i=0; i<m_species->size(); i++) {
        fracs[i] *= tot;
    }
}

// Returns the mole fraction of species i.
real Mixture::MoleFraction(unsigned int i) const
{
    if (i < m_species->size()) {
        return m_data[i];
    } else {
        return 0.0;
    }
}

// Returns the molar concentration of species i.
real Mixture::MolarConc(unsigned int i) const
{
    if (i < m_species->size()) {
        return m_data[i] * (*m_pdens);
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
        for (int j=0; j<m_species->size(); j++) {
            tot += m_data[j] * (*m_species)[j]->MolWt();
        }

        // Return mass fraction.
        return m_data[i] * (*m_species)[i]->MolWt() / tot;
    } else {
        return 0.0;
    }
}

// Sets the vector of species mole fractions.
void Mixture::SetFracs(const std::vector<real> &fracs)
{
    int i;
    real tot =0.0;

    // Set the mole fractions.
    for (i=0; i<m_species->size(); i++) {
        m_data[i] = fracs[i];
        tot += fracs[i];
    }

    // Ensure that the mole fractions are normalised.
    if (tot != 1.0) {
        for (i=0; i<m_species->size(); i++) {
            m_data[i] /= tot;
        }
    }
}

// Sets the species mole fractions from an array of values.
void Mixture::SetFracs(const Sprog::real fracs[], int n)
{
    // Check that the array of of sufficient length.
    if (n >= m_species->size()) {
        int i;
        real tot = 0.0;

        // Set the fractions.
        for (i=0; i<m_species->size(); i++) {
            m_data[i] = fracs[i];
            tot += fracs[i];
        }

        // Ensure that the mole fractions are normalised.
        if (tot != 1.0) {
            for (i=0; i<m_species->size(); i++) {
                m_data[i] /= tot;
            }
        }
    }
}

// Sets the species mole fractions using the supplied molar concentrations.
void Mixture::SetConcs(const std::vector<real> &concs)
{
    // Check that the concentration vector is of sufficient length.
    if (concs.size() >= m_species->size()) {
        int i;

        // Sum up the total concentration.
        *m_pdens = 0.0;
        for (i=0; i<m_species->size(); i++) {
            m_data[i] = concs[i];
            *m_pdens += concs[i];
        }

        // Convert values to mole fractions.
        real invdens = 1.0 / (*m_pdens);
        for (i=0; i<m_species->size(); i++) {
            m_data[i] *= invdens;
        }
    }
}

// Sets the species mole fractions using the supplied mass fractions.
void Mixture::SetMassFracs(const std::vector<real> &fracs)
{
    // Check that the mass fraction vector is of sufficient length.
    if (fracs.size() >= m_species->size()) {
        int i;
        real val = 0.0, tot = 0.0, totfrac = 0.0;

        // Check that the fractions are normalised by summing up
        // the total fractions, and dividing the values by this
        // sum in the next step.
        for (i=0; i<m_species->size(); i++) {
            totfrac += fracs[i];
        }

        // Convert to mole fractions:
        //   x = y / (wt * sum(y/wt))
        for (i=0; i<m_species->size(); i++) {
            val = fracs[i] / (totfrac * (*m_species)[i]->MolWt());
            m_data[i] = val;
            tot += val;
        }
        tot = 1.0 / tot;
        for (i=0; i<m_species->size(); i++) {
            m_data[i] *= tot;
        }

    }
}


// MIXTURE DENSITY.

// Returns the mixture molar density.
real Mixture::Density() const
{
    return *m_pdens;
}

// Returns the mixture mass density.
real Mixture::MassDensity() const
{
    int i;
    real rho = 0.0;

    // Calcualate mass density:
    //   rho_mass = rho_mole * sum(x * wt)
    for (i=0; i<m_species->size(); i++) {
        rho += m_data[i] * (*m_species)[i]->MolWt();
    }
    rho *= (*m_pdens);
    return rho;
}

// Sets the mixture molar density.
void Mixture::SetDensity(Sprog::real dens)
{
    *m_pdens = dens;
}

// Sets the molar density using the supplied mass density.
void Mixture::SetMassDensity(Sprog::real dens)
{
    int i;
    *m_pdens = 0.0;
    
    // Calcualate molar density:
    //   rho_mass = rho_mole * sum(x * wt)
    for (i=0; i<m_species->size(); i++) {
        *m_pdens += m_data[i] * (*m_species)[i]->MolWt();
    }
   *m_pdens = dens / (*m_pdens);
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
    m_data.resize(m_species->size()+2);
    m_pT = &m_data.at(m_species->size());
    m_pdens = &m_data.at(m_species->size()+1);
}


// RAW DATA.

real *const Mixture::RawData()
{
    return &(m_data[0]);
}
