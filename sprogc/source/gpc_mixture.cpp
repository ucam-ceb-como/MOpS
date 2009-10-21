/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Mixture class declared in the
    gpc_mixture.h header file.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#include "gpc_params.h"
#include "gpc_mixture.h"
#include "gpc_idealgas.h"
#include "gpc_transport_factory.h"
#include <vector>
#include <stdexcept>
//#include <bits/stl_vector.h>

using namespace Sprog;
using namespace Sprog::Thermo;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
Mixture::Mixture(void)
{
    m_data.clear();
    m_species = NULL;
}

// Default constructor (public, requires species list).
Mixture::Mixture(const SpeciesPtrVector &sp)
{
    SetSpecies(sp);
}

// Copy constructor.
Mixture::Mixture(const Mixture &copy)
{
    *this = copy;
}


// Stream-reading constructor.
Mixture::Mixture(std::istream &in, const SpeciesPtrVector &sp)
{
    Deserialize(in);
    SetSpecies(sp);
}

// Default destructor.
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
    }

    return *this;
}


// TEMPERATURE.

// Returns the mixture temperature.
real Mixture::Temperature() const
{
    return m_data[temperatureIndex()];
}


// Set the mixture temperature.
void Mixture::SetTemperature(Sprog::real T)
{
    m_data[temperatureIndex()] = T;
}


// CONCENTRATIONS/FRACTIONS.

// Returns the vector of mixture species mole fractions.
const fvector &Mixture::MoleFractions() const
{
    return m_data;
}

// Returns a vector of species concentrations.
void Mixture::GetConcs(fvector &concs) const
{
    // Resize output vector.
    const size_t numSpecies = m_species->size();
    concs.resize(numSpecies);

    // Loop over all mole fractions and convert to concentrations.
    const real density = m_data[densityIndex()];
    for (unsigned int i = 0; i != numSpecies; ++i) {
        concs[i] = m_data[i] * density;
    }
}

// Returns a vector of species mass fractions.
void Mixture::GetMassFractions(fvector &fracs) const
{
    // Clear output vector.
    fracs.resize(m_species->size());

    // Loop over all mole fractions and convert to mass fractions:
    //   y = x * wt / sum(x*wt)
    real val, tot = 0.0;
    for (unsigned int i=0; i!=m_species->size(); ++i) {
        val = m_data[i] * (*m_species)[i]->MolWt();
        fracs[i] = val;
        tot += val;
    }
    tot = 1.0 / tot;
    for (unsigned int i=0; i!=m_species->size(); ++i) {
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
        return m_data[i] * m_data[densityIndex()];
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
        for (unsigned int j=0; j!=m_species->size(); ++j) {
            tot += m_data[j] * (*m_species)[j]->MolWt();
        }

        // Return mass fraction.
        return m_data[i] * (*m_species)[i]->MolWt() / tot;
    } else {
        return 0.0;
    }
}

// Sets the vector of species mole fractions.
void Mixture::SetFracs(const fvector &fracs)
{
    real tot =0.0;

    // Set the mole fractions.
    for (unsigned int i=0; i!=m_species->size(); ++i) {
        m_data[i] = fracs[i];
        tot += fracs[i];
    }

    // Ensure that the mole fractions are normalised.
    if (tot != 1.0) {
        for (unsigned int i=0; i!=m_species->size(); ++i) {
            m_data[i] /= tot;
        }
    }
}

// Sets the species mole fractions from an array of values.
void Mixture::SetFracs(const Sprog::real fracs[], int n)
{
    // Check that the array of of sufficient length.
    if ((unsigned)n >= m_species->size()) {
        real tot = 0.0;

        // Set the fractions.
        for (unsigned int i=0; i!=m_species->size(); ++i) {
            m_data[i] = fracs[i];
            tot += fracs[i];
        }

        // Ensure that the mole fractions are normalised.
        if (tot != 1.0) {
            for (unsigned int i=0; i!=m_species->size(); ++i) {
                m_data[i] /= tot;
            }
        }
    }
}

// Sets the species mole fractions using the supplied molar concentrations.
void Mixture::SetConcs(const fvector &concs)
{
    // Check that the concentration vector is of sufficient length.
    if (concs.size() >= m_species->size()) {
        // Sum up the total concentration.
        m_data[densityIndex()] = 0.0;
        for (unsigned int i=0; i!=m_species->size(); ++i) {
            m_data[i] = concs[i];
            m_data[densityIndex()] += concs[i];
        }

        // Convert values to mole fractions.
        real invdens = 1.0 / m_data[densityIndex()];
        for (unsigned int i=0; i!=m_species->size(); ++i) {
            m_data[i] *= invdens;
        }
    }
}

// Sets the species mole fractions using the supplied mass fractions.
void Mixture::SetMassFracs(const fvector &fracs)
{
    // Check that the mass fraction vector is of sufficient length.
    if (fracs.size() >= m_species->size()) {
        real val = 0.0, tot = 0.0, totfrac = 0.0;

        // Check that the fractions are normalised by summing up
        // the total fractions, and dividing the values by this
        // sum in the next step.
        for (unsigned int i=0; i!=m_species->size(); ++i) {
            totfrac += fracs[i];
        }

        // Convert to mole fractions:
        //   x = y / (wt * sum(y/wt))
        for (unsigned int i=0; i!=m_species->size(); ++i) {
            val = fracs[i] / (totfrac * (*m_species)[i]->MolWt());
            m_data[i] = val;
            tot += val;
        }
        tot = 1.0 / tot;
        for (unsigned int i=0; i!=m_species->size(); ++i) {
            m_data[i] *= tot;
        }

    }
}

// Checks the vector of mole fractions for validity by settings all negative
// values to zero, and by normalising the values so that they sum
// to one.
void Mixture::Normalise()
{
    real xtot = 0.0;

    for (unsigned int i=0; i!=m_species->size(); ++i) {
        if (m_data[i] < 0.0) m_data[i] = 0.0;
        xtot += m_data[i];
    }

    if (xtot != 1.0) {
        for (unsigned int i=0; i!=m_species->size(); ++i) {
            m_data[i] /= xtot;
        }
    }
}


// MIXTURE DENSITY.

// Returns the mixture molar density.
real Mixture::Density() const
{
    return m_data[densityIndex()];
}

// Returns the mixture mass density.
real Mixture::MassDensity() const
{
    real rho = 0.0;

    // Calcualate mass density:
    //   rho_mass = rho_mole * sum(x * wt)
    for (unsigned int i=0; i!=m_species->size(); ++i) {
        rho += m_data[i] * (*m_species)[i]->MolWt();
    }
    rho *= m_data[densityIndex()];
    return rho;
}

// Sets the mixture molar density.
void Mixture::SetDensity(Sprog::real dens)
{
    m_data[densityIndex()] = dens;
}

// Sets the molar density using the supplied mass density.
void Mixture::SetMassDensity(Sprog::real dens)
{
    m_data[densityIndex()] = 0.0;
    
    // Calcualate molar density:
    //   rho_mass = rho_mole * sum(x * wt)
    for (unsigned int i=0; i!=m_species->size(); ++i) {
        m_data[densityIndex()] += m_data[i] * (*m_species)[i]->MolWt();
    }
   m_data[densityIndex()] = dens / m_data[densityIndex()];
}


// MIXTURE CONTEXT.

// Returns a pointer to the vector of species used to define the mixture.
const SpeciesPtrVector *const Mixture::Species() const
{
    return m_species;
}

// Sets the vector of species used to define the mixture.
void Mixture::SetSpecies(const Sprog::SpeciesPtrVector &sp)
{
    m_species = &sp;
    m_data.resize(m_species->size()+sNumNonSpeciesData);
}


// RAW DATA.

real *const Mixture::RawData()
{
    return &(m_data[0]);
}

const real *const Mixture::RawData() const
{
    return &(m_data[0]);
}

// READ/WRITE/COPY FUNCTIONS.

// Creates a copy of the mixture object.
Mixture *const Mixture::Clone() const
{
    return new Mixture(*this);
}

// Writes the mixture to a binary data stream.
void Mixture::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output the data vector size.
        unsigned int sz = m_data.size();
        out.write((char*)&sz, sizeof(sz));

        // Output all elements in the data vector.
        fvector::const_iterator i;
        for (i=m_data.begin(); i!=m_data.end(); i++) {
            out.write((char*)&(*i), sizeof(*i));
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sprog, Mixture::Serialize).");
    }
}

// Reads the mixture data from a binary data stream.
void Mixture::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                // Read the data vector size.
                unsigned int sz;
                in.read(reinterpret_cast<char*>(&sz), sizeof(sz));

                // Fill the data vector.
                real val;
                m_data.reserve(sz);
                for (unsigned int i=0; i<sz; i++) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_data.push_back(val);
                }

                // The mixture has no species associated it with right now.
                m_species = NULL;
                break;
            default:
                throw runtime_error("Mixture serialized version number "
                                    "is invalid (Sprog, Mixture::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, Mixture::Deserialize).");
    }
}

// Identifies the mixture type for serialisation.
Serial_MixtureType Mixture::SerialType() const
{
    return Serial_Mixture;
}

// returns the avg mol wt given the mass fractions added by Vinod
real Mixture::getAvgMolWt(Sprog::fvector &massFrac){
	real avgMolWt = 0.0;
	for(unsigned int i=0; i!= m_species->size(); i++)
		avgMolWt += massFrac[i]/(*m_species)[i]->MolWt();

	return 1.0/avgMolWt;
}

real Mixture::getAvgMolWt(){
    real avgMolWt = 0.0;
    vector<real> moleFrac = MoleFractions();
    for(unsigned int i=0; i!= m_species->size(); i++)
        avgMolWt += moleFrac[i]*(*m_species)[i]->MolWt();
    
    return avgMolWt;

}

// Following transport related routines are added by Vinod
// returns the mixture viscosity in Kg/m-s
real Mixture::getViscosity() const{

	Sprog::Transport::MixtureTransport mt;
	return mt.getViscosity(Temperature(),*this);
}
// returns the mixture thermal conductivity in J/m-s-K
real Mixture::getThermalConductivity(real pre) const{
	Sprog::Transport::MixtureTransport mt;
	return mt.getThermalConductivity(Temperature(),pre,*this);
}

const vector<real> Mixture::getMolarSpecificHeat(){

        vector<real> cpMols;
	Sprog::Thermo::IdealGas ig(*this->Species());
	ig.CalcCps(Temperature(),cpMols);
        return cpMols;
}

const vector<real> Mixture::getMolarEnthalpy(real T){
	vector<real> enthalpy;
	Sprog::Thermo::IdealGas ig(*this->Species());
	ig.CalcHs(T,enthalpy);
	return enthalpy;

}

const vector<real> Mixture::getMolarEnthalpy(){
    vector<real> enthalpy;
    Sprog::Thermo::IdealGas ig(*this->Species());
    ig.CalcHs(Temperature(),enthalpy);
    return enthalpy;
}

// returns the mixture specific heat capacity in J/Kg K
real Mixture::getSpecificHeatCapacity(Sprog::real T){
	real cp = 0.0;
	vector<real> cpMols, massFrac;
	Sprog::Thermo::IdealGas ig(*this->Species());
	ig.CalcCps(T,cpMols);
	GetMassFractions(massFrac);
	for(unsigned int i=0; i != m_species->size(); i++)
		cp += massFrac[i]*(cpMols[i])/(*m_species)[i]->MolWt();

	return cp;
}

//return the specific heat capacity for the given mixture in J/kg K
real Mixture::getSpecificHeatCapacity(){
	real cp = 0.0;
	vector<real> cpMols, massFrac;
	Sprog::Thermo::IdealGas ig(*this->Species());
	ig.CalcCps(Temperature(),cpMols);
	GetMassFractions(massFrac);
	for(unsigned int i=0; i != m_species->size(); i++)
		cp += massFrac[i]*(cpMols[i])/(*m_species)[i]->MolWt();

	return cp;
}

// returns the vector of mixture diffusion coefficient in m^2/s
const vector<real> Mixture::getMixtureDiffusionCoeff(const real pre) const{
	Sprog::Transport::MixtureTransport mt;
	return mt.getMixtureDiffusionCoeff(Temperature(),pre,*this);
}

