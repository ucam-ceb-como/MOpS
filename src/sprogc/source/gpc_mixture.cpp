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
#include "gpc_mech.h"
#include <vector>
#include <stdexcept>
#include <boost/serialization/vector.hpp>
//#include <bits/stl_vector.h>

using namespace Sprog;
using namespace Sprog::Thermo;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
Mixture::Mixture(void)
: m_vmodel(Sprog::iAir)
{
    m_data.clear();
    m_species = NULL;
}

// Default constructor (public, requires species list).
Mixture::Mixture(const SpeciesPtrVector &sp)
: m_vmodel(Sprog::iAir)
{
    SetSpecies(sp);
}

/*
Mixture::Mixture(const SpeciesPtrVector &sp, const int NumGasSp, const int NumSurfSp)
{
    SetSpecies(sp);
	gasSpeciesCount =  NumGasSp; 
	surfSpeciesCount = 	NumSurfSp; 
}
*/

// Copy constructor.
Mixture::Mixture(const Mixture &copy)
: m_vmodel(copy.m_vmodel)
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
		gasSpeciesCount = mix.gasSpeciesCount;
		surfSpeciesCount = mix.surfSpeciesCount; 
    }

    return *this;
}

// Returns the PAH formation rate. 
double Mixture::PAHFormationRate() const
{
    return m_data[PAHFormationIndex()];
}

// Set the PAH formation rate.
void Mixture::SetPAHFormationRate(double r)
{
    m_data[PAHFormationIndex()] = r;
}


// TEMPERATURE.

// Returns the mixture temperature.
double Mixture::Temperature() const
{
    return m_data[temperatureIndex()];
}


// Set the mixture temperature.
void Mixture::SetTemperature(double T)
{
    m_data[temperatureIndex()] = T;
}

// Set convective velocity
void Mixture::SetConvectiveVelocity(double u)
{
    m_data[ConvectiveVelocityIndex()] = u;
}

// Get convective velocity
double Mixture::GetConvectiveVelocity() const
{
    return m_data[ConvectiveVelocityIndex()];
}

// Set thermophoretic velocity
void Mixture::SetThermophoreticVelocity(double v)
{
    m_data[ThermophoreticVelocityIndex()] = v;
}

// Get thermophoretic velocity
double Mixture::GetThermophoreticVelocity() const
{
    return m_data[ThermophoreticVelocityIndex()];
}

// Set diffusion term
void Mixture::SetDiffusionTerm(double D)
{
    m_data[DiffusionTermIndex()] = D;
}

// Get diffusion term
double Mixture::GetDiffusionTerm() const
{
    return m_data[DiffusionTermIndex()];
}
// CONCENTRATIONS/FRACTIONS.

// Returns the vector of mixture species mole fractions. (NO NEED CHANGING)
const fvector &Mixture::MoleFractions() const
{
    return m_data;
}

// Returns a vector of species concentrations. (Modified) (C, Debugged)
void Mixture::GetConcs(fvector &concs) const
{
    // Resize output vector.
    const size_t numSpecies = m_species->size();
    concs.resize(numSpecies);

    // Loop over all mole fractions and convert to concentrations.
    const double density = m_data[densityIndex()];
	
    for (unsigned int i = 0; i != gasSpeciesCount; ++i) {
        concs[i] = m_data[i] * density;
    }

	// Concentration for surface species	
    for (unsigned int i =  gasSpeciesCount; i != m_species->size(); ++i) {
	
	string spName = (*m_species)[0]->Mechanism()->GetSpecies(i)->Name(); 
	string phName = (*m_species)[0]->Mechanism()->GetSpecies(spName)->PhaseName();
	double site_d =  (*m_species)[0]->Mechanism()->FindSiteDensity(phName);
	int sp_occ = (*m_species)[0]->Mechanism()->FindSiteOccup(spName);
	
        concs[i] = m_data[i] * site_d/sp_occ;
    }
}

// Returns a vector of species mass fractions. (Modified) (C, Not used in mops, sprogc)
void Mixture::GetMassFractions(fvector &fracs) const
{
    // Clear output vector.
    fracs.resize(m_species->size());

    // Loop over all mole fractions and convert to mass fractions:
    //   y = x * wt / sum(x*wt)
    double val, tot = 0.0;
    for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
        val = m_data[i] * (*m_species)[i]->MolWt();
        fracs[i] = val;
        tot += val;
    }
    tot = 1.0 / tot;
    for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
        fracs[i] *= tot;
    }
	
	
	
	std::vector<double> TOT; 
	
	for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
		val = 0.0;
		tot = 0.0;
		std::string phname =  (*m_species)[0]->Mechanism()->Phase(j)->Name();
		for (unsigned int i=gasSpeciesCount; i!=m_species->size(); ++i) {
			if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
			val = m_data[i] * (*m_species)[i]->MolWt();
			fracs[i] = val;
			tot += val;
			}
	
		}
	   TOT.push_back(tot); 
	}
	
	
	for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
		TOT[j] = 1.0/ TOT[j];
		std::string phname =  (*m_species)[0]->Mechanism()->Phase(j)->Name();
		for (unsigned int i=gasSpeciesCount; i!=m_species->size(); ++i) {
			if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
			fracs[i] *= TOT[j];
			}
		}
	   
	}
	
}

// Returns the mole fraction of species i. (UNMODIFIED) (C)
double Mixture::MoleFraction(unsigned int i) const
{
    if (i < m_species->size()) {
        return m_data[i];
    } else {
        return 0.0;
    }
}

// Returns the molar concentration of species i. (Modified by mm864, currently unused in mops-sprogc) (C)
double Mixture::MolarConc(unsigned int i) const
{
    if (i < gasSpeciesCount) {
		return m_data[i] * m_data[densityIndex()];
    } 
	else {
        return MolarSurfConc(i);
    }
}

// Returns the molar surface concentration of species i. (Modified by mm864) (C, currently unused in mops-sprogc)
double Mixture::MolarSurfConc(unsigned int i) const
{
    if ((i >= gasSpeciesCount) && (i < m_species->size()) ) {
      string phName = (*m_species)[i]->PhaseName();
      double siteDen = 0.0; 
      siteDen = (*m_species)[0]->Mechanism()->FindSiteDensity(phName);

      return m_data[i] * siteDen / ((*m_species)[i]->SiteOccupancy()); // m_data should contains site fraction !! 
    } else {
        return 0.0;
    }
}

// Returns the mass fraction of species i. (modified by mm864)  (C, currently unused in mops-sprogc)
double Mixture::MassFraction(unsigned int i) const
{
    if (i < gasSpeciesCount) {
        // Get total x * W.
        double tot = 0.0;
        for (unsigned int j=0; j!=gasSpeciesCount; ++j) {
            tot += m_data[j] * (*m_species)[j]->MolWt();
        }

        // Return mass fraction.
        return m_data[i] * (*m_species)[i]->MolWt() / tot;
    } 
	
	else if((i >= gasSpeciesCount) && (i < m_species->size()) ) {
	
	fvector massFrac;
	GetMassFractions(massFrac); 
	
	return massFrac[i];
	
	}
	
	else {
        return 0.0;
    }
}

// Sets the vector of species mole fractions. (modified by mm864) (C, debugged)
void Mixture::SetFracs(const fvector &fracs)
{
    double tot =0.0;

    // Set the mole fractions.
    for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
        m_data[i] = fracs[i];
        tot += fracs[i];
    }

    // Ensure that the mole fractions are normalised.
    if (tot != 1.0) {
        for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
            m_data[i] /= tot;
            
        }
    }
    
	
	// Surface species
	std::vector<double> Z; 

	for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
	double ztot = 0.0;
	std::string phname =  (*m_species)[0]->Mechanism()->Phase(j)->Name();
	cout << "phasename" << phname << endl;
		for (unsigned int i=gasSpeciesCount; i!=m_species->size(); ++i) {
			if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
			m_data[i] = fracs[i];
			ztot += fracs[i];
			}	
		}
	    Z.push_back(ztot); // NOTE Z will allocate Z = 0.0 for gas phase but it won't mess up the next stage below
	
	}

	// Ensure that the mole fractions are normalised. 
	for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
	
	if ((Z[j] != 1.0) && ( (*m_species)[0]->Mechanism()->Phase(j)->Name().compare("gas") != 0)) {
	
			std::string phname = (*m_species)[0]->Mechanism()->Phase(j)->Name();
	
			for (unsigned int i=gasSpeciesCount; i!=m_species->size(); ++i) {
				if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
				m_data[i] /= Z[j];
				}
			}
    }
	
    }

	
}

// Sets the species mole fractions from an array of values. (modified by mm864) (C, currently not used in mops-sprogc)
void Mixture::SetFracs(const double fracs[], int n)
{

	/*
    // Check that the array of of sufficient length.
    if ((unsigned)n >= gasSpeciesCount) {
        double tot = 0.0;

        // Set the fractions.
        for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
            m_data[i] = fracs[i];
            tot += fracs[i];
        }

        // Ensure that the mole fractions are normalised.
        if (tot != 1.0) {
            for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
                m_data[i] /= tot;
            }
        }
    }
	*/
		
	// Check that the array of of sufficient length.
    if ((unsigned)n >= m_species->size()) {
        
        // Set the fractions.
        
		// Surface and gas species
		std::vector<double> Z; 
		
		for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
		double fractot = 0.0;
		std::string phname =  (*m_species)[0]->Mechanism()->Phase(j)->Name();
		cout << "phasename" << phname << endl;
		for (unsigned int i=0; i!=m_species->size(); ++i) {
			if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
			m_data[i] = fracs[i];
			fractot += fracs[i];
			}	
		}
	    Z.push_back(fractot); 
	
	}
	
	// Ensure that the mole fractions are normalised.
		for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
	
		if ((Z[j] != 1.0)) {
	
			std::string phname = (*m_species)[0]->Mechanism()->Phase(j)->Name();
	
			for (unsigned int i=0; i!=m_species->size(); ++i) {
				if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
				m_data[i] /= Z[j];
				}
			}
		}
	
		}

		
		
    }
	
	
}

// Sets the species mole fractions using the supplied molar concentrations?? (modified by mm864 , currently not used in mops-sprogc)
void Mixture::SetConcs(const fvector &concs) // Is this correct?
{
    // Check that the concentration vector is of sufficient length.
    if (concs.size() >= m_species->size()) {
        // Sum up the total concentration.
		
		
        m_data[densityIndex()] = 0.0;
        for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
            m_data[i] = concs[i];
            m_data[densityIndex()] += concs[i];
        }

        // Convert values to mole fractions.
        double invdens = 1.0 / m_data[densityIndex()];
        for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
            m_data[i] *= invdens;
        }
		
	for (unsigned int i=gasSpeciesCount; i!=m_species->size(); ++i) {
		string phName = (*m_species)[i]->PhaseName();
		double siteDen = 0.0; 
		siteDen = (*m_species)[0]->Mechanism()->FindSiteDensity(phName);
		m_data[i] = concs [i] * ((*m_species)[i]->SiteOccupancy()) / siteDen; // m_data should contains site fraction !! 
	}
		

    }
}

// Sets the species mole fractions using the supplied mass fractions. (modified by mm864, currently not in used for mops-sprogc) (?) Can you have mass frac in surface
void Mixture::SetMassFracs(const fvector &fracs)
{
    // Check that the mass fraction vector is of sufficient length.
    if (fracs.size() >= m_species->size()) {
        double val = 0.0, tot = 0.0, totfrac = 0.0;
	
        // Check that the fractions are normalised by summing up
        // the total fractions, and dividing the values by this
        // sum in the next step.
        //std::cout << "Mass fractions:";
        for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
            totfrac += fracs[i];
            //std::cout << ' ' << fracs[i];
        }
        //std::cout << " tot " << totfrac << std::endl;

        // Convert to mole fractions:
        //   x = y / (wt * sum(y/wt))
        //std::cout << "Mole fractions (1):";
        for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
            val = fracs[i] / (totfrac * (*m_species)[i]->MolWt());
            m_data[i] = val;
            //std::cout << ' ' << val;
            tot += val;
        }
        //std::cout << " tot " << tot << std::endl;

        tot = 1.0 / tot;
        //std::cout << "Mole fractions (2):";
        for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
            m_data[i] *= tot;
            //std::cout << ' ' << m_data[i];
        }
        //std::cout << " tot " << 1.0 / tot << std::endl;

		// Surface Mole Fraction from mass fraction 
		
		
		std::vector<double> TOT; 
	
	for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
		val = 0.0, tot = 0.0, totfrac = 0.0;
		
		std::string phname =  (*m_species)[0]->Mechanism()->Phase(j)->Name();
		
		for (unsigned int i=gasSpeciesCount; i!=m_species->size(); ++i) {
			if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
			totfrac += fracs[i];
			val = fracs[i] / (totfrac * (*m_species)[i]->MolWt());
            m_data[i] = val;
			tot += val;
			}
	
		}
	   TOT.push_back(tot); 
	}
	
	
	for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
		TOT[j] = 1.0/ TOT[j];
		std::string phname =  (*m_species)[0]->Mechanism()->Phase(j)->Name();
		for (unsigned int i=gasSpeciesCount; i!=m_species->size(); ++i) {
			if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
			m_data[i] *= TOT[j];
			}
		}
	   
	}
	
		
    }
}

// Checks the vector of mole fractions for validity by settings all negative
// values to zero, and by normalising the values so that they sum
// to one. (Modified by mm864) (C, Debugged)
void Mixture::Normalise()
{
    double xtot = 0.0;

    for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
        if (m_data[i] < 0.0) m_data[i] = 0.0;
        xtot += m_data[i];
    }

    if (xtot != 1.0) {
        for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
            m_data[i] /= xtot;
        }
    }
	
	std::vector<double> Z; 
	
	for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
	double ztot = 0.0;
	std::string phname =  (*m_species)[0]->Mechanism()->Phase(j)->Name();
	
	for (unsigned int i=gasSpeciesCount; i!=m_species->size(); ++i) {
	if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
        if (m_data[i] < 0.0) m_data[i] = 0.0;
        ztot += m_data[i];
    }
	
	}
	Z.push_back(ztot); // NOTE Z will allocate Z = 0.0 for gas phase but it won't mess up the next stage below
	
	}
	
	for (unsigned int j = 0; j!=(*m_species)[0]->Mechanism()->PhaseCount(); ++j){
	
	if ((Z[j] != 1.0) && ( (*m_species)[0]->Mechanism()->Phase(j)->Name().compare("gas") != 0)) {
	
        std::string phname = (*m_species)[0]->Mechanism()->Phase(j)->Name();
	
		for (unsigned int i=gasSpeciesCount; i!=m_species->size(); ++i) {
		if (((*m_species)[i]->PhaseName()).compare(phname) == 0 ){
        m_data[i] /= Z[j];
        }
		}
    }
	
    }

	
}


// MIXTURE DENSITY. (UPDATED)

// Returns the mixture molar density. (C)
double Mixture::Density() const
{
    return m_data[densityIndex()];
}

// Returns the mixture mass density. (include gas phase only) (C, currently not in used for mops-sprogc)
double Mixture::MassDensity() const
{
    double rho = 0.0;

    // Calcualate mass density:
    //   rho_mass = rho_mole * sum(x * wt)
    for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
        rho += m_data[i] * (*m_species)[i]->MolWt();
    }
    rho *= m_data[densityIndex()];

    return rho;
}

// Sets the mixture molar density. (C, debugged)
void Mixture::SetDensity(double dens)
{
    m_data[densityIndex()] = dens;
}

// Sets the molar density using the supplied mass density. (THIS SHOULD BE: SET MOLAR DENSITY) (C, currently not in used for mops-sprogc)
void Mixture::SetMassDensity(double dens) // (include gas phase only)
{
    double sum = 0.0;

    // Calcualate molar density:
    //   rho_mass = rho_mole * sum(x * wt)
    for (unsigned int i=0; i!=gasSpeciesCount; ++i) {
        sum += m_data[i] * (*m_species)[i]->MolWt();
    }
   m_data[densityIndex()] = dens / sum;
}


// MIXTURE CONTEXT. (NO NEED CHANGING)

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

	gasSpeciesCount =0;


	for (unsigned int j = 0; j != m_species->size(); j++){
		if (((*m_species)[j]->PhaseName()).compare("gas") == 0){
	
			gasSpeciesCount++; 
		}
	}

	surfSpeciesCount = m_species->size() - gasSpeciesCount;  
}



// RAW DATA. (NO NEED CHANGING)

double *const Mixture::RawData()
{
    return &(m_data[0]);
}

const double *const Mixture::RawData() const
{
    return &(m_data[0]);
}

// READ/WRITE/COPY FUNCTIONS. (C)

// Creates a copy of the mixture object.
Mixture *const Mixture::Clone() const
{
    return new Mixture(*this);
}

// Identifies the mixture type for serialisation.
Serial_MixtureType Mixture::SerialType() const
{
    return Serial_Mixture;
}

// returns the avg mol wt given the mass fractions added by Vinod (modified by mm864, not in used for mops-sprogc)
double Mixture::getAvgMolWt(Sprog::fvector &massFrac) const{
	double avgMolWt = 0.0;
	for(unsigned int i=0; i!= gasSpeciesCount; i++)
		avgMolWt += massFrac[i]/(*m_species)[i]->MolWt();

	return 1.0/avgMolWt;
}

double Mixture::getAvgMolWt() const { // (modified by mm864, not in used for mops-sprogc)
    double avgMolWt = 0.0;
    vector<double> moleFrac = MoleFractions();
    for(unsigned int i=0; i!= gasSpeciesCount; i++) 
        avgMolWt += moleFrac[i]*(*m_species)[i]->MolWt();

    return avgMolWt;

}

/*!
 * Mean collision cross-section.  This is used, for example,
 * in mean path calculations, although for historical reasons
 * there may be some hard coded instances of related values in
 * Mopssuite.
 *
 * @return Mean collision cross-sectional area in \f$m^2\f$.
 */
double Mixture::getMeanCollisionSection() const { // Restricted to gas phase (mm864)
    double avgColDiam2 = 0.0;
    for (unsigned i = 0; i!= gasSpeciesCount; ++i) {
        avgColDiam2 += MoleFraction(i) * (*m_species)[i]->getCollisionDiameter() * (*m_species)[i]->getCollisionDiameter();
    }
    return avgColDiam2 * Sprog::PI;
}

// Following transport related routines are added by Vinod
// returns the mixture viscosity in Kg/m-s
double Mixture::getViscosity() const{

	Sprog::Transport::MixtureTransport mt;
	return mt.getViscosity(Temperature(),*this);
}


/*!
 * This prevents segmentation faults when accessing the Champan Enksog viscosity
 * model from Sweep.
 *
 * @param mix   Mixture for which to check integrity of transport data
 *
 * @exception   std::runtime_error      No transport data found for species
 */
void Mixture::checkForTransportData() const
{
    const SpeciesPtrVector *spv = Species();
    std::cout << "Checking for transport data." << std::endl;
    // Loop over all species and check they have transport data
    for (size_t k = 0; k != spv->size(); ++k) {
        if (!(*spv)[k]->hasTransportData()) {
            throw runtime_error("No transport data supplied for species " +
                    (*spv)[k]->Name() + ".");
        }
    }
}

// returns the mixture thermal conductivity in J/m-s-K
double Mixture::getThermalConductivity(double pre) const{
	Sprog::Transport::MixtureTransport mt;
	return mt.getThermalConductivity(Temperature(),pre,*this);
}

const vector<double> Mixture::getMolarSpecificHeat(){

        vector<double> cpMols;
	Sprog::Thermo::IdealGas ig(*this->Species());
	ig.CalcCps(Temperature(),cpMols);
        return cpMols;
}

const vector<double> Mixture::getMolarEnthalpy(double T){
	vector<double> enthalpy;
	Sprog::Thermo::IdealGas ig(*this->Species());
	ig.CalcHs(T,enthalpy);
	return enthalpy;

}

const vector<double> Mixture::getMolarEnthalpy(){
    vector<double> enthalpy;
    Sprog::Thermo::IdealGas ig(*this->Species());
    ig.CalcHs(Temperature(),enthalpy);
    return enthalpy;
}

// returns the mixture specific heat capacity in J/Kg K (modified by mm864)
double Mixture::getSpecificHeatCapacity(double T){
	double cp = 0.0;
	vector<double> cpMols, massFrac;
	Sprog::Thermo::IdealGas ig(*this->Species());
	ig.CalcCps(T,cpMols);
	GetMassFractions(massFrac);
	for(unsigned int i=0; i !=gasSpeciesCount; i++) // OVER GAS PHASE ONLY
		cp += massFrac[i]*(cpMols[i])/(*m_species)[i]->MolWt();

	return cp;
}

//return the specific heat capacity for the given mixture in J/kg K (modified by mm864)
double Mixture::getSpecificHeatCapacity(){
	double cp = 0.0;
	vector<double> cpMols, massFrac;
	Sprog::Thermo::IdealGas ig(*this->Species());
	ig.CalcCps(Temperature(),cpMols);
	GetMassFractions(massFrac);
	for(unsigned int i=0; i != gasSpeciesCount; i++) // OVER GAS PHASE ONLY
		cp += massFrac[i]*(cpMols[i])/(*m_species)[i]->MolWt();

	return cp;
}

// returns the vector of mixture diffusion coefficient in m^2/s
const vector<double> Mixture::getMixtureDiffusionCoeff(const double pre) const{
	Sprog::Transport::MixtureTransport mt;
	return mt.getMixtureDiffusionCoeff(Temperature(),pre,*this);
}

// Writes the mixture to a binary data stream.
void Mixture::Serialize(std::ostream &out) const
{

    /*unsigned int size = m_data.size();

    boost::archive::text_oarchive oa(out);
    oa << size << m_data;*/

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

        // Write the viscosity model
        sz = m_vmodel;
        out.write((char*)&sz, sizeof(sz));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sprog, Mixture::Serialize).");
    }
}

// Reads the mixture data from a binary data stream.
void Mixture::Deserialize(std::istream &in)
{

    /*unsigned int size = m_data.size();

    boost::archive::text_iarchive oa(in);
    oa >> size >> m_data;*/

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
                double val;
                m_data.reserve(sz);
                for (unsigned int i=0; i<sz; i++) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_data.push_back(val);
                }

                // The mixture has no species associated it with right now.
                m_species = NULL;

                // Read the viscosity model
                in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
                m_vmodel = (Sprog::ViscosityModel)sz;

                break;
            default:
                throw runtime_error("Mixture serialized version number "
                                    "is invalid (Sprog, Mixture::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, Mixture::Deserialize).");
    }
}
