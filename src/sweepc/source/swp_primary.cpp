/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Primary class declared in the
    swp_primary.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#include "swp_primary.h"
#include "swp_model_factory.h"
#include "swp_cell.h"
#include "swp_kmc_typedef.h"

#include <stdexcept>
#include <memory.h>

#include <boost/random/uniform_smallint.hpp>

using namespace Sweep;
using namespace Sweep::KMC_ARS;

using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

//! Default constructor (protected).
AggModels::Primary::Primary(void)
: m_pmodel(NULL), m_createt(0.0), m_time(0.0), m_diam(0.0), m_dcol(0.0), 
m_dmob(0.0), m_surf(0.0), m_vol(0.0), m_mass(0.0), m_numcarbon(0), m_frag(0), m_numOf6Rings(0), m_phaseterm(0.0)
{
}

// Initialising constructor.
AggModels::Primary::Primary(double time, const Sweep::ParticleModel &model)
{
    init();
    m_pmodel  = &model;
    
    // Set particle create and update times.
    m_createt = time;
    m_time    = time;

    // Resize component and tracker vectors to match number of
    // each in the particle model.
    m_comp.resize(model.ComponentCount(), 0.0);
    m_values.resize(model.TrackerCount(), 0.0);
}

// Copy constructor.
AggModels::Primary::Primary(const Primary &copy)
{
    init();
    *this = copy;
}

// Stream-reading constructor.
AggModels::Primary::Primary(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
AggModels::Primary::~Primary()
{
    releaseMem();
}


// OPERATOR OVERLOADS.

//! Assignment operator.
AggModels::Primary &AggModels::Primary::operator=(const Primary &rhs)
{
    if (this != &rhs) {
        // Check if the RHS uses the same particle model before copying
        // components and tracker values.
        if (rhs.m_pmodel == m_pmodel) {
            // Copy composition.
            memcpy(&m_comp[0], &rhs.m_comp[0], sizeof(double)*m_comp.size());        
            // Copy tracker variables.
            memcpy(&m_values[0], &rhs.m_values[0], sizeof(double)*m_values.size());
        } else {
            // Set the particle model.
            m_pmodel = rhs.m_pmodel;
            // Copy components.
            m_comp.assign(rhs.m_comp.begin(), rhs.m_comp.end());
            // Copy tracker variables.
            m_values.assign(rhs.m_values.begin(), rhs.m_values.end());
        }

        // Copy primary time.
        m_createt = rhs.m_createt;
        m_time    = rhs.m_time;

        //! Copy derived and basic properties.
        m_diam = rhs.m_diam;
        m_dcol = rhs.m_dcol;
        m_dmob = rhs.m_dmob;
        m_surf = rhs.m_surf;
        m_vol  = rhs.m_vol;
        m_mass = rhs.m_mass;
        m_numcarbon = rhs.m_numcarbon;
		m_numOf6Rings = rhs.m_numOf6Rings;
        m_frag = rhs.m_frag;
		m_phaseterm = rhs.m_phaseterm;
    }
    return *this;
}

// DEFINING PARTICLE MODEL.

// Returns the particle model used to create this primary.
const Sweep::ParticleModel *const AggModels::Primary::ParticleModel(void) const
{
    return m_pmodel;
}


// PARTICLE COMPOSITION.

//Returns named component
double AggModels::Primary::GetComponent(std::string name) const
{
	int index = m_pmodel->ComponentIndex(name);
	return Composition(index);
}

// Returns the composition vector.
const fvector &AggModels::Primary::Composition() const
{
	return m_comp;
}

// Returns the ith component value.  Returns 0.0 if i is invalid.
double AggModels::Primary::Composition(unsigned int i) const
{
	if (i < m_comp.size()) {
		return m_comp[i];
	} else {
		return 0.0;
	}
}

// Sets the composition vector.
void AggModels::Primary::SetComposition(const Sweep::fvector &comp)
{
	m_comp.assign(comp.begin(), comp.end());
}


// TRACKER VARIABLE VALUES.

// Returns the tracker value vector.
const fvector &AggModels::Primary::Values() const
{
	return m_values;
}

// Returns the ith tracker variable value.  Returns 0.0 if i is invalid.
double AggModels::Primary::Values(unsigned int i) const
{
	if (i < m_values.size()) {
		return m_values[i];
	} else {
		return 0.0;
	}
}

// Sets the values vector.
void AggModels::Primary::SetValues(const fvector &vals)
{
    m_values.assign(vals.begin(), vals.end());
}

// Sets the ith trackervalue.
void AggModels::Primary::SetValue(unsigned int i, double val)
{
	if (i < m_values.size()) {
		m_values[i] = val;
	}
}

// PHASE VARIABLES

// Get the mass of a specific phase
double AggModels::Primary::GetPhaseMass(int i) const{

	// Get component indices for phase 
	std::vector<unsigned int> indices = m_pmodel->GetPhaseComponents(i);

	// Loop over composition and calculate mass and volume.
	double mass = 0.0;
	for (std::vector<unsigned int>::const_iterator it = indices.begin(); it != indices.end(); ++it){
		mass += m_pmodel->Components((*it))->MolWt() * m_comp[(*it)] / NA;
	}

	return mass;
}

// PRIMARY CREATE TIME.

// Returns the primary create time.
double AggModels::Primary::CreateTime() const {return m_createt;}


// LAST UPDATE TIME.

// Returns the primary last update time.
double AggModels::Primary::LastUpdateTime() const {return m_time;}

// Sets the last update time of the primary.
void AggModels::Primary::SetTime(double t) {m_time = t;}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType AggModels::Primary::AggID(void) const {return AggModels::Spherical_ID;}

// BASIC DERIVED PARTICLE PROPERTIES.

//! Calculates the derived properties from the unique properties.
void AggModels::Primary::UpdateCache(void)
{
    double m = 0.0;
    int m_numcarbon_temp = 0;

    // Loop over composition and calculate mass and volume.
    m_mass = m_vol = 0.0;
    for (unsigned int i=0; i!=m_pmodel->ComponentCount(); ++i) {
        m = m_pmodel->Components(i)->MolWt() * m_comp[i] / NA;
        m_mass += m;
        if (m_pmodel->Components(i)->Density() > 0.0)
            m_vol  += m / m_pmodel->Components(i)->Density();
    }

    /**
     * Get the number of carbon atoms.
     * It is assumed that the 0th index of m_comp contains the number of carbon
     * atoms.
     */
    m_numcarbon = 0;
    m_numcarbon = m_comp[0];

    /**
     * For particles with one carbon atom, set the number of carbon atoms to 0.
     * m_numcarbon is used to calculate the rate of fragmentation. Particles
     * with an m_numcarbon of 1 cannot further fragment; therefore, to exclude
     * these particles from the rate calculation m_numcarbon is set to 0.  
     */
    if(m_numcarbon < 2) {
        m_numcarbon = 0;
    }

    /**
     * Indicate that a particle with one carbon atom cannot fragment.
     * m_frag is used in the uniform selection of particle for fragmentation.
     * By default, a particle is able to fragment.
     */
    m_frag = 0;
    if(m_numcarbon >= 2) {
        m_frag = 1;
    }

    // Calculate other properties (of sphere).
    m_diam = pow(6.0 * m_vol / PI, ONE_THIRD);
    m_dcol = m_diam;
    m_dmob = m_diam;
    m_surf = PI * m_diam * m_diam;
    m_numcarbon = 0;
    for (unsigned int i=0; i!=m_pmodel->ComponentCount(); ++i) {
        m_numcarbon_temp = m_comp[i];
        m_numcarbon += m_numcarbon_temp;
    }

	// phaseterm for kinetic titania phase transformation
	// Components An and Ru assumed, 0 if components don't exist
	double An = AggModels::Primary::GetComponent("An");
	double Ru = AggModels::Primary::GetComponent("Ru");
	m_phaseterm = pow( An, TWO_THIRDS) * pow((An + Ru), ONE_THIRD);
}

// Returns the particle equivalent sphere diameter.
double AggModels::Primary::SphDiameter(void) const {return m_diam;}

// Returns the collision diameter.
double AggModels::Primary::CollDiameter(void) const {return m_dcol;}

// Rethrns the mobility diameter.
double AggModels::Primary::MobDiameter(void) const {return m_dmob;}

// Returns the surface area.
double AggModels::Primary::SurfaceArea(void) const {return m_surf;}

// Returns the equivalent sphere surface area, based
// on the volume.
double AggModels::Primary::SphSurfaceArea(void) const
{
    return PI * m_diam * m_diam;
}

// Returns the volume.
double AggModels::Primary::Volume(void) const {return m_vol;}

// Returns the mass.
double AggModels::Primary::Mass(void) const {return m_mass;}

//! Returns the number of carbon atoms.
int AggModels::Primary::NumCarbon(void) const {return m_numcarbon;}

//! Returns fragmentation flag.
int AggModels::Primary::Frag(void) const {return m_frag;}
int AggModels::Primary::NumRings(void) const { return m_numOf6Rings; }

//! Returns the property with the given ID.
double AggModels::Primary::Property(const Sweep::PropID id) const
{
    switch (id) {
        case iDsph:      // Equivalent sphere diameter.
            return m_diam;
        case iDcol:   // Collision diameter.
            return m_dcol;
        case iDmob:   // Mobility diameter.
            return m_dmob;
        case iS:      // Surface area.
            return m_surf;
        case iASN:      // Surface area.
            return GetSites();
        case iV:      // Volume.
            return m_vol;
        case iM:      // Mass.
            return m_mass;
        case iNumCarbon:        //! Number of carbon atoms.
            return m_numcarbon;
        case iFrag:             //! Fragmentation flag.
            return m_frag;
        default:
            return 0.0;
    }
}


// BASIC DERIVED PROPERTY OVERWRITES.

// Sets the spherical particle diameter
void AggModels::Primary::SetSphDiameter(double diam) {m_diam = diam;}

// Sets the collision diameter of the particle.
void AggModels::Primary::SetCollDiameter(double dcol) {m_dcol = dcol;}

// Sets the mobility diameter.
void AggModels::Primary::SetMobDiameter(double dmob) {m_dmob = dmob;}

// Sets the surface area, subject to minimum spherical area condition.
void AggModels::Primary::SetSurfaceArea(double surf) {m_surf = surf;}

// Sets the volume.
void AggModels::Primary::SetVolume(double vol) {m_vol = vol;}

// Sets the mass.
void AggModels::Primary::SetMass(double m) {m_mass = m;}

//! Sets the number of carbon atoms.
void AggModels::Primary::SetNumCarbon(int numcarbon) {m_numcarbon = numcarbon;}

//! Sets fragmentation flag.
void AggModels::Primary::SetFrag(int frag) {m_frag = frag;}
/*!
 * Check that this primary is a physically valid particle.  This currently
 * means checking that the amount of each component is within a range specifed
 * for that component.
 *
 *@return   True iff primary particle is valid as a physical particle
 */
bool AggModels::Primary::IsValid() const {
    for(unsigned int i = 0; i < m_comp.size(); ++i) {
        // Check each component value, but stop as soon as an invalid value
        // is found
        if(!m_pmodel->Components(i)->IsValidValue(m_comp[i]))
            return false;
    }

    // Check one of the cached values, all the m_??? could be added here
    // if necessary.
    return (m_mass > 0);
}

// OPERATIONS.

// Adjusts the primary with the given composition and 
// tracker values changes n times.  If the particle cannot be adjust
// n times, then this function returns the number of times
// it was adjusted.
unsigned int AggModels::Primary::Adjust(const fvector &dcomp, const fvector &dvalues, rng_type &rng, unsigned int n)
{
	unsigned int i = 0;

	n = CalculateMaxAdjustments(dcomp, n);

	// Add the components.
	for (i=0; i!=min(m_comp.size(),dcomp.size()); ++i) {
		m_comp[i] += dcomp[i] * (double)n;
	}

	// Add the tracker values.
	for (i=0; i!=min(m_values.size(),dvalues.size()); ++i) {
		m_values[i] += dvalues[i] * (double)n;
	}

    // Update property cache.
	if (n > 0) {
		AggModels::Primary::UpdateCache();
	}

    return n;
}

/*!
 * Calculates the maximum number of LPDA adjustments allowed.
 *
 * This is relevant for Processes which remove components from a Primary,
 * particularly where multiple components are present in a particle.
 *
 * @param dcomp     The component changes for the Process
 * @param n         Number of times requested
 * @return          Actual number of times allowed
 */
unsigned int AggModels::Primary::CalculateMaxAdjustments(
        const fvector &dcomp,
        unsigned int n) const
{
    unsigned int i;

    // First check if there is a negative in dcomp
    bool neg_change(false);
    for (i=0; i!=dcomp.size(); ++i) {
        if (dcomp[i] < 0.0) neg_change = true;
    }

    // Calculate the actual number of events only if there is a negative
    // change in the component amount
    if (neg_change) {
        for (i=0; i!=min(m_comp.size(),dcomp.size()); ++i) {
            if (dcomp[i] < 0.0)
                n = std::min(n, (unsigned int) abs(m_comp[i]/dcomp[i]));
        }
    }

    return n;
}

/*!
* Calculates the maximum possible number of adjustments
*
* This is relevant for Processes which remove components from a Primary,
* particularly where multiple components are present in a particle.
*
* @param dcomp     The component changes for the Process
* @return          Number of times allowed
*/
unsigned int AggModels::Primary::CalculateMaxAdjustments(const fvector &dcomp) const
{
	unsigned int i;
	unsigned int n = 0;

	// First check if there is a negative in dcomp
	bool neg_change(false);
	for (i = 0; i != dcomp.size(); ++i) {
		if (dcomp[i] < 0.0) {
			neg_change = true;
			n = std::numeric_limits<unsigned int>::max();
		}
	}

	// Calculate the actual number of events only if there is a negative
	// change in the component amount
	if (neg_change) {
		for (i = 0; i != min(m_comp.size(), dcomp.size()); ++i) {
			if (dcomp[i] < 0.0)
				n = std::min(n, (unsigned int)abs(m_comp[i] / dcomp[i]));
		}
	}

	return n;
}

// Adjusts the particle n times for IntParticle reaction
unsigned int AggModels::Primary::AdjustIntPar(const fvector &dcomp, const fvector &dvalues, rng_type &rng, unsigned int n)
{
    n = AggModels::Primary::Adjust(dcomp, dvalues, rng, n);

    return n;
}

//Adjusts the particle n times for the phase transformation process
unsigned int AggModels::Primary::AdjustPhase(const fvector &dcomp,
                              const fvector &dvalues,
                              rng_type &rng,
							  unsigned int n)
{
	n = AggModels::Primary::Adjust(dcomp, dvalues, rng, n);
	
    return n;
}

//Thermodynamic phase transformation (melting model)
void AggModels::Primary::Melt(rng_type &rng, Cell &sys)
{
	
	fvector dcomp(m_pmodel->ComponentCount(), 0.0);
	std::vector<fvector> all_dcomp;
	fvector dvalues(m_pmodel->TrackerCount(), 0.0); //for now the transformation model doesn't adjust the tracker variables
	unsigned int n = 0;

	//check if liquid
	bool Liquid = ParticleModel()->MeltModel().IsLiquid(sys, *this);

	//has the particle melted?
	if (Liquid == true){
		// melt the particle
		ParticleModel()->MeltModel().MeltingCompositionChange(all_dcomp);
	}
	else{ //particle is solid
		//total composition excluding liquid phase
		double mass = 0.0;
		fvector phaseMass(m_pmodel->PhaseCount(), 0.0);
		
		//add all components excluding liquid components
		for (unsigned int i = 0; i < m_pmodel->PhaseCount(); i++){
			if (!m_pmodel->PhaseIsLiquid(i)){
				phaseMass[i] = GetPhaseMass(i);
				mass += phaseMass[i];
			}
		}

		if (mass > 0.0){
			//if solid phases exist then convert liquid to solid phase probabilistically
			//weighted by mass of existing solid phases
			boost::uniform_01<rng_type&, double> uniformGenerator(rng);
			double j = uniformGenerator() * mass; //generate random number
			for (unsigned int i = 0; i < m_pmodel->PhaseCount(); i++){

				if (j <= phaseMass[i]){ 
					//get change in composition	
					ParticleModel()->MeltModel().CompositionChange(i, all_dcomp);
					break;
				}
				else{
					j -= phaseMass[i];
				}
			}
		}
		else{
			//if everything is liquid then tranform based on melting model
			ParticleModel()->MeltModel().CompositionChange(sys, *this, all_dcomp); //get change in composition	
		}
	}

	//loop over all composition changes and transform components
	for (std::vector<fvector>::const_iterator it = all_dcomp.begin(); it != all_dcomp.end(); ++it){

		dcomp = *it;

		// transform components
		n = Primary::CalculateMaxAdjustments(dcomp); //transform everything: get max adjustment
		if (n > 0){ n = Primary::Adjust(dcomp, dvalues, rng, n); }; //adjust primary
	}
}

// Property for titania phase transformation model
double AggModels::Primary::GetPhaseTerm(void) const
{
	return m_phaseterm;	
}

/*!
 *  Combines this primary with another.
 *
 *  Note the very strange behaviour when the primaries do not
 *  have equal particle model pointers.  Users should also
 *  be very careful about not mixing different types that inherit
 *  from PrimaryParticle.
 *
 * \param[in]       rhs         Primary particle to add to current instance
 * \param[in,out]   rng         Random number generator
 *
 * \return      Reference to the current instance after rhs has been added
 */
AggModels::Primary &AggModels::Primary::Coagulate(const Primary &rhs, rng_type &rng)
{
    // Check if the RHS uses the same particle model.  If not, then
    // just use the assignment operator because you can't add apples 
    // and bananas!
    if (rhs.m_pmodel == m_pmodel) {
        // Add the components.
        for (unsigned int i=0; i!=min(m_comp.size(),rhs.m_comp.size()); ++i) {
            m_comp[i] += rhs.m_comp[i];
        }

        // Add the tracker values.
        for (unsigned int i=0; i!=min(m_values.size(),rhs.m_values.size()); ++i) {
            m_values[i] += rhs.m_values[i];
        }

        // Create time is the earliest time.
        m_createt = min(m_createt, rhs.m_createt);
        m_time    = min(m_time, rhs.m_time);
    } else {
        // Different particle models!
        *this = rhs;
        std::cerr << "Sweep::AggModels::Primary::Coagulate called for particles with models " << m_pmodel
                  << " and " << rhs.m_pmodel << std::endl;
    }
    AggModels::Primary::UpdateCache();
    return *this;
}

/*!
 *  Combines this primary with another.
 *
 *  Note the very strange behaviour when the primaries do not
 *  have equal particle model pointers.  Users should also
 *  be very careful about not mixing different types that inherit
 *  from PrimaryParticle.
 *
 * \param[in]       rhs         Primary particle to add to current instance
 * \param[in,out]   rng         Random number generator
 *
 * \return      Reference to the current instance after rhs has been added
 */
AggModels::Primary &AggModels::Primary::Fragment(const Primary &rhs, rng_type &rng)
{
    // Check if the RHS uses the same particle model.  If not, then
    // just use the assignment operator because you can't add apples 
    // and bananas!
    if (rhs.m_pmodel == m_pmodel) {
        // Add the components.
        for (unsigned int i=0; i!=min(m_comp.size(),rhs.m_comp.size()); ++i) {
            m_comp[i] += rhs.m_comp[i];
        }

        // Add the tracker values.
        for (unsigned int i=0; i!=min(m_values.size(),rhs.m_values.size()); ++i) {
            m_values[i] += rhs.m_values[i];
        }

        // Create time is the earliest time.
        m_createt = min(m_createt, rhs.m_createt);
        m_time    = min(m_time, rhs.m_time);
    } else {
        // Different particle models!
        *this = rhs;
        std::cerr << "Sweep::AggModels::Primary::Fragment called for particles with models " << m_pmodel
                  << " and " << rhs.m_pmodel << std::endl;
    }
    AggModels::Primary::UpdateCache();
    return *this;
}

// This routine sinters the Primary for the given length of
// time using the provided sintering model.
void AggModels::Primary::Sinter(double dt, Cell &sys,
                     const Processes::SinteringModel &model,
                     rng_type &rng,
                     double wt)
{
    // Spherical primaries don't sinter.
	
    return;
}

/*!
 * @brief       Gets the ratio of second to first component number
 *
 * A coarse estimation of the coverage fraction; as the amount of the second
 * component relative to the amount of the first. Returns 1.0 for single comp.
 *
 * @return      Ratio of number of second to first component
 */
double AggModels::Primary::GetCoverageFraction() const
{
    double val = 1.0;

    if (m_pmodel->ComponentCount() > 1 && m_comp[0] > 0.0) {
        val = m_comp[1]/m_comp[0];
    }

    return val;
}

// READ/WRITE/COPY.

// Returns a copy of the model data.
AggModels::Primary *const AggModels::Primary::Clone(void) const
{
    return new Primary(*this);
}

//! Writes the object to a binary stream.
void AggModels::Primary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write number of components.
        unsigned int n = (unsigned int)m_comp.size();
        out.write((char*)&n, sizeof(n));

        // Write components.
        double val = 0.0;
        for (unsigned int i=0; i!=n; ++i) {
            val = (double)m_comp[i];
            out.write((char*)&val, sizeof(val));
        }

        // Write number of tracker values.
        n = (unsigned int)m_values.size();
        out.write((char*)&n, sizeof(n));

        // Write values.
        for (unsigned int i=0; i!=n; ++i) {
            val = (double)m_values[i];
            out.write((char*)&val, sizeof(val));
        }

        // Write create time.
        val = (double)m_createt;
        out.write((char*)&val, sizeof(val));

        // Write last update time.
        val = (double)m_time;
        out.write((char*)&val, sizeof(val));

        // Write equivalent sphere diameter.
        val = (double)m_diam;
        out.write((char*)&val, sizeof(val));

        // Write collision diameter.
        val = (double)m_dcol;
        out.write((char*)&val, sizeof(val));

        // Write mobility diameter.
        val = (double)m_dmob;
        out.write((char*)&val, sizeof(val));

        // Write surface area.
        val = (double)m_surf;
        out.write((char*)&val, sizeof(val));

        // Write volume.
        val = (double)m_vol;
        out.write((char*)&val, sizeof(val));

        // Write mass.
        val = (double)m_mass;
        out.write((char*)&val, sizeof(val));

        //! Write number of carbon atoms.
        val = (int)m_numcarbon;
        out.write((char*)&val, sizeof(val));

        //! Write fragmentation flag.
        val = (int)m_frag;
        out.write((char*)&val, sizeof(val));

		// Write phaseterm.
        val = (double)m_phaseterm;
        out.write((char*)&val, sizeof(val));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, AggModels::Primary::Serialize).");
    }
}

//! Reads the object from a binary stream.
void AggModels::Primary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    releaseMem();
    m_pmodel = &model;

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double val = 0.0;

        switch (version) {
            case 0:
                // Read number of components.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read components.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_comp.push_back((double)val);
                }

                // Read number of tracker values.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read values.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_values.push_back((double)val);
                }

                // Read create time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_createt = (double)val;

                // Read last update time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_time = (double)val;

                // Read equivalent sphere diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_diam = (double)val;

                // Read collision diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_dcol = (double)val;

                // Read mobility diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_dmob = (double)val;

                // Read surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_surf = (double)val;

                // Read volume.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_vol = (double)val;

                // Read mass.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_mass = (double)val;

                //! Read number of carbon atoms.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_numcarbon = (int)val;

                //! Read fragmentation flag.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_frag = (int)val;

				// Read phaseterm.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_phaseterm = (double)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, AggModels::Primary::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, AggModels::Primary::Deserialize).");
    }
}


// DATA MANAGEMENT.

// Release all memory associated with object.
void AggModels::Primary::releaseMem(void)
{
    m_comp.clear();
    m_values.clear();
}

//! Initialisation routine.
void AggModels::Primary::init(void)
{
    m_pmodel = NULL;
    m_createt = 0.0;
    m_time = 0.0;
    m_diam = 0.0; // Equivalent spherical diameter.
    m_dcol = 0.0; // Collision diameter.
    m_dmob = 0.0; // Mobility diameter.
    m_surf = 0.0; // Surface area.
    m_vol = 0.0;  // Volume.
    m_mass = 0.0; // Mass.
    m_numcarbon = 0;
	m_numOf6Rings = 0;
    m_frag = 0;      //!< Fragmentation flag.
	m_phaseterm = 0.0; // Phase term.
    releaseMem();
}

//! Check whether the number of carbon atoms in the primary is equal to that of
//! the inception species.
int AggModels::Primary::InceptedPAH() const
{
    std::vector<ParticleModel::PostProcessStartingStr> str_list = ParticleModel()->InceptedPAH();
	for (std::vector<int>::size_type ii = 0; ii!=str_list.size(); ii++){
		ParticleModel::PostProcessStartingStr str = str_list[ii];
		//A1, A1CH3, A2, A4, A4CH3, R5A3, A5
		switch (str){
			case ParticleModel::A1:
				if (NumCarbon() == BENZENE_C)
					return 1;
				break;
			case ParticleModel::A1CH3:
				if (NumCarbon() == TOLUENE_C)
					return 1;
				break;
			case ParticleModel::A2:
				if (NumCarbon() == NAPHTHALENE_C)
					return 1;
				break;
			case ParticleModel::A4:
				if (NumCarbon() == PYRENE_C)
					return 1;
				break;
			case ParticleModel::A4CH3:
				if (NumCarbon() == METHYLPYRENE_C)
					return 1;
				break;
			case ParticleModel::R5A3:
				if (NumCarbon() == MPHENANTHRENER_C)
					return 1;
				break;
			case ParticleModel::A5:
				if (NumCarbon() == BENZOPYRENE_C)
					return 1;
				break;
			default:
				return 0;
		}
	}
	return 0;
}

//! Returns the frame position and orientation, and primary coordinates
//	Used by bintree_primary for particle tracking
void AggModels::Primary::GetFrameCoords(std::vector<fvector> &coords) const
{
	//unused
}

void AggModels::Primary::setTracking()
{
	//spherical particles don't have primaries
}

void AggModels::Primary::removeTracking()
{
	//spherical particles don't have primaries
}