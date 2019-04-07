/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  Merged with the SubParticle class by Robert Patterson, October 2012

  File purpose:
    Implementation of the Particle class declared in the
    swp_particle.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

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

#include "swp_particle.h"

#include "swp_particle_model.h"
#include "swp_particle_image.h"
#include "string_functions.h"

#include <cmath>
#include <stdexcept>
#include <cassert>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Particle::Particle(void)
: m_Position(0.0)
, m_PositionTime(0.0)
, m_StatWeight(1.0)
, m_primary(NULL)
, m_CoagCount(0)
, m_FragCount(0)
, m_createt(0.0)
, mLPDAtime(0.0)
{
}

// Initialising constructor.
Particle::Particle(double time, const Sweep::ParticleModel &model)
: m_Position(0.0)
, m_PositionTime(0.0)
, m_StatWeight(1.0)
, m_CoagCount(0)
, m_FragCount(0)
, m_createt(0.0)
, mLPDAtime(0.0)
{
    m_primary = new AggModels::Primary(time, model);
}

/*!
 * @param[in]	time	Create time
 * @param[in]	weight	Statistical weight of new particle
 * @param[in]	model	Particle model that interprets contents of particle
 */
Particle::Particle(double time, double weight, const Sweep::ParticleModel &model)
: m_Position(0.0)
, m_PositionTime(0.0)
, m_StatWeight(weight)
, m_CoagCount(0)
, m_FragCount(0)
, m_createt(0.0)
, mLPDAtime(0.0)
{
    m_primary = new AggModels::Primary(time, model);
}

// Initialising constructor (from Primary particle).
Particle::Particle(Sweep::AggModels::Primary &pri)
: m_Position(0.0)
, m_PositionTime(0.0)
, m_StatWeight(1.0)
, m_primary(&pri)
, m_CoagCount(0)
, m_FragCount(0)
{
    m_createt = pri.CreateTime();
    mLPDAtime = pri.CreateTime();
}

// Copy constructor.
Particle::Particle(const Sweep::Particle &copy)
: m_primary(NULL)
{
    // Use assignment operator.
    *this = copy;
}

/*!
 * @brief Read the particle from the binary stream. 
 *
 * This method must read the data in the same way in which it is written in Serialize.
 *
 * @param[in,out]	 in		             Input binary stream
 * @param[in]        model	             Particle model defining interpretation of particle data
 * @param[in,out]    duplicates          Addresses of PAHs for use when reading primary particles
 *
 * @exception		 invalid_argument    Stream not ready
 */
Particle::Particle(std::istream &in, const Sweep::ParticleModel &model, void *duplicates)
{
    if(in.good()) {
        m_primary = ModelFactory::ReadPrimary(in, model, duplicates);

        in.read(reinterpret_cast<char*>(&m_Position), sizeof(m_Position));
        in.read(reinterpret_cast<char*>(&m_PositionTime), sizeof(m_PositionTime));
        in.read(reinterpret_cast<char*>(&m_StatWeight), sizeof(m_StatWeight));
        in.read(reinterpret_cast<char*>(&m_CoagCount), sizeof(m_CoagCount));
        in.read(reinterpret_cast<char*>(&m_FragCount), sizeof(m_FragCount));
        in.read(reinterpret_cast<char*>(&m_createt), sizeof(m_createt));
        in.read(reinterpret_cast<char*>(&mLPDAtime), sizeof(mLPDAtime));
    }
    else {
        throw std::invalid_argument("Input stream not ready \
        (Sweep::Particle stream reading constructor");
    }
}

// Default destructor.
Particle::~Particle()
{
    // deleting a null pointer is not a problem - it just does nothing
    delete m_primary;
}

/*!
 * @param[in]		xml			XML node specifying the particle
 * @param[in]		model		Particle model that defines the interpretation of the particle data
 *
 * @return		Pointer to new particle constructed on the heap (caller must delete).
 *
 * @exception	runtime_error	Unrecognised component
 * @exception	runtime_error	Unrecognised tracker
 * @exception	runtime_error	Non-positive statistical weight
 */
Particle* Particle::createFromXMLNode(const CamXML::Element& xml, const Sweep::ParticleModel& model)
{
    // Read initial particle composition.
    vector<CamXML::Element*> subitems; 
    xml.GetChildren("component", subitems);
    fvector components(model.ComponentCount(), 0);
    
    for (vector<CamXML::Element*>::iterator j=subitems.begin(); j!=subitems.end(); ++j) {
        // Get component ID.
        string str = (*j)->GetAttributeValue("id");
        int id = model.ComponentIndex(str);

        if (id >= 0) {
            // Get component value (XML uses dx to match format for inception).
            str = (*j)->GetAttributeValue("dx");
            components[id] = Strings::cdble(str);
        } else {
            // Unknown component in mechanism.
            throw std::runtime_error(str + ": Component not found in mechanism \
                                       (Sweep, Particle::createFromXMLNode).");
        }
    }

    // Read initial tracker variable values.
    xml.GetChildren("track", subitems);
    fvector trackers(model.TrackerCount(), 0);
    
    for (vector<CamXML::Element*>::iterator j=subitems.begin(); j!=subitems.end(); j++) {
        // Get tracker ID.
        string str = (*j)->GetAttributeValue("id");
        int id = model.GetTrackerIndex(str);

        if (id >= 0) {
            // Get tracker value (XML uses dx to match format for inception).
            str = (*j)->GetAttributeValue("dx");
            trackers[id] = Strings::cdble(str);
        } else {
            // Unknown tracker variable in mechanism.
            throw std::runtime_error(str + ": Tracker variable not found in mechanism. \
                                       (Sweep, Particle::createFromXMLNode).");
        }
    }

    // Pointer to new particle
    // \TODO wrap in an auto_ptr for exception safety
    Particle* pNew = model.CreateParticle(0.0);

    // Read any statistical weight
    const CamXML::Element* const pWeightNode = xml.GetFirstChild("weight");
    if(pWeightNode != NULL) {
    	std::string str = pWeightNode->Data();
    	const double wt = Strings::cdble(str);
    	if(wt <= 0) {
    		throw std::runtime_error("Particle statistical weight must be >0, not " + str +
    				                 "(Sweep, Particle::createFromXMLNode).");
    	}
    	else
            pNew->setStatisticalWeight(wt);
    }
    
    // Read SurfVol specific parameters
    if (model.AggModel() == Sweep::AggModels::SurfVol_ID) {
        const CamXML::Element* const pSurfNode = xml.GetFirstChild("surf");
        if (pSurfNode != NULL) {
            std::string str = pSurfNode->Data();
            const double surf = Strings::cdble(str);
            if(surf <= 0) {
                throw std::runtime_error("Particle surface area must be >0, not " + str +
                                         "(Sweep, Particle::createFromXMLNode).");
            }
            else
                pNew->m_primary->SetSurfaceArea(surf);
        }
    }

    // Initialise the new particle.
    pNew->Primary()->SetComposition(components);
    pNew->Primary()->SetValues(trackers);
    pNew->UpdateCache();
    
    return pNew;
}


// Compound assignment (coagulation).
// OPERATOR OVERLOADING.

// Assignment operator.
Particle &Particle::operator=(const Sweep::Particle &rhs)
{
    if (this != &rhs) {
        // Copy primary.
        if (rhs.m_primary != NULL) {
            // Does not matter if m_primary is null
            delete m_primary;
            m_primary = rhs.m_primary->Clone();
        } else {
            delete m_primary;
            m_primary = NULL;
        }

        // Copy remaining data
        m_Position = rhs.m_Position;
        m_PositionTime = rhs.m_PositionTime;
        m_StatWeight = rhs.m_StatWeight;
        m_CoagCount = rhs.m_CoagCount;
        m_FragCount = rhs.m_FragCount;
        m_createt = rhs.m_createt;
        mLPDAtime = rhs.mLPDAtime;
    }
    return *this;
}

// PRIMARY PARTICLE CHILD.

/*!
 * Access the child primary particle which contains the physical details
 * of the particle
 *
 *@return   Pointer to primary particle or Null if none present
 */
Sweep::AggModels::Primary *const Particle::Primary()
{
    return m_primary;
}

/*!
 * Access the child primary particle which contains the physical details
 * of the particle
 *
 *@return   Pointer to primary particle or Null if none present
 */
const Sweep::AggModels::Primary *const Particle::Primary() const
{
    return m_primary;
}

/*!
 * Pass through to primary particle
 */
double Particle::SphDiameter() const
{
    return m_primary->SphDiameter();
}

/*!
 * Pass through to primary particle
 */
double Particle::CollDiameter() const
{
    return m_primary->CollDiameter();
}

/*!
 * Pass through to primary particle
 */
double Particle::MobDiameter() const
{
    return m_primary->MobDiameter();
}

/*!
 * Pass through to primary particle
 */
double Particle::SurfaceArea() const
{
    return m_primary->SurfaceArea();
}

/*!
 * Pass through to primary particle
 */
double Particle::SphSurfaceArea() const
{
    return m_primary->SphSurfaceArea();
}

/*!
 * Pass through to primary particle
 */
double Particle::Volume() const
{
    return m_primary->Volume();
}

/*!
 * Pass through to primary particle
 */
double Particle::Mass(void) const
{
    return m_primary->Mass();
}

//! Pass through to primary particle.
int Particle::NumCarbon(void) const
{
    return m_primary->NumCarbon();
}

//! Pass through to primary particle.
int Particle::Frag(void) const
{
    return m_primary->Frag();
}

int Particle::NumRings() const
{
	return m_primary->NumRings();
}
/*!
 * Pass through to primary particle
 */
double Particle::GetSites(void) const
{
    return m_primary->GetSites();
}

/*!
 * Pass through to primary particle
 */
double Particle::GetSintRate(void) const
{
    return m_primary->GetSintRate();
}

/*!
 * Pass through to primary particle
 */
double Particle::GetCoverageFraction(void) const
{
    return m_primary->GetCoverageFraction();
}

/*!
 * It is not clear to me why the number of primary particles is
 * not used as the number of sub-units (riap 01.10.2012).
 *
 * @param[in]   Reciprocal of number of sub-units in the aggregate
 *
 * @return      Geometric mean of sub-unit diameter
 */
double Particle::avgeomdiam(double oneovernumsubpart) const
{
    return std::pow(m_primary->Property(Sweep::iDsph),oneovernumsubpart);
}

//! Phase composition term for phase transformation process
double Particle::GetPhaseTerm(void) const
{
	return m_primary->GetPhaseTerm();
}

/*!
 * Provide an interface that allows run time specification of particle properties
 * for use in process rate calculations.  It is currently used for some surface
 * reactions.  Where possible, the use of specific accessors should be preferred.
 *
 *@param[in]    id      Symbolic index of property required
 *
 *@return       Requested particle property
 *
 *@pre          A valid primary particle is present (the data is not stored in SubParticle)
 */
double Particle::Property(PropID id) const
{
    switch (id) {
        case iDsph:      // Equivalent sphere diameter.
            return SphDiameter();
        case iDcol:   // Collision diameter.
            return CollDiameter();
        case iDmob:   // Mobility diameter.
            return MobDiameter();
        case iS:      // Surface area.
            return SurfaceArea();
        case iV:      // Volume.
            return Volume();
        case iM:      // Mass.
            return Mass();
        
        //! Number of carbon atoms.
        case iNumCarbon:
            return NumCarbon();

        //! Fragmentation flag.
        case iFrag:
            return Frag();

        // Collision rate properties:
        case iD2:
            return CollDiameter() * CollDiameter();
        case iD_1:
            return 1.0 / CollDiameter();
        case iD_2:
            return 1.0 / CollDiameter() / CollDiameter();
        case iM_1_2:
            return 1.0 / std::sqrt(Mass());
        case iD2_M_1_2:
            return CollDiameter() * CollDiameter() / std::sqrt(Mass());
        case iASN:
            return GetSites();
        case iSintRate:
            return GetSintRate();
        case iCoverage:
            return GetCoverageFraction();
        case iFS:
            throw std::logic_error("Free surface no longer supported (Particle::Property)");
            return 0.0;       
        case -1:
            // Special case property, used to select particles
            // uniformly.
            return 1.0;
		case iAn_2_3_comp:
			return GetPhaseTerm();
        default:
            throw std::logic_error("Unrecognised property requested (Particle::Property)");
            return 0.0;
    }
}


// PARTICLE COMPOSITION.

// Returns the composition vector.
const fvector &Particle::Composition() const
{
    return m_primary->Composition();
}

// Returns the ith component value.  Returns 0.0 if i is invalid.
double Particle::Composition(unsigned int i) const
{
    if (i < Composition().size()) {
        return Composition()[i];
    } else {
        return 0.0;
    }
}

// TRACKER VARIABLE VALUES.

/*!
 * Returns the tracker value vector.
 */
const fvector &Particle::Values() const
{
    return m_primary->Values();
}

// Returns the ith tracker variable value.  Returns 0.0 if i is invalid.
double Particle::Values(unsigned int i) const
{
    if (i < Values().size()) {
        return Values()[i];
    } else {
        return 0.0;
    }
}

/*!
 * Both position and the associated time must be updated together.  This is
 * because it makes no sense to specify a position without knowing when it
 * applies.
 *
 *@param[in]    x       New position of particle
 *@param[in]    t       Time at which new position is correct
 */
void Particle::setPositionAndTime(const double x, const double t) {
    m_Position = x;
    m_PositionTime = t;
}

/*!
 * Set last LPDA update time; function name is historical (and confusing).
 *
 *@param[in]    t   Time to which LPDA has just been performed
 */
void Particle::SetTime(double t)
{
    m_primary->SetTime(t);
    mLPDAtime = t;
}

/*!
 *  Recalculates the derived properties from the unique properties.
 *  This function moves down the tree of subparticles from the top to the bottom.
 */
void Particle::UpdateCache(void)
{
    // Get cache from primary particle.
    m_primary->UpdateCache();

    m_createt = m_primary->CreateTime();
    mLPDAtime    = m_primary->LastUpdateTime();
}

/*!
 *@return   Number of coagulation events since counter was reset
 */
unsigned int Particle::getCoagCount() const {
    return m_CoagCount;
}

/*!
 *@return   Number of coagulation events since counter was reset
 */
unsigned int Particle::getFragCount() const {
    return m_FragCount;
}

// PARTICLE ADJUSTMENT AND PROCESSES.

/*!
 * Apply the given composition and values changes n times
 *
 *@param[in]        n       Number of times to apply changes
 *
 *@return       Number of times changes actually applied
 */
unsigned int Particle::Adjust(const fvector &dcomp,
                              const fvector &dvalues,
                              rng_type &rng,
                              unsigned int n)
{
    unsigned int m = n;

    // This is a leaf-node sub-particle as it contains a
    // primary particle.  The adjustment is applied to
    // the primary.
    m = m_primary->Adjust(dcomp, dvalues, rng, n);

    // Where-ever the adjustment has been applied this sub-particle must
    // now update its cache.
    UpdateCache();

    return m;
}

/*!
 * Apply the given composition and values changes n times
 * while allowing for the special nature of IntParticleReaction.
 *
 *@param[in]        n       Number of times to apply changes
 *
 *@return       Number of times changes actually applied
 */
unsigned int Particle::AdjustIntPar(const fvector &dcomp,
                                    const fvector &dvalues,
                                    rng_type &rng,
                                    unsigned int n)
{
    unsigned int m = n;

    // This is a leaf-node sub-particle as it contains a
    // primary particle.  The adjustment is applied to
    // the primary.
    m = m_primary->AdjustIntPar(dcomp, dvalues, rng, n);

    // Where-ever the adjustment has been applied this sub-particle must
    // now update its cache.
    UpdateCache();

    return m;
}

/*!
 * Apply the given composition and values changes n times
 * for the phase transformation 
 *
 *@param[in]        n       Number of times to apply changes
 *
 *@return       Number of times changes actually applied
 */
unsigned int Particle::AdjustPhase(const fvector &dcomp,
                              const fvector &dvalues,
                              rng_type &rng,
                              unsigned int n)
{
    
	unsigned int m = n;

    // This is a leaf-node sub-particle as it contains a
    // primary particle.  The adjustment is applied to
    // the primary.
    n = m_primary->AdjustPhase(dcomp, dvalues, rng, n);

	int loops = 0; //csl37 - max number of loops set to 10
	while (n<m && loops < 20){  //csl37 - max number of loops set to 10
		n += m_primary->AdjustPhase(dcomp, dvalues, rng, m-n);
		loops ++;
	}

    // Where-ever the adjustment has been applied this sub-particle must
    // now update its cache.
    UpdateCache();

    return m;
}

/*!
 * Combines this particle with another.
 *
 * \param[in]       rhs         Particle to add to current instance
 * \param[in,out]   rng         Random number generator
 * \return      Reference to the current instance after rhs has been added
 */
Particle &Particle::Coagulate(const Particle &rhs, rng_type &rng)
{
    m_primary->Coagulate(*rhs.m_primary, rng);
    UpdateCache();

    return *this;
}

/*!
 * Combines this particle with another.
 *
 * \param[in]       rhs         Particle to add to current instance
 * \param[in,out]   rng         Random number generator
 * \return      Reference to the current instance after rhs has been added
 */
Particle &Particle::Fragment(const Particle &rhs, rng_type &rng)
{
    m_primary->Fragment(*rhs.m_primary, rng);
    UpdateCache();

    return *this;
}

/*!
 * Sinter the particle over the specified amount of time
 *
 *@param[in]        dt      Length of sintering time
 *@param[in]        sys     Cell containing the particle and providing details of the environment
 *@param[in]        model   Sintering model to define sintering rate
 *@param[in,out]    rng     Random number generator object
 *@param[in]        wt      Statistical weight of particle
 */
void Particle::Sinter(double dt, Cell &sys,
                      const Processes::SinteringModel &model,
                      rng_type &rng,
                      double wt)
{
    m_primary->Sinter(dt, sys, model, rng, wt);
}

// Creates a clone of the particle.
Particle *const Particle::Clone() const
{
    return new Particle(*this);
}

/*!
 * Perform checks on the internal data structure.  This is mainly for
 * testing and checking purposes; it should not be called from performance
 * critical sections.
 *
 * It was not quite clear in the old SubParticle code (now merged into this class)
 * to what extent this method should check the integrity
 * of the data structure and to what extent it should verify the physical
 * meaning.  I think the two purposes may have become a little mixed.
 * (riap2 24 Jan 2011)
 * 
 * @return	true iff internal structures pass tests
 */
bool Particle::IsValid() const {
	return (m_primary != NULL) && (m_primary->IsValid()) && (m_StatWeight > 0);
}

void Particle::writeParticlePOVRAY(std::ofstream &out) const
{
    // First construct the image
    Sweep::Imaging::ParticleImage image;
    image.Construct(*(this), *(Primary()->ParticleModel()));

    // Write it to file
    image.WritePOVRAY(out);

    // Now we can delete the image
}

/*!
 * @brief Writes the object to a binary stream
 *
 * @param[in,out]	 out		         Output binary stream
 * @param[in,out]    duplicates          Addresses of PAHs that have already been serialised
 * @exception		 invalid_argument    Output stream not ready
 *
 * @pre              IsValid             returns true indicating (amongst other things) that the primary particle is initialised
 */
void Particle::Serialize(std::ostream &out, void *duplicates) const
{
	assert(IsValid());

    if (out.good()) {
        ModelFactory::WritePrimary(*m_primary, out, duplicates);

        // Output the data members in this class
        out.write((char*)&m_Position, sizeof(m_Position));
        out.write((char*)&m_PositionTime, sizeof(m_PositionTime));
        out.write((char*)&m_StatWeight, sizeof(m_StatWeight));
        out.write((char*)&m_CoagCount, sizeof(m_CoagCount));
        out.write((char*)&m_FragCount, sizeof(m_FragCount));
        out.write(reinterpret_cast<const char*>(&m_createt), sizeof(m_createt));
        out.write(reinterpret_cast<const char*>(&mLPDAtime), sizeof(mLPDAtime));
    }
    else {
        throw std::invalid_argument("Output stream not ready \
                                    (Sweep, Particle::Serialize).");
    }
}

//csl37: get primary particle coords
void Particle::getParticleCoords(std::vector<fvector> &coords) const
{
	m_primary->GetPrimaryCoords(coords);
}

//csl37: initialise primary particle tracking for videos
void Particle::setTracking()
{
	m_primary->setTracking();
}

//csl37: remove primary tracking
void Particle::removeTracking()
{
	m_primary->removeTracking();}