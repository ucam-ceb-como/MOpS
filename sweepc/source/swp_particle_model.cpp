/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ParticleModel class declared in the
    swp_particle_model.h header file.

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

#include "swp_particle_model.h"
#include "swp_particle.h"
#include "swp_primary.h"
#include "swp_surfvol_primary.h"
#include "swp_pripart_primary.h"
#include "swp_model_factory.h"
#include "swp_actsites_type.h"
#include "swp_abf_model.h"
#include "swp_PAH_primary.h"
#include <stdexcept>
#include <cmath>
#include "swp_PAH_trajectory.h"

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ParticleModel::ParticleModel(void)
{
    init();
}

// Copy constructor.
ParticleModel::ParticleModel(const ParticleModel &copy)
{
	*this = copy;
}

// Default destructor.
ParticleModel::~ParticleModel(void)
{
    releaseMem();
}


// OPERATOR OVERLOADS.

// Assignment operator.
ParticleModel &ParticleModel::operator=(const ParticleModel &rhs)
{
    if (this != &rhs) {
        // Clear current ParticleModel from memory.
        releaseMem();

        // Copy species vector pointer.
        m_species = rhs.m_species;

        // Copy components.
        for (CompPtrVector::const_iterator i=rhs.m_components.begin();
             i!=rhs.m_components.end(); ++i) {
            m_components.push_back((*i)->Clone());
        }

        // Copy trackers.
        for (TrackPtrVector::const_iterator i=rhs.m_trackers.begin();
             i!=rhs.m_trackers.end(); ++i) {
            m_trackers.push_back((*i)->Clone());
        }

        // Copy particle sub-model info.
        for (SubModels::SubModelTypeSet::const_iterator i=rhs.m_submodels.begin();
            i!=rhs.m_submodels.end(); ++i) {
            m_submodels.insert(*i);
        }

        // Copy sub-particle tree flag.
        m_subpart_tree = rhs.m_subpart_tree;

        // Copy aggregation model.
        m_aggmodel = rhs.m_aggmodel;

        // Copy Knudsen drag formula parameters
        m_DragA = rhs.m_DragA;
        m_DragB = rhs.m_DragB;
        m_DragE = rhs.m_DragE;

        // Choice of transport expressions
        m_DragType = rhs.m_DragType;
        m_DiffusionType = rhs.m_DiffusionType;
        m_AdvectionType = rhs.m_AdvectionType;
    }
    return *this;
}


// CHEMICAL SPECIES.

// Returns the chemical species vector.
const Sprog::SpeciesPtrVector *const ParticleModel::Species(void) const
{
    return m_species;
}

// Sets the chemical species vector.
void ParticleModel::SetSpecies(const Sprog::SpeciesPtrVector &sp)
{
    m_species = &sp;
}


// COMPONENT DEFINITIONS.

// Returns the number of components in the ParticleModel.
unsigned int ParticleModel::ComponentCount(void) const
{
    return m_components.size();
}

// Returns the vector of particle components.
const CompPtrVector &ParticleModel::Components() const
{
    return m_components;
}

// Returns the component with the given index.
const Component *const ParticleModel::Components(unsigned int i) const
{
    // Remove index checking for optimization purpose
    return m_components[i];
    //if (i < m_components.size()) {
    //    return m_components[i];
    //} else {
    //    return NULL;
    //}
}

// Returns the index of the component with the 
// given name in the ParticleModel if found, otherwise 
// return negative.
int ParticleModel::ComponentIndex(const std::string &name) const
{
    for (unsigned int i=0; i!=m_components.size(); ++i) {
        if (name.compare(m_components[i]->Name())==0) {
            return i;
        }
    }
    return -1;
}

// Adds a component to the ParticleModel and returns the index
// of the component.
unsigned int ParticleModel::AddComponent(Component &comp)
{
    m_components.push_back(&comp);
    return m_components.size()-1;
}

// Overwrites the ith component with that given.  Previous
// component is deleted from memory.
void ParticleModel::ReplaceComponent(unsigned int i, Component &comp)
{
    if (i < m_components.size()) {
        delete m_components[i];
        m_components[i] = &comp;
    }
}

// Sets the particle components vector.
void ParticleModel::SetComponents(const CompPtrVector &comps)
{
    // Delete current components.
    CompPtrVector::iterator i;
    for (i=m_components.begin(); i!=m_components.end(); ++i) {
        delete (*i);
    }

    // Resize component vector and copy components.
    m_components.resize(comps.size());
    CompPtrVector::const_iterator ic;
    for (ic=comps.begin(); ic!=comps.end(); ++ic) {
        m_components.push_back((*ic)->Clone());
    }
}


// TRACKER VARIABLES.

// Returns the number of tracker variables.
unsigned int ParticleModel::TrackerCount(void) const 
{
    return m_trackers.size();
}

// Returns the vector of tracker variables.
const TrackPtrVector &ParticleModel::Trackers(void) const
{
    return m_trackers;
}

// Returns the ith tracker variable.
const Tracker *const ParticleModel::Trackers(unsigned int i) const
{
    if (i < m_trackers.size()) {
        return m_trackers[i];
    } else {
        return NULL;
    }
}

// Returns the index of the tracker variable with the given name 
// on success, otherwise returns negative.
int ParticleModel::GetTrackerIndex(const std::string &name) const
{
    for (unsigned int i=0; i!=m_trackers.size(); ++i) {
        if (name.compare(m_trackers[i]->Name())==0) {
            return i;
        }
    }
    return -1;
}

// Adds a tracker variable to the ParticleModel.
void ParticleModel::AddTracker(Tracker &track)
{
    m_trackers.push_back(&track);
}

// Replaces the tracker at the given index with the given tracker
// object.
void ParticleModel::ReplaceTracker(unsigned int i, Tracker &track)
{
    if (i < m_trackers.size()) {
        delete m_trackers[i];
        m_trackers[i] = &track;
    }
}

// Sets the vector of tracker variables.
void ParticleModel::SetTrackers(const TrackPtrVector &track)
{
    // Delete current trackers.
    TrackPtrVector::iterator i;
    for (i=m_trackers.begin(); i!=m_trackers.end(); ++i) {
        delete (*i);
    }

    // Resize tracker vector and copy tracker variables.
    m_trackers.resize(track.size());
    TrackPtrVector::const_iterator ic;
    for (ic=track.begin(); ic!=track.end(); ++i) {
        m_trackers.push_back((*ic)->Clone());
    }
}


// PARTICLE MODELS.

// Returns the set of particle model ID used by this ParticleModel
const SubModels::SubModelTypeSet &ParticleModel::SubModels(void) const
{
    return m_submodels;
}

// Returns true if the ParticleModel include the given model.
bool ParticleModel::ContainsSubModel(SubModels::SubModelType id) const
{
    return m_submodels.find(id) != m_submodels.end();
}

// Adds a model to the ParticleModel.  Any subsequent particles
// created with this ParticleModel will use this model.
void ParticleModel::AddSubModel(SubModels::SubModelType id)
{
    m_submodels.insert(id);

    // Initialise the model as necessary.
    switch(id) {
        case ActSites::ABFSites_ID:
            ActSites::ABFModel::Instance().Initialise(*this);
            break;
        default:
            break;
    }
}



// SUB-PARTICLE TREE.

// Returns true if the sub-particle tree is used in subsequently
// created particles, otherwise false.
bool ParticleModel::UseSubPartTree(void) const {return m_subpart_tree;}

// Enables the sub-particle tree in subsequently created particles.
void ParticleModel::EnableSubPartTree(void) {m_subpart_tree = true;}

// Disables the sub-particle tree in subsequently created particles.
void ParticleModel::DisableSubPartTree(void) {m_subpart_tree = false;}


// AGGREGATION MODEL TYPE.

// Returns the ID of the aggregation model used by dependent particles.
AggModels::AggModelType ParticleModel::AggModel(void) const {return m_aggmodel;}

// Sets the aggregation model ID for subsequently created particles.
void ParticleModel::SetAggModel(AggModels::AggModelType id) {m_aggmodel = id;}


// SINTERING MODEL.

// Returns the sintering model.
Processes::SinteringModel &ParticleModel::SintModel(void) const {return m_sint_model;}


// PARTICLE FUNCTIONS.

// Creates a new particle and sets it up with all the models
// required by the ParticleModel. A particle created here
// will have a single primary.
Sweep::Particle *const ParticleModel::CreateParticle(real time) const
{
    // Create new primary using the aggregation model currently
    // set in this model.
    Primary *pri = ModelFactory::CreatePrimary(m_aggmodel, time, *this);

    // Now create a particle from this primary.  This sets up the
    // sub-model cache within the particle automatically.
    Particle *part = new Particle(*pri);

    // Returns particle.
    return part;
}



void ParticleModel::LoadPAHProfile()
{
	m_PAHDatabase.LoadPAHProfiles();
}



//Collision Efficiency
double ParticleModel::CollisionEff(Particle *p1, Particle *p2) const
{	
	double redmass=0;
	int ncarbon1,ncarbon2;
	const AggModels::PAHPrimary *pah1 = NULL;
	const AggModels::PAHPrimary *pah2 = NULL;
	pah1 = dynamic_cast<AggModels::PAHPrimary*>(p1->Primary());
	pah2 = dynamic_cast<AggModels::PAHPrimary*>(p2->Primary());
	ncarbon2=(int)(1.0*pah2->NumCarbon());
	ncarbon1=(int)(1.0*pah1->NumCarbon());
	double x;
	double y;
	x=1e10*pah1->CollDiameter();
	y=1e10*pah2->CollDiameter();
	double reddiam;
	double ceffi;
	reddiam=min(x,y);
	redmass=12*min(ncarbon1,ncarbon2);
	ceffi=1/(1+exp(-2.0* (reddiam*reddiam*reddiam/(redmass)+pow(redmass/1100,6.0)-5.0)));
	if (pah1->NumPAH()>1 && pah2->NumPAH()>1)
		ceffi=1;
	return ceffi;


//  Used to read in the collision efficiency from a database instead from a formula
/*
	if (redmass<m_reduced_mass.at(0))
		return m_collision_eff.at(0);
	double cefflarger;
	double ceffsmaller;
	double redmasslarger;
	double redmasssmaller;
	int j;
	for (j=0;j<m_reduced_mass.size();j++)
	{
		if (redmass<m_reduced_mass.at(j))
		{
			redmasslarger=m_reduced_mass.at(j);
			cefflarger=m_collision_eff.at(j);
			break;
		}
		else
		{
			redmasssmaller=m_reduced_mass.at(j);
			ceffsmaller=m_collision_eff.at(j);
		}
		if (j==m_reduced_mass.size()-1) 
		{
			redmasslarger=m_reduced_mass.at(j);
			cefflarger=m_collision_eff.at(j);
		}
	}
	if (cefflarger==ceffsmaller) return cefflarger;
	double a=(cefflarger-ceffsmaller)/(redmasslarger-redmasssmaller);
	return ceffsmaller+(redmass-redmasssmaller)*a;*/
}

// READ/WRITE/COPY.

// Creates a copy of the ParticleModel.
ParticleModel *const ParticleModel::Clone(void) const
{
    return new ParticleModel(*this);
}

// Writes the object to a binary stream.
void ParticleModel::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write number of components.
        unsigned int n = (unsigned int)m_components.size();
        out.write((char*)&n, sizeof(n));

        // Write components.
        for (CompPtrVector::const_iterator i=m_components.begin(); 
             i!=m_components.end(); ++i) {
            (*i)->Serialize(out);
        }

        // Write number of trackers.
        n = (unsigned int)m_trackers.size();
        out.write((char*)&n, sizeof(n));

        // Write trackers.
        for (TrackPtrVector::const_iterator i=m_trackers.begin(); 
             i!=m_trackers.end(); ++i) {
            (*i)->Serialize(out);
        }

        // Write model count.
        n = (unsigned int)m_submodels.size();
        out.write((char*)&n, sizeof(n));

        // Write model set.
        for (SubModels::SubModelTypeSet::const_iterator i=m_submodels.begin(); i!=m_submodels.end(); ++i) {
            n = (unsigned int)(*i);
            out.write((char*)&n, sizeof(n));
        }

        // Write if sub-particle tree is enabled.
        if (m_subpart_tree) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write the aggregation model ID.
        n = (unsigned int)m_aggmodel;
        out.write((char*)&n, sizeof(n));

        // Write Knudsen drag parameters
        out.write(reinterpret_cast<const char *>(&m_DragA), sizeof(m_DragA));
        out.write(reinterpret_cast<const char *>(&m_DragB), sizeof(m_DragB));
        out.write(reinterpret_cast<const char *>(&m_DragE), sizeof(m_DragE));

        // Transport model choices
        out.write(reinterpret_cast<const char *>(&m_DragType), sizeof(m_DragType));
        out.write(reinterpret_cast<const char *>(&m_DiffusionType), sizeof(m_DiffusionType));
        out.write(reinterpret_cast<const char *>(&m_AdvectionType), sizeof(m_AdvectionType));

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleModel::Serialize).");
    }
}

// Reads the object from a binary stream.
void ParticleModel::Deserialize(std::istream &in)
{
    releaseMem();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n=0, id=0;

        switch (version) {
            case 0:
                // Read number of components.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read components.
                for (unsigned int i=0; i!=n; ++i) {
                    m_components.push_back(new Component(in));
                }

                // Read number of trackers.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read trackers.
                for (unsigned int i=0; i!=n; ++i) {
                    m_trackers.push_back(new Tracker(in));
                }

                // Read model count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read model set.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&id), sizeof(id));
                    m_submodels.insert((SubModels::SubModelType)id);
                }
                
                // Read if the sub-particle tree is enabled.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_subpart_tree = (n==1);

                // Read the aggregation model ID.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_aggmodel = (AggModels::AggModelType)n;

                // Read in Knudsen drag parameters
                in.read(reinterpret_cast<char*>(&m_DragA), sizeof(m_DragA));
                in.read(reinterpret_cast<char*>(&m_DragA), sizeof(m_DragA));
                in.read(reinterpret_cast<char*>(&m_DragA), sizeof(m_DragA));

                // Transport model choices
                in.read(reinterpret_cast<char*>(&m_DragType), sizeof(m_DragType));
                in.read(reinterpret_cast<char*>(&m_DiffusionType), sizeof(m_DiffusionType));
                in.read(reinterpret_cast<char*>(&m_AdvectionType), sizeof(m_AdvectionType));

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ParticleModel::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ParticleModel::Deserialize).");
    }
}


// MEMORY MANAGEMENT.

// Initialises the model to the default condition.  Used in
// constructors.
void ParticleModel::init(void)
{
    releaseMem();
    m_species      = NULL;
    m_subpart_tree = false;
    m_aggmodel     = AggModels::Spherical_ID;

    m_DragA = 0.0;
    m_DragB = 0.0;
    m_DragE = 0.0;

    // Not sure what to put as default for m_DragType etc
}

// Clears the current ParticleModel from memory.
void ParticleModel::releaseMem(void)
{
    // Clear species vector pointer.
    m_species = NULL;

    // Delete components.
    for (CompPtrVector::iterator i=m_components.begin(); 
         i!=m_components.end(); ++i) {
        delete *i;
    }
    m_components.clear();

    // Delete trackers.
    for (TrackPtrVector::iterator i=m_trackers.begin(); 
         i!=m_trackers.end();++i) {
        delete *i;
    }
    m_trackers.clear();

    // Clear sub-models.
    m_submodels.clear();

    // Set sub-particle tree and aggregation models to default values.
    m_subpart_tree = false;
    m_aggmodel     = AggModels::Spherical_ID;
}

/*!
 * The Knudsen correction to the Stokes formula for the drag coefficient depends
 * on three constants, see table I of Li & Wang, Phys. Rev. E 68, 061206 (2003).
 * The formula is
 * \f[
 *    \frac{6 \pi \mu R}{1 + Kn \left[A + B \exp\left(-E/Kn\right)\right]}
 * \f]
 *
 *@param[in]    A   Constant A from above table
 *@param[in]    B   Constant B from above table
 *@param[in]    E   Constant E from above table
 */
void ParticleModel::SetKnudsenDragConstants(const real A, const real B, const real E) {
    m_DragA = A;
    m_DragB = B;
    m_DragE = E;
}

/*!
 * The drag coefficient is calculated using the the Knudsen correction to the
 * Stokes formula for the drag coefficient, see table I of Li & Wang,
 * Phys. Rev. E 68, 061206 (2003).  Note that this article casts considerable
 * doubt on the general validity of this standard formula, which is implemented
 * here as a base case for further work.
 * The formula is
 * \f[
 *    \frac{6 \pi \mu R}{1 + Kn \left[A + B \exp\left(-E/Kn\right)\right]}.
 * \f]
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Drag coefficient
 */
real ParticleModel::KnudsenDragCoefficient(const Cell &sys, const Particle &sp) const {
    const real Kn = Sweep::KnudsenAir(sys.Temperature(), sys.Pressure(), sp.CollDiameter());

    // 3 * pi = 9.424777962 and note that diameter not radius is used below
    const real numerator = 9.424777962 * sp.CollDiameter() * Sweep::ViscosityAir(sys.Temperature());

    return numerator / (1 + Kn * (m_DragA + m_DragB * std::exp(-m_DragE / Kn)));
}

/*!
 * The drag coefficient for the free molecular regimes is calculated using a
 * reduced collision integral of 1 following the formula in table I of Li & Wang,
 * Phys. Rev. E 68, 061206 (2003).
 * The formula is
 * \f[
 *    \frac{8}{3}\sqrt{2 \pi m_{\mathrm{gas}} k T} N R^2,
 * \f]
 * where \f$ m_{\mathrm{gas}} = 0.028 \mathrm{kg mol}^{-1}\f$.
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Drag coefficient
 */
real ParticleModel::FreeMolDragCoefficient(const Cell &sys, const Particle &sp) const {
    return 5.355264342e-24 * std::sqrt(sys.Temperature()) * sys.Density() * sp.CollDiamSquared();
}

/*!
 * For testing purposes it is useful to have a drag coefficient that is
 * proprtional to temperature so that the diffusion constant calculated using
 * Einstein's relation \see{DiffusionCoefficient} is a constant independent of
 * particle properties and the gaseous environment.
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Drag coefficient
 */
real ParticleModel::TemperatureDragCoefficient(const Cell &sys, const Particle &sp) const {
    return m_DragA * sys.Temperature();
}

/*!
 * Calculate diffusion co-efficient using Einstein's relation
 * \f[
 *    D = \frac{k_B T}{k_d}.
 * \f]
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Diffusion coefficient
 */
real ParticleModel::EinsteinDiffusionCoefficient(const Cell &sys, const Particle &sp) const {
    switch(m_DragType) {
        case KnudsenDrag:
            return Sweep::KB * sys.Temperature() / KnudsenDragCoefficient(sys, sp);
            // Will not go any further because of return statement.
        case FreeMolDrag:
            return Sweep::KB * sys.Temperature() / FreeMolDragCoefficient(sys, sp);
        case TemperatureDrag:
            return Sweep::KB * sys.Temperature() / TemperatureDragCoefficient(sys, sp);
        default:
            throw std::runtime_error("Unrecognised drag type in Sweep::ParticleModel::EinsteinDiffusionCoefficient()");
    }

    return 0.0;
}

/*!
 * Calculate diffusion co-efficient using Einstein's relation
 * \f[
 *    D = \frac{k_B T}{k_d}.
 * \f]
 * In the flamelet case this is multiplied by a factor of
 * \f[
 *    \left(\frac{\partial Z}{\partial x}\right)^2.
 * \f]
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Diffusion coefficient
 */
real ParticleModel::DiffusionCoefficient(const Cell &sys, const Particle &sp) const {
    switch(m_DiffusionType) {
        case EinsteinDiffusion:
            return EinsteinDiffusionCoefficient(sys, sp);
            // Will not go any further because of return statement.
        case FlameletDiffusion:
            //@todo Multiply by dissipation rate
            return EinsteinDiffusionCoefficient(sys, sp)
                     * sys.GradientMixFrac() * sys.GradientMixFrac();
        default:
            throw std::runtime_error("Unrecognised diffusion type in Sweep::ParticleModel::DiffusionCoefficient()");
    }

    return 0.0;
}

/*!
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Advection velocity
 */
real ParticleModel::AdvectionVelocity(const Cell &sys, const Particle &sp) const {
    switch(m_AdvectionType) {
        case BulkAdvection:
            // Particle moves with the gas flow
            return sys.Velocity();

        case FlameletAdvection:
        {
            // Mass of soot per unit volume of gas [kg m^-3] is needed repeatedly
            const real sootMassDensity = sys.Particles().GetSum(ParticleCache::iM)
                                         / sys.SampleVolume();
            
            // Thermophoretic drift term
            real v = sys.GradientMixFrac() * sys.ThermoVelocity()
                     * sootMassDensity;

            // Difference in diffusion of soot and mixture fraction
            v += sys.LaplacianMixFrac() * sootMassDensity
                 * (sys.MixFracDiffCoeff() - EinsteinDiffusionCoefficient(sys, sp));

            // Another term I do not have a good way to calculate


            // Divide by density to get a velocity
            return v / sootMassDensity;
        }
        default:
            throw std::runtime_error("Unrecognised advectopm type in Sweep::ParticleModel::AdvectionVelocity()");
    }
    return 0.0;
}

/*!
 * Collision integral is calculated using the formula from
 * Z Li and H Wang, "Drag force, diffusion coefficient, and eletric mobility
 * of small particles. II. Application", Phys. Rev. E 68:061207 (2003)
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Collision integral interpolated between diffuse and specular limits
 */
real ParticleModel::Omega1_1_avg(const Cell &sys, const Particle &sp) const {
    // Knudsen number is used to interpolate between the specular and diffuse integrals
    const real knudsen = Sweep::KnudsenAir(sys.Temperature(), sys.Pressure(), sp.CollDiameter());

    // @todo Reduced collision diameter
    const real sigmaPrime = 1.0;

    //@todo Modified temperature to power -1/4
    const real TStar = 1.0;

    // Collision integrals calculated for pure specular and pure diffusion scattering
    const real specular = Omega1_1_spec(TStar, sigmaPrime);
    const real diffuse  = Omega1_1_diff(TStar, sigmaPrime);

    //@todo proper interpolation
    return (specular + diffuse) / 2.0;
}

/*!
 * Collision integral is calculated using the formula from
 * Z Li and H Wang, "Drag force, diffusion coefficient, and eletric mobility
 * of small particles. II. Application", Phys. Rev. E 68:061207 (2003)
 *
 *@param[in]    t_star_1_4      Modified temperature
 *@param[in]    sigma_prime     Reduced collision diameter
 *
 *@return   Collision integral calculated with diffuse scattering
 */
real ParticleModel::Omega1_1_diff(const real t_star_1_4, const real sigma_prime) const {
    // 1 + pi/8
    real integral = 1.3926;

    integral += sigma_prime * (1.072 + 2.078 * t_star_1_4
                               + 1.261 * t_star_1_4 * t_star_1_4);

    integral += sigma_prime * sigma_prime
                * (3.285 - 8.872 * t_star_1_4 + 5.225 * t_star_1_4 * t_star_1_4);

    return integral;
}

/*!
 * Collision integral is calculated using the formula from
 * Z Li and H Wang, "Drag force, diffusion coefficient, and eletric mobility
 * of small particles. II. Application", Phys. Rev. E 68:061207 (2003)
 *
 *@param[in]    t_star_1_4      Modified temperature
 *@param[in]    sigma_prime     Reduced collision diameter
 *
 *@return   Collision integral calculated with specular scattering
 */
real ParticleModel::Omega1_1_spec(const real t_star_1_4, const real sigma_prime) const {
    real integral = 1.0;

    integral += sigma_prime * (0.316 + 1.47 * t_star_1_4
                               + 0.476 * t_star_1_4 * t_star_1_4);

    integral += sigma_prime * sigma_prime
                * (1.53 - 85.013 * t_star_1_4 + 4.025 * t_star_1_4 * t_star_1_4);

    return integral;
}