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
#include "swp_PAH_primary.h"
#include "swp_surfvol_primary.h"
#include "swp_model_factory.h"
#include "swp_actsites_type.h"
#include "swp_abf_model.h"
#include "swp_PAH_primary.h"
#include <stdexcept>
#include <cmath>
#include <limits>

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
        m_ThermophoresisType = rhs.m_ThermophoresisType;

        m_MoleculeStories = rhs.m_MoleculeStories;

        // Secondary particles
        m_SecondaryParticles = rhs.m_SecondaryParticles;
        m_MinSecondaryMass = rhs.m_MinSecondaryMass;
        m_MaxSecondaryMass = rhs.m_MaxSecondaryMass;
        m_MinSecondaryCollDiam = rhs.m_MinSecondaryCollDiam;
        m_MaxSecondaryCollDiam = rhs.m_MaxSecondaryCollDiam;

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

/*!
 * @param[in]       time        Time at which particle is being created
 * @param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 *
 * Creates a new particle and sets it up with all the models
 * required by the ParticleModel. A particle created here
 * will have a single primary.  Ownership of the particle
 * (in particular responsibility for calling delete) is taken
 * by the caller.
 *
 */
Sweep::Particle *const ParticleModel::CreateParticle(const real time,
                                                     int (*rand_int)(int, int)) const
{
    // Create new primary using the aggregation model currently
    // set in this model.
    Primary *pri = ModelFactory::CreatePrimary(m_aggmodel, time, *this, rand_int);

    // Now create a particle from this primary.  This sets up the
    // sub-model cache within the particle automatically.
    Particle *part = new Particle(*pri);

    // No position information available
    part->setPositionAndTime(-1.0, -1.0);

    // Returns particle.
    return part;
}

/*!
 * @param[in]   time        Time at which particle is being created
 * @param[in]   position    Position at which particle is being created
 * @param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 *
 * Creates a new particle and sets it up with all the models
 * required by the ParticleModel. A particle created here
 * will have a single primary.  Ownership of the particle
 * (in particular responsibility for calling delete) is taken
 * by the caller.
 *
 */
Sweep::Particle *const ParticleModel::CreateParticle(const real time, const real position,
                                                     int (*rand_int)(int, int)) const
{
    // Create new primary using the aggregation model currently
    // set in this model.
    Primary *pri = ModelFactory::CreatePrimary(m_aggmodel, time, position, *this, rand_int);

    // Now create a particle from this primary.  This sets up the
    // sub-model cache within the particle automatically.
    Particle *part = new Particle(*pri);

    part->setPositionAndTime(position, time);

    // Returns particle.
    return part;
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
        out.write(reinterpret_cast<const char *>(&m_ThermophoresisType), sizeof(m_ThermophoresisType));

        // Secondary particle model details
        out.write(reinterpret_cast<const char *>(&m_SecondaryParticles), sizeof(m_SecondaryParticles));
        out.write(reinterpret_cast<const char *>(&m_MinSecondaryMass), sizeof(m_MinSecondaryMass));
        out.write(reinterpret_cast<const char *>(&m_MaxSecondaryMass), sizeof(m_MaxSecondaryMass));
        out.write(reinterpret_cast<const char *>(&m_MinSecondaryCollDiam), sizeof(m_MinSecondaryCollDiam));
        out.write(reinterpret_cast<const char *>(&m_MaxSecondaryCollDiam), sizeof(m_MaxSecondaryCollDiam));
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
                in.read(reinterpret_cast<char*>(&m_ThermophoresisType), sizeof(m_ThermophoresisType));

                // Secondary particvle details
                in.read(reinterpret_cast<char*>(&m_SecondaryParticles), sizeof(m_SecondaryParticles));
                in.read(reinterpret_cast<char*>(&m_MinSecondaryMass), sizeof(m_MinSecondaryMass));
                in.read(reinterpret_cast<char*>(&m_MaxSecondaryMass), sizeof(m_MaxSecondaryMass));
                in.read(reinterpret_cast<char*>(&m_MinSecondaryCollDiam), sizeof(m_MinSecondaryCollDiam));
                in.read(reinterpret_cast<char*>(&m_MaxSecondaryCollDiam), sizeof(m_MaxSecondaryCollDiam));

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
    m_ThermophoresisType = NoThermophoresis;

    m_MoleculeStories.clear();

    m_SecondaryParticles = false;
    m_MinSecondaryMass = 0.0;
    m_MaxSecondaryMass = std::numeric_limits<real>::max();
    m_MinSecondaryCollDiam = 0.0;
    m_MaxSecondaryCollDiam = std::numeric_limits<real>::max();
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

    m_SecondaryParticles = false;
    m_MinSecondaryMass = 0.0;
    m_MaxSecondaryMass = std::numeric_limits<real>::max();
    m_MinSecondaryCollDiam = 0.0;
    m_MaxSecondaryCollDiam = std::numeric_limits<real>::max();
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
    const real d = sp.CollDiameter();
    return 5.355264342e-24 * std::sqrt(sys.Temperature()) * sys.Density() * d * d;
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
            return EinsteinDiffusionCoefficient(sys, sp)
                     * sys.GradientMixFrac() * sys.GradientMixFrac();
        default:
            throw std::runtime_error("Unrecognised diffusion type in Sweep::ParticleModel::DiffusionCoefficient()");
    }

    return 0.0;
}

/*!
 * Calculate the advection term based on soot particle transport equations of
 * various kinds.  Information on neighbouring cells is used in the flamelet
 * cases, because there are terms in gradient in mixture fraction space.
 *
 *@param[in]    sys         System in which particle experiences drag
 *@param[in]    sp          Particle for which to calculate drag coefficient
 *@param[in]    neighbours  Pointers to neighbouring cells
 *@param[in]    geom        Information on layout of neighbouring cells
 *
 *@return       Advection velocity
 *@exception    std::runtime_error  Unrecognised advection type
 */
real ParticleModel::AdvectionVelocity(const Cell &sys, const Particle &sp,
                                      const std::vector<const Cell*> &neighbours,
                                      const Geometry::LocalGeometry1d &geom) const {
    switch(m_AdvectionType) {
        case BulkAdvection:
            // Particle moves with the gas flow, thermophoresis can be set to 0
            return sys.Velocity() + ThermophoreticVelocity(sys, sp);
        case FlameletAdvection:
        {
            // Mass of soot per unit volume of gas [kg m^-3] is needed repeatedly
            const real sootMassDensity = sys.Particles().GetSum(Sweep::iM)
                                         / sys.SampleVolume();

            // Thermophoretic drift term
            real v = sys.GradientMixFrac() * ThermophoreticVelocity(sys, sp)
                     * sootMassDensity;

            // Difference in diffusion of soot and mixture fraction
            v += sys.LaplacianMixFrac() * sootMassDensity
                 * (sys.MixFracDiffCoeff() - EinsteinDiffusionCoefficient(sys, sp));

            if((neighbours[0] != NULL) && (neighbours[1] != NULL))
            // Estimate two gradients in mixture fraction space using the
            // formula grad f ~ (f(z_{i+1}) - f(z_{i-1})) / (z_{i+1} - z_{i-1})
            // This will not work if one or both of the neghbouring cells
            // is missing (at the edge of a domain) and for simplicity
            // assume this term is 0 in this case.
            // @todo Use one sided derivatives at the boundary
            {
                // Distance between the two Z values
                const real dZ = geom.calcSpacing(Geometry::left) +
                                geom.calcSpacing(Geometry::right);

                // rho D_Z at z_{i-1}
                const real leftRhoDZ  = neighbours[0]->MassDensity() *
                                        neighbours[0]->MixFracDiffCoeff();
                // rho D_Z at z_{i+1}
                const real rightRhoDZ = neighbours[1]->MassDensity() *
                                        neighbours[1]->MixFracDiffCoeff();
                // Now calculate the gradient estimate for the produect
                // of gas mass density and mixture fraction diffusion
                // coefficient.
                const real gradRhoDZ = (rightRhoDZ - leftRhoDZ) / dZ;

                // Same process for product of soot mass per unit volume of gas
                // and particule diffusion coefficient of this particle.
                const real leftVal  = neighbours[0]->Particles().GetSum(Sweep::iM) /
                                      neighbours[0]->SampleVolume() *
                                      EinsteinDiffusionCoefficient(*neighbours[0], sp);
                const real rightVal = neighbours[1]->Particles().GetSum(Sweep::iM) /
                                      neighbours[1]->SampleVolume() *
                                      EinsteinDiffusionCoefficient(*neighbours[1], sp);
                // Finish the gradient calculation
                const real grad = (rightVal - leftVal) / dZ;


                // Put together this term in the equation and add it to the
                // velocity.
                v+= sys.GradientMixFrac() * sys.GradientMixFrac() *
                    (sootMassDensity / sys.MassDensity() * gradRhoDZ - grad);
            }


            // Divide by density to get a velocity
            return v / sootMassDensity;
        }
        default:
            throw std::runtime_error("Unrecognised advection type in Sweep::ParticleModel::AdvectionVelocity()");
    }
    return 0.0;
}

/*!
 * Formulae taken from Z Li and H Wang, "Thermophoretic force and velocity of
 * nanoparticles in the free molecule regime" Phys. Rev. E 70:021205 (2004).
 *
 *@param[in]    sys         System in which particle experiences drag
 *@param[in]    sp          Particle for which to calculate drag coefficient
 *
 *@return       Thermophoretic velocity
 *@exception    std::runtime_error  Unrecognised thermophoresis type
 */
real ParticleModel::ThermophoreticVelocity(const Cell &sys, const Particle &sp) const {

    // Declare outside the switch statement to keep the compiler happy
    real tempFactor = 0.0;

    switch(m_ThermophoresisType) {
        case WaldmannThermophoresis:
            tempFactor = sys.getThermalConductivity(sys.Pressure())
                         * sys.GradientTemperature() / sys.Pressure();

            // Equation 2 of the Li & Wang paper with phi = 0.9
            // 1 / (5 * (1 + pi * phi / 8)) == 0.14777
            return tempFactor * -0.14777;

        case NoThermophoresis:
            return 0.0;
        default:
            throw std::runtime_error("Unrecognised advection type in Sweep::ParticleModel::AdvectionVelocity()");
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
    // Reduced collision diameter
    const real sigmaPrime = collisionIntegralDiameter(sys, sp);

    // Modified temperature to power -1/4
    const real TStar = std::pow(collisionIntegralTemperature(sys, sp), -0.25);

    // Collision integrals calculated for pure specular and pure diffusion scattering
    const real specular = Omega1_1_spec(TStar, sigmaPrime);
    const real diffuse  = Omega1_1_diff(TStar, sigmaPrime);

    // Interpolation via the accomodation function
    const real phi = accomodationFunction(sys, sp);
    return specular * (1.0 - phi) + diffuse * phi;
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
real ParticleModel::Omega1_2_avg(const Cell &sys, const Particle &sp) const {
    // Reduced collision diameter
    const real sigmaPrime = collisionIntegralDiameter(sys, sp);

    // Modified temperature to power -1/4
    const real TStar = std::pow(collisionIntegralTemperature(sys, sp), -0.25);

    // Collision integrals calculated for pure specular and pure diffusion scattering
    const real specular = Omega1_2_spec(TStar, sigmaPrime);
    const real diffuse  = Omega1_2_diff(TStar, sigmaPrime);

    // Interpolation via the accomodation function
    const real phi = accomodationFunction(sys, sp);
    return specular * (1.0 - phi) + diffuse * phi;
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
                * (1.53 - 5.013 * t_star_1_4 + 4.025 * t_star_1_4 * t_star_1_4);

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
 *@return   Collision integral calculated with diffuse scattering
 */
real ParticleModel::Omega1_2_diff(const real t_star_1_4, const real sigma_prime) const {
    // 1 + 5 * pi/48
    real integral = 1.3272;

    integral += sigma_prime * (1.159 + 1.506 * t_star_1_4
                               + 1.204 * t_star_1_4 * t_star_1_4);

    integral += sigma_prime * sigma_prime
                * (3.028 - 7.719 * t_star_1_4 + 4.180 * t_star_1_4 * t_star_1_4);

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
real ParticleModel::Omega1_2_spec(const real t_star_1_4, const real sigma_prime) const {
    real integral = 1.0;

    integral += sigma_prime * (0.338 + 1.315 * t_star_1_4
                               + 0.412 * t_star_1_4 * t_star_1_4);

    integral += sigma_prime * sigma_prime
                * (1.503 - 4.654 * t_star_1_4 + 3.410 * t_star_1_4 * t_star_1_4);

    return integral;
}


/*!
 * Dimensionless collision diameter for use in formul\ae from
 * Z Li and H Wang, "Drag force, diffusion coefficient, and eletric mobility
 * of small particles. II. Application", Phys. Rev. E 68:061207 (2003)
 *
 *@param[in]    sys     Gaseous system in which particle is moving
 *@param[in]    sp      Particle
 *
 *@return       Lennard Jones collision radius divided by particle collision radius
 */
real ParticleModel::collisionIntegralDiameter(const Cell &sys, const Particle &sp) const {
    // Value of Lennard Jones radius (m) taken from Table IV of above cited paper
    // using combination rules to combine values for gaseous and solid phases.
    return 3.576e-10 / sp.CollDiameter();
}


/*!
 * Dimensionless temperature for use in formul\ae from
 * Z Li and H Wang, "Drag force, diffusion coefficient, and eletric mobility
 * of small particles. II. Application", Phys. Rev. E 68:061207 (2003)
 *
 *@param[in]    sys     Gaseous system in which particle is moving
 *@param[in]    sp      Particle
 *
 *@return       Temperature divided by Lennard Jones well depth and collision cross-section cubed
 */
real ParticleModel::collisionIntegralTemperature(const Cell &sys, const Particle &sp) const {
    // Value of well depth (K) taken from Table IV of above cited paper
    // using combination rules to combine values for gaseous and solid phases.
    const real wellDepth = 57.24;

    const real collisionCubed = std::pow(collisionIntegralDiameter(sys, sp), 3);

    // 1.107e-29 is volume occupied by 1 Carbon atom assuming bulk density 1800 kg m^-3 (value for soot)
    return 3.0 * Sprog::kB * sys.Temperature() * 1.107e-29 / (2.0 * Sprog::PI * wellDepth * collisionCubed);
}

/*!
 * Proportion of diffuse scattering result to use in formul\ae from
 * Z Li and H Wang, "Drag force, diffusion coefficient, and eletric mobility
 * of small particles. II. Application", Phys. Rev. E 68:061207 (2003).
 * This formula is most clearly presented in the follow on paper:
 * Z Li and H Wang, "Thermophoretic force and velocity of nanoparticles
 * in the free molecular regime", Phys, Rev. E 70:021205 (2004), see
 * equation 38.
 *
 *@param[in]    sys     Gaseous system in which particle is moving
 *@param[in]    sp      Particle
 *
 *@return       Accomodation function, \f$\phi\f$
 */
real ParticleModel::accomodationFunction(const Cell &sys, const Particle &sp) const {
    // Particle radius in nm divided by 2.5
    real switchTerm = sp.CollDiameter() / 5e-9;

    // High power to switch aggressively between regimes
    switchTerm = std::pow(switchTerm, 15);
    switchTerm += 1.0;

    switchTerm = 1.0 - 1.0 / switchTerm;

    // Knudsen number is used to interpolate between the specular and diffuse integrals
    const real knudsen = Sweep::KnudsenAir(sys.Temperature(), sys.Pressure(), sp.CollDiameter());

    return (1.0 + 0.9 * knudsen * switchTerm) / (1.0 + knudsen);
}

/*!
 * @param[in]   sp  Particle to check
 *
 * @return      True if sufficiently small and simple for a secondary population
 */
bool ParticleModel::isSecondary(const Particle &sp) const {
    if(m_SecondaryParticles){
        if((sp.Mass() > m_MinSecondaryMass) && (sp.CollDiameter() < m_MaxSecondaryCollDiam)) {
            if(m_aggmodel == AggModels::PAH_ID) {
                const AggModels::PAHPrimary* pri = dynamic_cast<const AggModels::PAHPrimary*>(sp.Primary());

                return (pri->Numprimary() == 1);
            }
            else
                return true;
        }
    }

    return false;

}

/*!
 * @param[in]   on_off      True iff some particles are to be treated as secondary particles
 */
void ParticleModel::setSecondary(bool on_off) {
    m_SecondaryParticles = on_off;
}
