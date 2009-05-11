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
// required by the ParticleModel.  A particle created here
// will have a single primary.
Particle *const ParticleModel::CreateParticle(real time) const
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



//Collision Efficiency
double ParticleModel::CollisionEff(Particle *p1, Particle *p2) const
{	double aeff=4e-7;
	double redmass=0;
	int ncarbon1,ncarbon2;
	const AggModels::PAHPrimary *pah1 = NULL;
	const AggModels::PAHPrimary *pah2 = NULL;
	pah1 = dynamic_cast<AggModels::PAHPrimary*>(p1->Primary());
	pah2 = dynamic_cast<AggModels::PAHPrimary*>(p2->Primary());
	ncarbon2=(int)(1.0*pah2->m_numcarbon);
	ncarbon1=(int)(1.0*pah1->m_numcarbon);
	double x;
	double y;
	x=1e10*pah1->CollDiameter();
	y=1e10*pah2->CollDiameter();
	double reddiam;
	double ceffi;
	reddiam=min(x,y);
	redmass=12*min(ncarbon1,ncarbon2);
	ceffi=1/(1+exp(-2.0* (reddiam*reddiam*reddiam/(redmass)+pow(redmass/1100,6.0)-5.0)));
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
