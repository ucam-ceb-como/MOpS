#include "swp_particledata.h"
#include "swp_model.h"
#include "swp_component.h"
#include "swp_coagmodel.h"
#include <vector>
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
ParticleData::ParticleData()
{
	m_components = NULL;
	m_trackers = NULL;
	m_createt = 0.0;
	m_time = 0.0;
	m_diam = 0.0;
	m_dcol = 0.0;
	m_dmob = 0.0;
	m_surf = 0.0;
	m_vol = 0.0;
	m_mass = 0.0;
	m_coag = new CoagModelData(*this);
}

// Default constructor (providing knowledge of the
// components and tracker variables).
ParticleData::ParticleData(const Sweep::CompPtrVector &components, 
						   const Sweep::TrackPtrVector &trackers)
{
	m_components = &components;
	m_trackers = &trackers;
	m_comp.assign(components.size(), 0.0);
	m_values.assign(trackers.size(), 0.0);
}

// Copy constructor.
ParticleData::ParticleData(const Sweep::ParticleData &copy)
{
	// Use assignment operator to define.
	*this = copy;
}

// Default destructor.
ParticleData::~ParticleData()
{
	// Must destruct all of the sub-models.
    delete m_coag;
    for (ModelMap::iterator i=m_models.begin(); i!=m_models.end(); ++i) {
        delete i->second;
    }
}


// OPERATOR OVERLOADS.

// Assignment operator.
ParticleData &ParticleData::operator=(const Sweep::ParticleData &rhs)
{
	// Assign components and tracker variable definitions.
	m_components = rhs.m_components;
	m_trackers = rhs.m_trackers;

	// Copy unique properties.
	m_comp.assign(rhs.m_comp.begin(), rhs.m_comp.end());
	m_values.assign(rhs.m_values.begin(), rhs.m_values.end());
	m_createt = rhs.m_createt;
	m_time = rhs.m_time;

	// Copy the coagulation model data.
	delete m_coag;
	m_coag = rhs.m_coag->Clone();
	m_coag->SetParent(*this);

	// Copy the data for other particle models.
	m_models.clear();
	for (ModelMap::const_iterator i=rhs.m_models.begin(); i!=rhs.m_models.end(); ++i) {
		m_models[i->first] = i->second->Clone();
	}

	// Copy the derived properties.
	m_diam = rhs.m_diam;
	m_dcol = rhs.m_dcol;
	m_dmob = rhs.m_dmob;
	m_surf = rhs.m_surf;
	m_vol  = rhs.m_vol;
	m_mass = rhs.m_mass;

	return *this;
}

// Addition-assignment operator.  Throws an exception
// if the particles have different component or tracker variable
// definitions.  This function is not used to coagulate particles, rather
// it is used to sum up particle properties in the ensemble binary tree.
ParticleData &ParticleData::operator+=(const Sweep::ParticleData &rhs)
{
	if ((m_components==rhs.m_components) && (m_trackers==rhs.m_trackers)) {
		unsigned int i = 0;

		// Add the components.
		for (i=0; i!=m_comp.size(); ++i) {
			m_comp[i] += rhs.m_comp[i];
		}

		// Add the tracker values.
		for (i=0; i!=m_values.size(); ++i) {
			m_values[i] += rhs.m_values[i];
		}

		// Now allow the particle models to deal with the additions.
        for (ModelMap::iterator j=m_models.begin(); j!=m_models.end(); ++j) {
            ModelMap::const_iterator k = rhs.m_models.find(j->first);
            if (k != rhs.m_models.end()) {
                *(j->second) += *(k->second);
            }
		}

		// Allow coagulation model to deal with the addition as well.
		*m_coag += *rhs.m_coag;

        // Sum cache variables.
	    m_diam += rhs.m_diam;
	    m_dcol += rhs.m_dcol;
	    m_dmob += rhs.m_dmob;
	    m_surf += rhs.m_surf;
	    m_vol  += rhs.m_vol;
	    m_mass += rhs.m_mass;
	} else {
		throw invalid_argument("Particles cannot have different component "
			                   "definitions (Sweep, ParticleData::operator+=).");
	}

	return *this;
}

// Addition operator.
const ParticleData ParticleData::operator+(const Sweep::ParticleData &rhs) const {
	// Use copy constructor and += operator to define.
	return ParticleData(*this) += rhs;
}


// DEFINING COMPONENTS.

// Returns the component vector.
const CompPtrVector *const ParticleData::Components() const
{
	return m_components;
}

// Sets the components vector.
void ParticleData::SetComponents(const Sweep::CompPtrVector &comps)
{
	m_components = &comps;
	m_comp.resize(comps.size(), 0.0);
}


// TRACKER VARIABLES.

// Returns the tracker variable vector.
const TrackPtrVector *const ParticleData::Trackers() const
{
	return m_trackers;
}

// Sets the tracker variable vector.
void ParticleData::SetTrackers(const Sweep::TrackPtrVector &track)
{
	m_trackers = &track;
	m_values.resize(track.size(), 0.0);
}


// PARTICLE COMPOSITION.

// Returns the composition vector.
const fvector &ParticleData::Composition() const 
{
	return m_comp;
}

// Returns the ith component value.  Returns 0.0 if i is invalid.
real ParticleData::Composition(unsigned int i) const
{
	if (i < m_comp.size()) {
		return m_comp[i];
	} else {
		return 0.0;
	}
}

// Sets the composition vector.
void ParticleData::SetComposition(const Sweep::fvector &comp)
{
	m_comp.assign(comp.begin(), comp.end());
}


// TRACKER VARIABLE VALUES.

// Returns the tracker value vector.
const fvector &ParticleData::Values() const
{
	return m_values;
}

// Returns the ith tracker variable value.  Returns 0.0 if i is invalid.
real ParticleData::Values(unsigned int i) const
{
	if (i < m_values.size()) {
		return m_values[i];
	} else {
		return 0.0;
	}
}


// PARTICLE CREATE TIME.

// Returns the particle create time.
real ParticleData::CreateTime() const
{
	return m_createt;
}


// LAST UPDATE TIME.

// Returns the particle last update time.
real ParticleData::LastUpdateTime() const
{
    return m_time;
}


// COAGULATION MODEL CACHE.

// Returns the coagulation model data.
CoagModelData *const ParticleData::CoagModelCache()
{
    return m_coag;
}

// Sets the coagulation model data.  The ParticleData object
// then takes control of the data for destruction purposes.
void ParticleData::SetCoagModelCache(Sweep::CoagModelData &data)
{
    delete m_coag;
    m_coag = &data;
    m_coag->SetParent(*this);
}


// MODEL CACHE.

// Returns the model data.
ModelMap &ParticleData::ModelCache()
{
    return m_models;
}

// Returns the data for the ith model.
IModelData *const ParticleData::ModelCache(ModelType id)
{
    if (id == CoagModel_ID) {
        return m_coag;
    } else {
        ModelMap::const_iterator i = m_models.find(id);
        if (i != m_models.end()) {
            return i->second;
        } else {
            return NULL;
        }
    }
}

// Add a model to the particle definition.
void ParticleData::AddModel(Sweep::IModelData *model)
{
    if (m_models.find(model->ID()) != m_models.end()) {
        delete m_models[model->ID()];
    }
    m_models[model->ID()] = model;
    model->SetParent(*this);
}


// BASIC DERIVED PARTICLE PROPERTIES.

// Returns the particle equivalent sphere diameter.
real ParticleData::SphDiameter(void) const {return m_diam;}

// Returns the collision diameter.
real ParticleData::CollDiameter(void) const {return m_dcol;}

// Rethrns the mobility diameter.
real ParticleData::MobDiameter(void) const {return m_dmob;}

// Returns the surface area.
real ParticleData::SurfaceArea(void) const {return m_surf;}

// Returns the volume.
real ParticleData::Volume(void) const {return m_vol;}

// Returns the mass.
real ParticleData::Mass(void) const {return m_mass;}

// Returns the property with the given ID.
real ParticleData::Property(PropertyID id) const
{
    switch (id) {
        case iCTime:  // Create time.
            return m_createt;
        case iLUTime: // Last update time.
            return m_time;
        case iD:      // Equivalent sphere diameter.
            return m_diam;
        case iDcol:   // Collision diameter.
            return m_dcol;
        case iDmob:   // Mobility diameter.
            return m_dmob;
        case iS:      // Surface area.
            return m_surf;
        case iV:      // Volume.
            return m_vol;
        case iM:      // Mass.
            return m_mass;
        default:
            return 0.0;
    }
}


// READ/WRITE/COPY.

// Returns a clone of the particle data.
ParticleData *const ParticleData::Clone(void) const
{
    return new ParticleData(*this);
}
