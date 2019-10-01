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

		// Copy phases
		for (PhasePtrVector::const_iterator i=rhs.m_phases.begin();
			i != rhs.m_phases.end(); ++i) {
			m_phases.push_back((*i)->Clone());
		}

        // Copy aggregation model.
        m_aggmodel = rhs.m_aggmodel;

        // Copy Knudsen drag formula parameters
        m_DragA = rhs.m_DragA;
        m_DragB = rhs.m_DragB;
        m_DragE = rhs.m_DragE;
        
        //! Collision efficiency model.
        colliParaA = rhs.colliParaA;
        colliParaB = rhs.colliParaB;
        colliParaC = rhs.colliParaC;
        m_inceptionThreshold = rhs.m_inceptionThreshold;
        m_condensationThreshold = rhs.m_condensationThreshold;
        m_mode = rhs.m_mode;

        // Choice of transport expressions
        m_DragType = rhs.m_DragType;
        m_DiffusionType = rhs.m_DiffusionType;
        m_AdvectionType = rhs.m_AdvectionType;
        m_ThermophoresisType = rhs.m_ThermophoresisType;

        // BinTree settings
        m_write_bintree         = rhs.m_write_bintree;
        m_bintree_coalthresh    = rhs.m_bintree_coalthresh;
        m_fract_dim             = rhs.m_fract_dim;

        // Indices for various mixture properties
         m_TemperatureGradientIndex = rhs.m_TemperatureGradientIndex;
         m_MixFracGradientIndex = rhs.m_MixFracGradientIndex;
         m_MixFracDiffusionIndex = rhs.m_MixFracDiffusionIndex;
         m_MixFracLaplacianIndex = rhs.m_MixFracLaplacianIndex;
         m_AvgMolWtIndex = rhs.m_AvgMolWtIndex;
         m_ThermalConductivityIndex = rhs.m_ThermalConductivityIndex;

         m_postprocessingType = rhs.m_postprocessingType;

         //! Tracking of primary separation or coordinates.
         m_trackPrimarySeparation = rhs.m_trackPrimarySeparation;
         m_trackPrimaryCoordinates = rhs.m_trackPrimaryCoordinates;
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

// PHASE

// Add phase to particle model
void ParticleModel::AddPhase(Phase &phase)
{
	m_phases.push_back(&phase);
}

// Returns the number of phases in the ParticleModel.
unsigned int ParticleModel::PhaseCount(void) const
{
	return m_phases.size();
}

// Returns if the phase in liquid
bool ParticleModel::PhaseIsLiquid(int i) const
{
	return m_phases[i]->GetLiquid();
}

// Returns the component indices
std::vector<unsigned int> ParticleModel::GetPhaseComponents(int i) const
{
	return m_phases[i]->GetComponents();
}

// Returns the index of the phase with the
// given name in the ParticleModel if found, otherwise
// return negative.
int ParticleModel::PhaseIndex(const std::string &name) const
{
	for (unsigned int i = 0; i != m_phases.size(); ++i) {
		if (name.compare(m_phases[i]->Name()) == 0) {
			return i;
		}
	}
	return -1;
}

// Returns the phase with the given index.
Phase *const ParticleModel::Phases(unsigned int i)
{
	return m_phases[i];
}

// AGGREGATION MODEL TYPE.

// Returns the ID of the aggregation model used by dependent particles.
AggModels::AggModelType ParticleModel::AggModel(void) const {return m_aggmodel;}

// Sets the aggregation model ID for subsequently created particles.
void ParticleModel::SetAggModel(AggModels::AggModelType id) {m_aggmodel = id;}


// SINTERING MODEL.

// Returns the sintering model.
Processes::SinteringModel &ParticleModel::SintModel(void) const {return m_sint_model;}

// Return titania phase transformation (melting) model
Processes::MeltingModel &ParticleModel::MeltModel(void) const { return m_melt_model; }

// PARTICLE FUNCTIONS.

/*!
 * @param[in]       time        Time at which particle is being created
 *
 * Creates a new particle and sets it up with all the models
 * required by the ParticleModel. A particle created here
 * will have a single primary.  Ownership of the particle
 * (in particular responsibility for calling delete) is taken
 * by the caller.
 *
 */
Sweep::Particle *const ParticleModel::CreateParticle(const double time) const
{
    // Create new primary using the aggregation model currently
    // set in this model.
    AggModels::Primary *pri = ModelFactory::CreatePrimary(m_aggmodel, time, *this);

    // Now create a particle from this primary.  This sets up the
    // sub-model cache within the particle automatically.
    Particle *part = new Particle(*pri);

    // No position information available
    part->setPositionAndTime(-1.0, -1.0);

    // Returns particle.
    return part;
}

/*!
 * @param[in]       time        Time at which particle is being created
 * @param[in]       position    Position at which particle is being created
 *
 * Creates a new particle and sets it up with all the models
 * required by the ParticleModel. A particle created here
 * will have a single primary.  Ownership of the particle
 * (in particular responsibility for calling delete) is taken
 * by the caller.
 *
 */
Sweep::Particle *const ParticleModel::CreateParticle(const double time, const double position) const
{
    // Create new primary using the aggregation model currently
    // set in this model.
    AggModels::Primary *pri = ModelFactory::CreatePrimary(m_aggmodel, time, position, *this);

    // Now create a particle from this primary.  This sets up the
    // sub-model cache within the particle automatically.
    Particle *part = new Particle(*pri);

    part->setPositionAndTime(position, time);

    // Returns particle.
    return part;
}

/*!
 * @brief Calculate collision efficiency based on the two selected particles.
 *
 * @param[in]    p1    Particle 1.
 * @param[in]    p2    Particle 2.
 */ 
double ParticleModel::CollisionEff(Particle *p1, Particle *p2) const
{
    int ncarbon1, ncarbon2;
    const AggModels::PAHPrimary *pah1 = NULL;
    const AggModels::PAHPrimary *pah2 = NULL;
    pah1 = dynamic_cast<AggModels::PAHPrimary*>(p1->Primary());
    pah2 = dynamic_cast<AggModels::PAHPrimary*>(p2->Primary()); 
    double ceffi;

    if (Mode() == "NONE" || Mode() == "") {
        double A = ColliParaA();
        double B = ColliParaB();
        double C = ColliParaC();
        double redmass=0;
        ncarbon2 = (int)(1.0 * pah2->NumCarbon());
        ncarbon1 = (int)(1.0 * pah1->NumCarbon());
        double x;
        double y;
        x = 1e10 * pah1->CollDiameter();
        y = 1e10 * pah2->CollDiameter();
        double reddiam;
        reddiam = min(x, y);
        redmass = 12 * min(ncarbon1, ncarbon2);
        ceffi = 1 / (1 + exp(-A * (reddiam * reddiam * reddiam / (redmass) + pow(redmass / B, 6.0) - C)));
        if (pah1->NumPAH() >1 && pah2->NumPAH() > 1)
            ceffi = 1;
        return ceffi;
    } else {
        //! Thresholds are based upon the total number of 6-member rings (excludes 5-member rings) in the PAH.
        int target_Rings_Inception = inceptionThreshold();
        int target_Rings_Condensation = condensationThreshold();

        int nRings1,nRings2;

        //! Inception: PAH + PAH = Particle.
        if (pah1->NumPAH() == 1 && pah2->NumPAH() == 1) {
            nRings2 = (int)(1.0 * pah2->NumRings());
            nRings1 = (int)(1.0 * pah1->NumRings());
            double redmass;
            if (Mode() == "MAX")
                redmass = max(nRings1, nRings2);
            else if (Mode() == "MIN")
                redmass = min(nRings1, nRings2);
            else if (Mode() == "COMBINED")
                redmass = nRings1 + nRings2;
            else if (Mode() == "REDUCED")
                redmass = nRings1 * nRings2 / (nRings1 + nRings2);
            else throw std::runtime_error("Mode of collision efficiency is modified by unknown process. Please check Sweep::ParticleModel::CollisionEff()."); 
            if (redmass >= target_Rings_Inception) ceffi = 1.0;
            else ceffi = 0;
        }

        //! Condensation: PAH + Particle (2 >= PAHs) = Particle.
        else if(pah1->NumPAH() > 1 && pah2->NumPAH() == 1){
            nRings2 = (int)(1.0*pah2->NumRings());
            if (nRings2 >= target_Rings_Condensation) ceffi = 1;
            else ceffi = 0;
        }
        else if(pah1->NumPAH() == 1 && pah2->NumPAH() > 1){
            nRings1 = (int)(1.0*pah1->NumRings());
            if (nRings1 >= target_Rings_Condensation) ceffi = 1;
            else ceffi = 0;
        }
        
        //! Coagulation event: Particle + Particle = Particle. All particles are able to coagulate.
        else {
            ceffi = 1;
        }

        return ceffi;
    }
}

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
    return ceffsmaller+(redmass-redmasssmaller)*a;
}*/

// READ/WRITE/COPY.

// Creates a copy of the ParticleModel.
ParticleModel *const ParticleModel::Clone(void) const
{
    return new ParticleModel(*this);
}

// Writes the object to a binary stream.
void ParticleModel::Serialize(std::ostream &out) const
{
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

		// Write number of phases
		n = (unsigned int)m_phases.size();
		out.write((char*)&n, sizeof(n));

		// Write phases
		for (PhasePtrVector::const_iterator i = m_phases.begin();
			i != m_phases.end(); ++i) {
			(*i)->Serialize(out);
		}

        // Write the aggregation model ID.
        n = (unsigned int)m_aggmodel;
        out.write((char*)&n, sizeof(n));

        // Write the mode.
        n = (unsigned int)m_mode.length();
        out.write((char*)&n, sizeof(n));
        out.write(m_mode.c_str(), n);

        //! Write the thresholds.
        n = (unsigned int)m_inceptionThreshold;
        out.write((char*)&n, sizeof(n));
        n = (unsigned int)m_condensationThreshold;
        out.write((char*)&n, sizeof(n));

        //write the postprocess species
        n = (unsigned int)m_InceptedPAH;
        out.write((char*)&n, sizeof(n));
        // Write whether binary trees should be serialised.
        bool flag(false);
        flag = m_write_bintree;
        out.write((char*)&flag, sizeof(flag));

        // Write the coalescence threshold
        double var(0.0);
        var = m_bintree_coalthresh;
        out.write((char*)&var, sizeof(var));

        // Write the fractal dimension
        var = m_fract_dim;
        out.write((char*)&var, sizeof(var));

        // Write the FM enhancement factor
        var = m_efm;
        out.write((char*)&var, sizeof(var));

        // Write Knudsen drag parameters
        out.write(reinterpret_cast<const char *>(&m_DragA), sizeof(m_DragA));
        out.write(reinterpret_cast<const char *>(&m_DragB), sizeof(m_DragB));
        out.write(reinterpret_cast<const char *>(&m_DragE), sizeof(m_DragE));

        // Transport model choices
        out.write(reinterpret_cast<const char *>(&m_DragType), sizeof(m_DragType));
        out.write(reinterpret_cast<const char *>(&m_DiffusionType), sizeof(m_DiffusionType));
        out.write(reinterpret_cast<const char *>(&m_AdvectionType), sizeof(m_AdvectionType));
        out.write(reinterpret_cast<const char *>(&m_ThermophoresisType), sizeof(m_ThermophoresisType));

        // Indices for various mixture properties
        out.write(reinterpret_cast<const char *>(&m_TemperatureGradientIndex), sizeof(m_TemperatureGradientIndex));
        out.write(reinterpret_cast<const char *>(&m_MixFracGradientIndex),     sizeof(m_MixFracGradientIndex));
        out.write(reinterpret_cast<const char *>(&m_MixFracDiffusionIndex),    sizeof(m_MixFracDiffusionIndex));
        out.write(reinterpret_cast<const char *>(&m_MixFracLaplacianIndex),    sizeof(m_MixFracLaplacianIndex));
        out.write(reinterpret_cast<const char *>(&m_AvgMolWtIndex),            sizeof(m_AvgMolWtIndex));
        out.write(reinterpret_cast<const char *>(&m_ThermalConductivityIndex), sizeof(m_ThermalConductivityIndex));

        out.write(reinterpret_cast<const char *>(&m_postprocessingType), sizeof(m_postprocessingType));

        //! Write whether the distance between the centres of primary particles
        //! is to be tracked.
        flag = m_trackPrimarySeparation;
        out.write((char*)&flag, sizeof(flag));

        //! Write whether the coordinates of primary particles are to be
        //! tracked.
        flag = m_trackPrimaryCoordinates;
        out.write((char*)&flag, sizeof(flag));

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleModel::Serialize).");
    }
}

// Reads the object from a binary stream.
void ParticleModel::Deserialize(std::istream &in)
{
    releaseMem();
    char *name = NULL;

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n=0;
        double var(0.0);
        bool flag(false);

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

				// Read number of phases.
				in.read(reinterpret_cast<char*>(&n), sizeof(n));

				// Read phases.
				for (unsigned int i = 0; i != n; ++i) {
					m_phases.push_back(new Phase(in));
				}

                // Read the aggregation model ID.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_aggmodel = (AggModels::AggModelType)n;

                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                name = new char[n];
                in.read(name, n);
                m_mode = string(name, n);
                delete [] name;

                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_inceptionThreshold = (int)n;

                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_condensationThreshold = (int)n;

                
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_InceptedPAH = (PostProcessStartingStr)n;
                // Read if binary trees should be read..
                in.read(reinterpret_cast<char*>(&flag), sizeof(flag));
                m_write_bintree = flag;

                // Read in the coalescence threshold
                in.read(reinterpret_cast<char*>(&var), sizeof(var));
                m_bintree_coalthresh    = var;

                // Read in the fractal dimension
                in.read(reinterpret_cast<char*>(&var), sizeof(var));
                m_fract_dim             = var;

                // Read in the FM enhancement factor
                in.read(reinterpret_cast<char*>(&var), sizeof(var));
                m_efm             = var;

                // Read in Knudsen drag parameters
                in.read(reinterpret_cast<char*>(&m_DragA), sizeof(m_DragA));
                in.read(reinterpret_cast<char*>(&m_DragB), sizeof(m_DragB));
                in.read(reinterpret_cast<char*>(&m_DragE), sizeof(m_DragE));

                // Transport model choices
                in.read(reinterpret_cast<char*>(&m_DragType), sizeof(m_DragType));
                in.read(reinterpret_cast<char*>(&m_DiffusionType), sizeof(m_DiffusionType));
                in.read(reinterpret_cast<char*>(&m_AdvectionType), sizeof(m_AdvectionType));
                in.read(reinterpret_cast<char*>(&m_ThermophoresisType), sizeof(m_ThermophoresisType));

                // Indices for various mixture properties
                in.read(reinterpret_cast<char*>(&m_TemperatureGradientIndex), sizeof(m_TemperatureGradientIndex));
                in.read(reinterpret_cast<char*>(&m_MixFracGradientIndex),     sizeof(m_MixFracGradientIndex));
                in.read(reinterpret_cast<char*>(&m_MixFracDiffusionIndex),    sizeof(m_MixFracDiffusionIndex));
                in.read(reinterpret_cast<char*>(&m_MixFracLaplacianIndex),    sizeof(m_MixFracLaplacianIndex));
                in.read(reinterpret_cast<char*>(&m_AvgMolWtIndex),            sizeof(m_AvgMolWtIndex));
                in.read(reinterpret_cast<char*>(&m_ThermalConductivityIndex), sizeof(m_ThermalConductivityIndex));

                in.read(reinterpret_cast<char*>(&m_postprocessingType), sizeof(m_postprocessingType));

                //! Read whether the distance between the centres of primary
                //! particles is to be tracked.
                in.read(reinterpret_cast<char*>(&flag), sizeof(flag));
                m_trackPrimarySeparation = flag;

                //! Read whether coordinates of primary particles are to be
                //! tracked.
                in.read(reinterpret_cast<char*>(&flag), sizeof(flag));
                m_trackPrimaryCoordinates = flag;
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
    m_aggmodel     = AggModels::Spherical_ID;

    // Parameters of the Knudsen Ansatz for drag, only used for the Knudsen drag model,
    // reasonable values are given in comments, but should be provided via the sweep.xml file.
    // m_DragA is also used for the temperature drag model, which has no physical meaning, but
    // may be useful for testing purposes.
    m_DragA = 0.0; //1.155;
    m_DragB = 0.0; //0.471;
    m_DragE = 0.0; //0.596;
    
    //! Initial CE model.
    colliParaA = 0.0; 
    colliParaB = 0.0;
    colliParaC = 0.0;
    m_inceptionThreshold = 0;
    m_condensationThreshold = 0;
    m_mode = "";

    // Not sure what to put as default for m_DragType etc
    m_ThermophoresisType = NoThermophoresis;

    // Default writing of binary trees is false.
    m_write_bintree = false;
    m_bintree_coalthresh    = 1.0;
    m_fract_dim             = 1.8;

    // Default FM enhancement factor is 2.2
    m_efm = 2.2;

    //! Postprocess based on the inception species concentration.
    m_postprocessingType = concentration;

    //! Primary particles are assumed to be in contact.
    m_trackPrimarySeparation = false;
    m_trackPrimaryCoordinates = false;
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


	// Delete phases.
	for (PhasePtrVector::iterator i = m_phases.begin();
		i != m_phases.end(); ++i) {
		delete *i;
	}
	m_phases.clear();

    // Set aggregation model to default value.
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
void ParticleModel::SetKnudsenDragConstants(const double A, const double B, const double E) {
    m_DragA = A;
    m_DragB = B;
    m_DragE = E;
}


//! SET PARAMETER FOR abhjeet's COLLISION EFFICIENCY MODEL
void Sweep::ParticleModel::SetCollisionEffPara(const double A, const double B, const double C) {
    colliParaA = A;
    colliParaB = B;
    colliParaC = C;
}

//! return parameters of abhjeet's collision efficiency model
double ParticleModel::ColliParaA() const {return colliParaA;}
double ParticleModel::ColliParaB() const {return colliParaB;}
double ParticleModel::ColliParaC() const {return colliParaC;}

//! Set the minimum number of 6-member rings (excludes 5-member rings) a PAH has to have to be able to incept.
void ParticleModel::setInceptionThreshold(const int target) {m_inceptionThreshold = target;}

//! Set the minimum number of 6-member rings (excludes 5-member rings) a PAH has to have to be able to condense onto a particle (2 or more PAHs).
void ParticleModel::setCondensationThreshold(const int target) {m_condensationThreshold = target;}

//! Return the minimum number of 6-member rings (excludes 5-member rings) a PAH has to have to be able to incept.
int ParticleModel::inceptionThreshold() const {return m_inceptionThreshold;}

//! Return the minimum number of 6-member rings (excludes 5-member rings) a PAH has to have to be able to condense onto a particle (2 or more PAHs).
int ParticleModel::condensationThreshold() const {return m_condensationThreshold;}

//! Set mode for collision efficiency model, currently 4 modes are supported, min, max, combined and reduced.
void ParticleModel::SetMode(const std::string &mode) {m_mode = mode;}

//! specify the incepting pah, currently only pyrene and benzene is supported for PAH-PP model, details see the kmc
void ParticleModel::SetInceptedPAH(const std::string &name) 
{
    if (!name.compare("A1"))
        m_InceptedPAH = A1;
    else if (!name.compare("A2"))
        m_InceptedPAH = A2;
    else if (!name.compare("A3"))
        throw std::runtime_error("A3 is not supported as InceptedPAH currently and please use A1, A2, or A4 (Sweep::ParticleModel::SetInceptedPAH())");
    else if (!name.compare("A4"))
            m_InceptedPAH = A4;
    else if (!name.compare("A5"))
            m_InceptedPAH = A5;
    else throw std::runtime_error("no information about the incepted PAH is available, only A1 A2 and A4 are supported now (Sweep::ParticleModel::SetInceptedPAH())");
}

//! return mode of collision efficency model
const std::string &ParticleModel::Mode() const {return m_mode;}

const ParticleModel::PostProcessStartingStr &ParticleModel::InceptedPAH() const {return m_InceptedPAH;}

//bool ParticleModel::IsPyreneInception() const
//{
//    if (m_InceptedPAH == "A4" || m_InceptedPAH == "pyrene")
//        return true;
//    else if (m_InceptedPAH == "A1" || m_InceptedPAH == "benzene")
//        return false;
//    else throw std::runtime_error("no information about the incepted PAH is available (Sweep::ParticleModel::IsPyreneInception())");
//}
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
double ParticleModel::KnudsenDragCoefficient(const Cell &sys, const Particle &sp) const {
    const double Kn = Sweep::KnudsenAir(sys.GasPhase().Temperature(), sys.GasPhase().Pressure(), sp.CollDiameter());

    // 3 * pi = 9.424777962 and note that diameter not radius is used below
    const double numerator = 9.424777962 * sp.CollDiameter() * sys.GasPhase().Viscosity();

//    std::cout << "Knudsen " << Kn << ", " << numerator
//              << ", MFP " << Sweep::MeanFreePathAir(sys.GasPhase().Temperature(), sys.GasPhase().Pressure())
//              << ", diam " << sp.CollDiameter() << '\n'
//              << ", T " << sys.GasPhase().Temperature() << ", " << sys.GasPhase().Pressure() << '\n';

    return numerator / (1 + Kn * (m_DragA + m_DragB * std::exp(-m_DragE / Kn)));
}

/*!
 * The drag coefficient for the free molecular regime is calculated using a
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
double ParticleModel::FreeMolDragCoefficient(const Cell &sys, const Particle &sp) const {
    const double d = sp.CollDiameter();
    // 4.818546232410640188 = sqrt(2 * pi * k * NA) * 2 / 3
    return 4.818546232410640188 * std::sqrt(sys.GasPhase().Temperature() * sys.GasPhase().PropertyValue(m_AvgMolWtIndex))
                                * sys.GasPhase().MolarDensity() * d * d;
}

/*!
 * For testing purposes it is useful to have a drag coefficient that is
 * proportional to temperature so that the diffusion constant calculated using
 * Einstein's relation \see{DiffusionCoefficient} is a constant independent of
 * particle properties and the gaseous environment.
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Drag coefficient
 */
double ParticleModel::TemperatureDragCoefficient(const Cell &sys, const Particle &sp) const {
    return m_DragA * sys.GasPhase().Temperature();
}

/*!
 * The drag coefficient for the free molecular regimes is calculated using a
 * reduced collision integral of 1 following equation (25) of Li & Wang,
 * Phys. Rev. E 68, 061207 (2003).
 * The formula is
 * \f[
 *    \frac{8}{3}\left(1 + \alpha^{-1.143}\right)^{-0.875}
 *    \sqrt{2 \pi m_{\mathrm{gas}} k T} N R^2
 *    \Omega^{(1,1)\ast}_{avg},
 * \f]
 * where for the molecular weight \f$ m_{\mathrm{gas}}\f$ of the gas molecules, a mixture average is used.
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Drag coefficient
 */
double ParticleModel::LiWangDragCoefficient(const Cell &sys, const Particle &sp) const {
    const double omega = Omega1_1_avg(sys, sp);

    // 3.247911159372846712e-24 = sqrt((2*pi)/(k*NA)) * 9 / 4 / NA
    double alpha = 3.24791159e-24 / sqrt(sys.GasPhase().Temperature() * sys.GasPhase().PropertyValue(m_AvgMolWtIndex));
    alpha *= sys.GasPhase().Viscosity() / sys.GasPhase().MolarDensity() / sp.CollDiameter();
    alpha /= omega;

    // 4.818546232410640188 = sqrt(2 * pi * k * NA) * 2 / 3
    double drag = 4.81854623 * std::sqrt(sys.GasPhase().Temperature() * sys.GasPhase().PropertyValue(m_AvgMolWtIndex))
                * sys.GasPhase().MolarDensity() * sp.CollDiameter() * sp.CollDiameter() * omega
                * pow(pow(1 + alpha, -1.143),-0.875);

    //std::cout << "Omega " << omega << ", Drag " << drag <<  '\n';

    return drag;
}

/*!
 * The drag coefficient for the free molecular regimes is calculated using a
 * reduced collision integral of 1 following equation (25) of Li & Wang,
 * Phys. Rev. E 68, 061207 (2003) and an alternative interpolation of my
 * (Patterson) own devising in order to avoid calls to pow().
 * The weights are
 * \f[
 *     \frac{1}{1 +\mathrm{Kn}}, \quad \frac{\mathrm{Kn}}{1 +\mathrm{Kn}},
 * \f]
 * where \f$\mathrm{Kn}\f$ is the Knudsen number.  The low Knudsen number
 * (continuum regime) limit is
 * \f[
 *     3 \pi \mu D
 * \f]
 * where \f$\mu\f$ is the fluid viscosity and \f$D\f$ the particle collision
 * diameter.  Interpolation is linear using the weights above and the large
 * Knudsen number (free molecular regime) limit is
 * \f[
 *    \frac{2}{3} \sqrt{2 \pi m_{\mathrm{r}} k T} N D^2
 *    \Omega^{(1,1)\ast}_{avg},
 * \f]
 * with the reduced mass given by
 * \f[
 *    m_{\mathrm{r}} = \frac{m_{\mathrm{p}}m_{\mathrm{g}}}{m_{\mathrm{p}} + m_{\mathrm{g}}}
 * \f]
 * where for the molecular weight \f$ m_{\mathrm{g}}\f$ of the gas molecules, a mixture average is used.
 * On the assumption that \f$ m_{\mathrm{r}} \gg m_{\mathrm{g}} \f$, a two term Taylor series
 * is used for the square root.
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *
 *@return       Drag coefficient
 */
double ParticleModel::LiWangPatDragCoefficient(const Cell &sys, const Particle &sp) const {
    const double omega = Omega1_1_avg(sys, sp);
    const double knudsen = Sweep::KnudsenAir(sys.GasPhase().Temperature(), sys.GasPhase().Pressure(),
                                           sp.CollDiameter());

    // 4.818546232410640188 = sqrt(2 * pi * k * NA) * 2 / 3
    // 1.2044286e+24 = 2 * NA
    double drag = 4.818546232410640188 * std::sqrt(sys.GasPhase().Temperature() * sys.GasPhase().PropertyValue(m_AvgMolWtIndex))
                * sys.GasPhase().MolarDensity() * sp.CollDiameter() * omega * knudsen
                * (1 - sys.GasPhase().PropertyValue(m_AvgMolWtIndex) / sp.Mass() / 1.2044283e+24);

    // 9.424777960769379348 = 3 * pi
    drag += 9.424777960769379348 * sys.GasPhase().Viscosity();

    drag *= sp.CollDiameter() / (1 + knudsen);

    return drag;
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
double ParticleModel::EinsteinDiffusionCoefficient(const Cell &sys, const Particle &sp) const {
    switch(m_DragType) {
        case KnudsenDrag:
            return Sweep::KB * sys.GasPhase().Temperature() / KnudsenDragCoefficient(sys, sp);
            // Will not go any further because of return statement.
        case FreeMolDrag:
            return Sweep::KB * sys.GasPhase().Temperature() / FreeMolDragCoefficient(sys, sp);
        case TemperatureDrag:
            return Sweep::KB * sys.GasPhase().Temperature() / TemperatureDragCoefficient(sys, sp);
        case LiWangDrag:
            return Sweep::KB * sys.GasPhase().Temperature() / LiWangDragCoefficient(sys, sp);
        case LiWangPatDrag:
            return Sweep::KB * sys.GasPhase().Temperature() / LiWangPatDragCoefficient(sys, sp);
        default:
            throw std::runtime_error("Unrecognised drag type in Sweep::ParticleModel::EinsteinDiffusionCoefficient()");
    }

    return 0.0;
}

/*!
 * Numerically calculate the gradient of the diffusion coefficient, \ref DiffusionCoefficient.
 *
 * Use the twosided approximation
 * \f[
 *     \nabla D(x) \approx \frac{D(x_{i+1} - D(x_{i-1}))}{x_{i+1} - x_{i-1}},
 * \f]
 * except at domain boundaries where 0 is assumed (\todo 1 sided approximation).
 *
 *@param[in]    sys     System in which particle experiences drag
 *@param[in]    sp      Particle for which to calculate drag coefficient
 *@param[in]    neighbours  Pointers to neighbouring cells
 *@param[in]    geom        Information on layout of neighbouring cells
 *
 *@return       Diffusion coefficient gradient
 */
double ParticleModel::GradDiffusionCoefficient(const Cell &sys, const Particle &sp,
                                             const std::vector<const Cell*> &neighbours,
                                             const Geometry::LocalGeometry1d &geom) const {
    if((neighbours[0] != NULL) && (neighbours[1] != NULL)) {
        const double dx = geom.calcSpacing(Geometry::left) +
                        geom.calcSpacing(Geometry::right);

        const double leftD  = DiffusionCoefficient(*(neighbours[0]), sp);
        const double rightD = DiffusionCoefficient(*(neighbours[1]), sp);

        return (rightD - leftD) / dx;
    }
    else
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
 *
 *@exception    std::runtime_error      Unrecognised diffusion type
 */
double ParticleModel::DiffusionCoefficient(const Cell &sys, const Particle &sp) const {
    switch(m_DiffusionType) {
        case EinsteinDiffusion:
            return EinsteinDiffusionCoefficient(sys, sp);
            // Will not go any further because of return statement.
        case FlameletDiffusion:
            return EinsteinDiffusionCoefficient(sys, sp)
                   * sys.GasPhase().PropertyValue(m_MixFracGradientIndex) * sys.GasPhase().PropertyValue(m_MixFracGradientIndex);
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
double ParticleModel::AdvectionVelocity(const Cell &sys, const Particle &sp,
                                      const std::vector<const Cell*> &neighbours,
                                      const Geometry::LocalGeometry1d &geom) const {
    switch(m_AdvectionType) {
        case BulkAdvection:
            // Particle moves with the gas flow, thermophoresis can be set to 0
            return sys.GasPhase().Velocity() + ThermophoreticVelocity(sys, sp);
        case FlameletAdvection:
        {
            // Mass of soot per unit volume of gas [kg m^-3] is needed repeatedly
            const double sootMassDensity = sys.Particles().GetSum(Sweep::iM)
                                         / sys.SampleVolume();

            // Thermophoretic drift term
            double v = sys.GasPhase().PropertyValue(m_MixFracGradientIndex) * ThermophoreticVelocity(sys, sp)
                     * sootMassDensity;

            // Difference in diffusion of soot and mixture fraction
            v += sys.GasPhase().PropertyValue(m_MixFracLaplacianIndex) * sootMassDensity
                 * (sys.GasPhase().PropertyValue(m_MixFracDiffusionIndex) - EinsteinDiffusionCoefficient(sys, sp));

            if((neighbours[0] != NULL) && (neighbours[1] != NULL))
            // Estimate two gradients in mixture fraction space using the
            // formula grad f ~ (f(z_{i+1}) - f(z_{i-1})) / (z_{i+1} - z_{i-1})
            // This will not work if one or both of the neghbouring cells
            // is missing (at the edge of a domain) and for simplicity
            // assume this term is 0 in this case.
            // @todo Use one sided derivatives at the boundary
            {
                // Distance between the two Z values
                const double dZ = geom.calcSpacing(Geometry::left) +
                                geom.calcSpacing(Geometry::right);

                // rho D_Z at z_{i-1}
                const double leftRhoDZ  = neighbours[0]->GasPhase().MassDensity() *
                                        neighbours[0]->GasPhase().PropertyValue(m_MixFracDiffusionIndex);
                // rho D_Z at z_{i+1}
                const double rightRhoDZ = neighbours[1]->GasPhase().MassDensity() *
                                        neighbours[1]->GasPhase().PropertyValue(m_MixFracDiffusionIndex);
                // Now calculate the gradient estimate for the product
                // of gas mass density and mixture fraction diffusion
                // coefficient.
                const double gradRhoDZ = (rightRhoDZ - leftRhoDZ) / dZ;

                // Same process for product of soot mass per unit volume of gas
                // and particle diffusion coefficient of this particle.
                const double leftVal  = neighbours[0]->Particles().GetSum(Sweep::iM) /
                                      neighbours[0]->SampleVolume() *
                                      EinsteinDiffusionCoefficient(*neighbours[0], sp);
                const double rightVal = neighbours[1]->Particles().GetSum(Sweep::iM) /
                                      neighbours[1]->SampleVolume() *
                                      EinsteinDiffusionCoefficient(*neighbours[1], sp);
                // Finish the gradient calculation
                const double grad = (rightVal - leftVal) / dZ;


                // Put together this term in the equation and add it to the
                // velocity.
                v+= sys.GasPhase().PropertyValue(m_MixFracGradientIndex) * sys.GasPhase().PropertyValue(m_MixFracGradientIndex) *
                    (sootMassDensity / sys.GasPhase().MassDensity() * gradRhoDZ - grad);
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
double ParticleModel::ThermophoreticVelocity(const Cell &sys, const Particle &sp) const {

    // Declare outside the switch statement to keep the compiler happy
    double tempFactor = 0.0;

    switch(m_ThermophoresisType) {
        case WaldmannThermophoresis:
            tempFactor = sys.GasPhase().PropertyValue(m_ThermalConductivityIndex);
            tempFactor *= sys.GasPhase().PropertyValue(m_TemperatureGradientIndex);
            tempFactor /= sys.GasPhase().Pressure();

            // Equation 2 of the Li & Wang paper with phi = 0.9
            // 1 / (5 * (1 + pi * phi / 8)) == 0.14777
            return tempFactor * -0.14777;

        case LiWangThermophoresis:
            // conductivity * temperature gradient / (NkT)
            tempFactor = sys.GasPhase().PropertyValue(m_ThermalConductivityIndex);
            tempFactor *= sys.GasPhase().PropertyValue(m_TemperatureGradientIndex);
            tempFactor /= sys.GasPhase().Pressure();

            // Now bring in the collision integrals
            return tempFactor * (1 - 6.0 * Omega1_2_avg(sys, sp) / 5.0 / Omega1_1_avg(sys, sp));

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
double ParticleModel::Omega1_1_avg(const Cell &sys, const Particle &sp) const {
    // Reduced collision diameter
    const double sigmaPrime = collisionIntegralDiameter(sys, sp);

    // Modified temperature to power -1/4
    const double TStar = std::pow(collisionIntegralTemperature(sys, sp), -0.25);

    // Collision integrals calculated for pure specular and pure diffusion scattering
    const double specular = Omega1_1_spec(TStar, sigmaPrime);
    const double diffuse  = Omega1_1_diff(TStar, sigmaPrime);
    const double Kn = Sweep::KnudsenAir(sys.GasPhase().Temperature(), sys.GasPhase().Pressure(), sp.CollDiameter());

    //std::cout << "Omega1_1_avg: " << specular << ", " << diffuse << ", " << Kn << '\n';

    // Build up the return value in stages, so that one can debug the process
    double omega = diffuse;
    omega += Kn * (0.9 * diffuse + 0.1 * specular);
    omega -= Kn * 0.9 * (diffuse - specular) / (1 + std::pow(sp.CollDiameter() / 5.0e-9, 15));
    omega /= (1 + Kn);

    return omega;
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
double ParticleModel::Omega1_2_avg(const Cell &sys, const Particle &sp) const {
    // Reduced collision diameter
    const double sigmaPrime = collisionIntegralDiameter(sys, sp);

    // Modified temperature to power -1/4
    const double TStar = std::pow(collisionIntegralTemperature(sys, sp), -0.25);

    // Collision integrals calculated for pure specular and pure diffusion scattering
    const double specular = Omega1_2_spec(TStar, sigmaPrime);
    const double diffuse  = Omega1_2_diff(TStar, sigmaPrime);
    const double Kn = Sweep::KnudsenAir(sys.GasPhase().Temperature(), sys.GasPhase().Pressure(), sp.CollDiameter());

    //std::cout << "Omega1_1_avg: " << specular << ", " << diffuse << ", " << Kn << '\n';

    // Build up the return value in stages, so that one can debug the process
    double omega = diffuse;
    omega += Kn * (0.9 * diffuse + 0.1 * specular);
    omega -= Kn * 0.9 * (diffuse - specular) / (1 + std::pow(sp.CollDiameter() / 5.0e-9, 15));
    omega /= (1 + Kn);

    return omega;
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
double ParticleModel::Omega1_1_diff(const double t_star_1_4, const double sigma_prime) const {
    // 1 + pi/8
    double integral = 1.3926;

    //std::cout << "Omega1_1_diff: " << t_star_1_4 << ", " << sigma_prime << '\n';

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
double ParticleModel::Omega1_1_spec(const double t_star_1_4, const double sigma_prime) const {
    double integral = 1.0;

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
double ParticleModel::Omega1_2_diff(const double t_star_1_4, const double sigma_prime) const {
    // 1 + 5 * pi/48
    double integral = 1.3272;

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
double ParticleModel::Omega1_2_spec(const double t_star_1_4, const double sigma_prime) const {
    double integral = 1.0;

    integral += sigma_prime * (0.338 + 1.315 * t_star_1_4
                               + 0.412 * t_star_1_4 * t_star_1_4);

    integral += sigma_prime * sigma_prime
                * (1.503 - 4.654 * t_star_1_4 + 3.410 * t_star_1_4 * t_star_1_4);

    return integral;
}


/*!
 * Dimensionless collision diameter for use in formula from
 * Z Li and H Wang, "Drag force, diffusion coefficient, and eletric mobility
 * of small particles. II. Application", Phys. Rev. E 68:061207 (2003)
 *
 *@param[in]    sys     Gaseous system in which particle is moving
 *@param[in]    sp      Particle
 *
 *@return       Lennard Jones collision radius divided by particle collision radius
 */
double ParticleModel::collisionIntegralDiameter(const Cell &sys, const Particle &sp) const {
    // Value of Lennard Jones radius (m) taken from Table IV of above cited paper
    // using combination rules to combine values for gaseous and solid phases.
    //@TODO this assumes Carbon (as protein!) and nitrogen - need to use double material properties
    return 7.152e-10 / sp.CollDiameter();
}


/*!
 * Dimensionless temperature for use in formula from
 * Z Li and H Wang, "Drag force, diffusion coefficient, and eletric mobility
 * of small particles. II. Application", Phys. Rev. E 68:061207 (2003)
 *
 *@param[in]    sys     Gaseous system in which particle is moving
 *@param[in]    sp      Particle
 *
 *@return       Temperature divided by Lennard Jones well depth and collision cross-section cubed
 */
double ParticleModel::collisionIntegralTemperature(const Cell &sys, const Particle &sp) const {
    // Value of well depth (K, ie divided by k) taken from Table IV of above cited paper
    // using combination rules to combine values for Carbon (as protein!) and N2 gas
    // @TODO This needs to be based on inputs
    const double wellDepth = 57.24;
    // Interaction diameter from same source as well depth: 3.576e-10 m
    const double interactionDiameterCubed = 4.5729086976e-29;

    // 1.107e-29 is volume occupied by 1 Carbon atom assuming bulk density 1800 kg m^-3 (value for soot)
    //@TODO this needs to be calculated for the material in use
    // 1.5 * 1.107e-29 / pi = 5.285535660081844307e-30
    return 5.285535660081844307e-30 * sys.GasPhase().Temperature() / (wellDepth * interactionDiameterCubed);
}

/*!
 * Proportion of diffuse scattering result to use in formula from
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
double ParticleModel::accomodationFunction(const Cell &sys, const Particle &sp) const {
    // Particle radius in nm divided by 2.5
    double switchTerm = sp.CollDiameter() / 5e-9;

    // High power to switch aggressively between regimes
    switchTerm = std::pow(switchTerm, 15);
    switchTerm += 1.0;

    switchTerm = 1.0 - 1.0 / switchTerm;

    // Knudsen number is used to interpolate between the specular and diffuse integrals
    const double knudsen = Sweep::KnudsenAir(sys.GasPhase().Temperature(), sys.GasPhase().Pressure(), sp.CollDiameter());

    return (1.0 + 0.9 * knudsen * switchTerm) / (1.0 + knudsen);
}
