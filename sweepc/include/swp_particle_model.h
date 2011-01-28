/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ParticleModel class defines the set of chemical species, particle
    components and particle tracker variables used to define particles and
    particle processes.  It is principally used as the base class of the
    Mechanism class, but also has use in creating and defining particles,
    which need not know of the particle processes.  This class also knows
    what sub-models are implemented in each particle, and provides routines
    for creating particles which are defined with these sub-models.

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

#ifndef SWEEP_PARTICLE_MODEL_H
#define SWEEP_PARTICLE_MODEL_H

#include "swp_params.h"
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_submodel_type.h"
#include "swp_aggmodel_type.h"
#include "swp_sintering_model.h"
#include "sprog.h"
#include "swp_molecule_evolution.h"

#include "local_geometry1d.h"

#include <vector>
#include <string>
#include <iostream>
namespace Sweep
{
// Forward declare the Particle class.
class Particle;

//! Class that defines the interpretation of the data in a particle
class ParticleModel
{
public:
	// Constructors.
    ParticleModel(void);                      // Default Constructor.
    ParticleModel(const ParticleModel &copy); // Copy-Constructor.
    ParticleModel(std::istream &in);          // Stream-reading constructor.

    // Destructor.
    ~ParticleModel(void);

    // Operators.
    ParticleModel &operator=(const ParticleModel &rhs);

	// CHEMICAL SPECIES.

    // Returns the chemical species vector.
    const Sprog::SpeciesPtrVector *const Species(void) const;

    // Sets the chemical species vector.
    void SetSpecies(const Sprog::SpeciesPtrVector &sp);


	// COMPONENT DEFINITIONS.

    // Returns the number of components in the mechanism.
    unsigned int ComponentCount(void) const;

    // Returns the vector of particle components.
    const CompPtrVector &Components() const;

    // Returns the component with the given index.
    const Component *const Components(unsigned int i) const;

    // Returns the index of the component with the
	// given name in the mechanism if found, otherwise
	// return negative.
    int ComponentIndex(const std::string &name) const;

    // Adds a component to the mechanism and returns the index
    // of the component.
    unsigned int AddComponent(Component &comp);

    // Overwrites the ith component with that given.  Previous
    // component is deleted from memory.
    void ReplaceComponent(
        unsigned int i, // Index of component to overwrite.
        Component &comp // New component.
        );

    // Sets the particle components vector.
    void SetComponents(const CompPtrVector &comps);


	// TRACKER VARIABLES.

    // Returns the number of tracker variables.
    unsigned int TrackerCount(void) const;

    // Returns the vector of tracker variables.
    const TrackPtrVector &Trackers(void) const;

    // Returns the ith tracker variable.
    const Tracker *const Trackers(unsigned int i) const;

    // Returns the index of the tracker variable with the given name
	// on success, otherwise returns negative.
    int GetTrackerIndex(const std::string &name) const;

    // Adds a tracker variable to the mechanism.
    void AddTracker(Tracker &track);

    // Replaces the tracker at the given index with the given tracker
    // object.
    void ReplaceTracker(
        unsigned int i, // Index of tracker variable to overwrite.
        Tracker &track  // New track variable.
        );

    // Sets the vector of tracker variables.
    void SetTrackers(const TrackPtrVector &track);


    // PARTICLE MODELS.

    // Returns the set of particle model ID used by this mechanism
    const SubModels::SubModelTypeSet &SubModels(void) const;

    // Returns true if the mechanism include the given model.
    bool ContainsSubModel(SubModels::SubModelType id) const;

    // Adds a model to the mechanism.  Any subsequent particles
    // created with this model will use that sub-model.
    void AddSubModel(SubModels::SubModelType id);


    // SUB-PARTICLE TREE.

    // Returns true if the sub-particle tree is used in subsequently
    // created particles, otherwise false.
    bool UseSubPartTree(void) const;

    // Enables the sub-particle tree in subsequently created particles.
    void EnableSubPartTree(void);

    // Disables the sub-particle tree in subsequently created particles.
    void DisableSubPartTree(void);


    // AGGREGATION MODEL TYPE.

    // Returns the ID of the aggregation model used by dependent particles.
    AggModels::AggModelType AggModel(void) const;

    // Sets the aggregation model ID for subsequently created particles.
    void SetAggModel(AggModels::AggModelType id);


    // SINTERING MODEL.

    // Returns the sintering model.
    Processes::SinteringModel &SintModel(void) const;


    // PARTICLE FUNCTIONS.

    //! Create a new particle on the heap
    Particle *const CreateParticle(const real time) const;

    //! Create a new particle on the heap
    Particle *const CreateParticle(const real time, const real position) const;


    // READ/WRITE/COPY.

    // Creates a copy of the mechanism.
    ParticleModel *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

    double CollisionEff(Particle *p1, Particle *p2) const;

    //! Initialise the Knudsem drag coefficient  calculation
    void SetKnudsenDragConstants(const real A, const real B, const real E);

    //! Calculate the drag coefficient for a particle using the Knudsen correction
    real KnudsenDragCoefficient(const Cell &sys, const Particle &sp) const;

    //! Calculate the drag coefficient for a particle using the Knudsen correction
    real FreeMolDragCoefficient(const Cell &sys, const Particle &sp) const;

    //! Calculate the drag coefficient for a particle as constant times temperature
    real TemperatureDragCoefficient(const Cell &sys, const Particle &sp) const;

    //! Calculate a particle diffusion coefficient
    real DiffusionCoefficient(const Cell &sys, const Particle &sp) const;

    //! Calculate the advection velocity
    real AdvectionVelocity(const Cell &sys, const Particle &sp,
                           const std::vector<const Cell*> &neighbours,
                           const Geometry::LocalGeometry1d &geom) const;

    //! Calculate the thermophoretic velocity
    real ThermophoreticVelocity(const Cell &sys, const Particle &sp) const;

    //! Select between possible drag models for use in diffusion coefficients
    enum DragType {
        //! Knudsen expression for drag coefficient
        KnudsenDrag,

        //! Free molecular regime drag coefficient
        FreeMolDrag,

        //! Drag coefficient proprtional to temperature
        TemperatureDrag,
    };

    //! Select between possible diffusion terms
    enum DiffusionType {
        //! Flamelet diffusion
        FlameletDiffusion,

        //! Ordinary physical diffusion according to Einstein's formula
        EinsteinDiffusion,
    };

    //! Select between possible advection terms
    enum AdvectionType {
        //! Move with bulk gas velocity in physical space
        BulkAdvection,

        //! Flamelet advection
        FlameletAdvection,
    };

    //! Select between possible thermophoresis terms
    enum ThermophoresisType {
        //! Traditional
        WaldmannThermophoresis,

        //! Thermphoresis turned off
        NoThermophoresis,

        // wang

    };

    //! Choose between drag models
    void SetDragType(const DragType& drag) {m_DragType = drag;}

    //! Choose the diffusion expression
    void setDiffusionType(const DiffusionType& diff) {m_DiffusionType = diff;}

    //! Choose the advection model
    void setAdvectionType(const AdvectionType& adv) {m_AdvectionType = adv;}

    //! Choose the thermophoresis model
    void setThermophoresisType(const ThermophoresisType& therm) {m_ThermophoresisType = therm;}

    //! Access molecule development histories
    const MoleculeEvolution::Database& getMoleculeStories() const {return m_MoleculeStories;}

    //! Set the contents of the database of molecule development histories
    void setMoleculeStories(const MoleculeEvolution::Database& db) {m_MoleculeStories = db;}

    //! See if a particle is eligible for the secondary population
    bool isSecondary(const Particle &sp) const;

    //! Control the use of secondary particles
    void setSecondary(bool on_off);

    //! Minimum mass of secondary particles
    real minSecondaryMass() const {return m_MinSecondaryMass;}

    //! Set minimum mass of secondary particles
    void setMinSecondaryMass(real min_val) {m_MinSecondaryMass = min_val;}

    //! Maximum mass of secondary particles
    real maxSecondaryMass() const {return m_MaxSecondaryMass;}

    //! Set maximum mass of secondary particles
    void setMaxSecondaryMass(real max_val) {m_MaxSecondaryMass = max_val;}

    //! Minimum collision diameter of secondary particles
    real minSecondaryCollDiam() const {return m_MinSecondaryCollDiam;}

    //! Set minimum collision diameter of secondary particles
    void setMinSecondaryCollDiam(real min_val) {m_MinSecondaryCollDiam = min_val;}

    //! Maximum collision diameter of secondary particles
    real maxSecondaryCollDiam() const {return m_MaxSecondaryCollDiam;}

    //! Set maximum collision diameter of secondary particles
    void setMaxSecondaryCollDiam(real max_val) {m_MaxSecondaryCollDiam = max_val;}

protected:
    // The species used to define the processes and the particles.
    const Sprog::SpeciesPtrVector *m_species;

    // Particle component and tracker variable definitions.
    CompPtrVector m_components; // The components used to build particles.
    TrackPtrVector m_trackers;  // Tracker variables.

    // Set of models which all particles produced by this model
    // must use.
    SubModels::SubModelTypeSet m_submodels;

    // Flag to enable/disable sub-particle tree in subsequently created
    // particles.
    bool m_subpart_tree;

    // The ID of the aggregation model used to create particles.
    AggModels::AggModelType m_aggmodel;

    // The sintering model.
    mutable Processes::SinteringModel m_sint_model;

    //! Calculate a phsyical particle diffusion coefficient from its drag
    real EinsteinDiffusionCoefficient(const Cell &sys, const Particle &sp) const;

    //==== Collision integrals ==========================

    //! Modified reduced diameter for use in collision integral calculations(should perhaps be input)
    real collisionIntegralDiameter(const Cell &sys, const Particle &sp) const;

    //! Dimensionless temperature for use in collision  integral calculations(should perhaps be input)
    real collisionIntegralTemperature(const Cell &sys, const Particle &sp) const;

    //! Blending function between diffuse and specular scattering for collision integrals
    real accomodationFunction(const Cell &sys, const Particle &sp) const;

    //! Knudsen based average of specular and diffusive 1,1 integrals
    real Omega1_1_avg(const Cell &sys, const Particle &sp) const;

    //! Diffuse scattering 1,1 collision integral approximation
    real Omega1_1_diff(const real t_star, const real sigma_prime) const;

    //! Specular scattering 1,1 collision integral approximation
    real Omega1_1_spec(const real t_star, const real sigma_prime) const;


    //! Knudsen based average of specular and diffusive 1,2 integrals
    real Omega1_2_avg(const Cell &sys, const Particle &sp) const;

    //! Diffuse scattering 1,2 collision integral approximation
    real Omega1_2_diff(const real t_star, const real sigma_prime) const;

    //! Specular scattering 1,2 collision integral approximation
    real Omega1_2_spec(const real t_star, const real sigma_prime) const;

    // MEMORY MANAGEMENT.

    // Initialises the model to the default condition.  Used in
    // constructors.
    void init(void);

    // Clears the model from memory.
    void releaseMem(void);

private:
    //! Constant A from Knudsen drag coefficient expression
    real m_DragA;

    //! Constant B from Knudsen drag coefficient expression
    real m_DragB;

    //! Constant E from Knudsen drag coefficient expression
    real m_DragE;

    //! Drag expression to use
    DragType m_DragType;

    //! Diffusion expression to use
    DiffusionType m_DiffusionType;

    //! Advection expression to use
    AdvectionType m_AdvectionType;

    //! Thermphoresis expression to use
    ThermophoresisType m_ThermophoresisType;

    //! Database of stories of the development of molecules that make up particles
    MoleculeEvolution::Database m_MoleculeStories;

    //! Is the secondary particle model in use
    bool m_SecondaryParticles;

    //! Minimum mass of a secondary particle / kg
    real m_MinSecondaryMass;

    //! Minimum mass of a secondary particle / kg
    real m_MaxSecondaryMass;

    //! Maximum collision diameter of a secondary particle / m
    real m_MinSecondaryCollDiam;

    //! Maximum collision diameter of a secondary particle / m
    real m_MaxSecondaryCollDiam;
};
} //namespace Sweep
#endif
