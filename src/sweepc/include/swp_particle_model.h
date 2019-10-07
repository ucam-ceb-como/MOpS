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
#include "swp_aggmodel_type.h"
#include "swp_sintering_model.h"
#include "swp_titania_melting_model.h"
#include "swp_environment_interface.h"
#include "swp_phase.h"

#include "local_geometry1d.h"

#include "gpc_species.h"

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

	// PHASES

	// Phases
	PhasePtrVector m_phases;

	// Add phase to particle model
	void AddPhase(Phase &phase);
	
	// Return component indices for a phase 
	std::vector<unsigned int> GetPhaseComponents(int i) const;

	// Returns the number of phases in the ParticleModel
	unsigned int PhaseCount(void) const;

	// Return if the phase is liquid
	bool PhaseIsLiquid(int i) const;

	// Return the phase index
	int PhaseIndex(const std::string &name) const;

	// Return phase with given index
	Phase *const Phases(unsigned int i);

    // AGGREGATION MODEL TYPE.

    // Returns the ID of the aggregation model used by dependent particles.
    AggModels::AggModelType AggModel(void) const;

    // Sets the aggregation model ID for subsequently created particles.
    void SetAggModel(AggModels::AggModelType id);


    // SINTERING MODEL.

    // Returns the sintering model.
    Processes::SinteringModel &SintModel(void) const;

	// MELTING MODEL

	// Return the titania phase melting model
	Processes::MeltingModel &MeltModel(void) const;
	
    // PARTICLE FUNCTIONS.

    //! Create a new particle on the heap
    Particle *const CreateParticle(const double time) const;

    //! Create a new particle on the heap
    Particle *const CreateParticle(const double time, const double position) const;


    // READ/WRITE/COPY.

    // Creates a copy of the mechanism.
    ParticleModel *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

    double CollisionEff(Particle *p1, Particle *p2) const;

    //! Initialise the Knudsem drag coefficient  calculation
    void SetKnudsenDragConstants(const double A, const double B, const double E);

    //! Particle diffusion coefficient
    double DiffusionCoefficient(const Cell &sys, const Particle &sp) const;

    //! Gradient of particle diffusion coefficient
    double GradDiffusionCoefficient(const Cell &sys, const Particle &sp,
                                  const std::vector<const Cell*> &neighbours,
                                  const Geometry::LocalGeometry1d &geom) const;

    //! Calculate the advection velocity
    double AdvectionVelocity(const Cell &sys, const Particle &sp,
                           const std::vector<const Cell*> &neighbours,
                           const Geometry::LocalGeometry1d &geom) const;

    //! Calculate the thermophoretic velocity
    double ThermophoreticVelocity(const Cell &sys, const Particle &sp) const;

    //! Select between possible drag models for use in diffusion coefficients
    enum DragType {
        //! Knudsen expression for drag coefficient
        KnudsenDrag,

        //! Free molecular regime drag coefficient
        FreeMolDrag,

        //! Drag coefficient proportional to temperature
        TemperatureDrag,

        //! Fit proposed in Li & Wang, Phys. Rev. E 68, 061207 (2003).
        LiWangDrag,

        //! Modified version of fit proposed in Li & Wang, Phys. Rev. E 68, 061207 (2003).
        LiWangPatDrag,
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

        //! Fit proposed in Li & Wang, Phys. Rev. E 70, 021205 (2004).
        LiWangThermophoresis,

    };

    //! Gas-phase trasfer species.
    enum PostProcessStartingStr {
        A1, A2, A4, A5
    };

    //! Postprocess based on the inception species concentration 
	//! (XA4 for soot, concentration otherwise) or based on the 
    //! molar rate of production by chemical reaction of the inception
    //! species per unit volume (wdotA4).
    enum postprocessingType {
        XA4,
        wdotA4,
		concentration,
    };

    //! Choose between drag models
    void SetDragType(const DragType& drag) {m_DragType = drag;}

    //! Choose the diffusion expression
    void setDiffusionType(const DiffusionType& diff) {m_DiffusionType = diff;}

    //! Choose the advection model
    void setAdvectionType(const AdvectionType& adv) {m_AdvectionType = adv;}

    //! Choose the thermophoresis model
    void setThermophoresisType(const ThermophoresisType& therm) {m_ThermophoresisType = therm;}

    double ColliParaA() const;
    double ColliParaB() const;
    double ColliParaC() const;

    void SetCollisionEffPara(double A, double B, double C);

    //! Set the minimum number of 6-member rings (excludes 5-member rings) a
    //! PAH has to have to be able to incept.
    void setInceptionThreshold(int target);

    //! Set the minimum number of 6-member rings (excludes 5-member rings) a
    //! PAH has to have to be able to condense onto a particle (2 or more
    //! PAHs).
    void setCondensationThreshold(int target);

    //! Return the minimum number of 6-member rings (excludes 5-member rings) a
    //! PAH has to have to be able to incept.
    int inceptionThreshold() const;

    //! Return the minimum number of 6-member rings (excludes 5-member rings) a
    //! PAH has to have to be able to condense onto a particle (2 or more
    //! PAHs).
    int condensationThreshold() const;

    void SetMode(const std::string &mode);
    const std::string &Mode() const;

    void SetInceptedPAH(const std::string &name);
    const PostProcessStartingStr &InceptedPAH() const;
    //bool IsPyreneInception() const;

    //! Activates writing of entire binary trees
    void SetWriteBinaryTrees(bool flag) {m_write_bintree = flag;}

    //! Return the binary tree writing status
    const bool WriteBinaryTrees() const {return m_write_bintree;}

    //! Sets the value of the binary tree coalescence threshold
    void SetBinTreeCoalThresh(double ct) {m_bintree_coalthresh = ct;}

    //! Return the coalescence threshold for bintree particle models.
    const double GetBinTreeCoalThresh() const {return m_bintree_coalthresh;}

    //! Sets the value of the binary tree coalescence threshold
    void SetFractDim(double df) {m_fract_dim = df;}

    //! Return the coalescence threshold for bintree particle models.
    const double GetFractDim() const {return m_fract_dim;}

    //! Sets the value of the free-molecular enhancement factor.
    void SetEnhancementFM(double efm) {m_efm = efm;}

    //! Returns the value of the free-molecular enhancement factor.
    const double GetEnhancementFM() const {return m_efm;}

    //! Activates tracking of the distance between the centres of primary particles.
    void setTrackPrimarySeparation(bool flag) {m_trackPrimarySeparation = flag;}

    //! Return the flag used to indicate whether to track the distance between the centres of primary particles.
    const bool getTrackPrimarySeparation() const {return m_trackPrimarySeparation;}

    //! Activates tracking of the coordinates of primary particles.
    void setTrackPrimaryCoordinates(bool flag) {m_trackPrimaryCoordinates = flag;}

    //! Return the flag used to indicate whether to track the coordinates of
    //! primary particles.
    const bool getTrackPrimaryCoordinates() const {return m_trackPrimaryCoordinates;}

    void setPostprocessingType(postprocessingType Type) {m_postprocessingType = Type;};

    postprocessingType Postprocessing(void) const {return m_postprocessingType;};

    //! Index for temperature gradient in gas phase interface
    void setTGradIndex(const EnvironmentInterface::PropertyIndex index) {m_TemperatureGradientIndex = index;}

    //! Index for mixture fraction gradient in gas phase interface
    void setMixFracGradIndex(const EnvironmentInterface::PropertyIndex index) {m_MixFracGradientIndex = index;}

    //! Index for mixture fraction diffusion coefficient in gas phase interface
    void setMixFracDiffIndex(const EnvironmentInterface::PropertyIndex index) {m_MixFracDiffusionIndex = index;}

    //! Index for mixture fraction Laplacian in gas phase interface
    void setMixFracLaplIndex(const EnvironmentInterface::PropertyIndex index) {m_MixFracLaplacianIndex = index;}

    //! Index for average molecular weight of mixture in gas phase interface
    void setAvgMolWtIndex(const EnvironmentInterface::PropertyIndex index) {m_AvgMolWtIndex = index;}

    //! Index for thermal conductivity of mixture in gas phase interface
    void setThermalConductIndex(const EnvironmentInterface::PropertyIndex index) {m_ThermalConductivityIndex = index;}

protected:
    // The species used to define the processes and the particles.
    const Sprog::SpeciesPtrVector *m_species;

    //! The value of each component is the amount of a particular kind of material in a particle
    CompPtrVector m_components;

    //! Tracked quantities are excluded from mass calculations
    TrackPtrVector m_trackers;

    //! The ID of the aggregation model used to create particles.
    AggModels::AggModelType m_aggmodel;

    //! The sintering model.
    mutable Processes::SinteringModel m_sint_model;

	//! Titania phase transformation (melting) model
	mutable Processes::MeltingModel m_melt_model;

    //! Calculate a phsyical particle diffusion coefficient from its drag
    double EinsteinDiffusionCoefficient(const Cell &sys, const Particle &sp) const;

    //! Calculate the drag coefficient for a particle using the Knudsen correction
    double KnudsenDragCoefficient(const Cell &sys, const Particle &sp) const;

    //! Calculate the drag coefficient for a particle in the free molecular regime
    double FreeMolDragCoefficient(const Cell &sys, const Particle &sp) const;

    //! Calculate the drag coefficient for a particle as constant times temperature
    double TemperatureDragCoefficient(const Cell &sys, const Particle &sp) const;

    //! Calculate the drag coefficient for a particle using the Li & Wang expressions
    double LiWangDragCoefficient(const Cell &sys, const Particle &sp) const;

    //! Calculate the drag coefficient for a particle similar to the Li & Wang expressions
    double LiWangPatDragCoefficient(const Cell &sys, const Particle &sp) const;

    //==== Collision integrals ==========================

    //! Modified reduced diameter for use in collision integral calculations(should perhaps be input)
    double collisionIntegralDiameter(const Cell &sys, const Particle &sp) const;

    //! Dimensionless temperature for use in collision  integral calculations(should perhaps be input)
    double collisionIntegralTemperature(const Cell &sys, const Particle &sp) const;

    //! Blending function between diffuse and specular scattering for collision integrals
    double accomodationFunction(const Cell &sys, const Particle &sp) const;

    //! Knudsen based average of specular and diffusive 1,1 integrals
    double Omega1_1_avg(const Cell &sys, const Particle &sp) const;

    //! Diffuse scattering 1,1 collision integral approximation
    double Omega1_1_diff(const double t_star, const double sigma_prime) const;

    //! Specular scattering 1,1 collision integral approximation
    double Omega1_1_spec(const double t_star, const double sigma_prime) const;


    //! Knudsen based average of specular and diffusive 1,2 integrals
    double Omega1_2_avg(const Cell &sys, const Particle &sp) const;

    //! Diffuse scattering 1,2 collision integral approximation
    double Omega1_2_diff(const double t_star, const double sigma_prime) const;

    //! Specular scattering 1,2 collision integral approximation
    double Omega1_2_spec(const double t_star, const double sigma_prime) const;

    // MEMORY MANAGEMENT.

    // Initialises the model to the default condition.  Used in
    // constructors.
    void init(void);

    // Clears the model from memory.
    void releaseMem(void);

private:
    //! Constant A from Knudsen drag coefficient expression
    double m_DragA;

    //! Constant B from Knudsen drag coefficient expression
    double m_DragB;

    //! Constant E from Knudsen drag coefficient expression
    double m_DragE;

    //! Drag expression to use
    DragType m_DragType;

    //! Diffusion expression to use
    DiffusionType m_DiffusionType;

    //! Advection expression to use
    AdvectionType m_AdvectionType;

    //! Thermphoresis expression to use
    ThermophoresisType m_ThermophoresisType;
    
    //! Three parameter for Abhjeet's collision efficiency model.
    double colliParaA, colliParaB, colliParaC;

    //! The minimum number of 6-member rings (excludes 5-member rings) a PAH
    //! has to have to be able to incept which depends on the particular mode
    //! of inception.
    int m_inceptionThreshold;

    //! The minimum number of 6-member rings (excludes 5-member rings) a PAH
    //! has to have to be able to condense onto a particle (2 or more PAHs).
    int m_condensationThreshold;

    //! Define three modes, collision efficience depends on the smaller, the bigger, the combined mass or reduced mass
    std::string m_mode;

    //! define a species that it will be transfered between gas and particle phase
    PostProcessStartingStr m_InceptedPAH;

    //! Add a flag to indicate reading/writing of full particle binary trees
    bool m_write_bintree;

    //! Coalescene threshold for multicomponent binary tree particle models
    double m_bintree_coalthresh;

    //! Fractal dimension (bintree particle model)
    double m_fract_dim;

    //! Free molecular enhancement factor (coag., cond. & incep.)
    double m_efm;

    //! Flag to indicate whether to track the distance between the centres of neighbouring primary particles.
    bool m_trackPrimarySeparation;

    //! Flag to indicate whether to track the coordinates of primary particles.
    bool m_trackPrimaryCoordinates;

    postprocessingType m_postprocessingType;

    //! Index for temperature gradient
    EnvironmentInterface::PropertyIndex m_TemperatureGradientIndex;

    //! Index for mixture fraction gradient
    EnvironmentInterface::PropertyIndex m_MixFracGradientIndex;

    //! Index for mixture fraction diffusion coefficient
    EnvironmentInterface::PropertyIndex m_MixFracDiffusionIndex;

    //! Index for mixture fraction Laplacian
    EnvironmentInterface::PropertyIndex m_MixFracLaplacianIndex;

    //! Index for average molecular weight of mixture
    EnvironmentInterface::PropertyIndex m_AvgMolWtIndex;

    //! Index for thermal conductivity of mixture
    EnvironmentInterface::PropertyIndex m_ThermalConductivityIndex;

};
} //namespace Sweep
#endif
