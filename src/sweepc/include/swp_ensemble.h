/*
  Author(s):      Matthew Celnik (msc37) and Peter Man (plwm2)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    A particle ensemble for Sweep.  The sweep particle ensemble uses a variable
    particle count algorithm and employs some scaling techniques to facilitate this.

    As most particle systems simulated involve an initial growth in particle count,
    followed by a coagulation phase in which particle count decreases (sometimes
    substantially) a particle doubling algorithm is used.  This algorithm copies
    one half of the population when the ensemble count falls below half capacity,
    thereby preventing the removal of all particles and hence statistically significant
    systematic error.

    The ensemble also uses a contraction algorithm if the ensemble requires more space for
    particles than is available.  This algorithm uniformly removes a particle from the
    ensemble to provide more space.

    To speed up particle selection the ensemble implements a binary tree.  Properties
    from the particle type are stored in the tree allowing for rapid selection by
    different properties.  The actual properties stored are not defined by the ensemble.

    The binary tree requires that the ensemble capacity must be a power of 2.

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

#ifndef SWEEP_ENSEMBLE_H
#define SWEEP_ENSEMBLE_H

#include "swp_params.h"
#include "swp_particle.h"
#include "swp_tree_weighted_cache.h"
#include "swp_tree_transcoag_weighted_cache.h"
#include "swp_property_indices.h"
#include "swp_gas_profile.h"
#include "swp_kmc_pah_structure.h"

#include "binary_tree.hpp"

#include <list>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace Sweep
{
class Mechanism;
class ParticleModel;
namespace KMC_ARS{
class CSV_data;
class KMCSimulator;}

/*!
 * \brief Manages the particle population inside a cell
 *
 *     A particle ensemble for Sweep.  The sweep particle ensemble uses a variable
    particle count algorithm and employs some scaling techniques to facilitate this.

    As most particle systems simulated involve an initial growth in particle count,
    followed by a coagulation phase in which particle count decreases (sometimes
    substantially) a particle doubling algorithm is used.  This algorithm copies
    one half of the population when the ensemble count falls below half capacity,
    thereby preventing the removal of all particles and hence statistically significant
    systematic error.

    The emsemble also uses a contraction algorithm if the ensemble requires more space for
    particles than is available.  This algorithm uniformly removes a particle from the
    ensemble to provide more space.

    To speed up particle selection the ensemble implements a binary tree.  Properties
    from the particle type are stored in the tree allowing for rapid selection by
    different properties.  The actual properties stored are not defined by the ensemble.

    The binary tree requires that the ensemble capacity must be a power of 2.

    \todo It should not be necessary to include "particle.h" in this header file
    because it exposes too much detail of the implementation of particle, when all
    we really need to know is what a vector of pointers to particles looks like.
 */
class Ensemble
{
public:
    // TYPEDEFS.

    //! The type of particle in the ensemble.
    typedef Particle particle_type;

    //! Iterator for the particles.
    typedef PartPtrVector::iterator iterator;

    //! Constant iterator for the particles.
    typedef PartPtrVector::const_iterator const_iterator;

    //! Particle value cache for specifying distributions on the particle list
    typedef Sweep::TreeTransCoagWeightedCache particle_cache_type;

    // Constructors.
    Ensemble(void); // Default constructor.
    Ensemble(                             // Initialising constructor (incl. particles).
        unsigned int count               //  - Capacity (max. number of particles).
        );
    Ensemble(const Ensemble &copy); // Copy constructor.
    Ensemble(                            // Stream-reading constructor.
        std::istream &in,                //   - Input stream.
        const Sweep::ParticleModel &mech //   - Mechanism used to define particles.
        );

    // Destructor.
    ~Ensemble(void);

    // Operators.
    Ensemble &operator=(const Ensemble &rhs);

    //! Overload of the << operator
    friend std::ostream& operator<<(
            std::ostream &os,
            const Sweep::Ensemble &e);


    // INITIALISATION.

    // Initialises the ensemble with the given capacity.
    void Initialise(
        unsigned int capacity             // Max. number of particles
        );

    //! Initialise with some particles, downsampling as necessary
    void SetParticles(std::list<Particle*>::iterator first, std::list<Particle*>::iterator last,
                      rng_type &rng);

    //! Empty the tree and pass on ownership of the particles
    PartPtrList TakeParticles();

    void SetDoubling(const bool val);

	unsigned int DoubleLimit();

	bool IsDoublingOn();

    // PARTICLE ADDITION AND REMOVAL.

    // Returns a pointer to the particle at index i.
    Particle *const At(unsigned int i);
    const Particle *const At(unsigned int i) const;

    // Adds the given particle to the ensemble.  Returns the new
    // particle's index in the ensemble.  The ensemble then takes
    // control of destruction of the particle.
    // For particle-number/particle (hybrid) model: the final two inputs 
    // flag events in which we cannot delete the particle with index i2 
    // if a contraction is triggered
    int Add(Particle &sp, rng_type &rng, int i2 = 0, bool hybrid_event_flag = false);

	//Find a particle that is a single PAH of a given structure
	int CheckforPAH(Sweep::KMC_ARS::PAHStructure &m_PAH, double t, int ind);

    //! Removes the particle at the given index from the ensemble.
    void Remove(
        unsigned int i, // Index of particle to remove.
        bool fdel=true  // Set true to delete particle from memory as well, otherwise false.
        );

    //! Removes invalid particles.
    void RemoveInvalids(void);

    //! Replaces the particle at the given index with the given particle.
    void Replace (
        unsigned int i, // Index of particle to replace.
        Particle &sp    // Particle to insert.
        );

    // Removes all particles from the ensemble.
    void Clear();


    // ITERATOR FUNCTIONS.

    // Returns an iterator to the first particle.
    iterator begin();

    // Returns const_iterator to the first particle.
    const_iterator begin() const;

    // Returns an iterator to the last particle.
    iterator end();

    // Returns const_iterator to the last particle.
    const_iterator end() const;


    // SELECTING PARTICLES.

    // Select a particle uniformly from the ensemble and returns
    // its index. Returns negative on failure.
    int Select(rng_type &rng) const;

    // Randomly selects a particle, weighted by the given particle
    // property index.  The particle properties are those stored in
    // the ParticleData type. Returns particle index on success, otherwise
    // negative.
    int Select(Sweep::PropID id, rng_type &rng) const;

    // Randomly selects a particle using a pre-chosen random number weighted
    // by the given particle property index. 
    // The particle properties are those stored in
    // the ParticleData type. Returns particle index on success, otherwise
    // negative.
    int Select_usingGivenRand(Sweep::PropID id, double rng_number, rng_type &rng) const;

    // ENSEMBLE CAPACITY AND PARTICLE COUNT.

    //! Returns the particle count.
    unsigned int Count(void) const;

    //!Returns the ensemble capacity.
    unsigned int Capacity(void) const;


    // SCALING AND PARTICLE DOUBLING.

    //! Returns the scaling factor due to particle operations.
    double Scaling() const;

    //! Resets the scaling parameters.
    void ResetScaling();

    //! Stops doubling algorithm.
    inline void FreezeDoubling();

    //! Restarts doubling if it was off, and checks if the ensemble should be doubled.
    inline void UnfreezeDoubling();


    // GET SUMS OF PROPERTIES.

    //! Returns the sums over all particles of all their cached properties.
    const particle_cache_type &GetSums(void) const;

    // Returns the sum of one particle property with the given index
    // from the binary tree.
    double GetSum(
        Sweep::PropID id // ID of the ParticleData property.
        ) const;

    //! Inform the ensemble that the particle at index i has been changed
    void Update(unsigned int i);

    //! Get alpha for the ensemble (ABF model)
    double Alpha(double T) const;

    // Hybrid particle-number/particle model functions
    // ===============================================
    // Set/get threshold size 
    void SetHybridThreshold(unsigned int threshold) { m_hybrid_threshold = threshold; }
    unsigned int GetHybridThreshold() const { return m_hybrid_threshold; }
    
    // Functions to update/reset the number count at a specific particle-number index
    unsigned int SetTotalParticleNumber();
    void ResetNumberAtIndex(unsigned int index);
    void UpdateNumberAtIndex(unsigned int index, int update);
    void UpdateTotalParticleNumber(int update) { m_total_number += update; }
    void UpdateTotalsWithIndex(unsigned int index, double change);
    void UpdateTotalsWithIndices(unsigned int i1, unsigned int i2);

    // Recalculate property sums to avoid accummulation of rounding errors
    void RecalcPNPropertySums();

    // Functions to initialise properties
    void InitialiseParticleNumberModel();
    void InitialiseDiameters(double molecularWeight, double density);
    int SetPNParticle(Particle &sp, unsigned int index);

    // Functions to get properties (at specified index)	
    double PropertyAtIndex(Sweep::PropID prop, unsigned int index) const;
    double GetPropertyTotal (Sweep::PropID prop) const;
    unsigned int GetTotalParticleNumber() const { return m_total_number; }
    unsigned int NumberAtIndex(unsigned int index) const;
    double Diameter2AtIndex(unsigned int index) const { return m_pn_diameters2[index]; }
    double DiameterAtIndex(unsigned int index) const { return m_pn_diameters[index]; }
    double MassAtIndex(unsigned int index) const { return m_pn_mass[index]; }
    Particle *const GetPNParticleAt(unsigned int index);
	
    // This could be a single function with a case statement
    // but it would be slower and some propIDs don't exist
    double GetTotalDiameter() const;
    double GetTotalDiameter2() const;
    double GetTotalDiameter_1() const;
    double GetTotalDiameter_2() const;
    double GetTotalDiameter3() const;
    double GetTotalDiameter2_mass_1_2() const;
    double GetTotalMass_1_2() const;
    double GetTotalMass() const;
    double GetTotalMass2() const;
    double GetTotalMass3() const;
    unsigned int GetTotalComponent() const;

    // Function to double totals when doubling is triggered
    void DoubleTotals();

    // ===============================================

    // READ/WRITE/COPY.

    //! Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,                // Input stream.
        const Sweep::ParticleModel &mech // Model used to define particles.
        );

    int NumOfInceptedPAH(int ID) const;// return the number of pyrene in current state.
    int IndexOfInceptedPAH(int ID) const; //move backwards.
    Sweep::KMC_ARS::KMCSimulator* Simulator();
    void SetSimulator(Sweep::GasProfile& gp);

    // modify the m_numofInceptedPAH according to processes,
    // there are two possible value for m_amount, 1 (increase by one ) and -1 (decrease by 1)
    //void SetNumOfInceptedPAH(int m_amount);
    //void SetNumOfInceptedPAH(int m_amount, Sweep::AggModels::Primary *m_primary);

	// PARTICLE TRACKING OUTPUT FOR VIDEOS

	//! Updates tracking after a coagulation event
	void UpdateTracking(int p_old, int p_merged);

	//! Returns a pointer to the given tracked particle.
	const Particle *const TrackedAt(unsigned int i) const;
	
	//! Set number of particle tracked for videos
	void SetParticleTrackingNumber(unsigned int val);

	//! Return number of particles currently being tracked for videos
	unsigned int TrackedParticleNumber() const;

	//! Initialise tracking of initial particle population 
	void InitialiseParticleTracking();

private:
    //! Vector of particles in the ensemble.
    PartPtrVector m_particles;
    Sweep::KMC_ARS::KMCSimulator *m_kmcsimulator;

	// PARTICLE TRACKING OUTPUT FOR VIDEOS

	//! (maximum) number of particles tracked for videos
	unsigned int m_tracked_number;

	//! Vector of pointers of tracked particles 
	PartPtrVector m_tracked_particles;    

// ENSEMBLE CAPACITY VARIABLES.
    unsigned int m_levels;   // Number of levels in the binary tree.
    unsigned int m_capacity; // The ensemble capacity (max. particle count).
    unsigned int m_halfcap;  // Half the ensemble capacity.
    unsigned int m_count;    // Number of particles currently in the ensemble.
    //unsigned int m_numofInceptedPAH;  // Number of starting PAH in the ensemble

    // ENSEMBLE SCALING VARIABLES.
    double m_contfactor;       // Contraction scaling factor, precalculated for speed.
    unsigned int m_ncont;    // Number of ensemble contractions that have occurred.
    double m_wtdcontfctr;    // Track change due to loss of particles (needed to account for unequal weights or PN/P model).
    bool m_contwarn;         // Has a contraction warning msg been printed?

    // DOUBLING ALGORITHM VARIABLES.
    unsigned int m_maxcount;   // The maximum particle count reached by the ensemble.
    unsigned int m_ndble;      // Number of ensemble doublings that have occurred.
    bool m_dbleactive;         // Is doubling active or not.  Condition is to reach cutoff (see below).
    unsigned int m_dblecutoff; // Min. particle count at which doubling is activated.
    unsigned int m_dblelimit;  // Particle count below which ensemble is doubled (if active).
    unsigned int m_dbleslack;  // Slack space at end of ensemble after doubling operation.
    bool m_dbleon;             // Allows user to manually switch off/on doubling.  Does not affect activation criterion.

    // Hybrid particle number model variables
    // ===============================================
    unsigned int m_hybrid_threshold; 
    unsigned int m_total_number;
    unsigned int m_total_component;
    double m_total_diameter;
    double m_total_diameter2;
    double m_total_diameter_1;
    double m_total_diameter_2;
    double m_total_diameter2_mass_1_2;
    double m_total_mass_1_2;
    double m_total_mass;
    double m_total_mass2;
    double m_total_mass3;
    double m_total_diameter3;
    std::vector<unsigned int> m_particle_numbers;
    std::vector<double> m_pn_diameters;
    std::vector<double> m_pn_diameters2;
    std::vector<double> m_pn_diameters_1;
    std::vector<double> m_pn_diameters_2;
    std::vector<double> m_pn_diameters2_mass_1_2;
    std::vector<double> m_pn_mass_1_2;
    std::vector<double> m_pn_mass;
    std::vector<double> m_pn_mass2;
    std::vector<double> m_pn_mass3;
    std::vector<double> m_pn_diameters3;
    PartPtrVector m_pn_particles;

    // ===============================================

    //! Reset the contents of the binary tree
    void rebuildTree();

    //! Type of weight tree for particle selection and property summation
    typedef Utils::BinaryTree<particle_cache_type, PartPtrVector::iterator> tree_type;

    //! Tree for inverting probability distributions on the particles and summing their properties
    tree_type m_tree;


    // MEMORY MANAGEMENT.

    //! Releases all memory resources used by the ensemble.
    void releaseMem(void);

    //! Sets the ensemble to its initial condition.  Used in constructors.
    void init(void);

    //! Performs the doubling algorithm.
    void dble();

    //! Empties the main population
    void ClearMain();

    //! Functor to extract weights from nodes of the new binary tree
    class WeightExtractor : public std::unary_function<const particle_cache_type&, double>
    {
    public:
        //! Set up an extractor for an indexed property
        WeightExtractor(const Sweep::PropID id);

        //! Extract a weight from a cache
        double operator()(const particle_cache_type& cache) const;

    private:
        //! Id of property to extract
        Sweep::PropID mId;

        //! Not possible to have an extractor of this kind without specifying an id
        WeightExtractor();
    };
}; // end of class Ensemble

} // end of namespace Sweep

// Include inline function definitions.
#include "swp_ensemble_inl.h"

#endif
