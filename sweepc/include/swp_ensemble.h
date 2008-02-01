/*
  Author(s):      Matthew Celnik (msc37) and Peter Man (plwm2)
  Project:        sweep (population balance solver)

  File purpose:
    A particle ensemble for Sweep.  The sweep particle ensemble uses a variable
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
*/

#ifndef SWEEP_ENSEMBLE_H
#define SWEEP_ENSEMBLE_H

#include "swpparams.h"
#include "swpparticle1d.h"
#include "swptreenode.h"
#include <cmath>

namespace Sweep
{
class Ensemble
{
public:
	/* typedefs */
    typedef Particle particle_type;
    typedef vector<Particle*>::iterator iterator;                // iterator through vector<Particle*>.
    typedef vector<Particle*>::const_iterator const_iterator;    // const_iterator through vector<Particle*>.

    vector<Particle*> m_particles;                               // List of particles in the ensemble.

	/* CONSTRUCTORS AND DESTRUCTOR */
	Ensemble(void);                                                  // Default constructor.
    Ensemble(const unsigned int count, const unsigned int nprop);    // Parameterised constructor.
    ~Ensemble(void);                                                 // Destructor.

    /* Initialises the ensemble with the given capacity. */
    void Initialise(const unsigned int capacity, const unsigned int nprop);  // The only function in the 
	                                                                         // parameterised constructor.
    /* Clears the ensemble from memory. */
    void Destroy(void);                                                      // The only function in the destructor.


	/* PARTICLE FUNCTIONS - ADD,REMOVE,CLEAR etc. */
    Particle *GetParticle(const unsigned int i) const;         // Returns a reference to the particle at index i.
    int AddParticle(Particle &sp);                             // Adds the given particle to the ensemble.
    void RemoveParticle(const unsigned int i);                 /* Removes the particle at the given index from the 
	                                                           ensemble. */
    void RemoveInvalids(void);                                 // Remove invalid particles.
    void ReplaceParticle(const unsigned int i, Particle &sp);  /* Replaces the particle at the given index with 
	                                                           the given particle. */
    void Clear();                                              // Removes all particles from the ensemble.

	/* ITERATOR FUNCTIONS */
    inline iterator begin() {return m_particles.begin();};              // Returns an iterator to the first particle.
    inline const_iterator begin() const {return m_particles.begin();};  // Returns const_iterator to the first particle.
    inline iterator end() {return m_particles.end();};                  // Returns an iterator to the last particle.
    inline const_iterator end() const {return m_particles.end();};      // Returns const_iterator to the last particle.

	/* SELECTING PARTICLES */
    int SelectParticle(void) const;              // Select a particle uniformly from the ensemble.
    int SelectParticle(const int wtid) const;    /* Randomly selects a particle, weighted by the given particle
                                                 property index.  The particle properties are those stored in
                                                 the binary tree. Returns particle index on success, otherwise
                                                 negative. */
	int SelectParticle(const Model * mod, const int mod_id) const;   /* Selects particle according to the particle property
	                                                                 specified by the given model and the given property id
																	 of the model */

	/* ENSEMBLE CAPACITY AND PARTICLE COUNT */
    inline unsigned int Count(void) const {return (unsigned int)m_particles.size();};  // Returns the particle count.
    inline unsigned int Capacity(void) const {return m_capacity;};                     // Returns the ensemble capacity.

	/* SCALING AND PARTICLE DOUBLING */
    real Scaling() const;                                      // Returns the scaling factor due to particle operations.
    void ResetScaling();                                       // Resets the scaling parameters.
    inline void FreezeDoubling() {m_dbleon = false;};          // Stops doubling algorithm.
    inline void UnfreezeDoubling() {m_dbleon=true; Double();}; /* Restarts doubling if it was off, and checks if the 
	                                                           ensemble should be doubled. */

	/* GET SUMS OF PROPERTIES */
    void GetSums(vector<real> &sums) const;  // Returns the sums of all properties in the binary tree.
    real GetSum(unsigned int i) const;       /* Returns the sum of one particle property with the given index 
	                                         from the binary tree. */

	/* UPDATE ENSEMBLE */
    void Update(const unsigned int i);    // Updates the ensemble tree from the given particle index.
    void Update(void);                    // Updates the ensemble tree completely.


protected: 
	
	/* ENSEMBLE CAPACITY VARIABLES */
    unsigned int m_levels;   // Number of levels in the binary tree.
    unsigned int m_capacity; // The ensemble capacity (max. particle count).
    unsigned int m_halfcap;  // Half the ensemble capacity.

    /* ENSEMBLE SCALING VARIABLES */
    real m_scale;            // Scaling factor due to internal operations (doubling etc.).
    real m_contfactor;       // Contraction scaling factor, precalculated for speed.
    unsigned int m_ncont;    // Number of ensemble contractions that have occurred.

	/* DOUBLING ALGORITHM VARIABLES */
    unsigned int m_ndble;      // Number of ensemble doublings that have occurred.
    bool m_dbleactive;         // Is doubling active or not.  Condition is to reach cutoff (see below).
    unsigned int m_dblecutoff; // Min. particle count at which doubling is activated.
    unsigned int m_dblelimit;  // Particle count below which ensemble is doubled (if active).
    unsigned int m_dbleslack;  // Slack space at end of ensemble after doubling operation.
    bool m_dbleon; // Allows user to manually switch off/on doubling.  Does not affect activation criterion.

	/* TREENODE AND TREE */
    struct TreeNode;
    vector<TreeNode> m_tree;   // The binary tree nodes.

	/* TREENODE ASSOCIATED FUNCTIONS AND VARIABLES */
    inline int TreeIndex(const unsigned int i) const // Returns the index of the lowest 
	{return m_halfcap + (i/2) - 1;};				 // binary tree node under which
                                                     // the given particle (by index) lies.

    /* Returns true if particle index is a left branch, false if it is a right one. */
    static inline bool IsLeftBranch(const unsigned int i) {return (unsigned int)fmod((double)i,2) == 0;}; 

    void AscendingRecalc(const unsigned int i);  // Recalculates a branch of the tree from the given node upwards.
    void Double();                               // Performs the doubling algorithm.
};

#include "swpensemble_inl.h"
};

#endif
