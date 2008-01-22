/*
  Author(s):      Matthew Celnik (msc37)
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
#include <cmath>

namespace Sweep
{
class Ensemble
{
public: 
    // The type of particle held in the ensemble.
    typedef DefaultParticle ptype;
    typedef vector<ptype*>::iterator iterator;
    typedef vector<ptype*>::const_iterator const_iterator;

    // Constructors.
    Ensemble(void); // Default constructor.
    Ensemble(const unsigned int count, const unsigned int nprop); // Initialising constructor.
    Ensemble(const Ensemble &ens); // Copy constructor.

    // Destructor.
    virtual ~Ensemble(void);

    // Initialisation and destruction.
    void Initialise(const unsigned int capacity, const unsigned int nprop); // Initialises the ensemble with the given capacity.
    void Destroy(void); // Clears the ensemble from memory.

    // STL style iterators.
    inline iterator begin();             // Returns an iterator to the first particle.
    inline const_iterator begin() const; // Returns a const iterator to the first particle.
    inline iterator end();               // Returns an iterator to the last particle.
    inline const_iterator end() const;   // Returns a const iterator to the last particle.

    // Particle access (get, add and replace).
    ptype *const GetParticle(const unsigned int i) const; // Returns a reference to the particle at index i.
    int AddParticle(ptype &sp); // Adds the given particle to the ensemble.
    void ReplaceParticle(const unsigned int i, ptype &sp); // Replaces the particle at the given index with the given particle.

    // Particle removal.
    void RemoveParticle(const unsigned int i); // Removes the particle at the given index from the ensemble.
    void RemoveInvalids(void);                 // Removes particles which are invalid from the ensemble.
    void Clear();                              // Removes all particles from the ensemble.

    // Particle selection.   
    int SelectParticle(void) const; // Select a particle uniformly from the ensemble.
    int SelectParticle(const int wtid) const; // Selects a particle randomly weighted by the given particle
                                              // property index.  The particle properties are those stored in
                                              // the binary tree. Returns particle index on success, otherwise
                                              // negative.

    // Ensemble properties.
    inline unsigned int Count(void) const;    // Returns the particle count.
    inline unsigned int Capacity(void) const; // Returns the ensemble capacity.

    // Scaling parameters due to ensemble resizing algorithms.
    real Scaling() const; // Returns the scaling factor due to particle operations.
    void ResetScaling();  // Resets the scaling parameters.

    // Particle doubling algorithm.
    inline void FreezeDoubling() {m_dbleon = false;}; // Stops doubling algorithm.
    inline void UnfreezeDoubling() {m_dbleon=true; Double();}; // Restarts doubling if it was off, and 
                                                               // checks if the ensemble should be doubled.

    // Summed particle properties.   
    void GetSums(vector<real> &sums) const; // Returns the sums of all properties in the binary tree.   
    real GetSum(unsigned int i) const;      // Returns the sum of one particle property with the
                                            // given index from the binary tree.

    // Ensemble update, to be called if a particle has changed.
    void Update(const unsigned int i); // Updates the ensemble tree from the given particle index.   
    void Update(void);                 // Updates the ensemble tree completely.

protected:
    // Ensemble capacity variables.
    unsigned int m_levels;   // Number of levels in the binary tree.
    unsigned int m_capacity; // The ensemble capacity (max. particle count).
    unsigned int m_halfcap;  // Half the ensemble capacity.
    // Ensemble scaling variables.
    real m_scale;            // Scaling factor due to internal operations (doubling etc.).
    real m_contfactor;       // Contraction scaling factor, precalculated for speed.
    unsigned int m_ncont;    // Number of ensemble contractions that have occured.
    // Doubling algorithm variables.
    unsigned int m_ndble;      // Number of ensemble doublings that have occured.
    bool m_dbleactive;         // Is doubling active or not.  Condition is reaching cutoff (see below).
    unsigned int m_dblecutoff; // Min. particle count at which doubling is activated.
    unsigned int m_dblelimit;  // Particle count below which ensemble is doubled (if active).
    unsigned int m_dbleslack;  // Slack space at end of ensemble after doubling operation.
    bool m_dbleon; // Allows user to manually switch off/on doubling.  Does not affect activation criterion.

    // List of particles in the ensemble.
    vector<DefaultParticle*> m_particles;

    // A node in the binary tree.
    struct NODE {
        vector<real> LeftSum, RightSum;
        NODE *Left;
        NODE *Right;
        NODE *Parent;
        NODE(void);
        NODE &operator=(const NODE &n);
        NODE &operator+=(const NODE &n);
        const NODE operator+(const NODE &n) const;
        void SetSize(const unsigned int size);
        void Clear();
    };

    vector<NODE> m_tree; // The binary tree nodes.
   
    // Binary tree.
    inline int TreeIndex(const unsigned int i) const; // Returns the index of the lowest binary tree node under which
                                                      // the given particle (by index) lies.
    static inline bool IsLeftBranch(const unsigned int i); // Returns true if particle index is a left branch, false if it is a right one.
    void AscendingRecalc(const unsigned int i); // Recalculates a branch of the tree from the given node upwards.

    // Ensemble scaling algorithms.
    void Double(); // Performs the doubling algorithm.
};

#include "swpensemble_inl.h"
};

#endif