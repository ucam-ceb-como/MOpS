/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    An ensemble for Sweep.
*/

#ifndef SWEEP_ENSEMBLE_H
#define SWEEP_ENSEMBLE_H

#include "swpparams.h"
#include "swpparticle1d.h"
#include <cmath>
#include <fstream>

namespace Sweep
{
class Ensemble
{
public:
    typedef DefaultParticle particle_type;
    typedef vector<DefaultParticle*>::iterator iterator;
    typedef vector<DefaultParticle*>::const_iterator const_iterator;
protected: // Ensemble capacity variables.
    unsigned int m_levels;   // Number of levels in the binary tree.
    unsigned int m_capacity; // The ensemble capacity (max. particle count).
    unsigned int m_halfcap;  // Half the ensemble capacity.
protected: // Ensemble scaling variables.
    real m_scale;            // Scaling factor due to internal operations (doubling etc.).
    real m_contfactor;       // Contraction scaling factor, precalculated for speed.
    unsigned int m_ncont;    // Number of ensemble contractions that have occured.
protected: // Doubling algorithm variables.
    unsigned int m_ndble;      // Number of ensemble doublings that have occured.
    bool m_dbleactive;         // Is doubling active or not.  Condition is reaching cutoff (see below).
    unsigned int m_dblecutoff; // Min. particle count at which doubling is activated.
    unsigned int m_dblelimit;  // Particle count below which ensemble is doubled (if active).
    unsigned int m_dbleslack;  // Slack space at end of ensemble after doubling operation.
    bool m_dbleon; // Allows user to manually switch off/on doubling.  Does not affect activation criterion.
public:
    /* List of particles in the ensemble. */
    vector<DefaultParticle*> m_particles;
protected:
    /* A node in the binary tree. */
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
    /* The binary tree nodes. */
    vector<NODE> m_tree;
public:
    Ensemble(void);
    Ensemble(const unsigned int count, const unsigned int nprop);
    ~Ensemble(void);
    /* Initialises the ensemble with the given capacity. */
    void Initialise(const unsigned int capacity, const unsigned int nprop);
    /* Clears the ensemble from memory. */
    void Destroy(void);
public:
    /* Returns a reference to the particle at index i. */
    DefaultParticle *GetParticle(const unsigned int i) const;
    /* Adds the given particle to the ensemble. */
    int AddParticle(DefaultParticle &sp);
    /* Removes the particle at the given index from the ensemble. */
    void RemoveParticle(const unsigned int i);
    void RemoveInvalids(void);
    /* Replaces the particle at the given index with the given particle. */
    void ReplaceParticle(const unsigned int i, DefaultParticle &sp);
    /* Removes all particles from the ensemble. */
    void Clear();
public:
    /* Returns an iterator to the first particle. */
    inline iterator begin() {return m_particles.begin();};
    inline const_iterator begin() const {return m_particles.begin();};
    inline iterator end() {return m_particles.end();};
    inline const_iterator end() const {return m_particles.end();};
public:
    /* Select a particle uniformly from the ensemble. */
    int SelectParticle(void) const;
    /* Selects a particle randomly weighted by the given particle
       property index.  The particle properties are those stored in
       the binary tree. Returns particle index on success, otherwise
       negative. */
    int SelectParticle(const int wtid) const;
public:
    /* Returns the particle count. */
    inline unsigned int Count(void) const {return (unsigned int)m_particles.size();};
    /* Returns the ensemble capacity. */
    inline unsigned int Capacity(void) const {return m_capacity;};
public:
    /* Returns the scaling factor due to particle operations. */
    real Scaling() const;
    /* Resets the scaling parameters. */
    void ResetScaling();
    /* Stops doubling algorithm. */
    inline void FreezeDoubling() {m_dbleon = false;};
    /* Restarts doubling if it was off, and checks if the ensemble should be doubled. */
    inline void UnfreezeDoubling() {m_dbleon=true; Double();};
public:
    /* Returns the sums of all properties in the binary tree. */
    void GetSums(vector<real> &sums) const;
    /* Returns the sum of one particle property with the given index from the binary tree. */
    real GetSum(unsigned int i) const;
public:
    /* Updates the ensemble tree from the given particle index. */
    void Update(const unsigned int i);
    /* Updates the ensemble tree completely. */
    void Update(void);
protected:
    /* Returns the index of the lowest binary tree node under which
       the given particle (by index) lies. */
    inline int TreeIndex(const unsigned int i) const {return m_halfcap + (i/2) - 1;};
    /* Returns true if particle index is a left branch, false if it is a right one. */
    static inline bool IsLeftBranch(const unsigned int i) {return (unsigned int)fmod((double)i,2) == 0;};
    /* Recalculates a branch of the tree from the given node upwards. */
    void AscendingRecalc(const unsigned int i);
    /* Performs the doubling algorithm. */
    void Double();
};
};

#endif