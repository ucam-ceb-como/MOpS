/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a 2D surface-volume particle used by sweep which includes
    a description of PAHs in the particle structure.  This is a rather restrictive
    definition of a stochastic particle and requires that certain conditions are
    met when defining the mechanism.  In particular the mechanism must define
    some value names, which will be used to track PAH structure.
*/

#ifndef SWEEP_PAHPARTICLE_H
#define SWEEP_PAHPARTICLE_H

#include <vector>
#include "swpsvparticle.h"

namespace Sweep
{
class PAHParticle : public Sweep::SVParticle
{
protected:
    // Addresses of PAH structural info in value array.
    real *m_ed, *m_zz, *m_r5, *m_ac, *m_hl;
public: // Default constructor and destructor.
    PAHParticle(void);
    ~PAHParticle(void);
    /* Initialises the particle with the given components and sets up PAH model. */
    virtual void Initialise(vector<Component*> &components, const unsigned int nvals);
public:
    /* Returns number of PAH edges in particle. */
    inline real PAHEdgeCount(void) const {return *m_ed;};
    /* Returns number of PAH zig-zags in particle. */
    inline real PAHZigZagCount(void) const {return *m_zz;};
    /* Returns number of PAH R5 rings in zig-zags in particle. */
    inline real PAHR5Count(void) const {return *m_r5;};
    /* Returns number of PAH armchairs in particle. */
    inline real PAHArmchairCount(void) const {return *m_ac;};
    /* Returns number of PAH holes in particle. */
    inline real PAHHoleCount(void) const {return *m_hl;};
public:
    /* Updates the particle with n R6s adding to edges.  Note if there
       are insufficient edges to perform this growth then n is truncated. */
    void PerformEdgeRingGrowth(unsigned int &n);
    /* Updates the particle with n rings adding to armchairs.  Note if there
       are insufficient armchairs to perform this growth then n is truncated. */
    void PerformArmchairRingGrowth(unsigned int &n);
    /* Updates the particle with n R5s adding to zig-zags.  Note if there
       are insufficient zig-zags to perform this process then n is truncated. */
    void PerformR5Addition(unsigned int &n);
    /* Updates the particle with n R5s converting to R6s at edges.  Note if there
       are insufficient R5s or edges to perform this process then n is truncated. */
    void PerformR5EdgeConversion(unsigned int &n);
    /* Updates the particle with n R5s converting to R6s at armchair.  Note if there
       are insufficient R5s or armchairs to perform this process then n is truncated. */
    void PerformR5ArmchairConversion(unsigned int &n);
    /* Updates the particle with n edge oxidations.  Note if there are insufficient
       edges to perform this process then n is truncated. */
    void PerformEdgeOxidation(unsigned int &n);
protected:
    /* Increments n PAH sites according to rules. */
    void IncrementSites(unsigned int n);
    /* Decrements n PAH sites according to rules. */
    void DecrementSites(unsigned int n);
};
}

#endif