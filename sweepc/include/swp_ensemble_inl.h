/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Inline function definitions for the Ensemble class.
*/

#ifndef SWEEP_ENSEMBLE_INL_H
#define SWEEP_ENSEMBLE_INL_H

#include "swp_ensemble.h"

// STL STYLE ITERATORS.

inline Sweep::Ensemble::iterator Ensemble::begin() {return m_particles.begin();};
inline Sweep::Ensemble::const_iterator Ensemble::begin() const {return m_particles.begin();};
inline Sweep::Ensemble::iterator Ensemble::end() {return m_particles.begin()+m_count;};
inline Sweep::Ensemble::const_iterator Ensemble::end() const {return m_particles.begin()+m_count;};

// ENSEMBLE PROPERTIES.

inline unsigned int Sweep::Ensemble::Count(void) const {return (unsigned int)m_count;};
inline unsigned int Sweep::Ensemble::Capacity(void) const {return m_capacity;};

// SCALING AND PARTICLE DOUBLING.

// Stops doubling algorithm.
inline void Sweep::Ensemble::FreezeDoubling() {m_dbleon = false;};

// Restarts doubling if it was off, and checks if the 
// ensemble should be doubled.
inline void Sweep::Ensemble::UnfreezeDoubling() {m_dbleon=true; dble();};


// BINARY TREE.

inline int Sweep::Ensemble::treeIndex(unsigned int i) const {
    return m_halfcap + (i/2) - 1;
};

inline bool Sweep::Ensemble::isLeftBranch(unsigned int i) {
    return (unsigned int)fmod((double)i,2) == 0;
};

#endif
