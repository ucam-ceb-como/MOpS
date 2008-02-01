/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Inline function definitions for the Ensemble class.
*/

#ifndef SWEEP_ENSEMBLE_INL_H
#define SWEEP_ENSEMBLE_INL_H

class Ensemble; // Forward declaration of Ensemble class.

// STL STYLE ITERATORS.

inline Ensemble::iterator begin() {return m_particles.begin();};
inline Ensemble::const_iterator begin() const {return m_particles.begin();};
inline Ensemble::iterator end() {return m_particles.end();};
inline Ensemble::const_iterator end() const {return m_particles.end();};

// ENSEMBLE PROPERTIES.

inline unsigned int Ensemble::Count(void) const {return (unsigned int)m_particles.size();};
inline unsigned int Ensemble::Capacity(void) const {return m_capacity;};

// BINARY TREE.

inline int Ensemble::TreeIndex(const unsigned int i) const {return m_halfcap + (i/2) - 1;};
static inline bool Ensemble::IsLeftBranch(const unsigned int i) {return (unsigned int)fmod((double)i,2) == 0;};


#endif
