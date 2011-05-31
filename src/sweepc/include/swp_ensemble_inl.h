/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Inline function definitions for the Ensemble class.

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

#ifndef SWEEP_ENSEMBLE_INL_H
#define SWEEP_ENSEMBLE_INL_H

#include <stdexcept>
#include <cassert>
#include <iostream>

#include "swp_ensemble.h"

// STL STYLE ITERATORS.
inline Sweep::Ensemble::iterator Sweep::Ensemble::begin() {
    return m_particles.begin();
};

inline Sweep::Ensemble::const_iterator Sweep::Ensemble::begin() const {
    return m_particles.begin();
};

inline Sweep::Ensemble::iterator Sweep::Ensemble::end() {
    return m_particles.begin()+m_count;
};

inline Sweep::Ensemble::const_iterator Sweep::Ensemble::end() const {
    return m_particles.begin()+m_count;
};


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


/**
 * Initialise the ensemble to hold particles of the type specified
 * by the model and containing the particular particles contained
 * in the range [particle_list_begin, particle_list_end).
 *@param[in]    particle_list_begin     Iterator to first in range of particle pointers to insert
 *@param[in]    particle_list_end       Iterator to one past end of range of particle pointers to insert
 */
template<class T> void Sweep::Ensemble::SetParticles(
    T particle_list_begin,
    T particle_list_end)
{
    // Clear any existing particles
    for(iterator it = m_particles.begin(); it != m_particles.end(); ++it) {
        delete *it;
    }
    // Read in the particles from the range
    m_particles.assign(particle_list_begin, particle_list_end);

    m_maxcount = m_count = m_particles.size();
    // By definition will not handle more particles than the ensemble capacity
    if(m_count > m_capacity)
    {
        // Would be nice to do something more sophisiticated here,
        // caller will need to clean up the particles if this
        // exception is thrown.
        throw std::runtime_error("Number of particles specified exceeds ensemble capacity (Sweep Ensemble::SetParticles)");
    }

    // Fill the particle list upto its full size, by filling any space with null pointers
    m_particles.resize(m_capacity, NULL);

    // Initialise scaling.
    m_ncont      = 0;
    m_scale      = 1.0;
    m_contfactor = (real)(m_capacity-1) / (real)(m_capacity);
    m_contwarn   = false;

    // Initialise doubling.
    m_ndble      = 0;
    m_dbleon     = true;
    m_dbleactive = false;
    m_dblecutoff = (3u * m_capacity) / 4u;

    // This is 2^(m_levels - 5) provided m_levels > 5, otherwise 1
    m_dbleslack  = 1u << ((m_levels > 5u) ? m_levels-5 : 0);
    m_dblelimit  = m_halfcap - m_dbleslack;

    // Build the tree with the weights for the new particles.
    rebuildTree();

    assert(m_tree.size() == m_count);
}


#endif
