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

#endif
