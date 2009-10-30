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

inline int Sweep::Ensemble::treeIndex(unsigned int i) const {
    return m_halfcap + (i/2) - 1;
};

/**
 * Particle at index i is under a left hand tree leaf iff i is even.
 *
 *@param[in]    i   Index of particle in particle list
 *
 *@return       True if the particle is to the left of its tree leaf
 */
inline bool Sweep::Ensemble::isLeftBranch(unsigned int i) {
    // i%2 will be 0 if 1 is even and 1 otherwise
    return !(i%2);
};

/**
 * Sum up all the tree nodes assuming that all the leaves have been
 * initialised.  This will set the values of tree elements in the
 * range [0, m_halfcap - 1) (note half open range).
 */
inline void Sweep::Ensemble::recalcAllNonLeaf() {

    // The top (root) node does not have a parent so handle separately
    for(size_t treeIndex = m_halfcap - 2; treeIndex > 0; --treeIndex)
    {
        // Define node as an alias for m_tree[treeIndex]
        TreeNode& node = m_tree[treeIndex];
        
        // Pointer to parent node
        node.Parent = &m_tree[(treeIndex - 1) / 2];
        
        // Node below and to the left
        TreeNode& left = m_tree[treeIndex * 2 + 1];
        node.Left = &left;
        node.LeftData = left.LeftData + left.RightData;
        
        // Node below and to the right
        TreeNode& right = m_tree[treeIndex * 2 + 2];
        node.Right = &right;
        node.RightData = right.LeftData + right.RightData;
    }
    
    // Root node has not parent, but is otherwise the same
    TreeNode& node = m_tree[0];
    node.Parent = NULL;

    TreeNode& left = m_tree[1];
    node.Left = &left;
    node.LeftData = left.LeftData + left.RightData;

    TreeNode& right = m_tree[2];
    node.Right = &right;
    node.RightData = right.LeftData + right.RightData;

    // Cache the overall sums
    m_sums = node.LeftData + node.RightData;

}

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
    // Reserve memory for ensemble.
    m_tree.resize(m_capacity-1, TreeNode(*m_model));
    
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
    
    // First particle sits under the left of tree node m_halfcap - 1
    size_t treeIndex = m_halfcap - 1;
    bool isLeft = true;

    PartPtrVector::const_iterator it = m_particles.begin();
    const PartPtrVector::const_iterator itEnd = m_particles.end();
    
    // Loop through the particles putting their data into the leaves of the tree
    // and initialising the pointers that connect tree nodes
    while(it != itEnd)
    {
        if(isLeft) {
            // Set the parent pointer the first time we get to this node
            // does not have to be done again when the right hand side
            // of the node is set.
            m_tree[treeIndex].Parent = &m_tree[(treeIndex - 1) / 2];

            // Put the particle data into the tree
            m_tree[treeIndex].LeftData = **it;
            
            // This is a leaf node so it does not have a left child
            m_tree[treeIndex].Left = NULL;
            
            // The next particle will be on the right of the current leaf (treeIndex not changed) 
            isLeft = false;
        }
        else {
            // Put the particle data into the tree
            m_tree[treeIndex].RightData = **it;
            
            // This is a leaf node so it does not have a right child
            m_tree[treeIndex].Right = NULL;
            
            // The next particle will be on the left of the next leaf
            ++treeIndex;
            isLeft = true;
        }

        // Move on to next particle in input list
        ++it;
    }
    
    if(!isLeft) {
        // Last inserted particle was to the left of m_tree[treeIndex] so need
        // to initialise the right hand part of m_tree[treeIndex]
        m_tree[treeIndex].RightData.Clear();
        m_tree[treeIndex].Right = NULL;
        ++treeIndex;
    }
    
    // Finally intialise any remaining leaves of the tree to show they have no particles
    while(treeIndex < m_capacity - 1)
    {
        TreeNode& node = m_tree[treeIndex];
        node.LeftData.Clear();
        node.Left = NULL;
        node.RightData.Clear();
        node.Right = NULL;
        node.Parent = &m_tree[(treeIndex - 1) / 2];
        ++treeIndex;
    }
    // All the leaves of the tree should now be initialised so sum
    // up the rest of the tree
    recalcAllNonLeaf();
    
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
    m_dblecutoff = (int)(3.0 * (real)m_capacity / 4.0);
    m_dblelimit  = m_halfcap - (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0));
    m_dbleslack  = (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0));
     
}

#endif
