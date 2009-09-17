/*
  Author(s):      Peter Man (plwm2)
  Project:        sweep (population balance solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This file contains the definition of the TreeNode struct.
	A binary tree is a vector of nodes.

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

#ifndef SWEEP_TREENODE_H
#define SWEEP_TREENODE_H

#include "swp_particle_cache.h"
#include <vector>

namespace Sweep
{
    // Forward declaration
    class ParticleModel;

//! Building block for binary tree
class TreeNode
{
friend class Ensemble;

public:
	// Constructors
    TreeNode(const ParticleModel &model); // Initialising constructor.
	TreeNode(const TreeNode &copy);       // Copy-constructor.

    // Destructor.
	~TreeNode(void);

	// Operators
	TreeNode & operator=(const TreeNode &rhs);
	TreeNode & operator+=(const TreeNode &rhs);
	const TreeNode operator+(const TreeNode &rhs) const;

    // Clears the node to a default empty state.
	void Clear(void);

private:
	// MEMBER VARIABLES.
    ParticleCache LeftData;  // Sum of the left child leaves.
    ParticleCache RightData; // Sum of the right child leaves.
    TreeNode *Left;          // Pointer to left child node.
    TreeNode *Right;         // Pointer to right child node.
    TreeNode *Parent;        // Pointer to parent node.

    // Default Constructor is private to prevent uninitialised nodes
    // being created.
    TreeNode(void);
};

typedef std::vector<TreeNode> Tree;
} //namespace Sweep

#endif
