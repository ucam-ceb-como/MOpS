/*
  Author(s):      Peter Man (plwm2)
  Project:        sweep (population balance solver).

  File purpose:
    This file contains the definition of the TreeNode struct.
	A binary tree is a vector of nodes.
*/

#ifndef SWEEP_TREENODE_H
#define SWEEP_TREENODE_H

#include "swp_particle_cache.h"
#include <vector>

namespace Sweep
{
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
};

#endif
