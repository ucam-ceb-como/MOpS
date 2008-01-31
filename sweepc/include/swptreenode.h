/*
  Author(s):      Peter Man (plwm2)
  Project:        sweep (population balance solver).

  File purpose:
    This file contains the definition of the TreeNode struct.
	A binary tree is a vector of nodes.
*/

#ifndef SWEEP_TREENODE_H
#define SWEEP_TREENODE_H

#include "swpparticlecache.h"

namespace Sweep
{
struct TreeNode
{
	// MEMBER VARIABLES.
	ParticleCache m_leftsum, m_rightsum;  // Sum of the left child leaves and the right child leaves respectively.
	TreeNode *m_left;      // Pointer to left child node.
	TreeNode *m_right;     // Pointer to right child node.
	TreeNode *m_parent;    // Pointer to parent node.

	// CONSTRUCTORS AND DESTRUCTOR.
	TreeNode(void);                 // Default Constructor.
	TreeNode(const TreeNode & tn);  // Copy-constructor.
	~TreeNode(void);                // Destructor.

	// OPERATOR OVERLOADING.
	TreeNode & operator=(const TreeNode &tn);            // Overloading Assignment operator.
	bool operator==(const TreeNode &tn) const;           // Overloading Comparison operator.
	bool operator!=(const TreeNode &tn) const;           // Overloading Inequality operator.
	TreeNode & operator+=(const TreeNode &tn);           // Overloading += operator.
	const TreeNode operator+(const TreeNode &tn) const;  // Overloading Addition operator.

	// OTHER MEMBER FUNCTIONS.
	void SetSize(const unsigned int size);  // Resize the m_leftsum and m_rightsum ParticleCaches to the given 'size'.
	void ClearNode(void);                   // Clear the node - set m_leftsum and m_rightsum ParticleCaches to 'zero'.
};
};

#endif