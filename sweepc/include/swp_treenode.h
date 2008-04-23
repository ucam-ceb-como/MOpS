/*
  Author(s):      Peter Man (plwm2)
  Project:        sweep (population balance solver).

  File purpose:
    This file contains the definition of the TreeNode struct.
	A binary tree is a vector of nodes.
*/

#ifndef SWEEP_TREENODE_H
#define SWEEP_TREENODE_H

#include "swp_particledata.h"
#include <vector>

namespace Sweep
{
class TreeNode
{
friend class Ensemble;

public:
	// CONSTRUCTORS AND DESTRUCTOR.
	TreeNode(const TreeNode &copy); // Copy-constructor.
    TreeNode(                       // Initialising constructor.
        const CompPtrVector &comp,  //  - Components.
        const TrackPtrVector &track //  - Tracker variables.
        );
	~TreeNode(void);                // Destructor.

	// OPERATOR OVERLOADING.
	TreeNode & operator=(const TreeNode &rhs);
	TreeNode & operator+=(const TreeNode &rhs);
	const TreeNode operator+(const TreeNode &rhs) const;

	// OTHER MEMBER FUNCTIONS.

    // Clear the node.
	void Clear(void);

private:
	// MEMBER VARIABLES.
	ParticleData LeftData;  // Sum of the left child leaves.
    ParticleData RightData; // Sum of right child leaves.
	TreeNode *Left;         // Pointer to left child node.
	TreeNode *Right;        // Pointer to right child node.
	TreeNode *Parent;       // Pointer to parent node.

    // Default Constructor is private to prevent uninitialised nodes
    // being created.
	TreeNode(void);
};

typedef std::vector<TreeNode> Tree;
};

#endif
