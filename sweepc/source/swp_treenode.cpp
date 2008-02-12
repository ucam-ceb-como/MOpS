#include "swp_treenode.h"

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default Constructor.
TreeNode::TreeNode(void)
: Left(NULL), Right(NULL), Parent(NULL)
{
}

// Copy-constructor.
TreeNode::TreeNode(const TreeNode & tn)
{
    // Use assignment operator.
	*this = tn;
}

// Destructor.
TreeNode::~TreeNode(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING

// Assignment operator.
TreeNode &TreeNode::operator=(const TreeNode &rhs)
{
    if (this != &rhs) {
        /*
        Left      = rhs.Left;
        Right     = rhs.Right;
        Parent    = rhs.Parent;
        */
        LeftData  = rhs.LeftData;
        RightData = rhs.RightData;
    }
    return *this;
}

// Compound assignment operator.
TreeNode & TreeNode::operator+=(const TreeNode &rhs)
{
	LeftData  += rhs.LeftData;
	RightData += rhs.RightData;
    return *this;
}

// Overloading Addition operator.
const TreeNode TreeNode::operator+(const TreeNode &rhs) const
{
    return TreeNode(*this) += rhs;
}


// OTHER MEMBER FUNCTIONS

// Clear the node - set m_leftsum and m_rightsum ParticleCaches to 'zero'.
void TreeNode::Clear(void)
{
	LeftData.Clear();
	RightData.Clear();
}

