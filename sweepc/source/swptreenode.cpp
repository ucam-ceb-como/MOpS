#include "swptreenode.h"

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default Constructor.
TreeNode::TreeNode(void)
{
    m_left = NULL;  m_right = NULL;  m_parent = NULL;
}

// Copy-constructor.
TreeNode::TreeNode(const TreeNode & tn)
{
	*this = tn;
}

// Destructor.
TreeNode::~TreeNode(void)
{
	m_left = NULL;  m_right = NULL;  m_parent = NULL;
}


// OPERATOR OVERLOADING

// Overloading Assignment operator.
TreeNode & TreeNode::operator=(const TreeNode &tn)
{
    if (this == &tn) return *this;
    m_leftsum = tn.m_leftsum;
    m_rightsum = tn.m_rightsum;
    return *this;
}

// Overloading Comparison operator.
bool TreeNode::operator==(const TreeNode &tn) const
{
	return (m_leftsum == tn.leftsum && m_rightsum == tn.rightsum);
}

// Overloading Inequality operator.
bool TreeNode::operator!=(const TreeNode &tn) const
{
	return !(*this == tn);
}

// Overloading += operator.
TreeNode & TreeNode::operator +=(const TreeNode &tn)
{
	m_leftsum += tn.m_leftsum;
	m_rightsum += tn.m_rightsum;
    return *this;
}

// Overloading Addition operator.
const TreeNode TreeNode::operator+(const TreeNode & tn) const
{
    TreeNode lhs = *this;
    lhs += tn;
    return lhs;
}


// OTHER MEMBER FUNCTIONS

// Resize the m_leftsum and m_rightsum ParticleCaches to the given 'size'.
void TreeNode::SetSize(const unsigned int size)
{
    m_leftsum.Resize(size);
	m_rightsum.Resize(size);
}

// Clear the node - set m_leftsum and m_rightsum ParticleCaches to 'zero'.
void TreeNode::ClearNode(void)
{
	m_leftsum.Clear();
	m_rightsum.Clear();
}

