/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the TreeNode class declared in the
    swp_treenode.h header file.

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

#include "swp_treenode.h"

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default Constructor (private).
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

// Initialising constructor.
TreeNode::TreeNode(const Sweep::ParticleModel &model)
: LeftData(0.0, model), RightData(0.0, model),
  Left(NULL), Right(NULL), Parent(NULL)
{
}

// Destructor.
TreeNode::~TreeNode(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING

// Assignment operator ignores the pointer members
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

