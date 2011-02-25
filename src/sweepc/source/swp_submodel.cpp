/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the SubModel class declared in the
    swp_submodel.h header file.

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

#include "swp_submodel.h"
#include "swp_primary.h"

using namespace Sweep;
using namespace Sweep::SubModels;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
SubModel::SubModel(void)
{
    m_parent = NULL;
}

// Initialising constructor.
SubModel::SubModel(Primary &parent)
{
    m_parent = &parent;
}

// Copy constructor (protected).
SubModel::SubModel(const SubModel &copy)
{
    *this = copy;
}

// Default destructor.
SubModel::~SubModel(void)
{
}


// OPERATOR OVERLOADS.

// Assignment operator.
SubModel &SubModel::operator=(const SubModel &rhs)
{
    if (this != &rhs) {
        // The parent primary particle is not set here.
    }
    return *this;
}


// PARENT PRIMARY PARTICLE.

// Returns the parent primary particle.
const Primary *const SubModel::Parent(void) const
{
    return m_parent;
}

// Sets the parent primary particle.
void SubModel::SetParent(Primary &parent)
{
    m_parent = &parent;
}
