/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the SubModelCache class declared in the
    swp_submodel_cache.h header file.

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

#include "swp_submodel_cache.h"

using namespace Sweep;
using namespace Sweep::SubModels;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
SubModels::SubModelCache::SubModelCache(void)
: m_parent(NULL)
{
}

// Default constructor (public).
SubModelCache::SubModelCache(Sweep::ParticleCache &parent)
{
    m_parent = &parent;
}

// Copy constructor.
SubModelCache::SubModelCache(const SubModelCache &copy)
    : m_parent(copy.m_parent)
{
}

// Default destructor.
SubModelCache::~SubModelCache()
{
    // Nothing special to destruct.
}

/*
// OPERATOR OVERLOADING.

// Assignment operator.
SubModelCache &SubModelCache::operator=(const SubModelCache &rhs)
{
    return *this;
}

// Compound assignment.
SubModelCache &SubModelCache::operator+=(const SubModelCache &rhs)
{
    return *this;
}
*/

// PARENT.

// Returns a pointer to the parent particle data.
ParticleCache *const SubModelCache::Parent(void) const
{
    return m_parent;
}

// Sets the parent particle data.
void SubModelCache::SetParent(ParticleCache &parent)
{
    m_parent = &parent;
}
