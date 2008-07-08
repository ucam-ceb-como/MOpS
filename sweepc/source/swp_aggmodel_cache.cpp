/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the AggModelCache class declared in the
    swp_aggmodel_cache.h header file.

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

#include "swp_aggmodel_cache.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
AggModelCache::AggModelCache(void)
{
    m_parent = NULL;
}

// Default constructor (public).
AggModelCache::AggModelCache(ParticleCache &parent)
{
    m_parent = &parent;
}

// Copy constructor.
AggModelCache::AggModelCache(const AggModelCache &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
AggModelCache::AggModelCache(std::istream &in, ParticleCache &parent)
{
    m_parent = &parent;
}

// Default destructor.
AggModelCache::~AggModelCache()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator (AggModelCache RHS).
AggModelCache &AggModelCache::operator=(const AggModelCache &rhs)
{
    return *this;
}

/*
// Assignment operator (Primary RHS).
AggModelCache &AggModelCache::operator=(const Primary &rhs)
{
    return (*this=rhs.TypedRef());
}
*/

// Compound assignment operator (AggModelCache RHS).
AggModelCache &AggModelCache::operator+=(const AggModelCache &rhs)
{
    return *this;
}

/*
// Compound assignment operator (Primary RHS).
AggModelCache &AggModelCache::operator+=(const Primary &rhs)
{
    return (*this+=rhs.TypedRef());
}
*/


// PARENT PARTICLE-CACHE.

// Returns a pointer to the parent particle data.
ParticleCache *const AggModelCache::Parent(void) const
{
    return m_parent;
}

// Sets the parent particle data.
void AggModelCache::SetParent(ParticleCache &parent)
{
    m_parent = &parent;
}
