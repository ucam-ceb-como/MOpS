/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The AggModelCache class is a base model cache for different primary
    particle aggregation models.

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

#ifndef SWEEP_AGGMODEL_CACHE_H
#define SWEEP_AGGMODEL_CACHE_H

#include "swp_params.h"
#include "swp_aggmodel_type.h"
#include <iostream>

namespace Sweep
{
// Forward declaration
class Primary;

namespace AggModels
{
/*!
 * \brief Interface for aggregation model specific properties
 */
class AggModelCache
{
public:
    // Destructors.
    virtual ~AggModelCache(void);

    // Operators.
    virtual AggModelCache& operator=(const Primary &rhs) = 0;
    virtual AggModelCache& operator+=(const Primary &rhs) = 0;

    // DATA MANAGEMENT.

    // Resets the cache to its "empty" condition.
    virtual void Clear(void) = 0;


    // READ/WRITE/COPY.

    // Creates a copy of the particle data object.
    virtual AggModelCache *const Clone(void) const = 0;

    // Returns the aggregation model type ID.  Required for serialization.
    virtual AggModelType ID(void) const = 0;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const = 0;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in     // Input stream.
        ) = 0;
};
}
}

#endif
