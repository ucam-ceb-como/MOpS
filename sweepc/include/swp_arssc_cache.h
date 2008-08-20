/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ARSSC_Cache class defines the additional data which is
    added to a particle to enable the ARS-SC model.

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

#ifndef SWEEP_ARSSC_CACHE_H
#define SWEEP_ARSSC_CACHE_H

#include "swp_params.h"
#include "swp_submodel.h"
#include "swp_submodel_cache.h"
#include "swp_submodel_type.h"
#include "swp_arssc_model.h"
#include <vector>
#include <map>
#include <iostream>

namespace Sweep
{
// Forward declare parent object.
class ParticleCache;

namespace SubModels
{
// Declare base class for storing sub-model data caches.
    class ARSSC_Cache : public SubModelCache
{
public:
    // Constructors.
    ARSSC_Cache(ParticleCache &parent);   // Default constructor.
    ARSSC_Cache(const ARSSC_Cache &copy); // Copy constructor.

    // Destructors.
    virtual ~ARSSC_Cache(void);

    // Assignment operators.
    virtual ARSSC_Cache &operator=(const ARSSC_Cache &rhs);
    virtual ARSSC_Cache &operator=(const ARSSC_Model &rhs);
    virtual ARSSC_Cache &operator=(const SubModelCache &rhs);
    virtual ARSSC_Cache &operator=(const SubModel &rhs);

    // Compound assignment operators.
    virtual ARSSC_Cache &operator+=(const ARSSC_Cache &rhs);
    virtual ARSSC_Cache &operator+=(const ARSSC_Model &rhs);
    virtual ARSSC_Cache &operator+=(const SubModelCache &rhs);
    virtual ARSSC_Cache &operator+=(const SubModel &rhs);


    // Resets the model data to the default state.
    virtual void Clear();


    // PROPERTIES.

    // Returns the property with the given ID.
    virtual real Property(unsigned int id) const;


    // READ/WRITE/COPY.

    // Returns a copy of the model data.
    virtual ARSSC_Cache *const Clone(void) const;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    virtual SubModelType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,     // Input stream.
        ParticleCache &parent // Parent object.
        );

protected:
    // Can't create a ARSSC_Cache object independently of a
    // parent ParticleCache.
    ARSSC_Cache(void);

private:
    // Principal and combined site counts.
    fvector m_sites;
};
};
};

#endif
