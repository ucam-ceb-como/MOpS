/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The PriPartCache class holds cached data for a pri-part primary particle.
    That is a primary which is described by the simple primary-particle
    aggregation model (West et al.).

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

#ifndef SWEEP_PRIPART_CACHE_H
#define SWEEP_PRIPART_CACHE_H

#include "swp_primary.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_surfvol_primary.h"
#include "swp_surfvol_cache.h"
#include <vector>

namespace Sweep
{
namespace AggModels
{
class PriPartPrimary;

class PriPartCache : public SurfVolCache
{
public:
    // Constructors.
    //    Note default constructor is protected to prevent creation
    //    without knowledge of parent.
    PriPartCache(ParticleCache &parent);    // Initialising constructor.
    PriPartCache(const PriPartCache &copy); // Copy constructor.
    PriPartCache(              // Stream-reading constructor.
        std::istream &in,      //  - Input stream.
        ParticleCache &parent  //  - Parent object.
        );

    // Destructors.
    virtual ~PriPartCache(void);

    // Assignment operators.
    virtual PriPartCache &operator=(const PriPartCache &rhs);
    virtual PriPartCache &operator=(const PriPartPrimary &rhs);
    virtual PriPartCache &operator=(const SurfVolCache &rhs);
    virtual PriPartCache &operator=(const SurfVolPrimary &rhs);
    virtual PriPartCache &operator=(const AggModelCache &rhs);
    virtual PriPartCache &operator=(const Primary &rhs);

    // Compound assignment operators.
    virtual PriPartCache &operator+=(const PriPartCache &rhs);
    virtual PriPartCache &operator+=(const PriPartPrimary &rhs);
    virtual PriPartCache &operator+=(const SurfVolCache &rhs);
    virtual PriPartCache &operator+=(const SurfVolPrimary &rhs);
    virtual PriPartCache &operator+=(const AggModelCache &rhs);
    virtual PriPartCache &operator+=(const Primary &rhs);


    // DATA MANAGEMENT.

    // Resets the model data to the default state.
    void Clear();


    // PROPERTIES.

    // Returns the number of primary particles.
    unsigned int Count(void) const;

    // Sets the number of primary particles.
    void SetCount(unsigned int n);

    // Returns the primary particle surface area.
    real PriSurfaceArea(void) const;

    // Sets the primary particle surface area.
    void SetPriSurfaceArea(real s);

    // Returns the average primary-particle diameter.
    real AvgPriDiameter(void) const;

    // Sets the average primary-particle diameter.
    void SetAvgPriDiameter(real d);


    // READ/WRITE/COPY.

    // Returns a copy of the model data.
    PriPartCache *const Clone(void) const;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    AggModelType ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,    // Input stream.  
        ParticleCache &cache // Parent.
        );

protected:
    // Can't create a PriPartCache object independently of a
    // parent.
    PriPartCache(void);

private:
    // Cached properties.
    unsigned int m_npri; // Number of primary particles.
    real m_prisurf;      // Primary particle surface area.
    real m_avgdiam;      // Avg. primary-particle diameter.
};
};
};

#endif
