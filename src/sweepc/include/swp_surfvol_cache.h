/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The SurfVolCache class is a specialisation of the AggModelCache
    class for holding cached data of surface-volume primary particles.

    In the surface-volume model primaries are described by the two coordinates
    of mass and surface area.  It is also necessary to know the equivalent
    spherical particle surface area (same mass), so that is stored in this
    class.

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

#ifndef SWEEP_SURFVOL_CACHE_H
#define SWEEP_SURFVOL_CACHE_H

#include "swp_params.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_primary.h"
#include <iostream>

namespace Sweep
{
namespace AggModels
{
// Forward declare SurfVolPrimary class.
class SurfVolPrimary;

/*!
 * \brief Implementation of the AggModelCache interface
 *
 * This class was probably intended for the old primary particle
 * model that has now been removed.  (riap2 25 Feb 2011)
 *
 * \deprecated
 */
class SurfVolCache : public AggModelCache
{
public:
    // Constructors.
    SurfVolCache();
    SurfVolCache(const SurfVolCache &copy); // Copy constructor.
    SurfVolCache(             // Stream-reading constructor.
        std::istream &in      //  - Input stream.
        );

    // Destructors.
    virtual ~SurfVolCache(void);

    // Assignment operators.
    virtual SurfVolCache &operator=(const SurfVolCache &rhs);
    virtual SurfVolCache &operator=(const SurfVolPrimary &rhs);
    virtual SurfVolCache &operator=(const AggModelCache &rhs);
    virtual SurfVolCache &operator=(const Primary &rhs);

    // Compound assignment operators.
    virtual SurfVolCache &operator+=(const SurfVolCache &rhs);
    virtual SurfVolCache &operator+=(const SurfVolPrimary &rhs);
    virtual SurfVolCache &operator+=(const AggModelCache &rhs);
    virtual SurfVolCache &operator+=(const Primary &rhs);


    // DATA MANAGEMENT.

    // Resets the model data to the default state.
    virtual void Clear();


    // AGGREGATION MODEL PARAMETERS.

    // Returns the equivalent spherical surface area.
    real SphSurfaceArea(void) const;

    // Sets the equivalent spherical surface area.
    void SetSphSurfaceArea(real surface);

    // Returns the actual surface area.
    real SurfaceArea(void) const;

    // Sets the actual surface area.
    void SetSurfaceArea(real surface);

    // Returns the primary particle count.
    unsigned int PP_Count(void) const;

    // Sets the primary particle count.
    void SetPP_Count(unsigned int n);

    // Returns the average primary particle diameter.
    real PP_Diameter(void) const;
    
    // Sets the average primary particle diameter.
    void SetPP_Diameter(real d);


    // READ/WRITE/COPY.

    // Returns a copy of the data.
    virtual SurfVolCache *const Clone(void) const;

    // Returns the model ID.
    virtual AggModelType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in     // Input stream.
        );

protected:
    // Surface-volume model particle properties.
    real m_sphsurf;     // Equivalent sphere surface area.
    real m_surf;        // Actual surface area.
    unsigned int m_ppn; // Primary-particle count.
    real m_ppd;         // Average primary-particle diameter.

};
}
}

#endif
