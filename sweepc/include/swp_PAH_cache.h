/*
  Author(s):      Markus Sander
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:


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

#ifndef SWEEP_PAH_CACHE_H
#define SWEEP_PAH_CACHE_H

#include "swp_params.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_primary.h"
#include <iostream>

namespace Sweep
{
namespace AggModels
{
// Forward declare PAHPrimary class.
class PAHPrimary;

class PAHCache : public AggModelCache
{
public:
    // Constructors.
    //    Note default constructor is protected to prevent creation
    //    without knowledge of parent.
    PAHCache(ParticleCache &parent);    // Initialising constructor.
    PAHCache(const PAHCache &copy); // Copy constructor.
    PAHCache(             // Stream-reading constructor.
        std::istream &in,     //  - Input stream.
        ParticleCache &parent //  - Parent object.
        );

    // Destructors.
    virtual ~PAHCache(void);

    // Assignment operators.
    virtual PAHCache &operator=(const PAHCache &rhs);
    virtual PAHCache &operator=(const PAHPrimary &rhs);
    virtual PAHCache &operator=(const AggModelCache &rhs);
    virtual PAHCache &operator=(const Primary &rhs);

    // Compound assignment operators.
    virtual PAHCache &operator+=(const PAHCache &rhs);
    virtual PAHCache &operator+=(const PAHPrimary &rhs);
    virtual PAHCache &operator+=(const AggModelCache &rhs);
    virtual PAHCache &operator+=(const Primary &rhs);


    // DATA MANAGEMENT.

    // Resets the model data to the default state.
    virtual void Clear();

	real NumPAH();
    // AGGREGATION MODEL PARAMETERS.





    // READ/WRITE/COPY.

    // Returns a correctly typed reference to the model cache object.
//    virtual PAHCache &TypedRef(void);
//    virtual const PAHCache &TypedRef(void) const;

    // Returns a copy of the data.
    virtual PAHCache *const Clone(void) const;

    // Returns the model ID.
    virtual AggModelType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,     // Input stream.
        ParticleCache &parent // Parent.
        );
	  real m_numPAH;
	  real testcachevar;

protected:
  

    // Can't create a PAHCache without knowledge
    // of the parent ParticleCache.
    PAHCache(void);

};
};
};

#endif
