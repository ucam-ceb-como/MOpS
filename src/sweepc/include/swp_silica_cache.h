/*
  Author(s):      Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Markus Sander.
  Extended by Shraddha Shekar (ss663) and William Menz (wjm34)

  File purpose:
	The SilicaCache class is a specialisation of the AggModelCache
	class for holding cached data of silica primary particles.

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

#ifndef SWEEP_SILICA_CACHE_H
#define SWEEP_SILICA_CACHE_H

#include "swp_params.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_primary.h"
#include <iostream>

namespace Sweep
{
namespace AggModels
{

// Forward declare the SilicaPrimary class.
class SilicaPrimary;

// Define the silicaCache class. Inherits from AggModelCache class.
class SilicaCache : public AggModelCache {
public:
    // Constructors.
    SilicaCache();    // Initialising constructor.
    SilicaCache(const SilicaCache &copy); // Copy constructor.
    SilicaCache(             // Stream-reading constructor.
        std::istream &in  //  - Input stream.
        );

    // Destructors.
    virtual ~SilicaCache(void);

    // Assignment operators.
    virtual SilicaCache &operator=(const SilicaCache &rhs);
    virtual SilicaCache &operator=(const SilicaPrimary &rhs);
    virtual SilicaCache &operator=(const AggModelCache &rhs);
    virtual SilicaCache &operator=(const Primary &rhs);

    // Compound assignment operators.
    virtual SilicaCache &operator+=(const SilicaCache &rhs);
    virtual SilicaCache &operator+=(const SilicaPrimary &rhs);
    virtual SilicaCache &operator+=(const AggModelCache &rhs);
    virtual SilicaCache &operator+=(const Primary &rhs);


    // DATA MANAGEMENT.

    // Resets the model data to the default state.
    virtual void Clear();

	int NumSi() const;
	int NumO() const;
	int NumOH() const;
    // AGGREGATION MODEL PARAMETERS.


    // Returns a copy of the data.
	virtual SilicaCache *const Clone(void) const;

    // Returns the model ID.
    virtual AggModelType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in     // Input stream.
        );

	  int m_numSi;
	  int m_numO;
	  int m_numOH;
	  int m_numprimary;
	  real m_sqrtLW;
	  real m_LdivW;
	  real m_primarydiam;
	  real m_fdim;
	  real m_Rg;
	  real m_avg_sinter;
	  real m_sintrate;
	  real m_sinttime;
	  real m_createtime;


}; // class SilicaCache

}  // namespace AggModels
}  // namespace Sweep

#endif
