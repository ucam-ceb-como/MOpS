/*
  Author(s):      Markus Sander (ms785)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2009 Markus Sander.

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

#ifndef SWEEP_PAH_PRIMARY_H
#define SWEEP_PAH_PRIMARY_H

#include "swp_params.h"
#include "swp_primary.h"
#include "swp_particle_model.h"
#include "swp_sintering_model.h"
#include "swp_aggmodel_type.h"
#include "swp_surfvol_cache.h"
#include "swp_PAH_cache.h"
#include <iostream>

namespace Sweep
{
namespace AggModels 
{
class PAHPrimary : public Primary
{
public:
	    // Constructors.
    //   Note:  Default constructor is protected to prevent a
    //          SurfVolPrimary being created without knowledge of the
    //          defining particle model.
    PAHPrimary(                       // Initialising constructor.
        real time,                        //  - Create time.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );
    PAHPrimary(const PAHPrimary &copy); // Copy constructor.
    PAHPrimary(                       // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );
	  
    // Destructors.
    virtual ~PAHPrimary(void);
      // READ/WRITE/COPY.

	PAHPrimary &operator=(const Primary &rhs);

    // Returns a copy of the primary.
    virtual PAHPrimary *const Clone(void) const;

    // Returns this object's instance.  This may seem rather circular, but
    // it has an important purpose for getting the correct object type reference
    // from a base class reference/pointer.
    virtual PAHPrimary &Instance();
    virtual const PAHPrimary &Instance() const;
    // Returns the collision diameter.
    real PAHCollDiameter(void) const;
	PAHPrimary &Coagulate(const Primary &rhs);
	void UpdateCache(void);
	void UpdateTime(double t);

	struct PAH {
        unsigned int m_numcarbon;
		unsigned int ID;
		double time_created;
		std::vector<double> time;
		std::vector<int> n_carbon_t;
    };
	int m_numcarbon;
	double m_PAHmass;
	double m_PAHCollDiameter;
	int m_numPAH;
	// Vector of PAHs.
    std::vector<PAH> m_PAH;

	AggModels::PAHCache *const CreateAggCache(ParticleCache &pcache) const;

	//serialize
	void Deserialize(std::istream &in, const Sweep::ParticleModel &model);
	void Serialize(std::ostream &out) const;

	AggModels::AggModelType AggID(void) const;

protected:

};
};
};
#endif
