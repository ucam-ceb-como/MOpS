/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Mixture class is the base class for gas-phase mixtures used
    by mops.  As mops only solves ideal gas systems, it inherits from
    the Sprog::IdealGas class.
    
    Additionally the Mixture class includes a description of a particle
    population which can be solved with Sweep.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#ifndef MOPS_MIXTURE_H
#define MOPS_MIXTURE_H

#include "mops_params.h"
#include "mops_mechanism.h"
#include "sprog.h"
#include "sweep.h"
#include "swp_gas_profile.h"
#include <iostream>

namespace Mops
{
class Mixture : public Sweep::Cell
{
public:
    // Constructors.
    Mixture(const Sweep::ParticleModel &model); // Default constructor.
    Mixture(                                    // Stream-reading constructor.
        std::istream &in,                       //   - Input stream.
        const Sweep::ParticleModel &model       //   - Mechanism which define the mixture.
        );

    // Destructors.
    ~Mixture(void); // Defaul destructors.

    // Operators.
    Mixture &operator=(const Mixture &rhs);


    // READ/WRITE/COPY.

    // Creates a clone of the mixture object.
    Mixture *const Clone() const;

    //!Returns the description of the gas-phase mixture as a sprog object
    const Sprog::Thermo::IdealGas &GasPhase(void) const;
    //!Returns the description of the gas-phase mixture  as a sprog object
    Sprog::Thermo::IdealGas &GasPhase(void);

protected:
    // As in Sprog, it is meaningless to define a mixture without knowledge
    // of the constituent species, therefore the default constructor is 
    // declared as protected.
    Mixture(void);
};
};

#endif
