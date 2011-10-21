/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ActiveSitesModel is a base abstract class which defines the interface
    for active site models.  These models are used by the ActiveSitesReaction
    which is a specialization of the SurfaceReaction class.  The concept
    is that particle surfaces have active "sites" with which gas-phase species
    may react.  The ActiveSitesModel uses gas-phase and particle properties
    to calculate the concentration of these active sites.

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

#ifndef SWEEP_ACTIVESITES_MODEL_H
#define SWEEP_ACTIVESITES_MODEL_H

#include "swp_params.h"
#include "swp_ensemble.h"
#include "swp_particle.h"
#include "swp_actsites_type.h"

//Forward declaration
namespace Sprog {
namespace Thermo {
    class IdealGas;
}
}

namespace Sweep
{
namespace ActSites
{
class ActSitesModel
{
public:
    // Calculates the active site density for the given gas-phase
    // and particle ensemble.
    virtual real SiteDensity(
        real t,                     // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Ensemble &particles   // Particle ensemble.
        ) const = 0;

    // Calculates the active site density for the given gas-phase
    // and particle.
    virtual real SiteDensity(
        real t,                     // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Particle &part        // Particle.
        ) const = 0;

    // Returns the active-sites model type.
    virtual ActSitesType ID(void) const = 0;
protected:
    // ActiveSitesModels are singletons.
    ActSitesModel(void); // Default constructor.
    ActSitesModel(const ActSitesModel &copy); // Copy constructor.
    virtual ~ActSitesModel(void); // Destructor.
    ActSitesModel &operator=(const ActSitesModel &rhs); // Assignment.
};
};
};

#endif
