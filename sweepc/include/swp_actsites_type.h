/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ModelType enumeration gives ID number for all the different
    possible particle models implemented in sweep.  This is used
    to identify models in ParticleData objects and to generate models
    using the ModelFactory class.

    Currently there is 1 active sites model implemented in sweepc.
    These are listed below with literature references.

    ACTIVE-SITES MODELS
    -------------------
    1.  ABF soot model (Appel et al., Combust. Flame, 121, 122-136, 2000).

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

#ifndef SWEEP_ACTSITES_TYPE_H
#define SWEEP_ACTSITES_TYPE_H

#include <set>

namespace Sweep
{
namespace ActSites
{
    // IMPORTANT:  DO NOT CHANGE THE VALUES HERE.  IT WILL INVALIDATE
    //             PREVIOUSLY GENERATED INPUT FILES.
    enum ActSitesType {
        // ACTIVE-SITE MODELS.
        // These models have no data hence cannot be generated with the
        // ModelFactory class.
        ActSites_ID = 20000, // Active sites model base.
        ABFSites_ID = 20004  // ABF active-sites model.
    };
    
    typedef std::set<ActSitesType> ActSitesTypeSet;
};
};

#endif
