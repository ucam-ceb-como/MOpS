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

    Currently there are 3 models implemented in sweepc, 2 particle models
    and 1 active-sites model.  These are listed below with literature references.

    PARTICLE MODELS
    ---------------
    1.  Surface-volume model (Patterson et al., Combust. Flame, 151, 160-172, 2007).
    2.  Primary-particle model (West et al., Ind. Eng. Chem. Res., 46, 6147-6156, 2007).

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

#ifndef SWEEP_SUBMODEL_TYPE_H
#define SWEEP_SUBMODEL_TYPE_H

#include <set>

namespace Sweep
{
    namespace SubModels
    {
        // IMPORTANT:  DO NOT CHANGE THE VALUES HERE.  IT WILL INVALIDATE
        //             PREVIOUSLY GENERATED INPUT FILES.
        enum SubModelType {
            BasicModel_ID = -1, // Refers to the basic ParticleData properties.

            // PARTICLE MODELS.
            //ARSSC_Model_ID = 1, // Aromatic-site site-counting model (Celnik et al., Combust. Flame, 2008, in press (PP51)).
            CNT_Model_ID   = 2, // Simple carbon nanotube model (Celnik et al., Carbon, 46(3), 422-433, 2008).

        };

        typedef std::set<SubModelType> SubModelTypeSet;
    };
};

#endif
