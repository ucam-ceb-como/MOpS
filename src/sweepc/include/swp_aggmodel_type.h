/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The AggModelType enumeration gives ID number for all the different
    possible aggregation models implemented in sweep for primary particles.
    Each model defines rules for how a primary particle is stored in memory
    and how they are modified by particle processes.

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

#ifndef SWEEP_AGGMODEL_TYPE_H
#define SWEEP_AGGMODEL_TYPE_H

namespace Sweep
{
    namespace AggModels
    {
        // IMPORTANT:  DO NOT CHANGE THE VALUES HERE.  IT WILL INVALIDATE
        //             PREVIOUSLY GENERATED INPUT FILES.
        enum AggModelType {
            Spherical_ID   = 10000, // Spherical particle model.
            SurfVol_ID     = 10001, // Surface-volume model (Patterson et al., Combust. Flame, 151, 160-172, 2007).
            //PriPartList_ID = 10002,  // Primary-particle list (West et al., Ind. Eng. Chem. Res., 46, 6147-6156, 2007). 
            //PAH_ID         = 10003,  // A particle storing multiple PAHs (Markus et al., Proceedings of the Combustion Institute, 33, 675-683, 2011)
			PAH_KMC_ID         = 10004, // coupled PAH-PP model (ms785) with KMC-ARS model (zz260 and ar447)
			Silica_ID		= 10005,  // Silica particle model (ms785 & ss663, integrated with git repository by wjm34)
            SurfVolHydrogen_ID = 10006, //! (Blanquart & Pitsch., Combust. Flame, 156, 1614-1626, 2009).
			BinTree_ID      = 10007,  // Generalised form of silica and PAH-PP models for multicomponent systems (wjm34)
			SurfVolSilica_ID = 10008, // Surfvol implementation of the silica model
			BinTreeSilica_ID = 10009, // Binary tree implementation of the silica model
			SurfVolCubic_ID = 10010, //! Surface volume approximation for cuboidal crystals
        };
    }
}
#endif
