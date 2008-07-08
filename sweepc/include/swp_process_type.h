/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ProcessType enumeration gives ID number for all the different
    possible processes implemented in sweep.  This is used
    to identify processes and to generate processes using the 
    ProcessFactory class.

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

#ifndef SWEEP_PROCESSTYPE_H
#define SWEEP_PROCESSTYPE_H

namespace Sweep
{
namespace Processes
{
    enum ProcessType {
        // Basic processes.
        Inception_ID=0,       // Inception process.
        Coagulation_ID=1,     // Coagulation process.
        SurfaceReaction_ID=2, // Surface reaction.
        Condensation_ID=3,    // Condensation process.
        ActSiteRxn_ID=4,      // Active-sites surface reaction.

        // ARSSC model processes.
        ARSSC_Inception_ID=1000  // ARS-SC inception process.
    };
};
};

#endif
