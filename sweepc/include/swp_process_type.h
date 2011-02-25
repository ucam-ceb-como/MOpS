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
        //! Form a soot particle with one PAH particle
        PAH_Inception_ID=0,         // Inception process for PAH inceptions
        //! Form a particle from two gas-phase molecules
        Dimer_Inception_ID=1,        // Inception process for dimer inceptions
        Death_ID=2,             // Death process.
        Birth_ID=3,             // Birth process.

        //! Form particles at a constant rate, independent of conditions
        Constant_Inception_ID=4,

        // Surface reactions.
        SurfaceReaction_ID=100, // Surface reaction.
        ActSiteRxn_ID=101,      // Active-sites surface reaction.

        // Condensation processes.
        Condensation_ID=200,    // Condensation process.

        // ARSSC model processes. (no longer used)
        ARSSC_Inception_ID=1000,   // ARS-SC inception process.
        ARSSC_Reaction_ID=1001,    // ARS-SC reaction process.
        ARSSC_Condensation_ID=1002, // ARS-SC condensation process.

        // Coagulation processes
        Transition_Coagulation_ID=10000,   // Transition regime coagulation
        Additive_Coagulation_ID=10001,     // Additive coagulation kernel
        //! Free molecular coagulation of secondary particles
        Secondary_FreeCoagulation_ID = 10002,
        //! Free molecular coagulation of secondary particles with particles from the main population
        Secondary_Primary_Coagulation_ID = 10003,

        // Transport processes
        //! Particle diffusion
        Diffusion_ID=100000,
        //! Particle advection
        Advection_ID=100001,

    };
}
}

#endif
