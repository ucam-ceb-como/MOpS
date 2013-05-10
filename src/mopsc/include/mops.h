/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This is the main include file for other projects which use
    mops.  All components of mops are included here.

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

//! This is the mops header file.
#ifndef MOPS_H
#define MOPS_H

#include "mops_reactor_type.h"
#include "mops_reactor.h"
#include "mops_psr.h"
#include "mops_reactor_factory.h"

#include "mops_params.h"
#include "mops_mixture.h"
#include "mops_mechanism.h"

#include "mops_reactor_network.h"
#include "mops_network_simulator.h"

#include "mops_settings_io.h"
#include "mops_timeinterval.h"

#include "mops_solver_factory.h"

#include "mops_flux_postprocessor.h"
#include "mops_gpc_sensitivity.h"

#include "swp_gas_profile.h"
#include "swp_flamesolver.h"

#endif
