/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This is the main include file for other projects which use
    mops.  All components of sprog are included here.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
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

#ifndef GPC_H
#define GPC_H

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// Global parameters.
#include "gpc_params.h"
#include "gpc_unit_systems.h"

// Chemical definitions.
#include "gpc_el_comp.h"
#include "gpc_element.h"
#include "gpc_species.h"

// Thermodynamics.
#include "gpc_thermo_params.h"
#include "gpc_thermo.h"
#include "gpc_mixture.h"
#include "gpc_gasphase.h"
#include "gpc_idealgas.h"
#include "gpc_mixture_type.h"
#include "gpc_mixture_factory.h"

// Kinetics.
#include "gpc_stoich.h"
#include "gpc_rate_params.h"
#include "gpc_reaction.h"
#include "gpc_reaction_set.h"

// Mechanisms.
#include "gpc_mech.h"
#include "gpc_mech_io.h"

#endif
