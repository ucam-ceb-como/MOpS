/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Header file to include in any program that uses sweep.  This file includes
    all header files in the sweep code.

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

#ifndef SWEEP_H
#define SWEEP_H

// Parameters.
#include "swp_params.h"
#include "swp_maths_functional.h"
#include "swp_coords.h"

// Driver.
#include "swp_solver.h"
#include "swp_mech_parser.h"
#include "swp_model_factory.h"

// General particle implementation.
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_particle_model.h"
#include "swp_primary.h"
#include "swp_particle.h"
#include "swp_PAH.h"

// Aggregation models.
#include "swp_aggmodel_type.h"

// KMC PAH growth model.
#include "swp_kmc_jump_process.h"
#include "swp_kmc_mech.h"
#include "swp_kmc_pah_structure.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_reaction.h"
#include "swp_kmc_simulator.h"
//#include "swp_kmc_starting_struct.h"
#include "swp_kmc_structure_comp.h"
#include "swp_kmc_typedef.h"
#include "swp_kmc_gaspoint.h"

// Particle system definition.
#include "swp_cell.h"

// Mechanism and process definitions.
#include "swp_mechanism.h"

// Statistics generation.
#include "swp_particle_stats.h"
#include "swp_ensemble_stats.h"
#include "swp_surfvol_stats.h"

#endif
