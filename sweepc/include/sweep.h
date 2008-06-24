/*
  Author(s):      Matthew Celnik (msc37)
  Project:        none

  File purpose:
    Header file to include in any program that uses sweep.  This file includes
    all header files in the sweep code.
*/

#ifndef SWEEP_H
#define SWEEP_H

// Parameters.
#include "swp_params.h"

// Driver.
#include "swp_solver.h"
#include "swp_mech_parser.h"
#include "swp_model_factory.h"

// General particle implementation.
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_particle_model.h"
#include "swp_primary.h"
#include "swp_particle_cache.h"
#include "swp_subparticle.h"
#include "swp_particle.h"

// Sub-models.
#include "swp_submodel_type.h"
#include "swp_submodel_cache.h"
#include "swp_submodel.h"

// Aggregation models.
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_surfvol_cache.h"
#include "swp_surfvol_primary.h"
#include "swp_pripart_cache.h"
#include "swp_pripart_primary.h"

// Particle system definition.
#include "swp_treenode.h"
#include "swp_ensemble.h"
#include "swp_cell.h"

// Mechanism and process definitions.
#include "swp_mechanism.h"
#include "swp_process_type.h"
#include "swp_process.h"
#include "swp_inception.h"
#include "swp_particle_process.h"
#include "swp_surface_reaction.h"
#include "swp_actsites_reaction.h"
#include "swp_condensation.h"
#include "swp_coagulation.h"
#include "swp_process_factory.h"

// Chemical models.
#include "swp_actsites_type.h"
#include "swp_actsites_model.h"
#include "swp_abf_model.h"

// Statistics generation.
#include "swp_particle_stats.h"
#include "swp_ensemble_stats.h"
#include "swp_surfvol_stats.h"
#include "swp_pripart_stats.h"
#include "swp_imager.h"
#include "swp_imgnode.h"

#endif
