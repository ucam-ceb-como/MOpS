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

// General particle implementation.
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_particledata.h"
#include "swp_particle.h"
#include "swp_modeldata.h"
#include "swp_model.h"

// Specific particle models.
#include "swp_modeltype.h"
#include "swp_coagmodeldata.h"
#include "swp_coagmodel.h"
#include "swp_surfvoldata.h"
#include "swp_surfvolmodel.h"
#include "swp_primary.h"
#include "swp_pripartdata.h"
#include "swp_pripartmodel.h"
#include "swp_modelfactory.h"

// Particle system definition.
#include "swp_treenode.h"
#include "swp_ensemble.h"
#include "swp_cell.h"

// Mechanism and process definitions.
#include "swp_mechanism.h"
#include "swp_process.h"
#include "swp_inception.h"
#include "swp_coagulation.h"
#include "swp_particleprocess.h"
#include "swp_surfacereaction.h"
#include "swp_condensation.h"
#include "swp_activesites_reaction.h"

// Chemical models.
#include "swp_activesites_model.h"
#include "swp_abfmodel.h"

#endif
