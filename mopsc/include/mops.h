/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    This is the main include file for other projects which use
    mops.  All components of mops are included here.
*/

#ifndef MOPS_H
#define MOPS_H

#include "mops_reactor_type.h"
#include "mops_reactor.h"
#include "mops_psr.h"
#include "mops_reactor_factory.h"

#include "mops_params.h"
#include "mops_mixture.h"
#include "mops_mechanism.h"

#include "mops_settings_io.h"
#include "mops_timeinterval.h"

#include "mops_solvertype.h"
#include "mops_solver.h"
#include "mops_strang_solver.h"
#include "mops_predcor_solver.h"

#include "swp_gas_profile.h"
#include "swp_flamesolver.h"
#endif
