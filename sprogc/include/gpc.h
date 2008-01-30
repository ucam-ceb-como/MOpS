#ifndef GPC_H
#define GPC_H

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
