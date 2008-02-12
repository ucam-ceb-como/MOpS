/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ModelFactory is a factory class for sweep particle models.  It
    provides routines for creating, reading and writing particle model
    data objects.
*/

#ifndef SWEEP_MODELFACTORY_H
#define SWEEP_MODELFACTORY_H

#include "swp_params.h"
#include "swp_modeltype.h"
#include "swp_modeldata.h"
#include "swp_particledata.h"

namespace Sweep
{
class ModelFactory
{
public:
    // Creates a new model data object of the given type.
    static IModelData *const CreateData(ModelType id, ParticleData &parent);
};
};

#endif
