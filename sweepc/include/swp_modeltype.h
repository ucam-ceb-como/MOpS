/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ModelType enumeration gives ID number for all the different
    possible particle models implemented in sweep.  This is used
    to identify models in ParticleData objects and to generate models
    using the ModelFactory class.
*/

#ifndef SWEEP_MODELTYPE_H
#define SWEEP_MODELTYPE_H

#include <set>

namespace Sweep
{
    enum ModelType {
        BasicModel_ID = -1, // Refers to the basic ParticleData properties.
        CoagModel_ID,       // Coagulation model (default).
        SVModel_ID,         // Surface-volume model.
        PriPartModel_ID,    // Primary particle model.
    };

    typedef std::set<ModelType> ModelTypeSet;
};

#endif
