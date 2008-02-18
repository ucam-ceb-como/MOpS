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
    // IMPORTANT:  DO NOT CHANGE THE VALUES HERE.  IT WILL INVALIDATE
    //             PREVIOUSLY GENERATED INPUT FILES.
    enum ModelType {
        BasicModel_ID=-1,  // Refers to the basic ParticleData properties.
        CoagModel_ID=0,    // Coagulation model (default).
        SVModel_ID=1,      // Surface-volume model.
        PriPartModel_ID=2, // Primary particle model.

        // These models have no data hence cannot be generated with the
        // ModelFactory class.
        ActSites_ID=3,     // Active sites model.
        ABFSites_ID=4      // ABF active-sites model.
    };

    typedef std::set<ModelType> ModelTypeSet;
};

#endif
