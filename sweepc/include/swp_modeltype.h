/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ModelType enumeration gives ID number for all the different
    possible particle models implemented in sweep.  This is used
    to identify models in ParticleData objects and to generate models
    using the ModelFactory class.

    Currently there are 3 models implemented in sweepc, 2 particle models
    and 1 active-sites model.  These are listed below with literature references.

    PARTICLE MODELS
    ---------------
    1.  Surface-volume model (Patterson et al., Combust. Flame, 151, 160-172, 2007).
    2.  Primary-particle model (West et al., Ind. Eng. Chem. Res., 46, 6147-6156, 2007).

    ACTIVE-SITES MODELS
    -------------------
    1.  ABF soot model (Appel et al., Combust. Flame, 121, 122-136, 2000).
*/

#ifndef SWEEP_MODELTYPE_H
#define SWEEP_MODELTYPE_H

#include <set>

namespace Sweep
{
    // IMPORTANT:  DO NOT CHANGE THE VALUES HERE.  IT WILL INVALIDATE
    //             PREVIOUSLY GENERATED INPUT FILES.
    enum ModelType {
        BasicModel_ID = -1, // Refers to the basic ParticleData properties.

        // PARTICLE MODELS.
        CoagModel_ID    = 10000, // Coagulation model (Default spherical particle model).
        SVModel_ID      = 10001, // Surface-volume model (Patterson et al., Combust. Flame, 151, 160-172, 2007).
        PriPartModel_ID = 10002, // Primary particle model (West et al., Ind. Eng. Chem. Res., 46, 6147-6156, 2007).

        // ACTIVE-SITE MODELS.
        // These models have no data hence cannot be generated with the
        // ModelFactory class.
        ActSites_ID = 20003, // Active sites model base.
        ABFSites_ID = 20004  // ABF active-sites model.
    };

    typedef std::set<ModelType> ModelTypeSet;
};

#endif
