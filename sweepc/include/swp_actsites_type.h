/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ModelType enumeration gives ID number for all the different
    possible particle models implemented in sweep.  This is used
    to identify models in ParticleData objects and to generate models
    using the ModelFactory class.

    Currently there is 1 active sites model implemented in sweepc.
    These are listed below with literature references.

    ACTIVE-SITES MODELS
    -------------------
    1.  ABF soot model (Appel et al., Combust. Flame, 121, 122-136, 2000).
*/

#ifndef SWEEP_ACTSITES_TYPE_H
#define SWEEP_ACTSITES_TYPE_H

namespace Sweep
{
namespace ActSites
{
    // IMPORTANT:  DO NOT CHANGE THE VALUES HERE.  IT WILL INVALIDATE
    //             PREVIOUSLY GENERATED INPUT FILES.
    enum ActSitesType {
        // ACTIVE-SITE MODELS.
        // These models have no data hence cannot be generated with the
        // ModelFactory class.
        ActSites_ID = 20000, // Active sites model base.
        ABFSites_ID = 20004  // ABF active-sites model.
    };
    
    typedef std::set<ActSitesType> ActSitesTypeSet;
};
};

#endif
