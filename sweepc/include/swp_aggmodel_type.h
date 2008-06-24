/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The AggModelType enumeration gives ID number for all the different
    possible aggregation models implemented in sweep for primary particles.
    Each model defines rules for how a primary particle is stored in memory
    and how they are modified by particle processes.
*/

#ifndef SWEEP_AGGMODEL_TYPE_H
#define SWEEP_AGGMODEL_TYPE_H

#include <set>

namespace Sweep
{
    namespace AggModels
    {
        // IMPORTANT:  DO NOT CHANGE THE VALUES HERE.  IT WILL INVALIDATE
        //             PREVIOUSLY GENERATED INPUT FILES.
        enum AggModelType {
            Spherical_ID   = 10000, // Spherical particle model.
            SurfVol_ID     = 10001, // Surface-volume model (Patterson et al., Combust. Flame, 151, 160-172, 2007).
            PriPartList_ID = 10002  // Primary-particle list (West et al., Ind. Eng. Chem. Res., 46, 6147-6156, 2007). 
        };

        typedef std::set<AggModelType> AggModelTypeSet;
    };
};
#endif
