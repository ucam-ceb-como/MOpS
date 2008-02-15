/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ProcessType enumeration gives ID number for all the different
    possible processes implemented in sweep.  This is used
    to identify processes and to generate processes using the 
    ProcessFactory class.
*/

#ifndef SWEEP_PROCESSTYPE_H
#define SWEEP_PROCESSTYPE_H

namespace Sweep
{
    enum ProcessType {
        Inception_ID,       // Inception process.
        Coagulation_ID,     // Coagulation process.
        SurfaceReaction_ID, // Surface reaction.
        Condensation_ID     // Condensation process.
    };
};

#endif
