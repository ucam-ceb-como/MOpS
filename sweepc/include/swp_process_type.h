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
namespace Processes
{
    enum ProcessType {
        // Basic processes.
        Inception_ID=0,       // Inception process.
        Coagulation_ID=1,     // Coagulation process.
        SurfaceReaction_ID=2, // Surface reaction.
        Condensation_ID=3,    // Condensation process.
        ActSiteRxn_ID=4,      // Active-sites surface reaction.

        // ARSSC model processes.
        ARSSC_Inception_ID=1000  // ARS-SC inception process.
    };
};
};

#endif
