/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    Definition of enumeration of different solver types.
*/

#ifndef MOPS_SOLVER_TYPE_H
#define MOPS_SOLVER_TYPE_H

namespace Mops
{
    enum SolverType {
        GPC,     // Gas-phase chemistry only, default.
        OpSplit, // Use simple operator splitting.
        Strang,  // Strang splitting.
        PredCor, // Split-Predictor---Split-Corrector.
        FlamePP, // Post-process a gas-phase profile (like sweep1).
        MoMIC,   // Method-of-moments for 1D particles.
    };
};

#endif
