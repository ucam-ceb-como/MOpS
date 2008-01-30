/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    Definition of enumeration of different reactor types.  Required
    for reactor serialisation.
*/

#ifndef MOPS_REACTOR_TYPE_H
#define MOPS_REACTOR_TYPE_H

namespace Mops
{
    enum Serial_ReactorType {Serial_Reactor, Serial_Batch, Serial_PSR, Serial_ShockTube};
};

#endif
