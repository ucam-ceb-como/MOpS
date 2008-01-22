/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains enumerations and routines to describe and convert
    unit systems.
*/

#ifndef GPC_UNIT_SYSTEMS_H
#define GPC_UNIT_SYSTEMS_H

namespace Sprog
{
    // Enumeration of different systems of units.
    enum UnitSystem {SI,CGS};

    // Enumeration of different concentration units.
    enum ConcUnits {MolarConcs, MoleFracs, MassFracs};
};

#endif
