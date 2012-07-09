/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The MixtureFactory class is a factory class for sprog mixture
    objects.  It contains routines for creating specific mixture types,
    in particular GasPhase objects.  As there are only a limited number
    of possible Mixture-derived classes, there is no requirement for an
    extensible factory class.  Instead the type IDs are maintained using
    an enumeration Serial_MixtureType which is defined in its own header.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef GPC_MIXTURE_FACTORY_H
#define GPC_MIXTURE_FACTORY_H

#include "gpc_params.h"
#include "gpc_mixture.h"
#include "gpc_gasphase.h"
#include "gpc_species.h"
#include <iostream>

namespace Sprog
{
namespace Thermo
{
class MixtureFactory
{
public:
    // MIXTURE CREATION.
    // Use these routines if a generic Mixture pointer is required.

    // Creates a new mixture object of the given type.
    static Mixture *const Create(
        Serial_MixtureType type,        // Type of mixture to create.
        const SpeciesPtrVector &species // List of species used to define the mixture.
        );

    // Reads a Mixture object from a binary stream.  The first thing
    // read from the stream is a MixtureType to define what type
    // of mixture to read from the stream.
    static Mixture *const Read(
        std::istream &in,               // Input stream from which to read mixture.
        const SpeciesPtrVector &species // List of species used to define the mixture.
        );


    // GAS-PHASE SPECIFIC CREATION.
    // Use these routines if the object is derived from GasPhase.

    // Creates a GasPhase object of the given type.  Throws an exception
    // if the specified type is not containing a gas-phase
    static GasPhase *const CreateGasPhase(
        Serial_MixtureType type,        // Type of mixture to create.
        const SpeciesPtrVector &species // List of species used to define the mixture.
        );

    // Reads a GasPhase object from a binary stream.  The first thing
    // read from the stream is a MixtureType to define what type
    // of gas phase to read from the stream.  Throws an exception if
    // the mixture type read is not a gas-phase
    static GasPhase *const ReadGasPhase(
        std::istream &in,               // Input stream from which to read mixture.
        const SpeciesPtrVector &species // List of species used to define the mixture.
        );


    // MIXTURE STREAM OUTPUT.

    static void Write(const Mixture &mix, std::ostream &out);
};
};
};

#endif
