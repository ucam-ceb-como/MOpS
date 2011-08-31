/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the MixtureFactory class declared in the
    gpc_mixture_factory.h header file.

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

#include "gpc_params.h"
#include "gpc_mixture_factory.h"
#include "gpc_mixture.h"
#include "gpc_gasphase.h"
#include "gpc_idealgas.h"
#include <vector>
#include <iostream>
#include <stdexcept>

using namespace Sprog;
using namespace Sprog::Thermo;
using namespace std;

// MIXTURE CREATION.

// Creates a mixture of the given type.
Mixture *const MixtureFactory::Create(Serial_MixtureType type,
                                      const Sprog::SpeciesPtrVector &species)
{
    Mixture *mix = NULL;

    switch (type) {
        case Serial_IdealGas:
            mix = new Sprog::Thermo::IdealGas(species);
            break;
        default:
            throw invalid_argument("Invalid mixture type "
                                   "(Sprog, MixtureFactory::Create).");
    }

    return mix;
}

// Reads a mixture from the given input stream.  The first item
// read from the file is the mixture type.
Mixture *const MixtureFactory::Read(std::istream &s,
                                    const Sprog::SpeciesPtrVector &species)
{
    if (s.good()) {
        Mixture *mix = NULL;

        // Read the mixture type from the input stream.
        unsigned int type;
        s.read((char*)&type, sizeof(type));

        // Read a mixture of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((Serial_MixtureType)type) {
            case Serial_IdealGas:
                mix = new IdealGas(s, species);
                break;
            default:
                throw runtime_error("Invalid mixture type read from "
                                    "input stream (Sprog, MixtureFactory::Read).");
        }

        return mix;
    } else {
        throw invalid_argument("Input stream not ready (Sprog, MixtureFactory::Read).");
    }
    return NULL;
}


// GAS-PHASE SPECIFIC CREATION.

// Creates a GasPhase of the given type.
GasPhase *const MixtureFactory::CreateGasPhase(Serial_MixtureType type,
                                               const Sprog::SpeciesPtrVector &species)
{
    GasPhase *gas = NULL;

    switch (type) {
        case Serial_IdealGas:
            gas = new IdealGas(species);
            break;
        default:
            throw invalid_argument("Invalid mixture type "
                                   "(Sprog, MixtureFactory::CreateGasPhase).");
    }

    return gas;
}

// Reads a GasPhase from the given input stream.  The first item
// read from the file is the mixture type.
GasPhase *const MixtureFactory::ReadGasPhase(std::istream &s,
                                             const Sprog::SpeciesPtrVector &species)
{
    if (s.good()) {
        GasPhase *gas = NULL;

        // Read the mixture type from the input stream.
        unsigned int type;
        s.read((char*)&type, sizeof(type));

        // Read a GasPhase of this particular type.
        switch ((Serial_MixtureType)type) {
            case Serial_IdealGas:
                gas = new IdealGas(s, species);
                break;
            default:
                throw runtime_error("Invalid mixture type read from "
                                    "input stream (Sprog, MixtureFactory::ReadGasPhase).");
        }

        return gas;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sprog, MixtureFactory::ReadGasPhase).");
    }
    return NULL;
}


// MIXTURE STREAM OUTPUT.

// Writes a Mixture object to an open output stream.
void MixtureFactory::Write(const Mixture &mix, std::ostream &out)
{
    if (out.good()) {
        // Write the Mixture Serial signature type to the stream.
        unsigned int type = (unsigned int)mix.SerialType();
        out.write((char*)&type, sizeof(type));

        // Serialize the mixture object.
        //mix.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sprog, MixtureFactory::Write).");
    }
}
