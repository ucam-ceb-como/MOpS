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
        out.write((char*)type, sizeof(type));

        // Serialize the mixture object.
        mix.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sprog, MixtureFactory::Write).");
    }
}
