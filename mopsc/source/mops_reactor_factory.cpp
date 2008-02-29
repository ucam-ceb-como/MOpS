#include "mops_reactor_factory.h"
#include "mops_reactor.h"
#include "mops_reactor_type.h"
#include "mops_psr.h"

#include <stdexcept>

using namespace Mops;
using namespace std;

// REACTOR CREATION.

Reactor *const ReactorFactory::Create(Mops::Serial_ReactorType type, 
                                      const Mops::Mechanism &mech)
{
    Reactor *r = NULL;

    switch(type) {
        case Serial_Reactor:
        case Serial_Batch:
            r = new Reactor(mech);
            break;
        case Serial_PSR:
            r = new PSR(mech);
            break;
        default:
            throw invalid_argument("Invalid reactor type "
                                   "(Mops, ReactorFactory::Create).");
    }

    return r;
}

Reactor *const ReactorFactory::Read(std::istream &in,
                                    const Mops::Mechanism &mech)
{
    if (in.good()) {
        Reactor *r = NULL;

        // Read the reactor type from the input stream.
        unsigned int type;
        in.read(reinterpret_cast<char*>(&type), sizeof(type));

        // Read a reactor of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((Serial_ReactorType)type) {
            case Serial_Reactor:
            case Serial_Batch:
                r = new Reactor(in, mech);
                break;
            case Serial_PSR:
                r = new PSR(in, mech);
                break;
            default:
                throw runtime_error("Invalid reactor type read from "
                                    "input stream (Mops, ReactorFactory::Read).");
        }

        return r;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Mops, ReactorFactory::Read).");
    }
}


// REACTOR STREAM OUTPUT.

void ReactorFactory::Write(const Mops::Reactor &r, std::ostream &out)
{
    if (out.good()) {
        // Write the Reactor Serial signature type to the stream.
        unsigned int type = (unsigned int)r.SerialType();
        out.write((char*)&type, sizeof(type));

        // Serialize the reactor object.
        r.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Mops, ReactorFactory::Write).");
    }
}
