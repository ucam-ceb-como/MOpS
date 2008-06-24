#include "swp_process_factory.h"
#include "swp_surface_reaction.h"
#include "swp_condensation.h"
#include "swp_actsites_reaction.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// PROCESS CREATION.

// Creates a new process data object of the given type.
Process *const ProcessFactory::Create(ProcessType id, const Sweep::Mechanism &mech)
{
    switch (id) {
        case Inception_ID:
            return new Inception(mech);
        case Coagulation_ID:
            return new Coagulation(mech);
        case SurfaceReaction_ID:
            return new SurfaceReaction(mech);
        case Condensation_ID:
            return new Condensation(mech);
        case ActSiteRxn_ID:
            return new ActSiteReaction(mech);
        default:
            throw invalid_argument("Invalid process ID (Sweep, "
                                   "ProcessFactory::Create).");
    }
}


// STREAM INPUT.

// Reads a process from a binary stream.  The first item read
// is the process ID which tells the ModelFactory what type
// of process to read.
Process *const ProcessFactory::Read(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        Process *proc = NULL;

        // Read the process type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a process of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((ProcessType)type) {
            case Inception_ID:
                proc = new Inception(in, mech);
                break;
            case Coagulation_ID:
                proc = new Coagulation(in, mech);
                break;
            case SurfaceReaction_ID:
                proc = new SurfaceReaction(in, mech);
                break;
            case Condensation_ID:
                proc = new Condensation(in, mech);
                break;
            case ActSiteRxn_ID:
                return new ActSiteReaction(in, mech);
            default:
                throw runtime_error("Invalid process type read from "
                                    "input stream (Sweep, ProcessFactory::Read).");
        }

        return proc;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ProcessFactory::Read).");
    }
}

// Reads an inception from a binary stream.  The first item read
// is the inception ID which tells the ModelFactory what type
// of inception to read.
Inception *const ProcessFactory::ReadInception(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        Inception *proc = NULL;

        // Read the process type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read an inception of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((ProcessType)type) {
            case Inception_ID:
                proc = new Inception(in, mech);
                break;
            default:
                throw runtime_error("Invalid inception type read from "
                                    "input stream (Sweep, "
                                    "ProcessFactory::ReadInception).");
        }

        return proc;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ProcessFactory::ReadInception).");
    }
}

// Reads a particle process from a binary stream.  The first item read
// is the process ID which tells the ModelFactory what type
// of process to read.
ParticleProcess *const ProcessFactory::ReadPartProcess(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        ParticleProcess *proc = NULL;

        // Read the process type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a process of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((ProcessType)type) {
            case SurfaceReaction_ID:
                proc = new SurfaceReaction(in, mech);
                break;
            case Condensation_ID:
                proc = new Condensation(in, mech);
                break;
            case ActSiteRxn_ID:
                return new ActSiteReaction(in, mech);
            default:
                throw runtime_error("Invalid inception type read from "
                                    "input stream (Sweep, "
                                    "ProcessFactory::ReadPartProcess).");
        }

        return proc;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ProcessFactory::ReadPartProcess).");
    }
}

// Reads an coagulation from a binary stream.  The first item read
// is the coagulation ID which tells the ModelFactory what type
// of coagulation to read.
Coagulation *const ProcessFactory::ReadCoag(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        Coagulation *proc = NULL;

        // Read the process type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read an inception of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((ProcessType)type) {
            case Coagulation_ID:
                proc = new Coagulation(in, mech);
                break;
            default:
                throw runtime_error("Invalid inception type read from "
                                    "input stream (Sweep, "
                                    "ProcessFactory::ReadCoag).");
        }

        return proc;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ProcessFactory::ReadInception).");
    }
}


// STREAM OUTPUT.

// Writes a process, along with its ID to an output stream.
void ProcessFactory::Write(const Process &proc, std::ostream &out)
{
    if (out.good()) {
        // Write the process Serial signature type to the stream.
        unsigned int type = (unsigned int)proc.ID();
        out.write((char*)&type, sizeof(type));

        // Serialize the process object.
        proc.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ProcessFactory::Write).");
    }
}
