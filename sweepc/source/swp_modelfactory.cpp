#include "swp_modelfactory.h"
#include "swp_coagmodeldata.h"
#include "swp_pointcontactdata.h"
#include "swp_pripartdata.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// MODEL DATA CREATION.

// Creates a new model data object of the given type.
IModelData *const ModelFactory::CreateData(ModelType id, ParticleData &parent)
{
    switch (id) {
        case CoagModel_ID:
            return new CoagModelData(parent);
        case SVModel_ID:
            return new PointContactData(parent);
        case PriPartModel_ID:
            return new PriPartModelData(parent);
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::CreateData).");
    }
}

// Reads a model from a binary stream.  The first item read
// is the model ID which tells the ModelFactory what type
// of model to read.
IModelData *const ModelFactory::Read(std::istream &in, 
                                     ParticleData &parent)
{
    if (in.good()) {
        IModelData *model = NULL;

        // Read the mixture type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a mixture of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((ModelType)type) {
            case CoagModel_ID:
                model = new CoagModelData(in, parent);
                break;
            case SVModel_ID:
                model = new PointContactData(in, parent);
                break;
            case PriPartModel_ID:
                model = new PriPartModelData(in, parent);
                break;
            default:
                throw runtime_error("Invalid model type read from "
                                    "input stream (Sweep, ModelFactory::Read).");
        }

        return model;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ModelFactory::Read).");
    }
}

// Reads a coagulation model from a binary stream.  The first item read
// is the model ID which tells the ModelFactory what type
// of coagulation model to read.
CoagModelData *const ModelFactory::ReadCoag(std::istream &in, 
                                            ParticleData &parent)
{
    if (in.good()) {
        CoagModelData *model = NULL;

        // Read the mixture type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a mixture of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((ModelType)type) {
            case CoagModel_ID:
                model = new CoagModelData(in, parent);
                break;
            case SVModel_ID:
                model = new PointContactData(in, parent);
                break;
            default:
                throw runtime_error("Invalid model type read from "
                                    "input stream (Sweep, ModelFactory::Read).");
        }

        return model;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ModelFactory::Read).");
    }
}


// STREAM OUTPUT.

// Writes a model, along with its ID to an output stream.
void ModelFactory::Write(const IModelData &model,  std::ostream &out)
{
    if (out.good()) {
        // Write the Mixture Serial signature type to the stream.
        unsigned int type = (unsigned int)model.ID();
        out.write((char*)type, sizeof(type));

        // Serialize the mixture object.
        model.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::Write).");
    }
}
