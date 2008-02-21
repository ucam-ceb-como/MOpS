#include "swp_modelfactory.h"
#include "swp_coagmodeldata.h"
#include "swp_pointcontactdata.h"
#include "swp_pripartdata.h"
#include "swp_particlestats.h"
#include "swp_coagmodel.h"
#include "swp_pointcontactmodel.h"
#include "swp_pripartmodel.h"
#include "swp_abfmodel.h"
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
        case ABFSites_ID:
            return NULL;
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

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
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

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
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

// Reads model stats from a binary stream.  The first item read
// is the model ID which tells the ModelFactory what type
// of model stats to read.
IModelStats *const ModelFactory::ReadStats(std::istream &in)
{
    if (in.good()) {
        IModelStats *stats = NULL;

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((ModelType)type) {
            case BasicModel_ID:
                stats = new ParticleStats(in);
                break;
            case CoagModel_ID:
            case SVModel_ID:
            case PriPartModel_ID:
            default:
                throw runtime_error("Invalid model type read from "
                                    "input stream (Sweep, ModelFactory::ReadStats).");
        }

        return stats;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ModelFactory::ReadStats).");
    }
}


// STREAM OUTPUT.

// Writes a model, along with its ID to an output stream.
void ModelFactory::Write(const IModelData &model,  std::ostream &out)
{
    if (out.good()) {
        // Write the model Serial signature type to the stream.
        unsigned int type = (unsigned int)model.ID();
        out.write((char*)type, sizeof(type));

        // Serialize the model object.
        model.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::Write).");
    }
}

// Writes a model stats object, along with its ID, to an output stream.
void ModelFactory::WriteStats(const IModelStats &stats, std::ostream &out)
{
    if (out.good()) {
        // Write the model Serial signature type to the stream.
        unsigned int type = (unsigned int)stats.ID();
        out.write((char*)type, sizeof(type));

        // Serialize the model object.
        stats.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::WriteStats).");
    }
}


// MODEL INSTANCE AQUISITION.

// Returns the instance of the model with the given ID.
IModel *const ModelFactory::GetModel(ModelType id)
{
    switch (id) {
        case CoagModel_ID:
            return &CoagModel::Instance();
        case SVModel_ID:
            return &PointContactModel::Instance();
        case PriPartModel_ID:
            return &PriPartModel::Instance();
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::GetModel).");
    }
}

// Returns the instance of the active-sites model with the given ID.
ActiveSitesModel *const ModelFactory::GetActSitesModel(ModelType id)
{
    switch (id) {
        case ABFSites_ID:
            return &ABFModel::Instance();
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::GetActSitesModel).");
    }
}
