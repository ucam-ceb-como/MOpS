#include "swp_modelfactory.h"
#include "swp_coagmodeldata.h"
#include "swp_pointcontactdata.h"
#include "swp_pripartdata.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

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
