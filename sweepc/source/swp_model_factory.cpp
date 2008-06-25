#include "swp_model_factory.h"
#include "swp_primary.h"
#include "swp_submodel.h"
#include "swp_submodel_type.h"
#include "swp_submodel_cache.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_surfvol_cache.h"
#include "swp_surfvol_primary.h"
#include "swp_pripart_cache.h"
#include "swp_pripart_primary.h"
#include "swp_particle_stats.h"
#include "swp_surfvol_stats.h"
#include "swp_pripart_stats.h"
#include "swp_abf_model.h"
#include <stdexcept>
#include <iostream>

using namespace Sweep;
using namespace std;

// PRIMARY PARTICLE CREATION.

// Creates a new primary particle of the given type.
Primary *const ModelFactory::CreatePrimary(AggModels::AggModelType id, 
                                           real time, const ParticleModel &model)
{
    switch (id) {
        case AggModels::SurfVol_ID:
            return new AggModels::SurfVolPrimary(time, model);
        case AggModels::PriPartList_ID:
            return new AggModels::PriPartPrimary(time, model);
        case AggModels::Spherical_ID:
            // Spherical primary model is default.
        default:
            return new Primary(time, model);
    }
}

// PRIMARY PARTICLE STREAM INPUT.

// Reads a primary particle from a binary stream.  First reads
// the primary type ID, in order to create a primary of the
// correct type.
Primary *const ModelFactory::ReadPrimary(std::istream &in, 
                                         const ParticleModel &model)
{
    if (in.good()) {
        Primary *pri = NULL;

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((AggModels::AggModelType)type) {
            case AggModels::Spherical_ID:
                pri = new Primary(in, model);
                break;
            case AggModels::SurfVol_ID:
                pri = new AggModels::SurfVolPrimary(in, model);
                break;
            case AggModels::PriPartList_ID:
                pri = new AggModels::PriPartPrimary(in, model);
                break;
            default:
                throw invalid_argument("Invalid model ID (Sweep, "
                                       "ModelFactory::ReadPrimary).");
        }

        return pri;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ModelFactory::ReadPrimary).");
    }
}

// PRIMARY PARTICLE STREAM OUTPUT.

// Writes a primary particle, along with its
// ID, to an output stream.
void ModelFactory::WritePrimary(const Primary &pri, std::ostream &out)
{
    if (out.good()) {
        // Write the model Serial signature type to the stream.
        unsigned int type = (unsigned int)pri.AggID();
        out.write((char*)&type, sizeof(type));

        // Serialize the model object.
        pri.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::WritePrimary).");
    }
}


// SUB-MODEL DATA CREATION.

// Creates a new sub-model data object of the given type.
SubModels::SubModel *const ModelFactory::Create(SubModels::SubModelType id, 
                                                Primary &parent)
{
    switch (id) {
        case SubModels::ARSSC_Model_ID:
        case SubModels::CNT_Model_ID:
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::Create).");
    }
}

// Creates a new sub-model cache object of the given type.
SubModels::SubModelCache *const ModelFactory::CreateCache(SubModels::SubModelType id, 
                                                          ParticleCache &parent)
{
    switch (id) {
        case SubModels::ARSSC_Model_ID:
        case SubModels::CNT_Model_ID:
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::CreateCache).");
    }
}


// SUB-MODEL STREAM INPUT.

// Reads a sub-model from a binary stream.  The first item read
// is the model ID which tells the ModelFactory what type
// of model to read.
SubModels::SubModel *const ModelFactory::Read(std::istream &in, 
                                              Primary &parent)
{
    if (in.good()) {
        SubModels::SubModel *model = NULL;

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((SubModels::SubModelType)type) {
            case SubModels::ARSSC_Model_ID:
            case SubModels::CNT_Model_ID:
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

// Reads a sub-model cache from a binary stream.  The first item read
// is the model ID which tells the ModelFactory what type
// of model to read.
SubModels::SubModelCache *const ModelFactory::ReadCache(std::istream &in, 
                                                        ParticleCache &parent)
{
    if (in.good()) {
        SubModels::SubModelCache *cache = NULL;

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((SubModels::SubModelType)type) {
            case SubModels::ARSSC_Model_ID:
            case SubModels::CNT_Model_ID:
            default:
                throw runtime_error("Invalid model type read from "
                                    "input stream (Sweep, ModelFactory::ReadCache).");
        }

        return cache;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ModelFactory::ReadCache).");
    }
}

// Reads model stats from a binary stream.  The first item read
// is the model ID which tells the ModelFactory what type
// of model stats to read.
Stats::IModelStats *const ModelFactory::ReadStats(std::istream &in, const ParticleModel &model)
{
    if (in.good()) {
        Stats::IModelStats *stats = NULL;

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((SubModels::SubModelType)type) {
            case SubModels::BasicModel_ID:
                stats = new Stats::ParticleStats(in, model);
                break;
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
void ModelFactory::Write(const SubModels::SubModel &model, std::ostream &out)
{
    if (out.good()) {
        // Write the model Serial signature type to the stream.
        unsigned int type = (unsigned int)model.ID();
        out.write((char*)&type, sizeof(type));

        // Serialize the model object.
        model.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::Write).");
    }
}

// Writes a sub-model cache, along with its ID to an output stream.
void ModelFactory::WriteCache(const SubModels::SubModelCache &cache,  
                              std::ostream &out)
{
    if (out.good()) {
        // Write the model Serial signature type to the stream.
        unsigned int type = (unsigned int)cache.ID();
        out.write((char*)&type, sizeof(type));

        // Serialize the model object.
        cache.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::Write).");
    }
}

// Writes a model stats object, along with its ID, to an output stream.
void ModelFactory::WriteStats(const Stats::IModelStats &stats, std::ostream &out)
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

// AGGREGATION MODEL CREATION.

// Creates a new aggregation model cache of the given type.
AggModels::AggModelCache *const ModelFactory::CreateAggCache(AggModels::AggModelType id, 
                                                             ParticleCache &parent)
{
    switch (id) {
        case AggModels::Spherical_ID:
            // There is no cache associated with the spherical
            // particle model.
            return NULL;
        case AggModels::SurfVol_ID:
            return new AggModels::SurfVolCache(parent);
        case AggModels::PriPartList_ID:
            return new AggModels::PriPartCache(parent);
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::CreateAggCache).");
    }
}

// Creates a new aggregation model stats object of the given type.
Stats::IModelStats *const ModelFactory::CreateAggStats(AggModels::AggModelType id, 
                                                       const ParticleModel &model)
{
    switch (id) {
        case AggModels::Spherical_ID:
            // There is no cache associated with the spherical
            // particle model.
            return NULL;
        case AggModels::SurfVol_ID:
            return new Stats::SurfVolStats();
        case AggModels::PriPartList_ID:
            return new Stats::PriPartStats();
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::CreateAggStats).");
    }
}


// AGGEGATION MODEL STREAM INPUT.

// Reads an aggregation model cache from a binary stream.  The first item read
// is the model ID which tells the ModelFactory what type
// of model to read.
AggModels::AggModelCache *const ModelFactory::ReadAggCache(std::istream &in,
                                                           ParticleCache &parent)
{
    if (in.good()) {
        AggModels::AggModelCache *model = NULL;

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((AggModels::AggModelType)type) {
            case AggModels::SurfVol_ID:
                model = new AggModels::SurfVolCache(in, parent);
                model->SetParent(parent);
                break;
            case AggModels::PriPartList_ID:
                model = new AggModels::PriPartCache(in, parent);
                model->SetParent(parent);
                break;
            default:
                throw invalid_argument("Invalid model ID (Sweep, "
                                       "ModelFactory::ReadAggCache).");
        }

        return model;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ModelFactory::ReadAggCache).");
    }
}

// Reads aggregation model stats from a binary stream.  The first 
// item read is the model ID which tells the ModelFactory what type
// of aggregation model stats to read.
Stats::IModelStats *const ModelFactory::ReadAggStats(std::istream &in, 
                                                     const ParticleModel &model)
{
    if (in.good()) {
        Stats::IModelStats *stats = NULL;

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((AggModels::AggModelType)type) {
            case AggModels::SurfVol_ID:
                stats = new Stats::SurfVolStats(in, model);
                break;
            case AggModels::PriPartList_ID:
                stats = new Stats::PriPartStats(in, model);
                break;
            default:
                throw invalid_argument("Invalid model ID (Sweep, "
                                       "ModelFactory::ReadAggStats).");
        }
        return stats;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ModelFactory::ReadAggStats).");
    }
}


// AGGREGATION MODEL STREAM OUTPUT.

// Writes an aggregation model cache, along with its ID to an output stream.
void ModelFactory::WriteCache(const AggModels::AggModelCache &cache,
                              std::ostream &out)
{
    if (out.good()) {
        // Write the model Serial signature type to the stream.
        unsigned int type = (unsigned int)cache.ID();
        out.write((char*)&type, sizeof(type));

        // Serialize the model object.
        cache.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::WriteCache[AggModelCache]).");
    }
}

// Writes an aggregation model stats object, along with 
// its ID, to an output stream.
void ModelFactory::WriteAggStats(const Stats::IModelStats &stats, std::ostream &out)
{
    if (out.good()) {
        // Write the model Serial signature type to the stream.
        unsigned int type = (unsigned int)stats.ID();
        out.write((char*)&type, sizeof(type));

        // Serialize the model object.
        stats.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::WriteAggStats).");
    }
}


    // ACTIVE-SITES MODEL INSTANCE AQUISITION.

// Returns the instance of the active-sites model with the given ID.
// Note that active-sites models are singleton classes.
ActSites::ActSitesModel *const ModelFactory::GetActSitesModel(ActSites::ActSitesType id)
{
    switch (id) {
        case ActSites::ABFSites_ID:
            return &ActSites::ABFModel::Instance();
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::GetActSitesModel).");
    }
}
