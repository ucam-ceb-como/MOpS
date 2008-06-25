/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ModelFactory is a factory class for sweep particle models.  It
    provides routines for creating, reading and writing particle model
    data objects.
*/

#ifndef SWEEP_MODELFACTORY_H
#define SWEEP_MODELFACTORY_H

#include "swp_params.h"
#include "swp_primary.h"
#include "swp_particle_model.h"
#include "swp_particle_cache.h"
#include "swp_submodel.h"
#include "swp_submodel_type.h"
#include "swp_submodel_cache.h"
#include "swp_model_stats.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_actsites_model.h"
#include "swp_actsites_type.h"
#include <iostream>

namespace Sweep
{
class ModelFactory
{
public:
    // PRIMARY PARTICLE CREATION.

    // Creates a new primary particle of the given type.
    static Primary *const CreatePrimary(
        AggModels::AggModelType id, // Model ID.
        real time,                  // Primary create time.
        const ParticleModel &model  // Defining particle model.
        );

    // PRIMARY PARTICLE STREAM INPUT.

    // Reads a primary particle from a binary stream.  First reads
    // the primary type ID, in order to create a primary of the
    // correct type.
    static Primary *const ReadPrimary(
        std::istream &in,          // Input stream.
        const ParticleModel &model // Defining particle model.
        );

    // PRIMARY PARTICLE STREAM OUTPUT.

    // Writes a primary particle, along with its
    // ID, to an output stream.
    static void WritePrimary(
        const Primary &pri, // Primary to write.
        std::ostream &out   // Output stream.
        );


    // SUB-MODEL CREATION.

    // Creates a new sub-model object of the given type.
    static SubModels::SubModel *const Create(
        SubModels::SubModelType id, // Model ID.
        Primary &parent             // Parent object.
        );

    // Creates a new sub-model cache object of the given type.
    static SubModels::SubModelCache *const CreateCache(
        SubModels::SubModelType id, // Model ID.
        ParticleCache &parent       // Parent object.
        );

    // Creates a new sub-model stats object of the given type.
    static Stats::IModelStats *const CreateStats(
        SubModels::SubModelType id, // Model ID.
        const ParticleModel &model  // Defining particle model.
        );


    // SUB-MODEL STREAM INPUT.

    // Reads a submodel from a binary stream.  The first item read
    // is the model ID which tells the ModelFactory what type
    // of sub-model to read.
    static SubModels::SubModel *const Read(
        std::istream &in, // Input stream.
        Primary &parent   // Parent object.
        );

    // Reads a sub-model cache from a binary stream.  The first item read
    // is the model ID which tells the ModelFactory what type
    // of model to read.
    static SubModels::SubModelCache *const ReadCache(
        std::istream &in,    // Input stream.
        ParticleCache &parent // Parent object.
        );

    // Reads sub-model stats from a binary stream.  The first item read
    // is the model ID which tells the ModelFactory what type
    // of model stats to read.
    static Stats::IModelStats *const ReadStats(
        std::istream &in,          // Input stream.
        const ParticleModel &model // Defining particle model.
        );


    // SUB-MODEL STREAM OUTPUT.

    // Writes a sub-model, along with its ID to an output stream.
    static void Write(
        const SubModels::SubModel &model, // Model to write.
        std::ostream &out                 // Output stream.
        );

    // Writes a sub-model cache, along with its ID to an output stream.
    static void WriteCache(
        const SubModels::SubModelCache &model, // Model to write.
        std::ostream &out                      // Output stream.
        );

    // Writes a sub-model stats object, along with its ID, to an output stream.
    static void WriteStats(
        const Stats::IModelStats &stats, // Stats to write.
        std::ostream &out                // Output stream.
        );


    // AGGREGATION MODEL CREATION.

    // Creates a new aggregation model cache of the given type.
    static AggModels::AggModelCache *const CreateAggCache(
        AggModels::AggModelType id, // Model ID.
        ParticleCache &parent       // Parent object.
        );

    // Creates a new aggregation model stats object of the given type.
    static Stats::IModelStats *const CreateAggStats(
        AggModels::AggModelType id, // Model ID.
        const ParticleModel &model  // Defining particle model.
        );


    // AGGEGATION MODEL STREAM INPUT.

    // Reads an aggregation model cache from a binary stream.  The first 
    // item read is the model ID which tells the ModelFactory what type
    // of model to read.
    static AggModels::AggModelCache *const ReadAggCache(
        std::istream &in,    // Input stream.
        ParticleCache &parent // Parent object.
        );

    // Reads aggregation model stats from a binary stream.  The first 
    // item read is the model ID which tells the ModelFactory what type
    // of aggregation model stats to read.
    static Stats::IModelStats *const ReadAggStats(
        std::istream &in,          // Input stream.
        const ParticleModel &model // Defining particle model.
        );


    // AGGREGATION MODEL STREAM OUTPUT.

    // Writes an aggregation model cache, along with its ID to an output stream.
    static void WriteCache(
        const AggModels::AggModelCache &model, // Model to write.
        std::ostream &out                      // Output stream.
        );

    // Writes an aggregation model stats object, along with 
    // its ID, to an output stream.
    static void WriteAggStats(
        const Stats::IModelStats &stats, // Stats to write.
        std::ostream &out                // Output stream.
        );


    // ACTIVE-SITES MODEL INSTANCE AQUISITION.

    // Returns the instance of the active-sites model with the given ID.
    // Note that active-sites models are singleton classes.
    static ActSites::ActSitesModel *const GetActSitesModel(
        ActSites::ActSitesType id // Model ID
        );
};
};

#endif
