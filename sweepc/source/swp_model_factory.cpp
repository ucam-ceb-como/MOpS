/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ModelFactory class declared in the
    swp_model_factory.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#include "swp_model_factory.h"
#include "swp_primary.h"
#include "swp_submodel.h"
#include "swp_submodel_type.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_surfvol_cache.h"
#include "swp_surfvol_primary.h"
#include "swp_PAH_primary.h"
#include "swp_PAH_cache.h"
#include "swp_particle_stats.h"
#include "swp_surfvol_stats.h"
#include "swp_PAH_stats.h"
#include "swp_abf_model.h"
#include <stdexcept>
#include <iostream>

using namespace Sweep;
using namespace std;

// PRIMARY PARTICLE CREATION.

/*!
 * @param[in]   id          Model for the aggregate structure of the particle
 * @param[in]   time        Time at which particle is being created
 * @param[in]   position    Position at which particle is being created
 * @param[in]   model       Model which defines the meaning of the particles
 * @param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 *
 * @return      Pointer to dynamically allocated primary (caller must delete)
 */
Primary *const ModelFactory::CreatePrimary(const AggModels::AggModelType id,
                                           const real time, const real position,
                                           const ParticleModel &model,
                                           int (*rand_int)(int, int))
{
    switch (id) {
        case AggModels::SurfVol_ID:
            return new AggModels::SurfVolPrimary(time, model);
        case AggModels::PAH_ID:
            return new AggModels::PAHPrimary(time, position, model, rand_int);
        case AggModels::Spherical_ID:
            // Spherical primary model is default.
        default:
            return new Primary(time, model);
    }
}

/*!
 * @param[in]   id          Model for the aggregate structure of the particle
 * @param[in]   time        Time at which particle is being created
 * @param[in]   model       Model which defines the meaning of the particles
 * @param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 *
 * @return      Pointer to dynamically allocated primary (caller must delete)
 */
Primary *const ModelFactory::CreatePrimary(const AggModels::AggModelType id,
                                           const real time, const ParticleModel &model,
                                           int (*rand_int)(int, int))
{
    switch (id) {
        case AggModels::SurfVol_ID:
            return new AggModels::SurfVolPrimary(time, model);
        case AggModels::PAH_ID:
            return new AggModels::PAHPrimary(time, model, rand_int);
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
            case AggModels::PAH_ID:
                pri = new AggModels::PAHPrimary(in, model);
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
        case SubModels::CNT_Model_ID:
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::Create).");
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
        case AggModels::PAH_ID:
            return new AggModels::PAHCache(parent);
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
        case AggModels::PAH_ID:
            return new Stats::PAHStats();	 	// ms785: postprocessing not yet implemented
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
            case AggModels::PAH_ID:
                model = new AggModels::PAHCache(in, parent);
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
            case AggModels::PAH_ID:
                stats = new Stats::PAHStats(in, model);
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
