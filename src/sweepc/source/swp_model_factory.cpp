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
#include "swp_aggmodel_type.h"
#include "swp_surfvol_primary.h"
#include "swp_surfvolhydrogen_primary.h"
#include "swp_surfvol_silica_primary.h"
#include "swp_surfvol_cubic_primary.h"
#include "swp_PAH_primary.h"
#include "swp_particle_stats.h"
#include "swp_surfvol_stats.h"
#include "swp_surfvolhydrogen_stats.h"
#include "swp_PAH_stats.h"
#include "swp_bintree_primary.h"
#include "swp_bintree_silica_primary.h"
#include "swp_bintree_stats.h"
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
 *
 * @return      Pointer to dynamically allocated primary (caller must delete)
 */
AggModels::Primary *const ModelFactory::CreatePrimary(const AggModels::AggModelType id,
                                                      const double time, const double position,
                                                      const ParticleModel &model)
{
    switch (id) {
        case AggModels::SurfVol_ID:
            return new AggModels::SurfVolPrimary(time, model);
        case AggModels::SurfVolHydrogen_ID:
            return new AggModels::SurfVolHydrogenPrimary(time, model);
        case AggModels::SurfVolSilica_ID:
            return new AggModels::SurfVolSilicaPrimary(time, model);
		case AggModels::PAH_KMC_ID:
            return new AggModels::PAHPrimary(time, position, model);
	    case AggModels::BinTree_ID:
            return new AggModels::BinTreePrimary(time, model);
        case AggModels::BinTreeSilica_ID:
            return new AggModels::BinTreeSilicaPrimary(time, model);
        case AggModels::SurfVolCubic_ID:
            return new AggModels::SurfVolCubicPrimary(time, model);
        case AggModels::Spherical_ID:
            // Spherical primary model is default.
        default:
            return new AggModels::Primary(time, model);
    }
}

/*!
 * @param[in]   id          Model for the aggregate structure of the particle
 * @param[in]   time        Time at which particle is being created
 * @param[in]   model       Model which defines the meaning of the particles
 *
 * @return      Pointer to dynamically allocated primary (caller must delete)
 */
AggModels::Primary *const ModelFactory::CreatePrimary(const AggModels::AggModelType id,
                                                      const double time, const ParticleModel &model)
{
    switch (id) {
        case AggModels::SurfVol_ID:
            return new AggModels::SurfVolPrimary(time, model);
        case AggModels::SurfVolHydrogen_ID:
            return new AggModels::SurfVolHydrogenPrimary(time, model);
        case AggModels::SurfVolSilica_ID:
            return new AggModels::SurfVolSilicaPrimary(time, model);
		case AggModels::PAH_KMC_ID:
            return new AggModels::PAHPrimary(time, model);
        case AggModels::BinTree_ID:
            return new AggModels::BinTreePrimary(time, model);
        case AggModels::BinTreeSilica_ID:
            return new AggModels::BinTreeSilicaPrimary(time, model);
        case AggModels::SurfVolCubic_ID:
            return new AggModels::SurfVolCubicPrimary(time, model);
        case AggModels::Spherical_ID:
            // Spherical primary model is default.
        default:
            return new AggModels::Primary(time, model);
    }
}
/*
 * @brief PRIMARY PARTICLE STREAM INPUT.
 *
 * Reads the primary particle from the binary stream. First reads the primary type ID, in order to create a primary of the correct type.
 *
 * @param[in,out]	 in		             Input binary stream
 * @param[in]        model	             Particle model defining interpretation of particle data
 * @param[in,out]    duplicates          Addresses of PAHs for use when reading primary particles
 *
 * @exception		 invalid_argument    Either invalid model type or stream is not ready
 */
AggModels::Primary *const ModelFactory::ReadPrimary(std::istream &in,
                                                    const ParticleModel &model, void *duplicates)
{
    if (in.good()) {
        AggModels::Primary *pri = NULL;

        // Read the model type from the input stream.
        unsigned int type;
        in.read((char*)&type, sizeof(type));

        // Read a model of this particular type.  This will throw
        // an exception if the type is invalid.
        switch ((AggModels::AggModelType)type) {
            case AggModels::Spherical_ID:
                pri = new AggModels::Primary(in, model);
                break;
            case AggModels::SurfVol_ID:
                pri = new AggModels::SurfVolPrimary(in, model);
                break;
            case AggModels::SurfVolHydrogen_ID:
                pri = new AggModels::SurfVolHydrogenPrimary(in, model);
                break;
            case AggModels::SurfVolSilica_ID:
                pri = new AggModels::SurfVolSilicaPrimary(in, model);
                break;
			case AggModels::PAH_KMC_ID:
            {
                // This is a dangerous conversion, but it avoids having to expose the type
                // of the duplicate information too widely
                AggModels::PAHPrimary::PahDeserialisationMap* pah_duplicates = reinterpret_cast<AggModels::PAHPrimary::PahDeserialisationMap*>(duplicates);
                pri = new AggModels::PAHPrimary(in, model, *pah_duplicates);
                break;
            }
            case AggModels::BinTree_ID:
                pri = new AggModels::BinTreePrimary(in, model);
                break;
            case AggModels::BinTreeSilica_ID:
                pri = new AggModels::BinTreeSilicaPrimary(in, model);
                break;
            case AggModels::SurfVolCubic_ID:
                pri = new AggModels::SurfVolCubicPrimary(in, model);
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

/*
 * @brief PRIMARY PARTICLE STREAM OUTPUT
 *
 * Writes a primary particle, along with its ID, to an output stream
 *
 * @param[in]        pri                 Primary particle object
 * @param[in,out]	 out		         Output binary stream
 * @param[in,out]    duplicates          Addresses of PAHs that have already been serialised
 *
 * @exception		 invalid_argument    Stream not ready
 */ 
void ModelFactory::WritePrimary(const AggModels::Primary &pri, std::ostream &out, void *duplicates)
{
    if (out.good()) {
        // Write the model Serial signature type to the stream.
        unsigned int type = (unsigned int)pri.AggID();
        out.write((char*)&type, sizeof(type));

        // Serialize the model object.
        const AggModels::PAHPrimary * pahPrimary = dynamic_cast<const AggModels::PAHPrimary*>(&pri);
        if(pahPrimary) {
            pahPrimary->Serialize(out, duplicates);
        }
        else
            pri.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::WritePrimary).");
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

        stats = new Stats::ParticleStats(in, model);

        return stats;
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ModelFactory::ReadStats).");
    }
}


// STREAM OUTPUT.

// Writes a model stats object, along with its ID, to an output stream.
void ModelFactory::WriteStats(const Stats::IModelStats &stats, std::ostream &out)
{
    if (out.good()) {
        // Write the model Serial signature type to the stream.
        unsigned int type = (unsigned int)stats.ID();
        out.write((char*)&type, sizeof(type));

        // Serialize the model object.
        stats.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ModelFactory::WriteStats).");
    }
}

// AGGREGATION MODEL CREATION.

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
        case AggModels::SurfVolHydrogen_ID:
            return new Stats::SurfVolHydrogenStats();
        case AggModels::SurfVolSilica_ID:
            return new Stats::SurfVolStats();
		case AggModels::PAH_KMC_ID:
            return new Stats::PAHStats();       // ms785: postprocessing not yet implemented
		case AggModels::BinTree_ID:
		    return new Stats::BinTreeStats();
        case AggModels::BinTreeSilica_ID:
	            return new Stats::BinTreeStats();
        case AggModels::SurfVolCubic_ID:
            return NULL;
        default:
            throw invalid_argument("Invalid model ID (Sweep, "
                                   "ModelFactory::CreateAggStats).");
    }
}


// AGGEGATION MODEL STREAM INPUT.

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
            case AggModels::SurfVolHydrogen_ID:
                stats = new Stats::SurfVolHydrogenStats(in, model);
                break;
            case AggModels::SurfVolSilica_ID:
                stats = new Stats::SurfVolStats(in, model);
                break;
			case AggModels::PAH_KMC_ID:
                stats = new Stats::PAHStats(in, model);
                break;
            case AggModels::BinTree_ID:
                stats = new Stats::BinTreeStats(in, model);
                break;
            case AggModels::BinTreeSilica_ID:
                stats = new Stats::BinTreeStats(in, model);
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
