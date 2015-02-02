/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ModelFactory is a factory class for sweep particle models.  It
    provides routines for creating, reading and writing particle model
    data objects.

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

#ifndef SWEEP_MODELFACTORY_H
#define SWEEP_MODELFACTORY_H

#include "swp_params.h"
#include "swp_primary.h"
#include "swp_particle_model.h"
#include "swp_model_stats.h"
#include "swp_aggmodel_type.h"

#include <iostream>

namespace Sweep
{
//! Factory class for particles, their statistics and sub-models
class ModelFactory
{
public:
    // PRIMARY PARTICLE CREATION.

    //! Creates a new primary particle of the given type.
    static AggModels::Primary *const CreatePrimary(
        const AggModels::AggModelType id, // Model ID.
        const double time,                  // Primary create time.
        const ParticleModel &model        // Defining particle model.
        );

    //! Creates a new primary particle of the given type.
    static AggModels::Primary *const CreatePrimary(
        const AggModels::AggModelType id, // Model ID.
        const double time,                  // Primary create time.
        const double position,
        const ParticleModel &model        // Defining particle model.
        );

    // PRIMARY PARTICLE STREAM INPUT.

    //! Reads a primary particle from a binary stream. First reads the primary type ID, in order to create a primary of the correct type.
    static AggModels::Primary *const ReadPrimary(
        std::istream &in,              // Input stream
        const ParticleModel &model,    // Defining particle model
        void *duplicates               // Information on duplicate PAHs
        );

    // PRIMARY PARTICLE STREAM OUTPUT.

    //! Writes a primary particle, along with its ID, to an output stream
    static void WritePrimary(
        const AggModels::Primary &pri,    // Primary to write
        std::ostream &out,                // Output stream
        void *duplicates                  // Information on duplicated PAHs
        );


    // Reads model stats from a binary stream.  The first item read
    // is the model ID which tells the ModelFactory what type
    // of model stats to read.
    static Stats::IModelStats *const ReadStats(
        std::istream &in,          // Input stream.
        const ParticleModel &model // Defining particle model.
        );


    // Writes a model stats object, along with its ID, to an output stream.
    static void WriteStats(
        const Stats::IModelStats &stats, // Stats to write.
        std::ostream &out                // Output stream.
        );


    // AGGREGATION MODEL CREATION.

    // Creates a new aggregation model stats object of the given type.
    static Stats::IModelStats *const CreateAggStats(
        AggModels::AggModelType id, // Model ID.
        const ParticleModel &model  // Defining particle model.
        );


    // AGGEGATION MODEL STREAM INPUT.

    // Reads aggregation model stats from a binary stream.  The first
    // item read is the model ID which tells the ModelFactory what type
    // of aggregation model stats to read.
    static Stats::IModelStats *const ReadAggStats(
        std::istream &in,          // Input stream.
        const ParticleModel &model // Defining particle model.
        );


    // AGGREGATION MODEL STREAM OUTPUT.

    // Writes an aggregation model stats object, along with
    // its ID, to an output stream.
    static void WriteAggStats(
        const Stats::IModelStats &stats, // Stats to write.
        std::ostream &out                // Output stream.
        );

};
}

#endif
