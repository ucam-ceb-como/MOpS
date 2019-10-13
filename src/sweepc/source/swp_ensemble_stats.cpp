/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the EnsembleStats class declared in the
    swp_ensemble_stats.h header file.

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

#include "swp_ensemble_stats.h"
#include "swp_model_factory.h"
#include "swp_mechanism.h"

#include <stdexcept>

using namespace Sweep;
using namespace std;
using namespace Sweep::Stats;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
EnsembleStats::EnsembleStats()
: m_basicstats(NULL)
{
}

// Default constructor (public).
EnsembleStats::EnsembleStats(const Sweep::ParticleModel &model)
{
    m_basicstats = new ParticleStats(model);
    m_aggstats   = ModelFactory::CreateAggStats(model.AggModel(), model);
}

// Copy constructor.
EnsembleStats::EnsembleStats(const EnsembleStats &copy)
{
    *this = copy;
}

// Default destructor.
EnsembleStats::~EnsembleStats()
{
    releaseMem();
}


// OPERATOR OVERLOADS.

// Assignment operator.
EnsembleStats &EnsembleStats::operator=(const EnsembleStats &rhs)
{
    if (this != &rhs) {
        releaseMem();

        // Copy basic stats.
        m_basicstats = new ParticleStats(*rhs.m_basicstats);

        // Copy aggregation model stats.
        if (rhs.m_aggstats != NULL) {
            m_aggstats = rhs.m_aggstats->Clone();
        } else {
            m_aggstats = NULL;
        }
    }
    return *this;
}


// IMPLEMENTATION.

// Returns the number of stats for this model.
unsigned int EnsembleStats::Count(void) const
{
    unsigned int n = m_basicstats->Count();
    if (m_aggstats!=NULL)
        n += m_aggstats->Count();

    return n;
}

// Calculates the model stats for a particle ensemble.
void EnsembleStats::Calculate(const Ensemble &e, double scale)
{
    m_basicstats->Calculate(e, scale);
    if (m_aggstats!=NULL)
        m_aggstats->Calculate(e, scale);
}

// Returns a vector containing the stats.
void EnsembleStats::Get(fvector &stats, unsigned int start) const
{
    // Get basic stats.
    m_basicstats->Get(stats, start);
    start += m_basicstats->Count();

    // Get aggregation model stats.
    if (m_aggstats!=NULL) {
        m_aggstats->Get(stats, start);
        start += m_aggstats->Count();
    }
}

// Returns a vector containing the stat names.
const std::vector<std::string> &EnsembleStats::Names(void) const
{
    static vector<string> names;
    Names(names);
    return names;
}

// Adds to a vector containing stat names.
void EnsembleStats::Names(std::vector<std::string> &names, unsigned int start) const
{
    // Get basic stats.
    m_basicstats->Names(names, start);
    start += m_basicstats->Count();

    // Get aggregation model stats.
    if (m_aggstats!=NULL) {
        m_aggstats->Names(names, start);
        start += m_aggstats->Count();
    }
}


// SUB-STATS INTERFACE.

// Returns the basic stats object.
const ParticleStats &EnsembleStats::BasicStats(void) const
{
    return *m_basicstats;
}


// PARTICLE SIZE LISTS.

// Returns the number of PSL output variables.
unsigned int EnsembleStats::PSL_Count(void) const
{
    unsigned int n = 1 + m_basicstats->PSL_Count();
    if (m_aggstats!=NULL)
        n += m_aggstats->PSL_Count();

    return n;
}

// Returns a vector of PSL variable names.
void EnsembleStats::PSL_Names(std::vector<std::string> &names, unsigned int start) const
{
    // First PSL variable is the particle weight.
    if (names.size() > start) {
        names[start] = "Weight (cm-3)";
    } else {
        names.resize(start+1);
        names[start] = "Weight (cm-3)";
    }
    ++start;

    // Get basic stats.
    m_basicstats->PSL_Names(names, start);
    start += m_basicstats->PSL_Count();

    // Get aggregation stats names.
    if (m_aggstats!=NULL) {
        m_aggstats->PSL_Names(names, start);
        start += m_aggstats->PSL_Count();
    }
}

// Returns the PSL entry for the given particle.  The particle weight
// is set to 1.0 because there is only one particle.
void EnsembleStats::PSL(const Sweep::Particle &sp, const Sweep::ParticleModel& model,
		                double time, fvector &psl, double scale)
{
    unsigned int start = 1;

    // First PSL variable is the particle weight.
    if (psl.size() > 0) {
        psl[0] = sp.getStatisticalWeight() * scale*1.0e-6;
    } else {
        psl.push_back(sp.getStatisticalWeight() * scale*1.0e-6); // m-3 to cm-3.
    }

    // Get basic stats.
    m_basicstats->PSL(sp, time, psl, start);
    start += m_basicstats->PSL_Count();

    // Get aggregation model PSL stats.
    if (m_aggstats!=NULL) {
        m_aggstats->PSL(sp, time, psl, start);
        start += m_aggstats->PSL_Count();
    }
}

// READ/WRITE/COPY.

// Writes the object to a binary stream.
void EnsembleStats::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize the basic stats.
        m_basicstats->Serialize(out);

        // Write if aggregation model stats exist.
        if (m_aggstats == NULL) {
            out.write((char*)&falseval, sizeof(falseval));
        } else {
            out.write((char*)&trueval, sizeof(trueval));

            // Serialise the aggregation model stats.
            ModelFactory::WriteAggStats(*m_aggstats, out);
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, EnsembleStats::Serialize).");
    }
}

// Reads the object from a binary stream.
void EnsembleStats::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    // TODO:  Need some way to set to a valid state here.
    releaseMem();

    const unsigned int trueval = 1;

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;

        switch (version) {
            case 0:
                // Deserialize the basic stats.
                m_basicstats = new ParticleStats(in, model);

                // Read if this ensemble stats has an aggregation
                // stats.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                if (n == trueval) {
                    // Deserialize the aggregation model stats.
                    m_aggstats = ModelFactory::ReadAggStats(in, model);
                } else {
                    m_aggstats = NULL;
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, EnsembleStats::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, EnsembleStats::Deserialize).");
    }
}

// Set statbound to this class and its children stats
void EnsembleStats::SetStatBoundary(const Sweep::Stats::IModelStats::StatBound &sb) {
    m_statbound = sb;
    m_basicstats->SetStatBoundary(sb);
    if (m_aggstats) m_aggstats->SetStatBoundary(sb);
}

// MEMORY MANAGEMENT.

// Clears all memory associated with the EnsembleStats object.
void EnsembleStats::releaseMem(void)
{
    delete m_basicstats;
    m_basicstats = NULL;

    delete m_aggstats;
    m_aggstats = NULL;
}

// Get primary particle details and connectivity
void EnsembleStats::PrintPrimary(const Sweep::Particle &sp, const Sweep::ParticleModel& model, vector<fvector> &nodes, std::vector<fvector> &primaries, int k) const
{
	if (m_aggstats != NULL){
		m_aggstats->PrintPrimary(sp,nodes, primaries, k);
	}
}
