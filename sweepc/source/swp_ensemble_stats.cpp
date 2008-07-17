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
#include "swp_pripart_stats.h"
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

        // Copy model stats.
        for (ModelStatsMap::const_iterator i=rhs.m_modelstats.begin(); 
             i!=rhs.m_modelstats.end(); ++i) {
            m_modelstats[i->first] = i->second->Clone();
        }
    }
    return *this;
}


// IMPLEMENTATION.

// Returns the number of stats for this model.
unsigned int EnsembleStats::Count(void) const
{
    unsigned int n = m_basicstats->Count();
    if (m_aggstats!=NULL) n += m_aggstats->Count();
    for (ModelStatsMap::const_iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        n += i->second->Count();
    }
    return n;
}

// Calculates the model stats for a particle ensemble.
void EnsembleStats::Calculate(const Ensemble &e, real scale)
{
    m_basicstats->Calculate(e, scale);
    if (m_aggstats!=NULL) m_aggstats->Calculate(e, scale);
    for (ModelStatsMap::iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        i->second->Calculate(e, scale);
    }
}

// Returns a vector containing the stats.
const fvector &EnsembleStats::Get(void) const
{
    static fvector stats;
    Get(stats);
    return stats;
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

    // Get stats from particle sub-models.
    for (ModelStatsMap::const_iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        i->second->Get(stats, start);
        start += i->second->Count();
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

    // Get stats from particle sub-models.
    for (ModelStatsMap::const_iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        i->second->Names(names, start);
        start += i->second->Count();
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
    if (m_aggstats!=NULL) n += m_aggstats->PSL_Count();
    for (ModelStatsMap::const_iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        n += i->second->PSL_Count();
    }
    return n;
}

// Returns a vector of PSL variable names.
void EnsembleStats::PSL_Names(std::vector<std::string> &names, unsigned int start) const
{
    // First PSL variable is the particle weight.
    if (names.size() > start) {
        names[start] = "Weight";
    } else {
        names.resize(start+1);
        names[start] = "Weight";
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

    // Get stats from particle sub-models.
    for (ModelStatsMap::const_iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        i->second->PSL_Names(names, start);
        start += i->second->PSL_Count();
    }
}

// Returns the particle size list (PSL) entry for particle i
// in the given ensemble.
void EnsembleStats::PSL(const Ensemble &ens, unsigned int i, 
                        real time, fvector &psl, real scale) const
{
    unsigned int start = 1;

    // First PSL variable is the particle weight.
    if (psl.size() > 0) {
        psl[0] = scale*1.0e-6;
    } else {
        psl.push_back(scale*1.0e-6); // m-3 to cm-3.
    }

    // Get basic stats.
    m_basicstats->PSL(ens, i, time, psl, start);
    start += m_basicstats->PSL_Count();

    // Get aggregation model PSL stats.
    if (m_aggstats!=NULL) {
        m_aggstats->PSL(ens, i, time, psl, start);
        start += m_aggstats->PSL_Count();
    }

    // Get stats from particle sub-models.
    for (ModelStatsMap::const_iterator j=m_modelstats.begin(); 
         j!=m_modelstats.end(); ++j) {
        j->second->PSL(ens, i, time, psl, start);
        start += j->second->PSL_Count();
    }
}

// Returns the PSL entry for the given particle.  The particle weight
// is set to 1.0 because there is only one particle.
void EnsembleStats::PSL(const Sweep::ParticleCache &sp, real time, fvector &psl)
{
    unsigned int start = 1;

    // First PSL variable is the particle weight.  This is set to
    // 1.0 for a single particle, because no information about the
    // ensemble is known.
    if (psl.size() > 0) {
        psl[0] = 1.0;
    } else {
        psl.push_back(1.0);
    }

    // Get basic stats.
    m_basicstats->PSL(sp, time, psl, start);
    start += m_basicstats->PSL_Count();

    // Get aggregation model PSL stats.
    if (m_aggstats!=NULL) {
        m_aggstats->PSL(sp, time, psl, start);
        start += m_aggstats->PSL_Count();
    }

    // Get stats from particle sub-models.
    for (ModelStatsMap::const_iterator j=m_modelstats.begin(); 
         j!=m_modelstats.end(); ++j) {
        j->second->PSL(sp, time, psl, start);
        start += j->second->PSL_Count();
    }
}


// PRIMARY-PARTICLE SIZE LISTS.

// Returns true if the current particle model allows
// primary-particle size lists to be generated.
bool EnsembleStats::GeneratesPPSL(void) const
{
    return (m_aggstats->ID() == (unsigned int)AggModels::PriPartList_ID);
}

// Returns the number of primary-PSL output variables.
unsigned int EnsembleStats::PPSL_Count(void) const
{
    if (m_aggstats->ID() == (unsigned int)AggModels::PriPartList_ID) {
        return static_cast<const PriPartStats*>(m_aggstats)->PPSL_Count();
    } else {
        return 0;
    }
}

// Returns a vector of primary-PSL variable names.
void EnsembleStats::PPSL_Names(std::vector<std::string> &names,
                               unsigned int start) const
{
    if (m_aggstats->ID() == (unsigned int)AggModels::PriPartList_ID) {
        static_cast<const PriPartStats*>(m_aggstats)->PPSL_Names(names, start);
    }
}

// Returns the primary-particle size list (PSL) entry for particle i
// in the given ensemble.
void EnsembleStats::PPSL(const Ensemble &ens, unsigned int i, real time,
                         vector<fvector> &ppsl, real scale) const
{
    if (m_aggstats->ID() == (unsigned int)AggModels::PriPartList_ID) {
        static_cast<const PriPartStats*>(m_aggstats)->PPSL(ens, i, time, ppsl);
    } else {
        ppsl.clear();
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


        // Output the number of sub-models.
        unsigned int n = (unsigned int)m_modelstats.size();
        out.write((char*)&n, sizeof(n));

        // Serialize the sub-model stats.
        for (ModelStatsMap::const_iterator i=m_modelstats.begin(); 
             i!=m_modelstats.end(); ++i) {
            ModelFactory::WriteStats(*i->second, out);
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

                // Read the number of sub-models.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Deserialize the sub-model stats.
                for (unsigned int i=0; i!=n; ++i) {
                    IModelStats *stats = ModelFactory::ReadStats(in, model);
                    m_modelstats[(SubModels::SubModelType)stats->ID()] = stats;
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


// MEMORY MANAGEMENT.

// Clears all memory associated with the EnsembleStats object.
void EnsembleStats::releaseMem(void)
{
    delete m_basicstats;
    m_basicstats = NULL;

    delete m_aggstats;
    m_aggstats = NULL;

    for (ModelStatsMap::iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        delete i->second;
    }
    m_modelstats.clear();
}
