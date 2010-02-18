/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the PAHStats class declared in the
    swp_PAH_stats.h header file.

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

#include "swp_PAH_stats.h"
#include "swp_aggmodel_type.h"
#include "swp_particle.h"
#include "swp_PAH_cache.h"
#include "swp_PAH_primary.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Stats;
using namespace std;

// STATIC CONST MEMBER VARIABLES.

const std::string PAHStats::m_statnames[PAHStats::STAT_COUNT] = {
    std::string("Avg. Number of PAHs"),
    std::string("Number of PAHs cm-3"),
    std::string("Surface real Part"),
    std::string("Avg. Mass real Part"),
    std::string("Avg. PAH real Part"),
    std::string("Avg. PAH Collision Diameter"),
	std::string("Avg. Number of Carbon Atoms"),
    std::string("Avg. Coalesc Threshold"),
    std::string("Num Primaries real Part"),

};

const IModelStats::StatType PAHStats::m_mask[PAHStats::STAT_COUNT] = {
    IModelStats::Avg,  // Avg. Number of PAHs.
    IModelStats::Sum,  // Number of PAHs.
    IModelStats::Sum,  // Surface real Part
    IModelStats::Avg,  // Avg. Mass real Part
    IModelStats::Avg,  // Avg. PAH real Part
    IModelStats::Avg,  // Avg. PAH Collision Diameter
	IModelStats::Avg,  // Avg. Number of Carbon atoms
    IModelStats::Avg,  // Avg. Coalesc Threshold
    IModelStats::Avg,  // Num Primaries real Part


};

const std::string PAHStats::m_const_pslnames[PAHStats::PSL_COUNT] = {
    std::string("Number of PAHs"),
    std::string("PAH Diameter"),
	std::string("Number of Carbon atoms"),
	std::string("Number primaries"),
	std::string("sqrt(LW)"),
	std::string("LdivW"),
	std::string("Avg. primary diameter"),
    std::string("Radius of gyration"),
    std::string("fdim"),

};



// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
PAHStats::PAHStats()
: m_stats(STAT_COUNT,0.0)
{
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        m_names.push_back(m_statnames[i]);
    }

    for (unsigned int i=0; i!=PSL_COUNT; ++i) {
        m_pslnames.push_back(m_const_pslnames[i]);
    }


}

// Copy constructor.
PAHStats::PAHStats(const Sweep::Stats::PAHStats &copy)
{
    *this = copy;
}

// Stream-reading constructor.
PAHStats::PAHStats(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
PAHStats::~PAHStats()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

PAHStats &PAHStats::operator=(const Sweep::Stats::PAHStats &rhs)
{
    if (this != &rhs) {
        m_stats.assign(rhs.m_stats.begin(), rhs.m_stats.end());
    }
    return *this;
}


// IMPLEMENTATION.

// Returns the number of basic particle stats.
unsigned int PAHStats::Count() const
{
    return STAT_COUNT;
}

// Calculates the model stats for a single particle.
void PAHStats::Calculate(const Particle &data)
{
    // Get surface-volume cache.
    const AggModels::PAHCache& cache =
        dynamic_cast<const AggModels::PAHCache&>(data.AggCache());

	m_stats[iNPAH]    = cache.m_numPAH;
    // Get stats.

}

// Calculates the model stats for a particle ensemble.
void PAHStats::Calculate(const Ensemble &e, real scale)
{
    // Empty the stats array.
    fill(m_stats.begin(), m_stats.end(), 0.0);

    // Loop over all particles, getting the stats from each.
    Ensemble::const_iterator ip;
    unsigned int n = 0;
    //particles with more then one PAH
    unsigned int nrealpart= 0;
    for (ip=e.begin(); ip!=e.end(); ++ip) {
        // Get surface-volume cache.
        const AggModels::PAHCache& cache =
            dynamic_cast<const AggModels::PAHCache&>((*ip)->AggCache());
		const AggModels::PAHPrimary *pah = NULL;
			pah = dynamic_cast<const AggModels::PAHPrimary*>((*(*ip)).Primary());
        real sz = cache.Parent()->Property(m_statbound.PID);
        // Check if the value of the property is within the stats bound
        if ((m_statbound.Lower < sz) && (sz < m_statbound.Upper) ) {
            // Sum stats from this particle.
			m_stats[iNPAH]    += cache.m_numPAH;
			m_stats[iPAHD]    += pah->PAHCollDiameter()*1e9;
			m_stats[iNCARB]	  += cache.m_numcarbon;
            m_stats[iNPAH+1]    += cache.m_numPAH;
            m_stats[iCOAL]    += cache.m_avg_coalesc;
			++n;
            if (cache.m_numPAH>1)
            {
                ++nrealpart;
                m_stats[iPARTSURF]+=(*ip)->SurfaceArea();
                m_stats[iNPRIM]+=cache.m_numprimary;
                m_stats[iPARTMASS]+=(*ip)->Primary()->Mass();
                m_stats[iNAVGPAH]+=cache.m_numPAH;
            }
        }
    }

    // Get the particle count.
    //real np    = (real)e.Count();
    real np    = (real) n;
    real invnp = (np>0) ? 1.0 / np : 0.0;

    // Scale the summed stats and calculate the averages.
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        if (m_mask[i] == Sum) {
            m_stats[i] *= (scale * 1.0e-6); // Convert scale from 1/m3 to 1/cm3.
        } else {
            if (i==iNPRIM || i==iPARTMASS || i==iNAVGPAH)
            {
                if (nrealpart>0)
                    m_stats[i] *= 1.0/nrealpart;
            }
            else
                m_stats[i] *= invnp;
        }
    }

}

// Returns a vector containing the stats.
const fvector &PAHStats::Get(void) const
{
    return m_stats;
}

// Returns a vector containing the stats.
void PAHStats::Get(fvector &stats, unsigned int start) const
{
    // Get an iterator to the first point of insertion in the
    // output stats array.
    fvector::iterator i;
    if (start < stats.size()) {
        i = stats.begin()+start;
    } else {
        i = stats.end();
    }

    // Add stats to output array until the end of that array is
    // reached or we have run out of stats to add.
    fvector::const_iterator j = m_stats.begin();
    for (; (i!=stats.end()) && (j!=m_stats.end()); ++i,++j) {
        *i = *j;
    }

    // If we still have stats to add, but the output array is full,
    // then we need to add elements to the end of the array.
    for (; j!=m_stats.end(); ++j) {
        stats.push_back(*j);
    }
}

// Returns a vector containing the stat names.
const std::vector<std::string> &PAHStats::Names(void) const
{
    return m_names;
}

// Adds to a vector containing stat names.
void PAHStats::Names(std::vector<std::string> &names,
                         unsigned int start) const
{
    // Get an iterator to the first point of insertion in the
    // output names array.
    std::vector<std::string>::iterator i;
    if (start < names.size()) {
        i = names.begin()+start;
    } else {
        i = names.end();
    }

    // Add stats to output array until the end of that array is
    // reached or we have run out of names to add.
    std::vector<std::string>::const_iterator j = m_names.begin();
    for (; (i!=names.end()) && (j!=m_names.end()); ++i,++j) {
        *i = *j;
    }

    // If we still have names to add, but the output array is full,
    // then we need to add elements to the end of the array.
    while(j!=m_names.end()) {
        names.push_back(*j);
        ++j;
    }
}




// PARTICLE SIZE LISTS.

// Returns the number of PSL output variables.
unsigned int PAHStats::PSL_Count(void) const
{
    return PSL_COUNT;
}

// Returns a vector of PSL variable names.
void PAHStats::PSL_Names(std::vector<std::string> &names,
                             unsigned int start) const
{
    // Get an iterator to the first point of insertion in the
    // output names array.
    std::vector<std::string>::iterator i;
    if (start < names.size()) {
        i = names.begin()+start;
    } else {
        i = names.end();
    }

    // Add stats to output array until the end of that array is
    // reached or we have run out of names to add.
    std::vector<std::string>::const_iterator j = m_pslnames.begin();
    for (; (i!=names.end()) && (j!=m_pslnames.end()); ++i,++j) {
        *i = *j;
    }

    // If we still have names to add, but the output array is full,
    // then we need to add elements to the end of the array.
    while(j!=m_pslnames.end()) {
        names.push_back(*j);
        ++j;
    }
}

// Returns the particle size list (PSL) entry for particle i
// in the given ensemble.
void PAHStats::PSL(const Ensemble &ens, unsigned int i,
                       real time, fvector &psl, unsigned int start) const
{
    // Get particle.
    const Sweep::Particle *const sp = ens.At(i);

    if (sp != NULL) {
        PSL(*sp, time, psl, start);
    } else {
        // Resize vector if too small.
        if (start+PSL_Count()-1 >= psl.size()) {
            psl.resize(start+PSL_Count(), 0.0);
        }
        // Clear data points.
        fill(psl.begin()+start, psl.begin()+start+PSL_Count()-1, 0.0);
    }
}

// Returns the PSL entry for the given particle.
void PAHStats::PSL(const Sweep::Particle &sp, real time,
                       fvector &psl, unsigned int start) const
{
    // Resize vector if too small.
    if (start+PSL_Count()-1 >= psl.size()) {
        psl.resize(start+PSL_Count(), 0.0);
    }

    // Get an iterator to the first point of insertion in the
    // output stats array.
    fvector::iterator j = psl.begin()+start-1;

    // Get surface-volume cache.
    const AggModels::PAHCache* cache =
        dynamic_cast<const AggModels::PAHCache*>(&sp.AggCache());

    // Get the PSL stats.
    if (cache != NULL) {
		*(++j) = (real)(cache->m_numPAH);
		*(++j) = (real)(cache->m_PAHDiameter)*1e9;			//convert to nm
		*(++j) = (real) (cache->m_numcarbon);
		*(++j) = (real) (cache->m_numprimary);
		*(++j) = (real) (cache->m_sqrtLW);
		*(++j) = (real) (cache->m_LdivW);
		*(++j) = (real) (cache->m_primarydiam)*1e9/(real)(cache->m_numprimary);//convert to nm
        *(++j) = (real) (cache->m_Rg)*1e9;
        *(++j) = (real) (cache->m_fdim);

    } else {
        fill (j+1, j+2, 0.0);
    }
}





// READ/WRITE/COPY.

// Creates a copy of the object.
PAHStats *const PAHStats::Clone(void) const
{
    return new PAHStats(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
unsigned int PAHStats::ID(void) const
{
    return (unsigned int)AggModels::PAH_ID;
}

// Writes the object to a binary stream.
void PAHStats::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output number of stats in vector.
        unsigned int n = (unsigned int)m_stats.size();
        out.write((char*)&n, sizeof(n));

        // Output stats.
        double v = 0.0;
        for (unsigned int i=0; i!=n; ++i) {
            v = m_stats[i];
            out.write((char*)&v, sizeof(v));
        }

        // Output number of stat names in vector.
        n = (unsigned int)m_names.size();
        out.write((char*)&n, sizeof(n));

        // Output names.
        for (unsigned int i=0; i!=n; ++i) {
            // Write the length of the stat name to the stream.
            unsigned int m = m_names[i].length();
            out.write((char*)&m, sizeof(m));

            // Write the stat name to the stream.
            out.write(m_names[i].c_str(), m);
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PAHStats::Serialize).");
    }
}

// Reads the object from a binary stream.
void PAHStats::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    // TODO:  Deserialize ParticleStats should reset to state with no components
    //        or tracker variables in the first instance.
    m_stats.clear();
    m_names.clear();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double val = 0.0;
        char *name = NULL;

        switch (version) {
            case 0:

                // Read number of stats in vector.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read stats.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_stats.push_back((real)val);
                }

                // Read number of stat names in vector.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read names.
                for (unsigned int i=0; i!=n; ++i) {
                    // Read the length of the name.
                    unsigned int m = 0;
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));

                    // Read the name.
                    name = new char[m];
                    in.read(name, m);
                    m_names.push_back(string(name, m));
                    delete [] name;
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, PAHStats::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PAHStats::Deserialize).");
    }
}
