/*
  Author(s):      Shraddha Shekar (ss663)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the silicaStats class declared in the
    swp_silica_stats.h header file.

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

#include "swp_silica_stats.h"
#include "swp_aggmodel_type.h"
#include "swp_particle.h"
#include "swp_silica_cache.h"
#include "swp_silica_primary.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Stats;
using namespace std;

// STATIC CONST MEMBER VARIABLES.

const std::string SilicaStats::m_statnames[SilicaStats::STAT_COUNT] = {
    std::string("Number of Si atoms (cm-3)"),
	std::string("Number of O atoms (cm-3)"),
	std::string("Number of OH atoms (cm-3)"),
    std::string("Avg. Number of Primaries per Particle (-)"),
	std::string("Avg. Primary Diameter (nm)"),
    std::string("Avg. Sintering Level (-)"),
	std::string("Avg. Particle Mass (kg)"),
	std::string("Avg. Sintering Rate (m2/s)"),
};

const IModelStats::StatType SilicaStats::m_mask[SilicaStats::STAT_COUNT] = {
    IModelStats::Sum,  // Number of Si atoms.
	IModelStats::Sum,  // Number of O atoms.
	IModelStats::Sum,  // Number of OH atoms.
    IModelStats::Avg,  // Avg.Number of primaries.
    IModelStats::Avg,  // Avg. Primary Particle diameter.
    IModelStats::Avg,  // Avg. Sintering level.
	IModelStats::Avg,  // Avg particle mass,
    IModelStats::Avg,  // Avg sint rate,

};

const std::string SilicaStats::m_const_pslnames[SilicaStats::PSL_COUNT] = {
    std::string("Number of Si atoms"),
	std::string("Number of O atoms"),
	std::string("Number of OH atoms"),
    std::string("Number of primaries"),
	std::string("Avg. primary diameter"),
	std::string("Avg. Sintering Level"),
	//std::string("Avg. Particle Mass"),

};



// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SilicaStats::SilicaStats()
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
SilicaStats::SilicaStats(const Sweep::Stats::SilicaStats &copy)
{
    *this = copy;
}

// Stream-reading constructor.
SilicaStats::SilicaStats(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
SilicaStats::~SilicaStats()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

SilicaStats &SilicaStats::operator=(const Sweep::Stats::SilicaStats &rhs)
{
    if (this != &rhs) {
        m_stats.assign(rhs.m_stats.begin(), rhs.m_stats.end());
    }
    return *this;
}


// IMPLEMENTATION.

// Returns the number of basic particle stats.
unsigned int SilicaStats::Count() const
{
    return STAT_COUNT;
}
/* TODO: remove? wjm34
// Calculates the model stats for a single particle.
void SilicaStats::Calculate(const Particle &data)
{
    // Get surface-volume cache.
    const AggModels::SilicaCache& cache = dynamic_cast<const AggModels::SilicaCache&>(data.AggCache());

	m_stats[iNSi] = cache.m_numSi;
	m_stats[iNO] = cache.m_numO;
	m_stats[iNOH] = cache.m_numOH;
	m_stats[iPRIMDIAM] = cache.m_primarydiam*1e9;

	m_stats[iCOAL] = cache.m_avg_sinter;
	m_stats[iNPRIM] = cache.m_numprimary;

} */

// Calculates the model stats for a particle ensemble.
void SilicaStats::Calculate(const Ensemble &e, real scale)
{
    // Empty the stats array.
    fill(m_stats.begin(), m_stats.end(), 0.0);

    // Calculate total weight
    real TotalWeight = e.Count()>0 ? e.GetSum(iW) : 0.0;
    real invTotalWeight = e.Count()>0 ? 1.0/e.GetSum(iW) : 0.0;

    // Loop over all particles, getting the stats from each.
    Ensemble::const_iterator ip;
    unsigned int n = 0;
    // particles with more then one silica
    unsigned int nrealpart= 0;
    for (ip=e.begin(); ip!=e.end(); ++ip) {

        // Get data from silica cache
        const AggModels::SilicaCache& cache =
            dynamic_cast<const AggModels::SilicaCache&>((*ip)->AggCache());

        real sz = (*ip)->Property(m_statbound.PID);
        real wt = (*ip)->getStatisticalWeight() * invTotalWeight;

        // Check if the value of the property is within the stats bound
        if ((m_statbound.Lower < sz) && (sz < m_statbound.Upper) ) {
            // Sum stats from this particle.
            m_stats[iNSi]   += (cache.m_numSi * wt);
            m_stats[iNO]    += (cache.m_numO  * wt);
            m_stats[iNOH]   += (cache.m_numOH  * wt);
            m_stats[iCOAL]    += (cache.m_avg_sinter  * wt);
            m_stats[iPRIMDIAM] += (cache.m_primarydiam * 1e9  * wt /cache.m_numprimary);
            m_stats[iSintRate] += cache.m_sintrate * wt;
            ++n;
            if (cache.m_numSi>1)
            {
                ++nrealpart;
                m_stats[iNPRIM]+=cache.m_numprimary  * wt;
                if((*ip)->Primary()!=NULL)
                {
                m_stats[iPARTMASS]+=(*ip)->Primary()->Mass() * wt;
                }
                else
                {
                    m_stats[iPARTMASS]+=0;
                }

            }
        }
    }
    // Scale the summed stats and calculate the averages.
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        if (m_mask[i] == Sum) {
            m_stats[i] *= (scale * 1.0e-6 * TotalWeight); // Convert scale from 1/m3 to 1/cm3.
        }
    // Don't need to scale by number of particles as this is included
    //     in the weighting scaling above.
    }

}

// Returns a vector containing the stats.
const fvector &SilicaStats::Get(void) const
{
    return m_stats;
}

// Returns a vector containing the stats.
void SilicaStats::Get(fvector &stats, unsigned int start) const
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
const std::vector<std::string> &SilicaStats::Names(void) const
{
    return m_names;
}

// Adds to a vector containing stat names.
void SilicaStats::Names(std::vector<std::string> &names,
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
unsigned int SilicaStats::PSL_Count(void) const
{
    return PSL_COUNT;
}

// Returns a vector of PSL variable names.
void SilicaStats::PSL_Names(std::vector<std::string> &names,
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

// Returns the PSL entry for the given particle.
void SilicaStats::PSL(const Sweep::Particle &sp, real time,
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
    const AggModels::SilicaCache* cache =
        dynamic_cast<const AggModels::SilicaCache*>(&sp.AggCache());

    // Get the PSL stats.
    if (cache != NULL) {
		*(++j) = (real)(cache->m_numSi)/(real)(cache->m_numprimary);
		*(++j) = (real)(cache->m_numO)/(real)(cache->m_numprimary);
		*(++j) = (real)(cache->m_numOH)/(real)(cache->m_numprimary);
		*(++j) = (real) (cache->m_numprimary);
		//*(++j) = (real) (cache->m_sqrtLW);
		//*(++j) = (real) (cache->m_LdivW);
		*(++j) = (real) (cache->m_primarydiam)*1e9/(real)(cache->m_numprimary);//convert to nm
		*(++j) = (real) (cache->m_avg_sinter);
        //*(++j) = (real) (cache->);
        //*(++j) = (real) (cache->m_Rg)*1e9;
        //*(++j) = (real) (cache->m_fdim);

    } else {
        fill (j+1, j+2, 0.0);
    }
}





// READ/WRITE/COPY.

// Creates a copy of the object.
SilicaStats *const SilicaStats::Clone(void) const
{
    return new SilicaStats(*this);
}


// Returns the model data type.  Used to identify different models
// and for serialisation.
unsigned int SilicaStats::ID(void) const
{
    return (unsigned int)AggModels::Silica_ID;
}

// Writes the object to a binary stream.
void SilicaStats::Serialize(std::ostream &out) const
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
                               "(Sweep, SilicaStats::Serialize).");
    }
}

// Reads the object from a binary stream.
void SilicaStats::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
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
                                    "(Sweep, SilicaStats::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SilicaStats::Deserialize).");
    }
}
