/*!
 * @file    swp_bintree_stats.cpp
 * @author  William J Menz
 * @brief   Implementation of statistics calculator for BinTreePrimary
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      Implementation of statistics calculator for BinTreePrimary.
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public License
 *      as published by the Free Software Foundation; either version 2
 *      of the License, or (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this program; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *      02111-1307, USA.
 *
 *   Contact:
 *      Prof Markus Kraft
 *      Dept of Chemical Engineering
 *      University of Cambridge
 *      New Museums Site
 *      Pembroke Street
 *      Cambridge
 *      CB2 3RA, UK
 *
 *      Email:       mk306@cam.ac.uk
 *      Website:     http://como.cheng.cam.ac.uk
*/

#include "swp_bintree_stats.h"
#include "swp_bintree_primary.h"
#include "swp_particle.h"
#include "swp_ensemble.h"

using namespace Sweep;
using namespace Sweep::Stats;
using namespace std;


const std::string BinTreeStats::m_statnames[BinTreeStats::STAT_COUNT] = {
    std::string("Avg. Number of Primaries per Particle (-)"),
    std::string("Avg. Primary Diameter (m)"),
    std::string("Avg. Sintering Level (-)"),
    std::string("Avg. Sintering Rate (m2/s)"),
    std::string("Avg. Sintering Time (s)")
};

const IModelStats::StatType BinTreeStats::m_mask[BinTreeStats::STAT_COUNT] = {
    IModelStats::Avg,  // Avg.Number of primaries.
    IModelStats::Avg,  // Avg. Primary Particle diameter.
    IModelStats::Avg,  // Avg. Sintering level.
    IModelStats::Avg,  // Avg. sint rate,
    IModelStats::Avg   // Avg. sint time,
};

const std::string BinTreeStats::m_const_pslnames[BinTreeStats::PSL_COUNT] = {
    std::string("Number of Primaries (-)"),
    std::string("Avg. Primary Diameter (nm)"),
    std::string("Avg. Sintering Level (-)"),
    std::string("Total Sintering Time (s)")
};

//! Default constructor.
BinTreeStats::BinTreeStats()
: m_stats(STAT_COUNT,0.0)
{
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        m_names.push_back(m_statnames[i]);
    }

    for (unsigned int i=0; i!=PSL_COUNT; ++i) {
        m_pslnames.push_back(m_const_pslnames[i]);
    }
}

/*!
 * @brief           Stream-reading constructor.
 *
 * @param in        Input binary stream
 * @param model     Particle model
 * @return          Initialised stats object
 */
BinTreeStats::BinTreeStats(std::istream &in,
        const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

//! Default destructor
BinTreeStats::~BinTreeStats() {
    // Auto-generated destructor stub
}

//! Equate operator overload
BinTreeStats &BinTreeStats::operator=(const BinTreeStats &rhs)
{
    if (this != &rhs) {
        m_stats.assign(rhs.m_stats.begin(), rhs.m_stats.end());
    }
    return *this;
}


/*!
 * @brief           Calculates the model stats for a particle ensemble.
 * @param e         Ensemble to do the stats for
 * @param scale     Scale factor
 */
void BinTreeStats::Calculate(const Ensemble &e, real scale)
{
    // Empty the stats array.
    fill(m_stats.begin(), m_stats.end(), 0.0);

    // Calculate total weight
    real TotalWeight = e.Count()>0 ? e.GetSum(iW) : 0.0;
    real invTotalWeight = e.Count()>0 ? 1.0/e.GetSum(iW) : 0.0;

    // Loop over all particles, getting the stats from each.
    Ensemble::const_iterator ip;
    unsigned int n = 0;

    for (ip=e.begin(); ip!=e.end(); ++ip) {

        const AggModels::BinTreePrimary * const prim =
                dynamic_cast<const AggModels::BinTreePrimary*>((*ip)->Primary());

        real sz = (*ip)->Property(m_statbound.PID);
        real wt = (*ip)->getStatisticalWeight() * invTotalWeight;

        // Check if the value of the property is within the stats bound
        if ((m_statbound.Lower < sz) && (sz < m_statbound.Upper) ) {
            // Sum stats from this particle.
            m_stats[iNPrim]     += prim->GetNumPrimary()  * wt;
            m_stats[iPrimDiam]  += prim->GetPrimaryDiam() * wt
                    / prim->GetNumPrimary();
            m_stats[iSintLevel] += prim->GetAvgSinterLevel() * wt;
            m_stats[iSintRate]  += prim->GetSintRate() * wt;
            m_stats[iSintTime]  += prim->GetSintTime() * wt;
            ++n;
        }
    }
    // Scale the summed stats and calculate the averages.
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        if (m_mask[i] == Sum) {
            // Convert scale from 1/m3 to 1/cm3.
            m_stats[i] *= (scale * 1.0e-6 * TotalWeight);
        }
    }

}

// Returns a vector containing the stats.
void BinTreeStats::Get(fvector &stats, unsigned int start) const
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

// Adds to a vector containing stat names.
void BinTreeStats::Names(std::vector<std::string> &names,
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


// Returns a vector of PSL variable names.
void BinTreeStats::PSL_Names(std::vector<std::string> &names,
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
void BinTreeStats::PSL(const Sweep::Particle &sp, real time,
                       fvector &psl, unsigned int start) const
{
    // Resize vector if too small.
    if (start+PSL_Count()-1 >= psl.size()) {
        psl.resize(start+PSL_Count(), 0.0);
    }

    // Get an iterator to the first point of insertion in the
    // output stats array.
    fvector::iterator j = psl.begin()+start-1;

    const AggModels::BinTreePrimary* const prim =
        dynamic_cast<const AggModels::BinTreePrimary *>(sp.Primary());

    // Get the PSL stats.
    if (prim != NULL) {
        *(++j) = (real) prim->GetNumPrimary();
        *(++j) = prim->GetPrimaryDiam() * 1.0e9 / (real)(prim->GetNumPrimary());
        *(++j) = prim->GetAvgSinterLevel();
        *(++j) = prim->GetSintTime();

    } else {
        fill (j+1, j+2, 0.0);
    }
}


//! Creates a copy of the object.
BinTreeStats *const BinTreeStats::Clone(void) const
{
    return new BinTreeStats(*this);
}



// Writes the object to a binary stream.
void BinTreeStats::Serialize(std::ostream &out) const
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
                               "(Sweep, BinTreeStats::Serialize).");
    }
}

// Reads the object from a binary stream.
void BinTreeStats::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
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
                                    "(Sweep, BinTreeStats::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, BinTreeStats::Deserialize).");
    }
}

