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
    std::string("Avg. Sintering Time (s)"),
    std::string("GStdev of Collision Diameter (-)"),
    std::string("GStdev of Avg. Primary Diameter (-)"),
    std::string("Mean GStdev of Primary Diameter (-)"),
    std::string("GMean of Collision Diameter (m)"),
    std::string("GMean of Avg. Primary Diameter (m)"),
	std::string("Avg. Radius of Gyration (m)")
};

const IModelStats::StatType BinTreeStats::m_mask[BinTreeStats::STAT_COUNT] = {
    IModelStats::Avg,  // Avg.Number of primaries.
    IModelStats::Avg,  // Avg. Primary Particle diameter.
    IModelStats::Avg,  // Avg. Sintering level.
    IModelStats::Avg,  // Avg. sint rate,
    IModelStats::Avg,  // Avg. sint time,
    IModelStats::Avg,  // Gstdev of mean dcol
    IModelStats::Avg,  // Gstdev of mean dpri
    IModelStats::Avg,  // Mean gstdev of dpri
    IModelStats::Avg,  // Gmean of dcol
    IModelStats::Avg,  // Gmean of avg dpri
	IModelStats::Avg   // Radius of gyration
};

const std::string BinTreeStats::m_const_pslnames[BinTreeStats::PSL_COUNT] = {
    std::string("Number of Primaries (-)"),
    std::string("Avg. Primary Diameter (nm)"),
    std::string("Avg. Sintering Level (-)"),
    std::string("Total Sintering Time (s)"),
    std::string("Arithmetic Stdev of Primary Diameter (nm)"),
    std::string("Geometric Mean of Primary Diameter (nm)"),
    std::string("Geometric Stdev of Primary Diameter (-)"),
    std::string("Radius of Gyration (nm)")
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
void BinTreeStats::Calculate(const Ensemble &e, double scale)
{
    // Empty the stats array.
    fill(m_stats.begin(), m_stats.end(), 0.0);

    // Calculate total weight
    double TotalWeight = (e.Count() + e.GetTotalParticleNumber())>0 ? (e.GetSum(iW) + e.GetTotalParticleNumber()) : 0.0;
    double invTotalWeight = (e.Count() + e.GetTotalParticleNumber())>0 ? 1.0/(e.GetSum(iW) + e.GetTotalParticleNumber()) : 0.0;

    // Loop over all particles, getting the stats from each.
    Ensemble::const_iterator ip;
    unsigned int n = 0;

    // Create some vectors to store diameter lists
    // Keep dcol in [0], dpri in [1] of d, theses are pushed into diams.
    std::vector<fvector> diams;
    fvector weights;
    fvector d;

    for (ip=e.begin(); ip!=e.end(); ++ip) {

        const AggModels::BinTreePrimary * const prim =
                dynamic_cast<const AggModels::BinTreePrimary*>((*ip)->Primary());

        double sz = (*ip)->Property(m_statbound.PID);
        double wt = (*ip)->getStatisticalWeight() * invTotalWeight;

        //! Check if the value of the property is within the stats bound.
        if ((m_statbound.Lower < sz) && (sz < m_statbound.Upper) ) {
            // Sum stats from this particle.
            m_stats[iNPrim]     += prim->GetNumPrimary()  * wt;
            m_stats[iPrimDiam]  += prim->GetPrimaryDiam() * wt
                    / (double) prim->GetNumPrimary();
            m_stats[iSintLevel] += prim->GetAvgSinterLevel() * wt;
            m_stats[iSintRate]  += prim->GetSintRate() * wt;
            m_stats[iSintTime]  += prim->GetSintTime() * wt;
            m_stats[iGStdevMean]+= prim->GetPrimaryGStdDev() * wt;
            m_stats[iRg]         += prim->GetRadiusOfGyration() * wt;

            // Collect the collision and primary diameters
            d.push_back(prim->CollDiameter());
            d.push_back(prim->GetPrimaryDiam() / (double) prim->GetNumPrimary());
            diams.push_back(d);
            weights.push_back(wt);
            d.clear();

            ++n;
        }
    }

    // Add contributions from the hybrid particle model
    if (e.GetTotalParticleNumber() != 0)
    {
        double wt = e.GetTotalParticleNumber() * invTotalWeight;
        double dpri = e.GetTotalDiameter() * invTotalWeight;
		
        m_stats[iNPrim] += wt;
        m_stats[iPrimDiam] += dpri;
        m_stats[iSintLevel] += wt;
        m_stats[iSintRate] += 0;
        m_stats[iSintTime] += 0;
        m_stats[iGStdevMean] += wt; 
    }

    // Now get the geometric standard devs, using [0] for dcol, [1] for dpri
    // Default to 1.0 GSTDEV (Can't have GSTDEV=0)
    fvector gstdevs, gmeans;
    GetGeometricStdev(2u, diams, weights, gmeans, gstdevs);
    m_stats[iCollGStdev] = gstdevs[0];
    m_stats[iPrimGStdev] = gstdevs[1];
    m_stats[iCollGMean] = gmeans[0];
    m_stats[iPrimGMean] = gmeans[1];

    // Scale the summed stats and calculate the averages.
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        if (m_mask[i] == Sum) {
            // Convert scale from 1/m3 to 1/cm3.
            m_stats[i] *= (scale * 1.0e-6 * TotalWeight);
        }
    }

}


/*!
 * Gets the geometric standard deviation of a vector of (vector of)
 * diameters. Currently does so for collision and primary diameters.
 * Also assumes here that the weights sum up to 1.0.
 *
 * @param num           Number of diameter types to calculate gstdev for
 * @param diams         Vector of length number of particles, containing
 *                          a fvector of diameters for that particle.
 * @param weights       Vector of length number of particles storing weights
 * @param gmeans        Output vector for gmeans
 * @param gstdevs       Output vector for gstdevs
 */
void BinTreeStats::GetGeometricStdev(
        const unsigned int num,
        std::vector<fvector> diams,
        fvector weights,
        fvector &gmeans,
        fvector &gstdevs) const {

    // Some checks first
    if (diams.size() != weights.size())
        throw std::runtime_error("Failed getting weights and diameters "
                "in BinTreeStats::GetGeometricStdev()");
    if (diams.size() < size_t(1u)) {
        // Don't calculate if there aren't enough particles.
        gmeans.resize(num, 0.0);
        gstdevs.resize(num, 1.0);
    } else {
        unsigned int i(0u); // Iterate number of particles
        unsigned int j(0u); // Iterate diameter types

        // Then we must calculate the geometric means
        gmeans.resize(num, 1.0);

        for (i = 0; i != diams.size(); i++) {
            // Loop over diameter types
            for (j = 0; j != diams[i].size(); j++) {
                gmeans[j] *= pow(diams[i].at(j), weights[i]);
            }
        }

        // Now we can get the geometric stdevs
        gstdevs.resize(num, 0.0);
        double dev(0.0);
        for (i = 0; i != diams.size(); i++) {
            // Loop over diameter types
            for (j = 0; j != num; j++) {
                dev = log(diams[i].at(j) / gmeans[j]);
                gstdevs[j] += weights[i] * dev * dev;
            }
        }

        // Just need to do a bit more to the sums...
        for (j = 0; j != num; j++) {
            gstdevs[j] = exp(sqrt(gstdevs[j]));
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
void BinTreeStats::PSL(const Sweep::Particle &sp, double time,
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
        *(++j) = (double) prim->GetNumPrimary();
        *(++j) = prim->GetPrimaryDiam() * 1.0e9 / (double)(prim->GetNumPrimary());
        *(++j) = prim->GetAvgSinterLevel();
        *(++j) = prim->GetSintTime();
        *(++j) = 1.0e9 * prim->GetPrimaryAStdDev();
        *(++j) = 1.0e9 * prim->GetPrimaryGMean();
        *(++j) = prim->GetPrimaryGStdDev();
        *(++j) = prim->GetRadiusOfGyration() * 1.0e9;

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
                    m_stats.push_back((double)val);
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

// Return primary particle details and connectivity
void BinTreeStats::PrintPrimary(const Sweep::Particle &sp, std::vector<fvector> &nodes, std::vector<fvector> &primaries, int k) const
{
    const AggModels::BinTreePrimary* const prim =
        dynamic_cast<const AggModels::BinTreePrimary *>(sp.Primary());

    if (prim != NULL) {
        prim->PrintPrimary(nodes, primaries, k);
    }
}
