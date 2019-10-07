/*
  Author(s):      Robert I A Patterson Robert.Patterson@wias-berlin.de
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2012 Robert I A Patterson

  File purpose:
    Implementation of the SurfVolHydrogenStats class declared in the
    swp_surfvol_stats.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#include "swp_surfvolhydrogen_stats.h"
#include "swp_surfvolhydrogen_primary.h"
#include "swp_aggmodel_type.h"
#include "swp_particle.h"

#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Stats;
using namespace std;

// STATIC CONST MEMBER VARIABLES.

const std::string SurfVolHydrogenStats::m_statnames[SurfVolHydrogenStats::STAT_COUNT] = {
    std::string("Equiv. Sphere Surface Area (cm2/cm3)"),
    std::string("Avg. Equiv. Sphere Surface Area (cm2)"),
    std::string("Est. Primary Particle Count"),
    std::string("Avg. Est. Primary Particle Count"),
    std::string("Avg. Est. Primary Particle Diameter (nm)")
};

const IModelStats::StatType SurfVolHydrogenStats::m_mask[SurfVolHydrogenStats::STAT_COUNT] = {
    IModelStats::Sum,  // Equiv. sphere surface area.
    IModelStats::Avg,  // Avg. equiv. sphere surface area.
    IModelStats::Sum,  // Est. primary particle count.
    IModelStats::Avg,  // Avg. est. primary particle count.
    IModelStats::Avg   // Avg. est. primary particle diameter.
};

const std::string SurfVolHydrogenStats::m_const_pslnames[SurfVolHydrogenStats::PSL_COUNT] = {
    std::string("Equiv. Sphere Surface Area (cm2)"),
    std::string("Est. Primary Count"),
    std::string("Est. Avg. Primary Diameter (nm)")
};


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SurfVolHydrogenStats::SurfVolHydrogenStats()
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
SurfVolHydrogenStats::SurfVolHydrogenStats(const Sweep::Stats::SurfVolHydrogenStats &copy)
{
    *this = copy;
}

// Stream-reading constructor.
SurfVolHydrogenStats::SurfVolHydrogenStats(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
SurfVolHydrogenStats::~SurfVolHydrogenStats()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

SurfVolHydrogenStats &SurfVolHydrogenStats::operator=(const Sweep::Stats::SurfVolHydrogenStats &rhs)
{
    if (this != &rhs) {
        m_stats.assign(rhs.m_stats.begin(), rhs.m_stats.end());
    }
    return *this;
}


// IMPLEMENTATION.

// Returns the number of basic particle stats.
unsigned int SurfVolHydrogenStats::Count() const
{
    return STAT_COUNT;
}

// Calculates the model stats for a particle ensemble.
void SurfVolHydrogenStats::Calculate(const Ensemble &e, double scale)
{
    fill(m_stats.begin(), m_stats.end(), 0.0);

    // Loop over all particles, getting the stats from each.
    Ensemble::const_iterator ip;

    for (ip=e.begin(); ip!=e.end(); ++ip) {
        // Get surface-volume cache.
        const AggModels::SurfVolHydrogenPrimary * const primary =
            dynamic_cast<const AggModels::SurfVolHydrogenPrimary *>((*ip)->Primary());
        double sz = (*ip)->Property(m_statbound.PID);
        double wt = (*ip)->getStatisticalWeight();

        // Check if the value of the property is within the stats bound
        if ((m_statbound.Lower < sz) && (sz < m_statbound.Upper) ) {
            // Sum stats from this particle.
            m_stats[iS]      += (*ip)->SphSurfaceArea() * 1.0e4 * wt; // Convert from m2 to cm2.
            m_stats[iS+1]    += (*ip)->SphSurfaceArea() * 1.0e4 * wt; // Convert from m2 to cm2.
            m_stats[iPPN]    += primary->PP_Count() * wt;
            m_stats[iPPN+1]  += primary->PP_Count() * wt;
            m_stats[iPPD]    += primary->PP_Diameter() * 1.0e9 * wt; // Convert from m to nm.
        }
    }

    // Calculate total weight
    double invTotalWeight = e.Count()>0 ? 1.0/e.GetSum(iW) : 0.0;

    // Scale the summed stats and calculate the averages.
    for (unsigned int i=1; i!=STAT_COUNT; ++i) {
        if (m_mask[i] == Sum) {
            m_stats[i] *= (scale * 1.0e-6); // Convert scale from 1/m3 to 1/cm3.
        } else {
            m_stats[i] *= invTotalWeight;
        }
    }
}

// Returns a vector containing the stats.
const fvector &SurfVolHydrogenStats::Get(void) const
{
    return m_stats;
}

// Returns a vector containing the stats.
void SurfVolHydrogenStats::Get(fvector &stats, unsigned int start) const
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
const std::vector<std::string> &SurfVolHydrogenStats::Names(void) const
{
    return m_names;
}

// Adds to a vector containing stat names.
void SurfVolHydrogenStats::Names(std::vector<std::string> &names,
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


// AVAILABLE BASIC STATS.

// Returns the total equivalent-sphere surface area.
double SurfVolHydrogenStats::SphSurfaceArea(void) const {return m_stats[iS];}

// Returns the average equivalent-sphere surface area.
double SurfVolHydrogenStats::AvgSphSurfaceArea(void) const {return m_stats[iS+1];}

// Returns the total primary particle count.
double SurfVolHydrogenStats::PriPartCount(void) const {return m_stats[iPPN];}

// Returns the average primary particle count.
double SurfVolHydrogenStats::AvgPriPartCount(void) const {return m_stats[iPPN+1];}

// Returns the average primary particle diameter.
double SurfVolHydrogenStats::AvgPriPartDiameter(void) const {return m_stats[iPPD];}


// PARTICLE SIZE LISTS.

// Returns the number of PSL output variables.
unsigned int SurfVolHydrogenStats::PSL_Count(void) const
{
    return PSL_COUNT;
}

// Returns a vector of PSL variable names.
void SurfVolHydrogenStats::PSL_Names(std::vector<std::string> &names,
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
void SurfVolHydrogenStats::PSL(const Sweep::Particle &sp, double time,
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
    const AggModels::SurfVolHydrogenPrimary* const primary =
        dynamic_cast<const AggModels::SurfVolHydrogenPrimary *>(sp.Primary());

    // Get the PSL stats.
    if (primary != NULL) {
        *(++j) = sp.SphSurfaceArea() * 1.0e4; // m2 to cm2.
        *(++j) = primary->PP_Count();
        *(++j) = primary->PP_Diameter() * 1.0e9; // m to nm.
    } else {
        fill(j+1, j+3, 0.0);
    }
}

// READ/WRITE/COPY.

// Creates a copy of the object.
SurfVolHydrogenStats *const SurfVolHydrogenStats::Clone(void) const
{
    return new SurfVolHydrogenStats(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
unsigned int SurfVolHydrogenStats::ID(void) const
{
    return (unsigned int)AggModels::SurfVol_ID;
}

// Writes the object to a binary stream.
void SurfVolHydrogenStats::Serialize(std::ostream &out) const
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
                               "(Sweep, SurfVolHydrogenStats::Serialize).");
    }
}

// Reads the object from a binary stream.
void SurfVolHydrogenStats::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
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
                                    "(Sweep, SurfVolHydrogenStats::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SurfVolHydrogenStats::Deserialize).");
    }
}

// Get primary particle details and connectivity (only used for PAHPrimary and BintreePrimary)
void SurfVolHydrogenStats::PrintPrimary(const Sweep::Particle &sp, std::vector<fvector> &nodes, std::vector<fvector> &primaries, int k) const
{
}
