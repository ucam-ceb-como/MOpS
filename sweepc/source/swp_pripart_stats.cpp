#include "swp_pripart_stats.h"
#include "swp_aggmodel_type.h"
#include "swp_particle.h"
#include "swp_pripart_cache.h"
#include "swp_pripart_primary.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Stats;
using namespace std;

// STATIC CONST MEMBER VARIABLES.

const std::string PriPartStats::m_statnames[PriPartStats::STAT_COUNT] = {
    std::string("Primary Particle Count"),
    std::string("Avg. Primary Particle Count"),
    std::string("Avg. Primary Particle Diameter (nm)"),
    std::string("Equiv. Sphere Surface Area (cm2/cm3)"),
    std::string("Avg. Equiv. Sphere Surface Area (cm2)"),
    std::string("Primary Surface Area (cm2/cm3)"),
    std::string("Avg. Primary Surface Area (cm2)")
};

const IModelStats::StatType PriPartStats::m_mask[PriPartStats::STAT_COUNT] = {
    IModelStats::Sum,  // Primary particle count.
    IModelStats::Avg,  // Avg. primary particle count.
    IModelStats::Avg,  // Avg. primary particle diameter.
    IModelStats::Sum,  // Equiv. sphere surface area.
    IModelStats::Avg,  // Avg. equiv. sphere surface area.
    IModelStats::Sum,  // Primary surface area.
    IModelStats::Avg   // Avg. primary surface area.
};

const std::string PriPartStats::m_const_pslnames[PriPartStats::PSL_COUNT] = {
    std::string("Primary Count"),
    std::string("Avg. Primary Diameter (nm)"),
    std::string("Equiv. Sphere Surface Area (cm2)"),
    std::string("Primary Surface Area (cm2)")
};

const std::string PriPartStats::m_const_ppslnames[PriPartStats::PPSL_COUNT] = {
    std::string("Monomer Count"),
    std::string("Diameter (nm)"),
    std::string("Surface Area (cm2)"),
    std::string("Volume (cm3)"),
    std::string("Mass (g)")
};

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
PriPartStats::PriPartStats()
: m_stats(STAT_COUNT,0.0)
{
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        m_names.push_back(m_statnames[i]);
    }

    for (unsigned int i=0; i!=PSL_COUNT; ++i) {
        m_pslnames.push_back(m_const_pslnames[i]);
    }

    for (unsigned int i=0; i!=PPSL_COUNT; ++i) {
        m_ppslnames.push_back(m_const_ppslnames[i]);
    }
}

// Copy constructor.
PriPartStats::PriPartStats(const Sweep::Stats::PriPartStats &copy)
{
    *this = copy;
}

// Stream-reading constructor.
PriPartStats::PriPartStats(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
PriPartStats::~PriPartStats()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

PriPartStats &PriPartStats::operator=(const Sweep::Stats::PriPartStats &rhs)
{
    if (this != &rhs) {
        m_stats.assign(rhs.m_stats.begin(), rhs.m_stats.end());
    }
    return *this;
}


// IMPLEMENTATION.

// Returns the number of basic particle stats.
unsigned int PriPartStats::Count() const 
{
    return STAT_COUNT;
}

// Calculates the model stats for a single particle.
void PriPartStats::Calculate(const ParticleCache &data)
{
    // Get surface-volume cache.
    const AggModels::PriPartCache* cache = 
        dynamic_cast<const AggModels::PriPartCache*>(data.AggCache());

    // Get stats.
    m_stats[iPPN]    = cache->Count();
    m_stats[iPPN+1]  = m_stats[iPPN];
    m_stats[iPPD]    = cache->AvgPriDiameter() * 1.0e9; // Convert from m to nm.
    m_stats[iS]      = cache->SphSurfaceArea() * 1.0e4; // Convert from m2 to cm2.
    m_stats[iS+1]    = m_stats[iS];
    m_stats[iPPS]    = cache->PriSurfaceArea() * 1.0e4; // Convert from m2 to cm2.
    m_stats[iPPS+1]  = m_stats[iS];
}

// Calculates the model stats for a particle ensemble.
void PriPartStats::Calculate(const Ensemble &e, real scale)
{
    // Empty the stats array.
    fill(m_stats.begin(), m_stats.end(), 0.0);

    // Loop over all particles, getting the stats from each.
    Ensemble::const_iterator ip;
    for (ip=e.begin(); ip!=e.end(); ++ip) {
        // Get surface-volume cache.
        const AggModels::PriPartCache* cache = 
            dynamic_cast<const AggModels::PriPartCache*>((*ip)->AggCache());

        // Sum stats from this particle.
        m_stats[iPPN]    += cache->Count();
        m_stats[iPPN+1]  += cache->Count();
        m_stats[iPPD]    += cache->AvgPriDiameter() * cache->Count() * 1e9; // Convert from m to nm.
        m_stats[iS]      += cache->SphSurfaceArea() * 1.0e4; // Convert from m2 to cm2.
        m_stats[iS+1]    += cache->SphSurfaceArea() * 1.0e4; // Convert from m2 to cm2.
        m_stats[iS]      += cache->PriSurfaceArea() * 1.0e4; // Convert from m2 to cm2.
        m_stats[iS+1]    += cache->PriSurfaceArea() * 1.0e4; // Convert from m2 to cm2.
    }
    
    // Get the particle count.
    real np    = (real)e.Count();
    real invnp = (np>0) ? 1.0 / np : 0.0;

    // Scale the summed stats and calculate the averages.
    for (unsigned int i=1; i!=STAT_COUNT; ++i) {
        if (m_mask[i] == Sum) {
            m_stats[i] *= (scale * 1.0e-6); // Convert scale from 1/m3 to 1/cm3.
        } else {
            m_stats[i] *= invnp;
        }
    }

    // Calculate the correct primary particle diameter, which is averaged
    // over the primaries, not the aggregates.
    if (m_stats[iPPN] > 0.0) {
        m_stats[iPPD] /= (invnp * m_stats[iPPN]);
    } else {
        m_stats[iPPD] = 0.0;
    }
}

// Returns a vector containing the stats.
const fvector &PriPartStats::Get(void) const
{
    return m_stats;
}

// Returns a vector containing the stats.
void PriPartStats::Get(fvector &stats, unsigned int start) const
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
const std::vector<std::string> &PriPartStats::Names(void) const
{
    return m_names;
}

// Adds to a vector containing stat names.
void PriPartStats::Names(std::vector<std::string> &names,
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

// Returns the total primary particle count.
real PriPartStats::PriPartCount(void) const {return m_stats[iPPN];}

// Returns the average primary particle count.
real PriPartStats::AvgPriPartCount(void) const {return m_stats[iPPN+1];}

// Returns the average primary particle diameter.
real PriPartStats::AvgPriPartDiameter(void) const {return m_stats[iPPD];}

// Returns the total equivalent-sphere surface area.
real PriPartStats::SphSurfaceArea(void) const {return m_stats[iS];}

// Returns the average equivalent-sphere surface area.
real PriPartStats::AvgSphSurfaceArea(void) const {return m_stats[iS+1];}

// Returns the total primary-particle surface area.
real PriPartStats::PriSurfaceArea(void) const {return m_stats[iPPS];}

// Returns the average primary-particle surface area.
real PriPartStats::AvgPriSurfaceArea(void) const {return m_stats[iPPS+1];}

// PARTICLE SIZE LISTS.

// Returns the number of PSL output variables.
unsigned int PriPartStats::PSL_Count(void) const
{
    return PSL_COUNT;
}

// Returns a vector of PSL variable names.
void PriPartStats::PSL_Names(std::vector<std::string> &names,
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
void PriPartStats::PSL(const Ensemble &ens, unsigned int i, 
                       real time, fvector &psl, unsigned int start) const
{
    // Get particle.
    const Sweep::ParticleCache *const sp = ens.At(i);
    
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

// Returns the PSL entry for the given particle cache.
void PriPartStats::PSL(const Sweep::ParticleCache &sp, real time,
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
    const AggModels::PriPartCache* cache = 
        dynamic_cast<const AggModels::PriPartCache*>(sp.AggCache());

    // Get the PSL stats.
    if (cache != NULL) {
        *(++j) = (real)(cache->Count());
        *(++j) = cache->AvgPriDiameter() * 1.0e9; // m to nm.
        *(++j) = cache->SphSurfaceArea() * 1.0e4; // m2 to cm2.
        *(++j) = cache->PriSurfaceArea() * 1.0e4; // m2 to cm2.
    } else {
        fill (j+1, j+4, 0.0);
    }
}


// PRIMARY-PARTICLE SIZE LIST.

// Returns the number of PPSL output variables.
unsigned int PriPartStats::PPSL_Count(void) const
{
    return PPSL_COUNT;
}

// Returns a vector of PPSL (primary-particle size list) variable names.
void PriPartStats::PPSL_Names(std::vector<std::string> &names,
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
    std::vector<std::string>::const_iterator j = m_ppslnames.begin();
    for (; (i!=names.end()) && (j!=m_ppslnames.end()); ++i,++j) {
        *i = *j;
    }

    // If we still have names to add, but the output array is full,
    // then we need to add elements to the end of the array.
    while(j!=m_ppslnames.end()) {
        names.push_back(*j);
        ++j;
    }
}

// Returns the primary-particle size list (PPSL) entry for particle i
// in the given ensemble.
void PriPartStats::PPSL(const Ensemble &ens, unsigned int i, 
                        real time, std::vector<fvector> &ppsl, 
                        unsigned int start) const
{
    // It only makes sense to call this function if the sub-particle
    // tree is not in use.
    if (ens.ParticleModel() == NULL) return;
    if (ens.ParticleModel()->UseSubPartTree()) return;

    // Cannot get the primary-particle list if this particle
    // does not have a primary (a double-check for sub-particle tree
    // because this statement can only be true if the sub-particle
    // tree is enabled).
    if (ens.At(i)->Primary() == NULL) return;

    // Get primary-particle data.
    const AggModels::PriPartPrimary* pp = 
        dynamic_cast<const AggModels::PriPartPrimary*>(ens.At(i)->Primary());

    // Resize vector of vectors to fit all primaries.
    ppsl.resize(pp->PriCount());

    // Loop over all primary particles.
    for (unsigned int k=0; k<pp->PriCount(); ++k) {
        // Get vector for this primary particle.
        vector<fvector>::reference psl = ppsl.at(k);

        // Resize vector if too small.
        if (start+PSL_Count()-1 >= psl.size()) {
            psl.resize(start+PPSL_Count(), 0.0);
        }

        // Get an iterator to the first point of insertion in the
        // primary output stats array.
        fvector::iterator j = psl.begin()+start-1;

        // Get the PPSL stats.
        *(++j) = (real)pp->MonomerCount(k);
        *(++j) = pp->PriDiameter(k);
        *(++j) = pp->PriSurface(k);
        *(++j) = pp->PriVolume(k);
        *(++j) = pp->PriMass(k);
    }
}


// READ/WRITE/COPY.

// Creates a copy of the object.
PriPartStats *const PriPartStats::Clone(void) const
{
    return new PriPartStats(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
unsigned int PriPartStats::ID(void) const 
{
    return (unsigned int)AggModels::PriPartList_ID;
}

// Writes the object to a binary stream.
void PriPartStats::Serialize(std::ostream &out) const
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
                               "(Sweep, PriPartStats::Serialize).");
    }
}

// Reads the object from a binary stream.
void PriPartStats::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
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
                                    "(Sweep, PriPartStats::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PriPartStats::Deserialize).");
    }
}
