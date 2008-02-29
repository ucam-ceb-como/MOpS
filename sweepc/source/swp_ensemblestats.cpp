#include "swp_ensemblestats.h"
#include "swp_modelfactory.h"
#include "swp_mechanism.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
EnsembleStats::EnsembleStats()
: m_basicstats(NULL)
{
}

// Default constructor (public).
EnsembleStats::EnsembleStats(const Mechanism &mech)
{
    m_basicstats = new ParticleStats(mech.Components(), mech.Trackers());
}

// Copy constructor.
EnsembleStats::EnsembleStats(const Sweep::EnsembleStats &copy)
{
    *this = copy;
}

// Default destructor.
EnsembleStats::~EnsembleStats()
{
    delete m_basicstats;
    for (ModelStatsMap::iterator i=m_modelstats.begin(); i!=m_modelstats.end(); ++i) {
        delete i->second;
    }
}


// OPERATOR OVERLOADS.

// Assignment operator.
EnsembleStats &EnsembleStats::operator=(const Sweep::EnsembleStats &rhs)
{
    if (this != &rhs) {
        releaseMem();

        // Copy basic stats.
        m_basicstats = new ParticleStats(*rhs.m_basicstats);

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

    // Get stats from particle sub-models.
    start += m_basicstats->Count();
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

    // Get stats from particle sub-models.
    start += m_basicstats->Count();
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
    unsigned int n = m_basicstats->PSL_Count();
    for (ModelStatsMap::const_iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        n += i->second->PSL_Count();
    }
    return n;
}

// Returns a vector of PSL variable names.
void EnsembleStats::PSL_Names(std::vector<std::string> &names) const
{
    unsigned int start = 1;

    // First PSL variable is the particle weight.
    if (names.size() > 0) {
        names[0] = "Weight";
    } else {
        names.push_back("Weight");
    }

    // Get basic stats.
    m_basicstats->PSL_Names(names, start);

    // Get stats from particle sub-models.
    start += m_basicstats->Count();
    for (ModelStatsMap::const_iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        i->second->PSL_Names(names, start);
        start += i->second->Count();
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

    // Get stats from particle sub-models.
    start += m_basicstats->Count();
    for (ModelStatsMap::const_iterator j=m_modelstats.begin(); 
         j!=m_modelstats.end(); ++j) {
        j->second->PSL(ens, i, time, psl, start);
        start += j->second->Count();
    }
}


// READ/WRITE/COPY.

// Writes the object to a binary stream.
void EnsembleStats::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize the basic stats.
        m_basicstats->Serialize(out);

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
void EnsembleStats::Deserialize(std::istream &in)
{
    // TODO:  Need some way to set to a valid state here.
    releaseMem();

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
                m_basicstats = new ParticleStats(in);

                // Read the number of sub-models.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Deserialize the sub-model stats.
                for (unsigned int i=0; i!=n; ++i) {
                    IModelStats *stats = ModelFactory::ReadStats(in);
                    m_modelstats[stats->ID()] = stats;
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

    for (ModelStatsMap::iterator i=m_modelstats.begin(); 
         i!=m_modelstats.end(); ++i) {
        delete i->second;
    }
    m_modelstats.clear();
}
