/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ParticleStats class is a specialization of the IModelStats
    interface which produces stats from the basic properties of the
    ParticleData object.  It does not consider any particle sub-models.

    The stats stored in this class are:

    1.      Particle count.
    2.      Number density (M0).
    3.      Average particle diameter.
    4.      Average particle collision diameter.
    5.      Average particle mobility diameter.
    6/7.    Total and average particlesurface area.
    8/9.    Total and average particle volume.
    10/11.  Total and average particle mass.

    Additionally the total and average values of all particle components
    and tracker values are appended to the end of the array.  For this 
    reason it is necessary to initialise the class with knowledge of
    the components and tracker variables.
*/

#ifndef SWEEP_PARTICLESTATS_H
#define SWEEP_PARTICLESTATS_H

#include "swp_params.h"
#include "swp_modelstats.h"
#include "swp_particledata.h"
#include "swp_ensemble.h"
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_modeltype.h"
#include <vector>
#include <string>
#include <iostream>

namespace Sweep
{
class ParticleStats : public IModelStats
{
public:
    // Constructors.
    ParticleStats( // Default constructor.
        const CompPtrVector &comp,  // Components.
        const TrackPtrVector &track // Tracker variables.
        );
    ParticleStats(const ParticleStats &copy); // Copy constructor.
    ParticleStats(std::istream &in);          // Stream-reading constructor.

    // Destructor.
    ~ParticleStats(void);

    // Operators.
    ParticleStats &operator=(const ParticleStats &rhs);


    // IMPLEMENTATION.

    // Returns the number of stats for this model.
    unsigned int Count(void) const;

    // Returns the stats mask which informs whether stats should
    // be summed, averaged or both.
//    const std::vector<StatType> &Mask(void) const;

    // Calculates the model stats for a single particle.
    void Calculate(const ParticleData &data);

    // Calculates the model stats for a particle ensemble.
    void Calculate(
        const Ensemble &e, // Ensemble from which to get stats.
        real scale = 1.0   // Scaling factor to unit volume (summed stats).
        );

    // Returns a vector containing the stats.
    const fvector &Get(void) const;

    // Returns a vector containing the stats.
    void Get(
        fvector &stats,        // Output vector.
        unsigned int start = 0 // Optional start index for the first stat.
        ) const;

    // Returns a vector containing the stat names.
    const std::vector<std::string> &Names(void) const;

    // Adds to a vector containing stat names.
    void Names(
        std::vector<std::string> &names, // Output vector.
        unsigned int start = 0           // Optional start index for the first stat.
        ) const;


    // AVAILABLE BASIC STATS.

    // Returns the particle count.
    real PCount(void) const;

    // Returns the number density.
    real M0(void) const;

    // Returns the avg. equiv. sphere diameter.
    real AvgDiam(void) const;

    // Returns the avg. collision diameter.
    real AvgCollDiam(void) const;

    // Returns the avg. mobility diameter.
    real AvgMobDiam(void) const;

    // Returns the total surface area.
    real SurfaceArea(void) const;

    // Returns the avg. surface area.
    real AvgSurfaceArea(void) const;

    // Returns the total volume.
    real Fv(void) const;

    // Returns the average volume.
    real AvgVolume(void) const;

    // Returns the total mass.
    real Mass(void) const;

    // Returns the average mass.
    real AvgMass(void) const;


    // PARTICLE SIZE LISTS.

    // Returns the number of PSL output variables.
    unsigned int PSL_Count(void) const;

    // Returns a vector of PSL variable names.
    void PSL_Names(
        std::vector<std::string> &names, // Vector in which to return names.
        unsigned int start = 0 // Optional start index for the first name.
        ) const;

    // Returns the particle size list (PSL) entry for particle i
    // in the given ensemble.
    void PSL(
        const Ensemble &ens,   // Ensemble from which to get properties.
        unsigned int i,        // Index of particle in ensemble to get.
        real time,             // The current time.
        fvector &psl,          // Output vector.
        unsigned int start = 0 // Optional start index for the first variable.
        ) const;


    // READ/WRITE/COPY.

    // Creates a copy of the object.
    ParticleStats *const Clone(void) const;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    ModelType ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

private:
    // Stats count and indices.
    static const unsigned int STAT_COUNT = 11;
    enum StatIndices {iNP=0, iM0=1, iD=2, iDcol=3, 
                      iDmob=4, iS=5, iV=7, iM=9};

    // PSL count and indices.
    static const unsigned int PSL_COUNT = 7;

    // Component and tracker counts.
    unsigned int m_ncomp, m_ntrack;

    // The stats.
    fvector m_stats;

    // The stat names.
    static const std::string m_statnames[STAT_COUNT];
    std::vector<std::string> m_names;

    // The stat mask.
    static const IModelStats::StatType m_mask[STAT_COUNT];

    // The PSL names.
    static const std::string m_const_pslnames[PSL_COUNT];
    std::vector<std::string> m_pslnames;

    // Cannot instantiate this class without knowing how many 
    // components and tracker variables there are for
    // each particle.
    ParticleStats(void); // Default constructor.
};
};

#endif
