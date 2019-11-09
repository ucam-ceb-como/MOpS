/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ParticleStats class is a specialization of the IModelStats
    interface which produces stats from the basic properties of the
    ParticleData object.  It does not consider any particle sub-models.

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

#ifndef SWEEP_PARTICLESTATS_H
#define SWEEP_PARTICLESTATS_H

#include "swp_params.h"
#include "swp_model_stats.h"
#include "swp_ensemble.h"
#include "swp_particle_model.h"

#include <vector>
#include <string>
#include <iostream>

namespace Sweep
{
namespace Stats
{
/*!
The stats stored in this class are:

1.      Particle count.
2.      Number density (M0).
3.      Average particle diameter.
4.      Average particle collision diameter.
5.      Average particle mobility diameter.
6/7.    Total and average particle surface area.
8/9.    Total and average particle volume.
10/11.  Total and average particle mass.
12.     Total of squared masses
13.     Total of cubed masses
14      Average number of coagulation events experienced

Additionally the total and average values of all particle components
and tracker values are appended to the end of the array.  For this
reason it is necessary to initialise the class with knowledge of
the components and tracker variables.
*/
class ParticleStats : public IModelStats
{
public:
    // Constructors.
    ParticleStats(const Sweep::ParticleModel &model); // Default constructor.
    ParticleStats(const ParticleStats &copy); // Copy constructor.
    ParticleStats(                        // Stream-reading constructor.
        std::istream &in,                 //  - Input stream
        const Sweep::ParticleModel &model //  - Defining particle model.
        );

    // Destructor.
    ~ParticleStats(void);

    // Operators.
    ParticleStats &operator=(const ParticleStats &rhs);


    // IMPLEMENTATION.

    // Returns the number of stats for this model.
    unsigned int Count(void) const;

    // Calculates the model stats for a particle ensemble.
    void Calculate(
        const Ensemble &e, // Ensemble from which to get stats.
        double scale         // Scaling factor to unit volume (summed stats).
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

    //! Returns the particle count.
    double PCount(void) const;

    //! Returns the number density.
    double M0(void) const;

    // Returns the avg. equiv. sphere diameter.
    double AvgDiam(void) const;

    // Returns the avg. collision diameter.
    double AvgCollDiam(void) const;

    // Returns the avg. mobility diameter.
    double AvgMobDiam(void) const;

    // Returns the total surface area.
    double SurfaceArea(void) const;

    // Returns the avg. surface area.
    double AvgSurfaceArea(void) const;

    // Returns the total volume.
    double Fv(void) const;

    // Returns the average volume.
    double AvgVolume(void) const;

    // Returns the total mass.
    double Mass(void) const;

    // Returns the average mass.
    double AvgMass(void) const;


    // PARTICLE SIZE LISTS.

    // Returns the number of PSL output variables.
    unsigned int PSL_Count(void) const;

    // Returns a vector of PSL variable names.
    void PSL_Names(
        std::vector<std::string> &names, // Vector in which to return names.
        unsigned int start = 0 // Optional start index for the first name.
        ) const;

    //! Build the PSL entry for the given particle.
    void PSL(
        const Sweep::Particle &sp,      // Particle from which to get PSL data.
        double time,                      // Current flow time (used to calculate particle age).
        fvector &psl,                   // Output vector.
        unsigned int start = 0          // Optional start index for the first variable.
        ) const;

	// Get primary particle details and connectivity (only used for PAHPrimary and BintreePrimary)
	void PrintPrimary(const Sweep::Particle &sp, std::vector<fvector> &nodes, std::vector<fvector> &primaries, int k) const;

    // READ/WRITE/COPY.

    // Creates a copy of the object.
    ParticleStats *const Clone(void) const;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    unsigned int ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,                 // Input stream
        const Sweep::ParticleModel &model // Defining particle model.
        );

private:
    // Stats count and indices.
    static const unsigned int STAT_COUNT = 15;
    enum StatIndices {iNP=0, iM0=1, iD=2, iDcol=3,
                      iDmob=4, iS=5, iV=7, iM=9,
                      iM2 = 11, iM3 = 12,
                      iCoag=13, iMaxCoag=14};

    // PSL count and indices.
    static const unsigned int PSL_COUNT = 9;

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
}

#endif
