/*
  Author(s):      Markus Sander (ms785)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The PAHStats class is a specialization of the IModelStats
    interface which produces stats from the basic properties of the
    simple primary-particle aggregation model.

    The stats stored in this class are:


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

#ifndef SWEEP_PAH_STATS_H
#define SWEEP_PAH_STATS_H

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
class PAHStats : public IModelStats
{
public:
    // Constructors.
    PAHStats(void);  // Default constructor.
    PAHStats(const PAHStats &copy); // Copy constructor.
    PAHStats(                         // Stream-reading constructor.
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

    // Destructor.
    ~PAHStats(void);

    // Operators.
    PAHStats &operator=(const PAHStats &rhs);


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

	//! Get primary particle details and connectivity
	void PrintPrimary(const Sweep::Particle &sp, std::vector<fvector> &nodes, std::vector<fvector> &primaries, int k) const;

    // READ/WRITE/COPY.

    // Creates a copy of the object.
    PAHStats *const Clone(void) const;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    unsigned int ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

private:
    // Stats count and indices.
    static const unsigned int STAT_COUNT = 12; //modified by hdy, used to be 12
    enum StatIndices {
        iNPAH=0,      // m_numPAH
        iPARTSURF=2,  // m_surf, the surface area for primary partilcle
        iPARTMASS=3,  // m_mass, the mass for primary partilcle
        iNAVGPAH=4,   // Avg. PAH double Part
        iPAHD=5,      // m_PAHCollDiameter
        iNCARB=6,     // m_numcarbon
        iNHYDROGEN=7, // m_numH
        iNEDGEC=8,    // m_numOfEdgeC
        iNRINGS=9,    // m_numOfRings
        iSINT=10,     // m_avg_sinter, if primary coordinates are tracked 
        iNPRIM=11,    // m_numprimary
        iSQRTLW=12,   // sqrt(LW)
        iLDIVW=13,    // LdivW
        iavgdim=14,   // m_primarydiam/m_numprimary, Avg. primary diameter
        irgyr=15,     // m_Rg, Radius of gyration
        ifdim=16,     // m_fdim
		//iCOAL=17      // m_avg_coalesc, if primary coordinates are not tracked
    };

    // PSL count and indices.
    static const unsigned int PSL_COUNT  = 17;
    static const unsigned int PPSL_COUNT = 0;

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

};
};
};

#endif
