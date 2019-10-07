/*!
 * @file    swp_bintree_stats.h
 * @author  William J Menz
 * @brief   Calculates statistics for BinTreePrimary model
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      Calculates statistics for BinTreePrimary model
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

#ifndef SWP_BINTREE_STATS_H_
#define SWP_BINTREE_STATS_H_

#include "swp_model_stats.h"
#include "swp_particle_model.h"

#include <vector>
#include <string>
#include <iostream>

namespace Sweep {

namespace Stats {
class BinTreeStats : public IModelStats {
public:
    //! Default constructor
    BinTreeStats();

    //! Copy constructor
    BinTreeStats(const BinTreeStats &copy) {*this = copy;}

    //! Stream-reading constructor
    BinTreeStats(std::istream &in, const Sweep::ParticleModel &model);

    //! Default destructor
    virtual ~BinTreeStats();

    //! Equate operator
    BinTreeStats &operator=(const BinTreeStats &rhs);


    // IMPLEMENTATION

    //! Returns the number of stats for this model.
    unsigned int Count() const {return STAT_COUNT;}

    //! Calculate the stats for the given ensemble
    void Calculate(const Ensemble &e, double scale);

    //! Returns a vector containing the stats.
    const fvector &Get() const {return m_stats;}

    //! Adds to a vector containing the stats from given point
    void Get(fvector &stats, unsigned int start = 0) const;

    //! Returns a vector containing the stat names.
    const std::vector<std::string> &Names(void) const {return m_names;}

    //! Adds to a vector containing stat from given point
    void Names(std::vector<std::string> &names, unsigned int start = 0) const;

    //! Returns the number of PSL output variables.
    unsigned int PSL_Count() const {return PSL_COUNT;}

    //! Adds to a vector of PSL variable names from given point
    void PSL_Names(std::vector<std::string> &names,
            unsigned int start = 0) const;

    //! Build the PSL entry for the given particle.
    void PSL(
        const Sweep::Particle &sp,      // Particle from which to get PSL data.
        double time,                      // Current flow time (used to calculate particle age).
        fvector &psl,                   // Output vector.
        unsigned int start = 0          // Optional start index for the first variable.
        ) const;

	//! Return primary particle details and connectivity
	void PrintPrimary(const Sweep::Particle &sp, std::vector<fvector> &nodes, std::vector<fvector> &primaries, int k) const;

    // READ/WRITE/COPY
    //! Creates a copy of the object.
    BinTreeStats *const Clone(void) const;

    //! Returns the model data type
    unsigned int ID() const {return (unsigned int)AggModels::BinTree_ID;}

    //! Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

private:
    //! Get geometric standard deviation of a list of diameters
    void GetGeometricStdev(
            const unsigned int num,
            std::vector<fvector> diams,
            fvector weights,
            fvector &gmeans,
            fvector &gstdevs
            ) const;

    //! Number of statistics entries.
    static const unsigned int STAT_COUNT = 11;

    //! Indices for stats headers.
    enum StatIndices {iNPrim=0, iPrimDiam=1, iSintLevel=2, iSintRate=3,
        iSintTime=4, iCollGStdev=5, iPrimGStdev=6, iGStdevMean=7, iCollGMean=8,
        iPrimGMean=9, iRg=10};

    //! The stats
    fvector m_stats;

    // The stat names.
    static const std::string m_statnames[STAT_COUNT];
    std::vector<std::string> m_names;

    //! The stat mask.
    static const IModelStats::StatType m_mask[STAT_COUNT];

    //! PSL count.
    static const unsigned int PSL_COUNT = 8;

    // The PSL names.
    static const std::string m_const_pslnames[PSL_COUNT];
    std::vector<std::string> m_pslnames;

};

} // Stats
} // Sweep

#endif /* SWP_BINTREE_STATS_H_ */
