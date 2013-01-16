/*!
 * @file    swp_bintree_silica_primary.h
 * @author  William J Menz
 * @brief   Silica particle model using generic binary tree structure
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      The silica model was originally published in Shekar et al (2012)
 *      J. Aerosol Sci. 44, 93-98. It was merged into the Git repository in
 *      August 2011 as its own class directly derived from Primary. Due to a
 *      large amount of code duplication and a number of errors in the
 *      original code, it was substantially modified in Oct 2012 to derive
 *      from a generic binary tree class.
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

#ifndef SWP_BINTREE_SILICA_PRIMARY_H_
#define SWP_BINTREE_SILICA_PRIMARY_H_

#include "swp_bintree_primary.h"

namespace Sweep {

namespace AggModels {

class BinTreeSilicaPrimary: public Sweep::AggModels::BinTreePrimary {
public:

    //! Build a new primary with one molecule
    BinTreeSilicaPrimary(
        const double time,
        const Sweep::ParticleModel &model
    );

    //! Copy constructor
    BinTreeSilicaPrimary(const BinTreeSilicaPrimary &copy);

    //! Stream-reading constructor
    BinTreeSilicaPrimary(
        std::istream &in,
        const Sweep::ParticleModel &model
    );

    //! Destructor
    ~BinTreeSilicaPrimary();

    //! Returns a copy of the primary.
    BinTreeSilicaPrimary *const Clone(void) const;

    //! Primary assignment operator
    BinTreeSilicaPrimary &operator=(const Primary &rhs);

    //! Returns the aggregation model which this primary describes
    AggModels::AggModelType AggID(void) const;

    //! Overload of parent class sinter
    void Sinter(double dt,
            Cell &sys,
            const Processes::SinteringModel &model,
            rng_type &rng,
            double wt);

    //! Return the number of hydroxide sites
    double GetSites() const { return GetComponent("hydroxide"); }

private:
    //! Default constructor (private)
    BinTreeSilicaPrimary();

    //! Sinter a node for time dt
    void SinterNode(double dt,
            Cell &sys,
            const Processes::SinteringModel &model,
            rng_type &rng,
            double wt);
};

}

}

#endif /* SWP_BINTREE_SILICA_PRIMARY_H_ */
