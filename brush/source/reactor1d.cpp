/*!
 * \file   reactor1d.cpp
 * \author Robert I A Patterson
 *
 * \brief 1d system of touching mops reactors
 *
 *  Copyright (C) 2009 Robert I A Patterson.

 Licence:
    This file is part of "brush".

    brush is free software; you can redistribute it and/or
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
    Prof Markus Kraft
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

#include "reactor1d.h"

#include "reset_chemistry.h"

#include <cassert>

/*!
 * Create an object with a collection of empty reactors, all sharing
 * the same mechanism instance.  A copy of the mechanism is stored in
 * the new Reactor1d instance and it is this copy that is shared
 * between the reactors in the cells.
 *
 *\param[in]    geom                    Geometry of new 1d reactor
 *\param[in]    mech                    Mechanism to use in all the reactors
 *\param[in]    max_particle_counts     Maximum number of computatioal particles per cell
 *\param[in]    max_m0s                 Maximum particle concentration
 */
Brush::Reactor1d::Reactor1d(const Geometry::Geometry1d &geom, const Mops::Mechanism &mech,
                            const Utils::LinearInterpolator<real, real> &max_particle_counts,
                            const Utils::LinearInterpolator<real, real> &max_m0s)
    : mMech(mech)
    , mGeometry(geom)
    , mReactors(geom.numCells(), Mops::Reactor(mMech))
{
    // Initialise the particle ensembles in each cell
    for(unsigned int i = 0; i < mReactors.size(); ++i) {
        // Create a reference to the reactor for syntactic convenience
        Mops::Reactor &r = mReactors[i];
        r.Fill(*(new Mops::Mixture(mMech.ParticleMech())));

        // Interpolate particle count and max m0 to the cell centre
        real posn = getCellCentre(i);
        real maxM0 = max_m0s.interpolate(posn);
        const unsigned int particleCount = static_cast<unsigned int>(max_particle_counts.interpolate(posn) + 0.5);

        r.Mixture()->Particles().Initialise(particleCount, mMech.ParticleMech());
        r.Mixture()->Reset(maxM0);

        assert(r.Mixture() != NULL);
    }
}


/*!
 * Mechanism will be copied, so mechanism pointers must be updated
 */
Brush::Reactor1d::Reactor1d(const Reactor1d &src)
    : mMech(src.mMech)
    , mGeometry(src.mGeometry)
    , mReactors(src.mReactors) {

    for(unsigned int i = 0; i < mReactors.size(); ++i) {
        mReactors[i].SetMech(mMech);
    }


}

/*!
 * Use the apply method on the supplied object to overwrite the chemical
 * contents of each cell, while leaving the particles untouched.
 *
 * The option to ignore the use of gaseous species in particle events is provided
 * for use in cases where gas phase profiles have been precalculated by an
 * independent method that included an approximation for soot formation.
 *
 *\param[in]    chem_source     Class containing data to set
 *\param[in]    fixed_chem      Ignore use of gaseous species in particle events
 */
void Brush::Reactor1d::ReplaceChemistry(const ResetChemistry& chem_source,
                                        const bool fixed_chem) {
    for(size_t i = 0; i < getNumCells(); ++i) {
        chem_source.apply(getCellCentre(i), mReactors[i]);
        mReactors[i].Mixture()->SetFixedChem(fixed_chem);
    }

}

/*!
 * Set the time on each of the sub-reactors contained
 * within this instance.
 *
 *\param[in]        t       Time to set
 */
void Brush::Reactor1d::setTime(const real t) {
    for(unsigned int i = 0; i < mReactors.size(); ++i) {
        mReactors[i].SetTime(t);
    }

}
