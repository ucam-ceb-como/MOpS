/*!
 * \file   swp_fixed_mixture.cpp
 * \author Robert I A Patterson
 *
 * \brief  Very simple chemical mixture
 *
 *  Copyright (C) 2012 Robert I A Patterson.
 *

 Licence:
    This file is part of "sweep".

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
#include "swp_fixed_mixture.h"

#include "gpc_mech.h"

#include <sstream>
#include <stdexcept>


/*!
 *\param[in]    data    Values of species concentrations and other quantities
 *\param[in]    species Details of the species to which the species concentrations refer
 *
 * The order of the entries in data is important and must be as follows (all SI units)
 * - 0 Skipped (intended for position)
 * - 1 Temperature
 * - 2 Mass density
 * - 3 Velocity
 * - 4 Overall molar concentration/density (mol/m3)
 * - 5 Pressure
 * - 6 Viscosity
 * - 7+ Species concentrations (mol/m3) in the same order as the species in mech
 *
 *\exception    std::runtime_error  Input data does not match number of species plus other items
 */
Sweep::FixedMixture::FixedMixture(const fvector& data,
                                  const Sprog::SpeciesPtrVector &species) {
    if(data.size() != (7 + species.size())) {
        std::ostringstream msg;
        msg << "Found " << data.size() << " data items, but mechanism has "
            << species.size() << " and 7 further items are also required "
            << "(Sweep::FixedMixture::FixedMixture)";
        throw std::runtime_error(msg.str());
    }

    // next item of data
    fvector::const_iterator dataIt = data.begin();

    // first item should be position, which is not stored
    ++dataIt;

    // now copy the data we need, for the order see the FixedChemistry case in
    // Brush::ResetChemistry::ResetChemistry.
    mTemperature = *dataIt++;
    mMassDensity = *dataIt++;
    mVelocity = *dataIt++;
    mMolarDensity = *dataIt++;
    mPressure = *dataIt++;
    mViscosity = *dataIt++;

    mSpecies.assign(dataIt, data.end());
}

/*!
 *@return   Pointer to a deep copy
 */
Sweep::EnvironmentInterface* Sweep::FixedMixture::Clone() const {
    return new FixedMixture(*this);
}

/*!
 * @param[in]    index    Index of species in chemical mechanism
 *
 * @return        Concentration of species
 */
double Sweep::FixedMixture::SpeciesConcentration(const SpeciesIndex index) const {
    return mSpecies[index];
}

/*!
 * @param[in]    index    Index of property
 *
 * @return       Value associated with specified property
 */
double Sweep::FixedMixture::PropertyValue(const PropertyIndex index) const {
    throw std::runtime_error("PropertyValue not yet supported in Sweep::FixedMixture::PropertyValue");
    return 0.0;
}


