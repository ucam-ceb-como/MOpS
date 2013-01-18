/*!
 * \file   swp_sprog_idealgas.cpp
 * \author Robert I A Patterson
 *  Copyright (C) 2012 Robert I A Patterson
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief Implement the Sweep::SprogIdealGasWrapper around a Sprog::IdealGas
 * 

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

#include "swp_sprog_idealgas_wrapper.h"

#include <stdexcept>
#include <sstream>


/*
 * These constants should be used to initialise processes and other classes that need to
 * refer to the respective properties.  They should not be used directly as an argument
 * in calls to PropertyValue, because this would defeat the abstraction that this
 * class is meant to provide.
 *
 * Developers are free to change the numerical values at any time.
 */
const Sweep::EnvironmentInterface::PropertyIndex Sweep::SprogIdealGasWrapper::sTemperatureGradientIndex = 1000;
const Sweep::EnvironmentInterface::PropertyIndex Sweep::SprogIdealGasWrapper::sMixFracGradientIndex = 1001;
const Sweep::EnvironmentInterface::PropertyIndex Sweep::SprogIdealGasWrapper::sMixFracLaplacianIndex = 1002;
const Sweep::EnvironmentInterface::PropertyIndex Sweep::SprogIdealGasWrapper::sMixFracDiffusionIndex = 1003;
const Sweep::EnvironmentInterface::PropertyIndex Sweep::SprogIdealGasWrapper::sAvgMolWtIndex = 1004;
const Sweep::EnvironmentInterface::PropertyIndex Sweep::SprogIdealGasWrapper::sThermalConductivityIndex = 1005;
const Sweep::EnvironmentInterface::PropertyIndex Sweep::SprogIdealGasWrapper::sAlphaIndex = 1006;
const Sweep::EnvironmentInterface::PropertyIndex Sweep::SprogIdealGasWrapper::sPAHFormationIndex = 1007;

/*!
 * @param[in]   model       Particle model containing gas phase mechanism
 */
Sweep::SprogIdealGasWrapper::SprogIdealGasWrapper(const Sprog::SpeciesPtrVector& species)
: mGas(species)
{}

/*!
 *@return   Pointer to a deep copy
 */
Sweep::EnvironmentInterface* Sweep::SprogIdealGasWrapper::Clone() const {
    return new SprogIdealGasWrapper(*this);
}

/*!
 * @return    Pressure in \f$\mathrm{Pa}\f$
 */
double Sweep::SprogIdealGasWrapper::Pressure() const {
    return mGas.Pressure();
}

/*!
 * @return      Absolute temperature of mixture
 */
double Sweep::SprogIdealGasWrapper::Temperature() const {
    return mGas.Temperature();
}

/*!
 *
 */
double Sweep::SprogIdealGasWrapper::Velocity() const {
    return mGas.Velocity();
}

/*!
 *@return Viscosity in \f$\mathrm{Pa.s}\f$
 */
double Sweep::SprogIdealGasWrapper::Viscosity() const {
    return mGas.getViscosity();
}

/*!
 *@return Mass Density in \f$\mathrm{kg}\,\mathrm{m}^{-3}\f$
 */
double Sweep::SprogIdealGasWrapper::MassDensity() const {
    return mGas.MassDensity();
}

/*!
 * Historically in Sprog the quantity returned by this function
 * was known as "Density".
 *
 *@return total concentration in \f$\mathrm{kg}\,\mathrm{m}^{-3}\f$
 */
double Sweep::SprogIdealGasWrapper::MolarDensity() const {
    // See note about historic naming in sprog
    return mGas.Density();
}

/*!
 *@param[in]    index   Identifier for species
 *
 *@return       Concentration of species in \f$\mathrm{mol}\,\mathrm{m}^{-3}\f$
 */
double Sweep::SprogIdealGasWrapper::SpeciesConcentration(const SpeciesIndex index) const {
    return mGas.MolarConc(index);
}


/*!
 *@param[in]    index   Identifier for property
 *
 *@return      Value of specified property
 *
 *@exception    std::runtime_error      Unrecognised property index
 */
double Sweep::SprogIdealGasWrapper::PropertyValue(const PropertyIndex index) const {
    switch(index) {
    case sTemperatureGradientIndex:
        return mGas.GradientTemperature();
    case sMixFracGradientIndex:
        return mGas.MixFracDiffCoeff();
    case sMixFracLaplacianIndex:
        return mGas.LaplacianMixFrac();
    case sAvgMolWtIndex:
        return mGas.getAvgMolWt();
    case sThermalConductivityIndex:
        return mGas.getThermalConductivity(mGas.Pressure());
    case sAlphaIndex:
        return mGas.Alpha();
    case sPAHFormationIndex:
        return mGas.PAHFormationRate();
    default:
        std::ostringstream msg;
        msg << "Unrecognised property index " << index << " in Sweep::SprogIdealGasWrapper::PropertyValue";
        throw std::runtime_error(msg.str());
    }
    // Avoid compiler warnings by providing a return value
    return std::numeric_limits<double>::quiet_NaN();
}

/*!
 * @param[in,out]   out     Output stream
 */
void Sweep::SprogIdealGasWrapper::Serialize(std::ostream &out) const {
    mGas.Serialize(out);
}

/*!
 * @param[in,out]   in      Input stream
 * @param[in]       species Gas phase species that help define the meaning of the data
 */
void Sweep::SprogIdealGasWrapper::Deserialize(std::istream &in, const Sprog::SpeciesPtrVector &species) {
    mGas.Deserialize(in);
    mGas.SetSpecies(species);
}
