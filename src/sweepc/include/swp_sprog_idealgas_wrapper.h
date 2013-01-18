/*!
 * \file   swp_sprog_idealgas_wrapper.h
 * \author Robert I A Patterson
 *  Copyright (C) 2012 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Interface to a Mops class describing an ideas gas
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

#ifndef SWEEP_SPROG_IDEALGAS_WRAPPER_H
#define SWEEP_SPROG_IDEALGAS_WRAPPER_H

#include "swp_environment_interface.h"

#include "gpc_idealgas.h"

#include "swp_params.h"

#include <iostream>

namespace Sweep {

/*!
 * \brief Wrap Sprog::IdealGas into the form expected by Sweep
 */
class SprogIdealGasWrapper : public EnvironmentInterface {
public:
    //! Construct from a mechanism
    SprogIdealGasWrapper(const Sprog::SpeciesPtrVector& species);

    //! Virtual equivalent of copy constructor
    EnvironmentInterface* Clone() const;

    //! Pressure in
    double Pressure() const;

    //! Temperature in K
    double Temperature() const;

    //! Velocity in
    double Velocity() const;

    //! Viscosity in Pa.s
    double Viscosity() const;

    //! Density in \f$\mathrm{kg}\,\mathrm{m}^{-3}\f$
    double MassDensity() const;

    //! Total concentration in \f$\mathrm{mol}\,\mathrm{m}^{-3}\f$
    double MolarDensity() const;

    //! Concentration of species in \f$\mathrm{mol}\,\mathrm{m}^{-3}\f$
    double SpeciesConcentration(const SpeciesIndex index) const;

    //! Value of some mixture property that is not a species concentration
    double PropertyValue(const PropertyIndex index) const;

    //! Access an older interface provided by sprog
    const Sprog::Thermo::IdealGas * Implementation() const {return &mGas;}

    //! Access an older style interface provided by sprog
    Sprog::Thermo::IdealGas * Implementation() {return &mGas;}

    //! Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                          // Input stream.
        const Sprog::SpeciesPtrVector &species     // Details of gas phase species
        );

    //! Index for the temperature gradient in \f$\mathrm{K}\,\mathrm{m}^{-1}\f$
    static const PropertyIndex sTemperatureGradientIndex;

    //! Index for the gradient of the flamelet mixture fraction
    static const PropertyIndex sMixFracGradientIndex;

    //! Index for the Laplacian of the flamelet mixture fraction
    static const PropertyIndex sMixFracLaplacianIndex;

    //! Index for the diffusion coefficient of the flamelet mixture fraction
    static const PropertyIndex sMixFracDiffusionIndex;

    //! Index for the average molecular weight in \f$\mathrm{kg}\,\mathrm{mol}^{-1}\f$
    static const PropertyIndex sAvgMolWtIndex;

    //! Index for the thermal conductivity of the mixture
    static const PropertyIndex sThermalConductivityIndex;

    //! Index for the active sites factor
    static const PropertyIndex sAlphaIndex;

    //! Index for the PAH formation rate
    static const PropertyIndex sPAHFormationIndex;

private:
    //! Cannot construct without knowledge of species
    SprogIdealGasWrapper();

    //! This is the data containing object that is being wrapped
    Sprog::Thermo::IdealGas mGas;

}; //class SprogIdealGasWrapper

} //namespace Sweep
#endif
