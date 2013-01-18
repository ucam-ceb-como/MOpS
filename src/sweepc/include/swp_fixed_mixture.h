/*!
 * \file   swp_fixed_mixture.h
 * \author Robert I A Patterson
 *
 * \brief Very simple chemical mixture
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
#ifndef SWEEP_FIXED_MIXTURE_H
#define	SWEEP_FIXED_MIXTURE_H

#include "gpc_species.h"

#include "swp_params.h"
#include "swp_environment_interface.h"

#include <vector>

// Forward declaration
namespace Sprog {
    class Mechanism;
}

namespace Sweep {

//! Hold a fixed mixture without renormalisation or equations of state..
/*!
 * This is a very simple implementation of the Sweep::EnvironmentInterface.  It has
 * no equation of state and no guarantee that the data satisifies the ideal gas law.
 * The read-only data is initialised in the constructor.
 */
class FixedMixture : public EnvironmentInterface{
public:
    //! Read in the data needed for the given mechanism
    FixedMixture(const fvector& data, const Sprog::SpeciesPtrVector &species);

    //==== Methods implementing the EnvironmentInterface interface ====//
    //! Virtual equivalent of copy constructor
    virtual EnvironmentInterface* Clone() const;

    //! Pressure in Pa
    virtual double Pressure() const {return mPressure;}

    //! Temperature in K
    virtual double Temperature() const {return mTemperature;}

    //! Velocity in
    virtual double Velocity() const {return mVelocity;}

    //! Viscosity in Pa.s
    virtual double Viscosity() const {return mViscosity;}

    //! Density in \f$\mathrm{kg}\,\mathrm{m}^{-3}\f$
    virtual double MassDensity() const {return mMassDensity;}

    //! Total concentration in \f$\mathrm{mol}\,\mathrm{m}^{-3}\f$
    virtual double MolarDensity() const {return mMolarDensity;}

    //! Concentration of species in \f$\mathrm{mol}\,\mathrm{m}^{-3}\f$
    virtual double SpeciesConcentration(const SpeciesIndex index) const;

    //! Value of some mixture property that is not a species concentration

    virtual double PropertyValue(const PropertyIndex index) const;
    //===== End of EnvironmentInterface implementation ================//

private:
    double mMassDensity; //! Mixture density (kg/m3)
    double mMolarDensity; //! Concentration of molecules (mol/m3)
    double mPressure; //! Pressure (Pa)
    double mTemperature; //! Temperature (K)
    double mVelocity; //! Velocity (m/s)
    double mViscosity; //! Viscosity (kg/m/s = Pa s)

    fvector mSpecies; //! Species concentrations (mol/m3)
};

}

#endif	/* SWEEP_FIXED_MIXTURE_H */

