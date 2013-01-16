/*!
 * \file   swp_environment_interface.h
 * \author Robert I A Patterson
 *  Copyright (C) 2012 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Interface for the fluid containing the particles
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

#ifndef SWEEP_ENVIRONMENT_INTERFACE_H
#define SWEEP_ENVIRONMENT_INTERFACE_H

#include "swp_params.h"

#include <string>

namespace Sweep {

/*!
 * \brief Specify interface for accessing properties of the fluid from within sweep
 *
 * This pure interface should provide a way to wrap all kinds of implementations of
 * fluid phases that may be used in calling code and provide the basic information
 * that sweep needs.  No data should be stored in this class.
 *
 * All named quantities are returned in SI units, concentrations in mol/m^3.
 * Support is also provided for additional properties that may depend on the
 * mixture and the problem under consideration, for example the ABF alpha
 * parameter.
 */
class EnvironmentInterface {
public:
    //! Virtual destructor, to call through to implementation
    virtual ~EnvironmentInterface() {}

    //! Virtual equivalent of copy constructor
    virtual EnvironmentInterface* Clone() const = 0;

    //! Indices to identify species
    typedef unsigned SpeciesIndex;

    //! Indices to identify non-species property
    typedef unsigned PropertyIndex;

    //! Pressure in
    virtual double Pressure() const = 0;

    //! Temperature in K
    virtual double Temperature() const = 0;

    //! Velocity in
    virtual double Velocity() const = 0;

    //! Viscosity in Pa.s
    virtual double Viscosity() const = 0;

    //! Density in \f$\mathrm{kg}\,\mathrm{m}^{-3}\f$
    virtual double MassDensity() const = 0;

    //! Total concentration in \f$\mathrm{mol}\,\mathrm{m}^{-3}\f$
    virtual double MolarDensity() const = 0;

    //! Concentration of species in \f$\mathrm{mol}\,\mathrm{m}^{-3}\f$
    virtual double SpeciesConcentration(const SpeciesIndex index) const = 0;

    //! Value of some mixture property that is not a species concentration
    virtual double PropertyValue(const PropertyIndex index) const = 0;

}; //class EnvironmentInterface

} //namespace Sweep

#endif
