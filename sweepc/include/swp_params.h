/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Defines parameters and typedefs used by sweep.

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

#ifndef SWEEP_PARAMS_H
#define SWEEP_PARAMS_H

#include "sprog.h"
#include <cmath>
#include <vector>

namespace Sweep
{
    typedef Sprog::real real;
    typedef Sprog::fvector fvector;

    const real PI         = Sprog::PI;
    const real ONE_THIRD  = Sprog::ONE_THIRD;
    const real TWO_THIRDS = Sprog::TWO_THIRDS;

    // Avogadro's constant.
    const real NA = Sprog::NA;

    // Gas constant.
    const real R     = Sprog::R;     // J/molK.
    const real R_CGS = Sprog::R_CGS; // erg/molK.
    const real RCAL  = Sprog::RCAL;  // kcal/molK.

    // Boltzmann constant (source = NIST website, physics.nist.gov).
    // Error = 2.4e-29 J/K.
    const real KB     = 1.3806504e-23; // J/K.
    const real KB_CGS = 1.3806504e-16; // erg/K.

    // Constant term in calculation on Knudsen number.  Actual
    // Kn = KNUDSEN_K * T [Kelvin] / (P [bar] * d [cm]).
//    const real KNUDSEN_K = 4.74151636e-8;  // bar.cm/K. TODO:  Check units here for Knudsen number!
//    const real KNUDSEN_K = 4.74289918e-5; // Pa.m/K (for gas d = 0.362 nm).

    // Free-molecular coagulation kernel parameters.
    const real CFM     = 4.65695224e-12; // Sqrt(KB * PI / 2) = sqrt(J/K).
    const real CFM_CGS = 1.47265760e-08; // = sqrt(erg/K).
    const real CFMMAJ  = 1.4178;		 //ATTENTION deleted by ms785    this value can only be used if the particles are spherical
  //	const real CFMMAJ  = 2;
	
    // Slip-flow coagulation kernel parameters.
    const real CSF     = 9.2046667e-24; // = KB * 2/3 (J/K).
    const real CSF_CGS = 9.2046667e-17; // = KB * 2/3 (erg/K).

    // Returns the viscosity of air at 1 atm 
    // for the given temperature.
    inline real ViscosityAir(real T)
    {
        // Ref: Kazakov & Frenklach (1998), Combust. Flame 114, 484-501.
        return 1.458e-6 * T * sqrt(T) / (T + 110.4); // kg/m.s.
    };

    // Returns the mean free path in air assuming ideal gas
    // law and a collision cross-section of molecules of 0.362nm.
    // Returns units are metres (m).
    inline real MeanFreePathAir (
        real T, // Temperature (K).
        real P) // Pressure (Pa).
    {
        //static const real dgas = 3.62e-10; // = 0.362 nm.
        //return (T/P) * R /(sqrt(2.0) * PI * NA * dgas * dgas);
        static const real k = 2.37075818e-5; // Pa.m/K (from sweep1 & sweep2).
//        static const real k = 2.37144959e-5; // Pa.m/K (for gas d = 0.362 nm).
        return k * T / P;
    };

    // Returns the mean free path of an ideal gas given
    // the temperature, pressure and collision cross section
    // of the gas molecules.
    inline real MeanFreePath(
        real T, // Temperature (K).
        real P, // Pressure (Pa).
        real d) // Collision cross-section (m).
    {
        return R * T / (P * sqrt(2.0) * PI * NA * d * d);
    };

    // Returns the Knudsen number for a particle of the given diameter
    // in air at the given temperature and pressure.
    inline real KnudsenAir(
        real T, // Temperature (K).
        real P, // Pressure (Pa).
        real d) // Particle collision diameter (m).
    {
        return 2.0 * MeanFreePathAir(T, P) / d;
    }
};

#endif
