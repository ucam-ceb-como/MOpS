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

#include "gpc_params.h"
#include <cmath>
#include <vector>
#include <boost/random/mersenne_twister.hpp>

namespace Sweep
{
    
    typedef Sprog::fvector fvector;

    //! Type of random number generator to use throughout sweep
    typedef boost::mt19937 rng_type;

    const double PI         = Sprog::PI;
    const double ONE_THIRD  = Sprog::ONE_THIRD;
    const double TWO_THIRDS = Sprog::TWO_THIRDS;
    const double TWO_ONE_THIRD = 0.793700526;      // 2^(-1/3)

    // Avogadro's constant.
    const double NA = Sprog::NA;

    // Gas constant.
    const double R     = Sprog::R;     // J/molK.
    const double R_CGS = Sprog::R_CGS; // erg/molK.
    const double RCAL  = Sprog::RCAL;  // kcal/molK.

    // Boltzmann constant (source = NIST website, physics.nist.gov).
    // Error = 2.4e-29 J/K.
    const double KB     = 1.3806504e-23; // J/K.
    const double KB_CGS = 1.3806504e-16; // erg/K.

    // Constant term in calculation on Knudsen number.  Actual
    // Kn = KNUDSEN_K * T [Kelvin] / (P [bar] * d [cm]).
//    const double KNUDSEN_K = 4.74151636e-8;  // bar.cm/K. TODO:  Check units here for Knudsen number!
//    const double KNUDSEN_K = 4.74289918e-5; // Pa.m/K (for gas d = 0.362 nm).

    // Free-molecular coagulation kernel parameters.
    const double CFM     = 4.65695224e-12; // Sqrt(KB * PI / 2) = sqrt(J/K).
    const double CFM_CGS = 1.47265760e-08; // = sqrt(erg/K).
	const double CFMMAJ  = 2;  //ms785    1.41 can only be used if the particles are spherical
   // const double CFMMAJ  = 1.4178;	

    // Slip-flow coagulation kernel parameters.
    const double CSF     = 9.2046667e-24; // = KB * 2/3 (J/K).
    const double CSF_CGS = 9.2046667e-17; // = KB * 2/3 (erg/K).

    // Returns the viscosity of air at 1 atm 
    // for the given temperature.
    inline double ViscosityAir(double T)
    {
        // Ref: Kazakov & Frenklach (1998), Combust. Flame 114, 484-501.
        return 1.458e-6 * T * sqrt(T) / (T + 110.4); // kg/m.s.
    };

    // Returns the mean free path in air assuming ideal gas
    // law and a collision cross-section of molecules of 0.362nm.
    // Returns units are metres (m).
    inline double MeanFreePathAir (
        double T, // Temperature (K).
        double P) // Pressure (Pa).
    {
        //static const double dgas = 3.62e-10; // = 0.362 nm.
        //return (T/P) * R /(sqrt(2.0) * PI * NA * dgas * dgas);
        static const double k = 2.37075818e-5; // Pa.m/K (from sweep1 & sweep2).
//        static const double k = 2.37144959e-5; // Pa.m/K (for gas d = 0.362 nm).
        return k * T / P;
    };

    // Returns the mean free path of an ideal gas given
    // the temperature, pressure and collision cross section
    // of the gas molecules.
    inline double MeanFreePath(
        double T, // Temperature (K).
        double P, // Pressure (Pa).
        double d) // Collision cross-section (m).
    {
        return R * T / (P * sqrt(2.0) * PI * NA * d * d);
    };

    // Returns the Knudsen number for a particle of the given diameter
    // in air at the given temperature and pressure.
    inline double KnudsenAir(
        double T, // Temperature (K).
        double P, // Pressure (Pa).
        double d) // Particle collision diameter (m).
    {
        return 2.0 * MeanFreePathAir(T, P) / d;
    }
};

#endif
