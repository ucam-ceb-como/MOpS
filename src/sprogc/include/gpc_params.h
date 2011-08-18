/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Defines parameters and typedefs used by sprog.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
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

#ifndef GPC_PARAMS_H
#define GPC_PARAMS_H

#include <vector>

namespace Sprog
{
    // COMMON TYPEDEFS.

    // Real numbers in Sprog, so they can be easily changed.
    typedef double real;

    // Real number STL vector.
    typedef std::vector<real> fvector;


    // MATHEMATICAL CONSTANTS.

    const real PI = 3.1415926535897932384626433832795;
    const real ONE_THIRD  = 3.3333334e-01;
    const real TWO_THIRDS = 6.6666667e-01;

    // PHYSICAL CONSTANTS.
	// These numbers are added by Vinod

	const real kB	= 1.3806504e-23; // Boltzmann constant m^2 kg s^-2 K^-1
	const real EPSILON0 = 8.854187816e-12;
	// conversion factors
	const real Debye__ = 3.33564e-30; //
	const real Angstroem__ = 1.0e-10;
    // Avogadro's number (source = NIST website, physics.nist.gov).
    // Error = 3.0e16 /mol.
    const real NA    = 6.02214179e23; // 1/mol.

    // Gas constant (source = NIST website, physics.nist.gov).
    // Error = 1.5e-5 J/mol/K.
    const real R     = 8.314472e0; // J/mol/K   (SI).
	//const real R = kB*NA;
    const real R_CGS = 8.314472e7; // ergs/mol/K (CGS).
    const real RCAL  = 1.9872e-3;  // kcal/mol/K (calories).


} //namespace Sprog

#endif
