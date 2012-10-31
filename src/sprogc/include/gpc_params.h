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


    // double number STL vector.
    typedef std::vector<double> fvector;


    // MATHEMATICAL CONSTANTS.

    const double PI = 3.1415926535897932384626433832795;
    const double ONE_THIRD  = 3.3333334e-01;
    const double TWO_THIRDS = 6.6666667e-01;

    // PHYSICAL CONSTANTS.
	// These numbers are added by Vinod

	const double kB	= 1.3806504e-23; // Boltzmann constant m^2 kg s^-2 K^-1
	const double EPSILON0 = 8.854187816e-12;
	// conversion factors
	const double Debye__ = 3.33564e-30; //
	const double Angstroem__ = 1.0e-10;
    // Avogadro's number (source = NIST website, physics.nist.gov).
    // Error = 3.0e16 /mol.
    const double NA    = 6.02214179e23; // 1/mol.

    // Gas constant (source = NIST website, physics.nist.gov).
    // Error = 1.5e-5 J/mol/K.
    const double R     = 8.314472e0; // J/mol/K   (SI).
	//const double R = kB*NA;
    const double R_CGS = 8.314472e7; // ergs/mol/K (CGS).
    const double RCAL  = 1.9872e-3;  // kcal/mol/K (calories).

    //! Enum of the possible calculations of viscosity
    enum ViscosityModel {
        iAir, iChapmanEnskog, iArgon, iHydrogen
    };

} //namespace Sprog

#endif
