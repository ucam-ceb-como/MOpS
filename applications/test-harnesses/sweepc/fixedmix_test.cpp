/*!
 * \file   fixedmixture_test.cpp
 * \author Robert I A Patterson
 *
 * \brief  Test harness for Sweep::FixedChemistry
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

#include "gpc_mech.h"
#include "gpc_mech_io.h"

#include "reset_chemistry.h"

#include "swp_fixed_mixture.h"

#include <iostream>
#include <limits>

int main() {
    // Set up a gas phase mechanism to use in the test
    Sprog::Mechanism mech;
    Sprog::IO::MechanismParser::ReadChemkin("chem.inp", mech, "therm.dat", 1, "tran.dat");

    // Read gas phase data from file
    Brush::ResetChemistry chemInputData("profile.dat", Brush::ResetChemistry::FixedChemistry, mech, 2);
    // Interpolate the data we need
    Sweep::fvector data(chemInputData.interpolateData(1.55e-4));

    // Now set up the FixedMixture object to test
    Sweep::FixedMixture fixedMixTestObject(data, mech);

    // tolerance for the numerical comparisons
    const double tol = 10.0 * std::numeric_limits<Sweep::real>::epsilon();

    // Now check the individual quantities
    {
        const double T = fixedMixTestObject.Temperature();
        const double TCorrect = 1e3;
        if(abs((T - TCorrect) / TCorrect) > tol) {
            std::cerr << "Got temperature of " << T << " when " << TCorrect << " expected\n";
            return 1;
        }
    }

    {
        const double m = fixedMixTestObject.MassDensity();
        const double mCorrect = 3.367623161e-01;
        if(abs((m - mCorrect) / mCorrect) > tol) {
            std::cerr << "Got mass density of " << m << " when " << mCorrect << " expected\n";
            return 1;
        }
    }

    {
        const double u = fixedMixTestObject.Velocity();
        const double uCorrect = 16.0;
        if(abs((u - uCorrect) / uCorrect) > tol) {
            std::cerr << "Got velocity of " << u << " when " << uCorrect << " expected\n";
            return 1;
        }
    }

    {
        const double n = fixedMixTestObject.MolarDensity();
        const double nCorrect = 7.274e-08;
        if(abs((n - nCorrect) / nCorrect) > tol) {
            std::cerr << "Got molar density of " << n << " when " << nCorrect << " expected\n";
            return 1;
        }
    }

    {
        const double P = fixedMixTestObject.Pressure();
        const double PCorrect = 1.0e+5;
        if(abs((P - PCorrect) / PCorrect) > tol) {
            std::cerr << "Got pressure of " << P << " when " << PCorrect << " expected\n";
            return 1;
        }
    }

    {
        const double mu = fixedMixTestObject.Viscosity();
        const double muCorrect = 9.99992786e-01;
        if(abs((mu - muCorrect) / muCorrect) > tol) {
            std::cerr << "Got viscosity of " << mu << " when " << muCorrect << " expected\n";
            return 1;
        }
    }

    // Now check the concentration of OH
    {
        const int iOH = mech.FindSpecies("OH");
        const double OH = fixedMixTestObject.SpeciesConcentration(iOH);
        const double OHCorrect = 4.1e-2;
        if(abs((OH - OHCorrect) / OHCorrect) > tol) {
            std::cerr << "Got [OH] (" << iOH << " as " << OH << " when " << OHCorrect << " expected\n";
            return 1;
        }
    }

    return 0;
}