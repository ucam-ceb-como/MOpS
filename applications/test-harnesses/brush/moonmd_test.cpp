/*!
 *\file moonmd_test.cpp
 *\author Robert I A Patterson
 *
 *\brief Test harness for the moonmd interface
 *
 *  Copyright (C) 2011 Robert I A Patterson.
 *

 Licence:
    This file is part of "brush".

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
 *
 *
 *\main Brush is a program for the coupled simulation of laminar reacting flow
 * and particle population development in 1d systems.  Moonmd is ...
 */

#include <cassert>
#include <string>
#include <vector>
#include <iostream>

#include "moonmd_interface.h"

#include "reactor1d.h"
#include "mt19937.h"

//! Run brush
int main(int argc, char* argv[])
{
    // reset the random number generator
    Sweep::init_genrand(789);

    std::cout << "Creating reactor\n";
    const std::string chemfile("chem.inp");
    const std::string thermfile("therm.dat");
    const std::string settfile("brush.xml");
    const std::string swpfile("sweep.xml");
    const std::string partsolnfile("partsoln.xml");
    const std::string chemsolnfile("chemsoln.dat");
    const double gridNodes[41] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                                 0.10 ,0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
                                 0.20 ,0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
                                 0.30 ,0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
                                 0.4};

    // Check we can create a reactor
    Brush::MooNMDInterface::particle_reactor_pointer pReac = NULL;
    pReac = Brush::MooNMDInterface::InitialiseBrush(chemfile, thermfile, settfile, swpfile,
                                                    partsolnfile, chemsolnfile, 5, gridNodes);
    assert(pReac != NULL);
    std::cout << "Reactor created successfully\n";

    //==================== File to store the moments of the particle distribution
    std::ofstream momentsFile("moonmd1-moments.csv");

    // Write the column headings to the file, with commas between names, but not after the last name
    momentsFile << "t,x";
    typedef std::vector<std::string> string_vector;
    {
        const string_vector &colNames = Sweep::Stats::EnsembleStats(pReac->getMechanism().ParticleMech()).Names();
        const string_vector::const_iterator itEnd = colNames.end();
        for(string_vector::const_iterator it = colNames.begin(); it != itEnd; ++it) {
            momentsFile << ',' << *it;
        }
    }
    momentsFile << '\n';

    // Now see if the reactor is useable
    std::cout << "Testing reactor\n";
    const double solnGrid[3] = {0.05, 0.3, 0.5};
    const double temperature[3] = {300.0, 300.0, 300.0};
    const double massConcentration[6] = {999.0, 999.0, 999.0, 1.0, 1.0, 1.0};
    const double velocity[3] = {2.0, 2.0, 2.0};
    double energySource[3] = {-9.9, -8.9, -7.9};
    double massConcSource[6] = {-19.9, -18.9, -17.9, -29.9, -28.9, -27.9};
    pReac = Brush::MooNMDInterface::RunParticlePhase(*pReac, 2.0, 3, 2, solnGrid, temperature, velocity,
                                                     massConcentration, energySource, massConcSource, momentsFile);


    delete pReac;
}
