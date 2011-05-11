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

//! Run brush
int main(int argc, char* argv[])
{
    const std::string chemfile("chem.inp");
    const std::string thermfile("therm.dat");
    const std::string settfile("brush.xml");
    const std::string swpfile("sweep.xml");
    const std::string partsolnfile("partsoln.xml");
    std::vector<double> gridNodes(5);
    gridNodes[0] = 0.0;
    gridNodes[1] = 0.1;
    gridNodes[2] = 0.2;
    gridNodes[3] = 0.3;
    gridNodes[4] = 0.4;

    Brush::MooNMDInterface::particle_reactor_pointer pReac = NULL;
    pReac = Brush::MooNMDInterface::InitialiseBrush(chemfile, thermfile, settfile, swpfile, partsolnfile, gridNodes);
    assert(pReac != NULL);

    // Return 0 (ie FALSE) if the reactor has the correct number of cells
    return (pReac->getNumCells() + 1 != gridNodes.size());
}
