/*!
  * \Author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_gaspoint.cpp
  *
  * \brief        Implementation file for swp_kmc_gaspoint.h.
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Holds the data structure which contains information on the molecular structure
    of PAH accounted by the kMC model.

  Licence:
    This file is part of "sweep".

    Sweep is free software; you can redistribute it and/or
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

#include "swp_kmc_gaspoint.h"

using namespace Sweep;
using namespace Sweep::KMC_ARS;
using namespace std;

// Constructors and Destructor
//! Default Constructor
KMCGasPoint::KMCGasPoint() 
:Time(0),
    T(1),
    H2(2),
    H(3),
    O2(4),
    OH(5),
    C2H2(6),
    C2H6(7),
    C6H6(8),
    H2O(9),
    CH4(10),
    CO(11),
    CO2(12),
    P(13),
    None(14),
    total(15)
{
    
    initData();
}
//! Copy Constructor
KMCGasPoint::KMCGasPoint(KMCGasPoint &gp) 
: Time(0),
    T(1),
    H2(2),
    H(3),
    O2(4),
    OH(5),
    C2H2(6),
    C2H6(7),
    C6H6(8),
    H2O(9),
    CH4(10),
    CO(11),
    CO2(12),
    P(13),
    None(14),
    total(15)
{
    m_data = gp.m_data;
}
// Default Destructor
KMCGasPoint::~KMCGasPoint() {
}
//! Initialise data point
void KMCGasPoint::initData() {
    m_data.clear();
    // set all to zero
    for(int i=0; i!=total-1; i++) {
        m_data.push_back(0);
    }
    // set 1 for no variables
    m_data.push_back(1);
}