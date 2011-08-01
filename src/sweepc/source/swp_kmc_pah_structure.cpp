/*!
  * \Author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_pah_structure.cpp
  *
  * \brief        Implementation file for swp_kmc_pah_structure.h.
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

#include "swp_kmc_pah_structure.h"
#include "swp_kmc_structure_comp.h"
#include "swp_kmc_typedef.h"
#include "swp_kmc_jump_process.h"
#include "swp_kmc_pah_process.h"
#include "rng.h"
#include <iostream>
#include <fstream>
#include <list>
#include <cmath>
#include <map>

using namespace Sweep;
using namespace Sweep::KMC_ARS;
using namespace std;



// Constructors and Destructor
//! Default Constructor
PAHStructure::PAHStructure() {
	//NULLC = m_carbonList.insert(Carbon()).first;
    m_cfirst = NULLC;
    m_clast = NULLC;
    m_rings = 0;
    m_counts = intpair(2, 0);
    m_parent = NULL;
}

//! Copy Constructor (private member)
PAHStructure::PAHStructure(const PAHStructure& copy){
}
//! Default Destructor
PAHStructure::~PAHStructure() {
    delete m_cfirst;
    delete m_clast;
    PAHProcess pp(*this);
	pp.clearStructure();
    m_siteMap.clear();
    m_siteList.clear();
}

void PAHStructure::setParent(Sweep::AggModels::PAH* parent) {
    m_parent = parent;
}


//! Overloaded Operators

bool PAHStructure::operator==(PAHStructure &rhs) const
{   return this->m_cpositions==rhs.m_cpositions&&
           this->m_counts==rhs.m_counts;
}

bool PAHStructure::operator!=(PAHStructure &rhs) const{
    return !(*this==rhs);
}

int PAHStructure::numofC(){
    return m_counts.first;
} 
void PAHStructure::initialise(StartingStructure ss) {
    PAHProcess p(*this);
    p.initialise(ss);
}
PAHStructure*  PAHStructure::Clone() {
    PAHProcess p(*this);
    return p.clonePAH();
}

bool PAHStructure::havebridgeC(){
    PAHProcess p(*this);
    return p.havebridgeC();
}