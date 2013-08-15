/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_jump_process.cpp
  *
  * \brief        Implementation file for swp_kmc_jump_process.h
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    The JumpProcess class contains the rate constants for a jump process related to
    several Reaction classes as described by Celnik (2008) and Raj (2009) for the 
    calculation of reaction rates of growth processes with intermediates assumed as steady 
    states.

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

#include "swp_kmc_jump_process.h"
//#include "swp_kmc_mech.h"
#include <vector>

using namespace std;
using namespace Sweep;
using namespace Sweep::KMC_ARS;

//! Default Constructor
JumpProcess::JumpProcess():
	m_sType(),
	m_rxnvector0p0267(),
	m_rxnvector0p12(),
	m_rxnvector1(),
	m_r(),
	m_rate(0.0),
	m_name(""),
	m_ID(0)
{}

//! Copy Constructor
JumpProcess::JumpProcess(const JumpProcess& p):
	m_sType(),
	m_rxnvector0p0267(),
	m_rxnvector0p12(),
	m_rxnvector1(),
	m_r(),
	m_rate(p.m_rate),
	m_name(p.m_name),
	m_ID(p.m_ID) {
    m_sType = p.m_sType;
    m_rxnvector0p0267 = p.m_rxnvector0p0267;
    m_rxnvector0p12 = p.m_rxnvector0p12;
    m_rxnvector1 = p.m_rxnvector1;
}
//! Destructor
JumpProcess::~JumpProcess() {
    m_rxnvector0p0267.clear();
    m_rxnvector0p12.clear();
    m_rxnvector1.clear();
}

// Set Processes
//! Loads elementary reaction details (defined in derived classes), kmcSiteType and StructureProc related to it
void JumpProcess::initialise() {
    cout<<"....Base Class initialise called....\n\n";
}
//! Adds reaction
void JumpProcess::addReaction(std::vector<Sweep::KMC_ARS::Reaction>& rxnv, const Sweep::KMC_ARS::Reaction& rxn) {
    rxnv.push_back(rxn);
}
//! Calculate rates of each elementary reaction
void JumpProcess::calculateElemRxnRate(std::vector<Sweep::KMC_ARS::Reaction>& rxnv, const KMCGasPoint& gp/*, const double t_now*/) {

	m_r.resize(rxnv.size(), 0.0);
	for (size_t i(0); i != rxnv.size(); ++i) {
		m_r[i] = rxnv[i].getRate(gp);
	}
}
//! Calculates jump process rates and store (for Pressures 0.0267, 0.12 & 1 atm; defined in derived classes)
double JumpProcess::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/){
    cout<<"....Base Class setRate0p0267 called....\n\n";
    return 0.00;
}
double JumpProcess::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/){
    cout<<"....Base Class setRate0p12 called....\n\n";
    return 0.00;
}
double JumpProcess::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/){
    cout<<"....Base Class setRate1 called....\n\n";
    return 0.00;
}

// Read Process
//! Jump rate
double JumpProcess::getRate() const {
    return m_rate;
}
//! Gets site type associated with the jump process
kmcSiteType JumpProcess::getSiteType() const {
    return m_sType;
}
//! Returns name of process
std::string JumpProcess::getName() const {
    return m_name;
}
//! Returns process ID
int JumpProcess::getID() const {
    return m_ID;
}
//! Returns references to elementary reactions vectors
std::vector<Sweep::KMC_ARS::Reaction>& JumpProcess::getVec0p0267() {
    return m_rxnvector0p0267;
}
std::vector<Sweep::KMC_ARS::Reaction>& JumpProcess::getVec0p12() {
    return m_rxnvector0p12;
}
std::vector<Sweep::KMC_ARS::Reaction>& JumpProcess::getVec1() {
    return m_rxnvector1;
}
