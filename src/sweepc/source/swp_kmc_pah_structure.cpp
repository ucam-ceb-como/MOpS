/*!
  * \author     Zakwan Zainuddin (zz260)
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
#include "string_functions.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <cmath>
#include <map>

using namespace Sweep;
using namespace Sweep::KMC_ARS;
using namespace std;
using namespace Strings;


// Constructors and Destructor
//! Default Constructor
PAHStructure::PAHStructure() {
    //NULLC = m_carbonList.insert(Carbon()).first;
    m_cfirst = NULLC;
    m_clast = NULLC;
    m_rings = 0;
    m_rings5 = 0;
    m_counts = intpair(0, 0);
    m_parent = NULL;
}

//! Copy Constructor (private member)
PAHStructure::PAHStructure(const PAHStructure& copy){
   *this = copy;
}

PAHStructure &PAHStructure::operator=(const PAHStructure &rhs)
{
    if (this != &rhs) {
        m_siteList = rhs.m_siteList;
        // the below variables are not interested now, and they has been commented out to save computation effort. If user wanna to actitave them, please be careful, particularly m_carbonList, since each element of m_carbonList will be destroied in the deconstructor. If you do wanna to copy the m_carbonList, make sure you locate new space for the objects instead of "simply copy". dc516
        //m_cfirst   = rhs.m_cfirst;
        //m_clast    = rhs.m_clast;
        //m_siteMap  = rhs.m_siteMap;
        //m_rings5   = rhs.m_rings5;
        //m_rings    = rhs.m_rings;
        //m_counts   = rhs.m_counts;
        //m_carbonList = rhs.m_carbonList;
        //m_cpositions = rhs.m_cpositions;
    }
    return *this;
}
//! Default Destructor
PAHStructure::~PAHStructure() {
    //delete m_cfirst;
    //delete m_clast;
    PAHProcess pp(*this);
    clear();
    m_siteMap.clear();
    m_siteList.clear();
}

//! Remove all data and release any memory
void PAHStructure::clear() {
    for(set<Cpointer>::iterator i=m_carbonList.begin(); i!=m_carbonList.end(); i++)
        delete *i;
    // clear all data
    m_carbonList.clear();
    m_siteMap.clear();
    m_siteList.clear();
    m_cfirst = NULL;
    m_clast = NULL;
    m_counts.first = 0;
    m_counts.second = 0;
    m_cpositions.clear();
    m_rings = 0;
    m_rings5 = 0;
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

int PAHStructure::numofC() const{
    return m_counts.first;
} 

int PAHStructure::numofH() const{
    return m_counts.second;
} 

int PAHStructure::numofRings() const{
    return m_rings;
}

int PAHStructure::numofRings5() const{
    return m_rings5;
}

int PAHStructure::numofEdgeC() const{
    // the m_cpositions stores the coordinates of PAH, which means the num of edge C equals the size of m_cpositions
    return m_cpositions.size();
}

int PAHStructure::numofSite() const
{
    return m_siteList.size();
}
void PAHStructure::setnumofC(int val)
{
    m_counts.first=val;
}

void PAHStructure::setnumofH(int val)
{
    m_counts.second=val;
}

void PAHStructure::setnumofRings(int val)
{
    m_rings=val;
}

void PAHStructure::setnumofRings5(int val)
{
    m_rings5=val;
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

void PAHStructure::saveDOTperLoop(int PAH_ID, int i)
{
    PAHProcess p(*this);
    string filename = "KMC_DEBUG/ID_";
    filename.append(Strings::cstr(PAH_ID));
    filename.append("_");
    filename.append(Strings::cstr(i));
    filename.append(".dot");
    p.saveDOT(filename);
}


// currently the serialization is incomplete, and some info is lost during this process
// for instance, m_carbonList, m_siteList, m_siteMap, m_cfirst, m_clast
void PAHStructure::Serialize(std::ostream &out) const
{
    int val=0;

    // output info for PAHProcess::createPAH().
    val=numofRings();
    out.write((char*)&(val), sizeof(val));

    val=numofRings5();
    out.write((char*)&(val), sizeof(val));

    PAHStructure m_copy (*this);
    PAHProcess p(m_copy);
    std::string m_SiteName = p.SiteString(',');

    val = (unsigned int)m_SiteName.length();
    out.write((char*)&val, sizeof(val));
    out.write(m_SiteName.c_str(), val);
}

void PAHStructure::Deserialize(std::istream &in)
{
    int val = 0;
    char *name = NULL;

    in.read(reinterpret_cast<char*>(&val), sizeof(val));
    int temp_numofRings = val;

    in.read(reinterpret_cast<char*>(&val), sizeof(val));
    int temp_numofRings5 = val;

    in.read(reinterpret_cast<char*>(&val), sizeof(val));
    name = new char[val];
    in.read(name, val);
    std::string m_SiteName = string(name, val);
    delete [] name;

    PAHProcess p(*this);
    p.initialise(m_SiteName, temp_numofRings, temp_numofRings5);
}

void PAHStructure::WriteCposition(std::ostream &out) const
{
	double val = 0.0;

	std::set<cpair>::iterator itEnd = m_cpositions.end();
	for (std::set<cpair>::iterator it = m_cpositions.begin(); it != itEnd; ++it)
	{
		val = (*it).first;
		out.write((char*)&val, sizeof(val));
		val = (*it).second;
		out.write((char*)&val, sizeof(val));
	}
}

std::list<Site> PAHStructure::GetSiteList() const {
	return m_siteList;
}

std::map<kmcSiteType, svector> PAHStructure::GetSiteMap() const {
	return m_siteMap;
}
// the size for m_cpositions is required obviously, otherwise, the codes will not know when to stop
void PAHStructure::ReadCposition(std::istream &in, const int size)
{
    double val=0.0;
    m_cpositions.clear();
    int m_first=0;
    int m_second=0;

    cpair position;
    for (int i=0; i!=size;++i)
    {
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_first = (int)val;
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_second = (int)val;

        position = make_pair(m_first, m_second);
        m_cpositions.insert(m_cpositions.end(), position);
    }
}


