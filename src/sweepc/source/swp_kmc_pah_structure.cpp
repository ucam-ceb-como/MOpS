/*!
  * \author     Zakwan Zainuddin (zz260) && Gustavo Leon (gl413)
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
	m_methyl_counts = 0;
    m_rings = 0;
	m_rings5_Lone = 0;
	m_rings5_Embedded = 0;
	m_rings7_Lone = 0;
	m_rings7_Embedded = 0;
    m_counts = intpair(0, 0);
	m_optimised = false;
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
		//m_rings5_Lone = rhs.m_rings5_Lone;
		//m_rings5_Embedded = rhs.m_rings5_Embedded;
		//m_rings7_Lone = rhs.m_rings7_Lone;
		//m_rings7_Embedded = rhs.m_rings7_Embedded;
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
	m_methyl_counts = 0;
    m_cpositions.clear();
	m_rings5_Lone = 0;
	m_rings5_Embedded = 0;
	m_rings7_Lone = 0;
	m_rings7_Embedded = 0;
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

int PAHStructure::numofCH3() const {
    return m_methyl_counts;
} 

int PAHStructure::numofRings() const{
    return m_rings;
}

int PAHStructure::numofLoneRings5() const{
	return m_rings5_Lone;
}

int PAHStructure::numofEmbeddedRings5() const{
	return m_rings5_Embedded;
}

int PAHStructure::numofLoneRings7() const{
	return m_rings7_Lone;
}

int PAHStructure::numofEmbeddedRings7() const{
	return m_rings7_Embedded;
}

int PAHStructure::numofEdgeC() const{
    // the m_cpositions stores the coordinates of PAH, which means the num of edge C equals the size of m_cpositions
    return m_cpositions.size();
}

int PAHStructure::numofSite() const
{
    return m_siteList.size();
}

int PAHStructure::numofBridges() const
{
	int num = 0;
	for(Ccontainer::const_iterator i=m_carbonList.begin(); i!=m_carbonList.end(); i++) {
        if((*i)->bridge) num++;
    }
    num /= 2;
    return num;
}

void PAHStructure::setnumofC(int val)
{
    m_counts.first=val;
}

void PAHStructure::setnumofH(int val)
{
    m_counts.second=val;
}

void PAHStructure::setnumofCH3(int val)
{
    m_methyl_counts=val;
}

void PAHStructure::setnumofRings(int val)
{
    m_rings=val;
}

void PAHStructure::setnumofLoneRings5(int val)
{
	m_rings5_Lone = val;
}

void PAHStructure::setnumofEmbeddedRings5(int val)
{
	m_rings5_Embedded = val;
}

void PAHStructure::setnumofLoneRings7(int val)
{
	m_rings7_Lone = val;
}

void PAHStructure::setnumofEmbeddedRings7(int val)
{
	m_rings7_Embedded = val;
}

void PAHStructure::initialise(StartingStructure ss) {
    PAHProcess p(*this);
    p.initialise(ss);
}

void PAHStructure::initialise_fromfile() {
    PAHProcess p(*this);
    p.initialise_fromfile();
}

PAHStructure*  PAHStructure::Clone() {
    PAHProcess p(*this);
    return p.clonePAH();
}

bool PAHStructure::havebridgeC(){
    PAHProcess p(*this);
    return p.havebridgeC();
}

bool PAHStructure::hasCurvedSubunits(){
	if (numofRings() >= 8 || numofEmbeddedRings5() >= 1 || numofEmbeddedRings7() >= 1) return true;
	else {
		std::list<Site> site_list =  GetSiteList();
		for(Spointer i=site_list.begin(); i!= site_list.end(); i++) {
			if ( (int)i->type >= 501 && (int)i->type <= 504 && (int)i->type%10 >=3) return true;
			else if ( (int)i->type >= 600 ) return true;
		}
		return false;
	}
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
	coordtype val_pos = 0.0;

    // output info for PAHProcess::createPAH().
    val=numofRings();
    out.write((char*)&(val), sizeof(val));

	val = numofLoneRings5();
	out.write((char*)&(val), sizeof(val));

	val = numofEmbeddedRings5();
	out.write((char*)&(val), sizeof(val));

	val = numofLoneRings7();
	out.write((char*)&(val), sizeof(val));

	val = numofEmbeddedRings7();
	out.write((char*)&(val), sizeof(val));
	
	val = numofC();
	out.write((char*)&(val), sizeof(val));
	
	val = numofH();
	out.write((char*)&(val), sizeof(val));
	
	val = numofCH3();
	out.write((char*)&(val), sizeof(val));
	
	val = numofSite();
	out.write((char*)&(val), sizeof(val));
	
	//Serializes m_siteList
	for(std::list<Site>::const_iterator it = m_siteList.begin(); it != m_siteList.end(); ++it){
		val = (int)(*it).type;
		out.write((char*)&(val), sizeof(val));
		val = (int)(*it).comb;
		out.write((char*)&(val), sizeof(val));
		val_pos = std::get<0>((*it).C1->coords);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<1>((*it).C1->coords);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<2>((*it).C1->coords);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<0>((*it).C2->coords);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<1>((*it).C2->coords);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<2>((*it).C2->coords);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
	}

	//Serializes m_siteMap
	val = m_siteMap.size();
	out.write((char*)&(val), sizeof(val));

	std::map<kmcSiteType, svector>::const_iterator map_it;

	for(map_it=m_siteMap.begin(); map_it!=m_siteMap.end(); map_it++){
		val = (int)map_it->first;
		out.write((char*)&(val), sizeof(val));

		std::vector<Spointer> site_vector = map_it->second;
		val = site_vector.size();
		out.write((char*)&(val), sizeof(val));


		for(int site_it=0; site_it<site_vector.size(); site_it++){
			Spointer S1 = site_vector[site_it];
			val = (int)S1->type;
			out.write((char*)&(val), sizeof(val));
			val = (int)S1->comb;
			out.write((char*)&(val), sizeof(val));
			val_pos = std::get<0>(S1->C1->coords);
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
			val_pos = std::get<1>(S1->C1->coords);
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
			val_pos = std::get<2>(S1->C1->coords);
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
			val_pos = std::get<0>(S1->C2->coords);
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
			val_pos = std::get<1>(S1->C2->coords);
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
			val_pos = std::get<2>(S1->C2->coords);
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		}
	}

	//Old method that stored the site list as a string
    //PAHStructure m_copy (*this);
    //PAHProcess p(m_copy);
    //std::string m_SiteName = p.SiteString(',');

    //val = (unsigned int)m_SiteName.length();
    //out.write((char*)&val, sizeof(val));
    //out.write(m_SiteName.c_str(), val);
	
	//Save number of bridges
	val = numofBridges();
	out.write((char*)&(val), sizeof(val));
	
	//Save edge carbon coordinates
	val = numofEdgeC();
	out.write((char*)&(val), sizeof(val));
	
	Cpointer Cnow = m_cfirst;
	Cpointer Cprev = m_cfirst;
	
	//If a carbon is bridged it is stored twice. This way we ensure to recover the pointers correctly.
	int values_to_write = numofBridges()*2 + numofEdgeC();
	int values_written = 0;
	do {
		val_pos = std::get<0>(Cnow->coords);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<1>(Cnow->coords);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<2>(Cnow->coords);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		
		if (Cnow-> bridge) val = 1;
		else val = 0;
		out.write((char*)&val, sizeof(val));
		if (Cnow->A != 'C'){
			if (Cnow->A == 'H') val = 1;
			if (Cnow->A == 'M') val = 2;
			out.write((char*)&val, sizeof(val));
			val_pos = std::get<0>(Cnow->growth_vector);
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
			val_pos = std::get<1>(Cnow->growth_vector);
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
			val_pos = std::get<2>(Cnow->growth_vector);
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		}
		else {
			val = 0;
			out.write((char*)&val, sizeof(val));
			val_pos = 0.0;
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
			out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		}
		
		if (Cnow-> bridge && Cprev != Cnow->C3){
			Cprev = Cnow;
			Cnow = Cnow->C3;
		}
		else{
			Cprev = Cnow;
			Cnow = Cnow->C2;
		}
		values_written++;
	}while (values_written!=values_to_write);
	
	//Save internal coordinates
	val = m_InternalCarbons.size();
	out.write((char*)&(val), sizeof(val));
	for(std::list<cpair>::const_iterator it = m_InternalCarbons.begin(); it != m_InternalCarbons.end(); ++it){
		val_pos = std::get<0>(*it);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<1>(*it);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<2>(*it);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
	}
	
	//Save R5 locations
	val = m_R5loc.size();
	out.write((char*)&(val), sizeof(val));
	for(std::list<cpair>::const_iterator it = m_R5loc.begin(); it != m_R5loc.end(); ++it){
		val_pos = std::get<0>(*it);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<1>(*it);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<2>(*it);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
	}
	
	//Save R7 locations
	val = m_R7loc.size();
	out.write((char*)&(val), sizeof(val));
	for(std::list<cpair>::const_iterator it = m_R7loc.begin(); it != m_R7loc.end(); ++it){
		val_pos = std::get<0>(*it);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<1>(*it);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
		val_pos = std::get<2>(*it);
		out.write(reinterpret_cast<const char *>(&val_pos), sizeof(val_pos));
	}
	
	
}

void PAHStructure::Deserialize(std::istream &in)
{
    int val = 0;
	coordtype val_pos = 0.0;
    char *name = NULL;

    in.read(reinterpret_cast<char*>(&val), sizeof(val));
    int temp_numofRings = val;

	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofLoneRings5 = val;

	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofEmbeddedRings5 = val;

	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofLoneRings7 = val;

	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofEmbeddedRings7 = val;
	
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofC = val;
	
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofH = val;
	
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofCH3 = val;
	
	//Old method that stored sites as a list
    //in.read(reinterpret_cast<char*>(&val), sizeof(val));
    //name = new char[val];
    //in.read(name, val);
    //std::string m_SiteName = string(name, val);
	
	//delete [] name;
	
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numSite = val;
	
	std::vector<std::tuple<int, int, cpair, cpair>> temp_site_vector;
	for(int i=0; i!= temp_numSite; i++){
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		int site_type = val;
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		int comb_site_type = val;
		
		coordtype x1, y1, z1, x2, y2, z2;
		std::tuple<coordtype, coordtype, coordtype> s_C1, s_C2;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		x1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		y1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		z1 = val_pos;
		s_C1 = std::make_tuple(x1, y1, z1);
		
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		x2 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		y2 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		z2 = val_pos;
		s_C2 = std::make_tuple(x2, y2, z2);
		std::tuple<int, int, cpair, cpair> temp_site = std::make_tuple(site_type, comb_site_type, s_C1, s_C2);
		
		temp_site_vector.push_back(temp_site);
	}

	//Deserializes m_siteMap
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_mapSize = val;

	std::map<int, std::vector<std::tuple<int, int, cpair, cpair>>> temp_map;

	for(int map_it = 0; map_it!=temp_mapSize; map_it++){
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		int temp_sitetype = val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		int temp_sitetypesize = val;

		std::vector<std::tuple<int, int, cpair, cpair>> temp_sitemap_vector;

		for(int site_it=0; site_it!=temp_sitetypesize; site_it++){
			in.read(reinterpret_cast<char*>(&val), sizeof(val));
			int site_type = val;
			in.read(reinterpret_cast<char*>(&val), sizeof(val));
			int comb_site_type = val;
			
			coordtype x1, y1, z1, x2, y2, z2;
			std::tuple<coordtype, coordtype, coordtype> s_C1, s_C2;
			in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
			x1 = val_pos;
			in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
			y1 = val_pos;
			in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
			z1 = val_pos;
			s_C1 = std::make_tuple(x1, y1, z1);
			
			in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
			x2 = val_pos;
			in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
			y2 = val_pos;
			in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
			z2 = val_pos;
			s_C2 = std::make_tuple(x2, y2, z2);
			std::tuple<int, int, cpair, cpair> temp_site = std::make_tuple(site_type, comb_site_type, s_C1, s_C2);
			
			temp_sitemap_vector.push_back(temp_site);
		}
		temp_map[temp_sitetype] = temp_sitemap_vector;
	}
	
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofBridges = val;
	
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofEdgeC = val;
	
	std::vector<std::tuple<cpair, int, int, cpair>> edgeCarbons;
	int c_counter = 0;
	do{
		double x1, y1, z1;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		x1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		y1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		z1 = val_pos;
		cpair Ccoords = std::make_tuple(x1, y1, z1);
		
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		int bridge_int = val;
		
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		int Hatom_int = val;
		
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		x1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		y1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		z1 = val_pos;
		cpair Cgrovec = std::make_tuple(x1, y1, z1);
		
		std::tuple<cpair, int, int, cpair> C_read = std::make_tuple(Ccoords, bridge_int, Hatom_int, Cgrovec);
		edgeCarbons.push_back(C_read);
		c_counter += 1;
	}while (c_counter < temp_numofEdgeC + temp_numofBridges * 2);
	
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofInternalC = val;
		
	std::list<cpair> temp_internalcoords;
	for (int i =0; i!= temp_numofInternalC; i++){
		double x1, y1, z1;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		x1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		y1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		z1 = val_pos;
		cpair Ccoords = std::make_tuple(x1, y1, z1);
		temp_internalcoords.push_back(Ccoords);
	}
	
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofR5loc = val;
		
	std::list<cpair> temp_R5loc;
	for (int i =0; i!= temp_numofR5loc; i++){
		double x1, y1, z1;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		x1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		y1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		z1 = val_pos;
		cpair R5coords = std::make_tuple(x1, y1, z1);
		temp_R5loc.push_back(R5coords);
	}
	
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	int temp_numofR7loc = val;
		
	std::list<cpair> temp_R7loc;
	for (int i =0; i!= temp_numofR7loc; i++){
		double x1, y1, z1;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		x1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		y1 = val_pos;
		in.read(reinterpret_cast<char*>(&val_pos), sizeof(val_pos));
		z1 = val_pos;
		cpair R7coords = std::make_tuple(x1, y1, z1);
		temp_R7loc.push_back(R7coords);
	}
	
    PAHProcess p(*this);
	//Previous method
	//p.initialise_sitelist_string(m_SiteName, temp_numofRings, temp_numofLoneRings5, temp_numofEmbeddedRings5, temp_numofLoneRings7, temp_numofEmbeddedRings7, temp_numofC, temp_numofH, temp_numofCH3, temp_internalcoords);
	p.initialise(temp_site_vector, temp_map, temp_numofRings, temp_numofLoneRings5, temp_numofEmbeddedRings5, temp_numofLoneRings7, temp_numofEmbeddedRings7, edgeCarbons, temp_internalcoords, temp_R5loc, temp_R7loc);
}

void PAHStructure::WriteCposition(std::ostream &out) const
{
	double val = 0.0;

	std::set<cpair>::iterator itEnd = m_cpositions.end();
	for (std::set<cpair>::iterator it = m_cpositions.begin(); it != itEnd; ++it)
	{
		val = std::get<0>(*it);
		out.write((char*)&val, sizeof(val));
		val = std::get<1>(*it);
		out.write((char*)&val, sizeof(val));
		val = std::get<2>(*it);
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
    double m_first=0;
	double m_second = 0;
	double m_third = 0;

    cpair position;
    for (int i=0; i!=size;++i)
    {
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_first = (int)val;
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_second = (int)val;
		in.read(reinterpret_cast<char*>(&val), sizeof(val));
		m_third = (int)val;

        position = std::make_tuple(m_first, m_second, m_third);
        m_cpositions.insert(m_cpositions.end(), position);
    }
}


