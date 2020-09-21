/*!
  * \author     Zakwan Zainuddin (zz260) && Gustavo Leon (gl413)
  * \file       swp_kmc_pah_process.cpp
  *
  * \brief        Implementation file for swp_kmc_pah_process.h.
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Defines the data structure which holds information of PAH molecules

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

#include "swp_kmc_pah_process.h"
#include "swp_kmc_jump_process.h"
#include "string_functions.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "choose_index.hpp"

#include <openbabel/babelconfig.h>
#include <openbabel/plugin.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h> //Used to export to xyz
#include <openbabel/forcefield.h>
#include <openbabel/ring.h>
#include <openbabel/math/vector3.h>

using namespace Sweep;
using namespace Sweep::KMC_ARS;
using namespace std;

static std::vector<kmcSiteType> PHsites = vectPHsites();
// Constructors and Destructor
//! Default Constructor
PAHProcess::PAHProcess() {
    m_pah = NULL;
    m_rates_save = false;
}
//! Overloaded constructor; specifying PAH structure to perform processes on
PAHProcess::PAHProcess(PAHStructure& pah) {
    m_pah = &pah;
    m_rates_save = false;
}
//! Copy Constructor
PAHProcess::PAHProcess(const PAHProcess &pahp) {
    m_pah = pahp.m_pah;
    m_rates_save = pahp.m_rates_save;
}
//! Default Destructor
PAHProcess::~PAHProcess() {
}

//! Sets target PAH structure to perform processes on
void PAHProcess::setPAH(PAHStructure& pah) {
    m_pah = &pah;
}

PAHStructure* PAHProcess::returnPAH(){
    return this->m_pah;
}
//! Returns a copy of PAH structure
PAHStructure* PAHProcess::clonePAH() const {
    PAHStructure* temp = new PAHStructure();
    PAHProcess p(*temp);
    //std::vector<kmcSiteType> sites = SiteVector();
	//std::vector<int> carbon_per_site = SiteIntVector();
	//Spointer first_site = SiteList().begin();
	//cpair first_c_pos = first_site->C1->coords;
	
	std::vector<std::tuple<int, int, cpair, cpair>> sites = SiteVector_clone();
	std::map<int, std::vector<std::tuple<int, int, cpair, cpair>>> site_map = SiteMap_clone();
	std::vector<std::tuple<cpair, int, int, cpair>> edgeCarbons = EdgeCarbonVector_clone();
	p.initialise(sites, site_map, m_pah->m_rings, m_pah->m_rings5_Lone, m_pah->m_rings5_Embedded, m_pah->m_rings7_Lone, m_pah->m_rings7_Embedded, edgeCarbons, m_pah->m_InternalCarbons, m_pah->m_R5loc, m_pah->m_R7loc);
	//p.createPAH(sites, carbon_per_site, m_pah->m_rings, m_pah->m_rings5_Lone, m_pah->m_rings5_Embedded, m_pah->m_rings7_Lone, m_pah->m_rings7_Embedded, m_pah->m_carbonList, m_pah->m_InternalCarbons, m_pah->m_R5loc, m_pah->m_R7loc, first_c_pos);
	//p.createPAH(sites, m_pah->m_rings, m_pah->m_rings5_Lone, m_pah->m_rings5_Embedded, m_pah->m_rings7_Lone, m_pah->m_rings7_Embedded, m_pah->m_InternalCarbons);
    return temp;
}
// Public Read Processes
//! Returns C and H counts
intpair PAHProcess::getCHCount() const {
    return m_pah->m_counts;
}
/*!
 * @brief Get site count.
 *
 * @param[in]    st    Site type.
 *
 * @return The number of the specified site type.
 */ 
unsigned int PAHProcess::getSiteCount(const kmcSiteType& st) const {
    if(m_rates_save) {
        if(st==benz) return 5;
		else if(st==FE2) return 2;
        else return 1;
    }
    if(st==benz) {
        unsigned int sum=0;
        for(int i=0; i!=(int)PHsites.size(); i++) {
			for(int j=0; j<(int)m_pah->m_siteMap[PHsites[i]].size(); j++) {
				if (m_pah->m_siteMap[PHsites[i]][j]->C1->A == 'H' && m_pah->m_siteMap[PHsites[i]][j]->C2->A == 'H' ){
					sum ++;
				}
			}
        }
        return sum;
    }
	if(st==Methyl){
		int num = numberOfMethyl();
		return num;
	}
    if(st==FE3) {
        if(getCHCount().first == 6) return 0;
        return (unsigned int) (m_pah->m_siteMap[st].size());
    }
	unsigned int site_number = m_pah->m_siteMap[st].size();
	if (site_number >= 1) {
		for(int i=0; i<(int)m_pah->m_siteMap[st].size(); i++) {
			if(m_pah->m_siteMap[st][i]->C1->A == 'M' || m_pah->m_siteMap[st][i]->C2->A == 'M') site_number--;
		}
	}
	return site_number;
}
/*
//Gets site count and then substracts the number of sites containing methyl.
unsigned int PAHProcess::getSiteCount(const kmcSiteType& st) const {
    if(m_rates_save) {
        if(st==benz) return 5;
        else return 1;
    }
    if(st==benz) {
        unsigned int sum=0;
        for(int i=0; i!=(int)PHsites.size(); i++) {
            sum += (unsigned int) m_pah->m_siteMap[PHsites[i]].size();
        }
        return sum;
    }
    if(st==FE3) {
        if(getCHCount().first == 6) return 0;
        return (unsigned int) (m_pah->m_siteMap[st].size());
    }
    return (unsigned int) m_pah->m_siteMap[st].size();
}*/


//! Get Ring Counts
std::tuple <int, int, int> PAHProcess::getRingsCount() const {
	std::tuple <int, int, int> rings_tuple;
	if (m_rates_save) {
		rings_tuple = std::make_tuple(1,1,1);
	}
	else{
		rings_tuple = std::make_tuple(m_pah->m_rings, m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded, m_pah->m_rings7_Lone + m_pah->m_rings7_Embedded);
	}
	return rings_tuple;
}

//! Get R5 embedded count
int PAHProcess::getR5EmbeddedCount() const{
	int R5_embedded = m_pah->m_rings5_Embedded;
	return R5_embedded;
}

//! Get R7 embedded count
int PAHProcess::getR7EmbeddedCount() const{
	int R7_embedded = m_pah->m_rings7_Embedded;
	return R7_embedded;
}
//! Get number of bridges
int PAHProcess::numberOfBridges() const {
    int num = 0;
    for(Ccontainer::const_iterator i=m_pah->m_carbonList.begin(); i!=m_pah->m_carbonList.end(); i++) {
        if((*i)->bridge) num++;
    }
    num /= 2;
    return num;
}
//! Get number of methyl moeities
int PAHProcess::numberOfMethyl() const {
    int num = 0;
    for(Ccontainer::const_iterator i=m_pah->m_carbonList.begin(); i!=m_pah->m_carbonList.end(); i++) {
        if((*i)->A=='M') num++;
    }
    return num;
}
//! Is probably curved (has embedded or partially embedded R5s, R7s, more than 12 rings)
bool PAHProcess::HasCurvedMoeities() const {
	if (std::get<0>(getRingsCount()) >= 12 || getR5EmbeddedCount() >= 1 || getR7EmbeddedCount() >= 1 || getSiteCount(SPIRAL) >= 1) return true;
	else {
		for(Spointer i=m_pah->m_siteList.begin(); i!= m_pah->m_siteList.end(); i++) {
			if ( (int)i->type >= 1000 && (int)i->type <= 1005) return true;
			if ( (int)i->type >= 2100) return true;
		}
		return false;
	}
}

//! Print Structure in console
void PAHProcess::printStruct() const{
    // iterator set to first C in structure
    Cpointer prev = m_pah->m_cfirst;
    Cpointer now;
	bool bridged_first=false;
	if (prev->bridge) {
		now = m_pah->m_cfirst->C3;
		bridged_first=true;
	}
	else now = m_pah->m_cfirst->C2;
    // counts the number of C at edge
    unsigned int count = 0;
    // a string to display if current C is a bridge
    std::string btxt;
    angletype angle;
	std::cout << "-----------------------------" << std::endl;
	std::cout << "Printing detailed structure." << std::endl;
    do {
        count++;
        // checks if bridge
        if(prev->bridge) {
            btxt = " bridge";
            angle = prev->bondAngle2;
        } else {
            btxt = "       ";
            angle = prev->bondAngle1;
        }
        // outputs text in the form of "X: (C-H or C-C)    bond angle with next C"
		std::cout << "coords (" << std::get<0>(prev->coords) << "," << std::get<1>(prev->coords) << "," << std::get<2>(prev->coords) << ")" << btxt << ": C-" << prev->A << "\t" << angle << std::endl;
        if (bridged_first){
			if (now == m_pah->m_cfirst->C2 && prev == m_pah->m_cfirst){
				bridged_first = false;
			}
		}
		// moves iterator to next C
        Cpointer oldnow = now;
        now = moveCPointer(prev, now);
        prev = oldnow;
    } while 
		(!(count != 1 && std::get<0>(prev->coords) == std::get<0>(m_pah->m_cfirst->coords)
		&& std::get<1>(prev->coords) == std::get<1>(m_pah->m_cfirst->coords)
		&& std::get<2>(prev->coords) == std::get<2>(m_pah->m_cfirst->coords)
		&& !bridged_first));
    // displays C and H counts
    std::cout << "C Count: " << m_pah->m_counts.first << std::endl;
    std::cout << "H Count: " << m_pah->m_counts.second << std::endl;
	std::cout << "-----------------------------" << std::endl;
}

bool PAHProcess::havebridgeC(){
    for(Ccontainer::const_iterator i=m_pah->m_carbonList.begin(); i!=m_pah->m_carbonList.end(); i++) {
        if((*i)->bridge) return true;
    }
    return false;
}
//! Print Structure in console, with arrow pointing at current C
void PAHProcess::printStruct(Cpointer c) const{
    // iterator set to first C in structure
    Cpointer prev = m_pah->m_cfirst;
    Cpointer now = m_pah->m_cfirst->C2;
    // counts the number of C at edge
    unsigned int count = 0;
    // a string to display if current C is a bridge
    std::string btxt;
    angletype angle;
    do {
        count++;
        // checks if bridge
        if(prev->bridge) {
            btxt = " bridge";
            angle = prev->bondAngle2;
        } else {
            btxt = "       ";
            angle = prev->bondAngle1;
        }
        // outputs text in the form of "X: (C-H or C-C)    bond angle with next C"
        std::cout << prev << btxt << ": C-" << prev->A << "\t" << angle;
        if(prev == c) std::cout << "\t<---";
        std::cout << '\n';
        // moves iterator to next C
        Cpointer oldnow = now;
        now = moveCPointer(prev, now);
        prev = oldnow;
    } while (prev != m_pah->m_cfirst);
    // displays C and H counts
    std::cout << "C Count: " << m_pah->m_counts.first << '\n';
    std::cout << "H Count: " << m_pah->m_counts.second << '\n';
}
//! Print Sites in console
void PAHProcess::printSites() const{
    Spointer i;
    std::string st;
    cout << "*******************\n";
    cout << "Sites List:\n_____\n";
    // displays total site count
    cout << "Total Site Count: " << m_pah->m_siteList.size() << "\t Total Edge Carbons: " << m_pah->m_carbonList.size() << "\t\t\t Total Internal Carbons: " << m_pah->m_InternalCarbons.size()<<'\n';
	cout << "Total R6 rings: " << m_pah->m_rings << "\t Total R5 rings: " << m_pah->m_rings5_Lone << " lone + "<< m_pah->m_rings5_Embedded << " embedded \t Total R7 rings: " << m_pah->m_rings7_Embedded << '\n';
    for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
        // convert site type into string
        st = kmcSiteName(i->type);
        // displays site type
        cout << st << '\n';
	}
    cout << "********************\n";
}
//! Print Sites in console, with an arrow pointing at site stt
void PAHProcess::printSites(Spointer& stt) const{
    Spointer i;
    std::string st;
    string pointer = "  <-  ";
    cout << "*******************\n";
    cout << "Sites List:\n_____\n";
    // displays total site count
    cout << "Total Site Count: " << m_pah->m_siteList.size() << "\t Total Edge Carbons: " << m_pah->m_carbonList.size() << "\t\t\t Total Internal Carbons: " << m_pah->m_InternalCarbons.size()<<'\n';
	cout << "Total R6 rings: " << m_pah->m_rings << "\t Total R5 rings: " << m_pah->m_rings5_Lone << " lone + "<< m_pah->m_rings5_Embedded << " embedded \t Total R7 rings: " << m_pah->m_rings7_Embedded << '\n';
    for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
        // convert site type into string
        st = kmcSiteName(i->type);
        st.append("\t");
        st.append(kmcSiteName(i->comb));
        // checks if site i equivalent to stt, if yes put an arrow
        if(i == stt) st.append(pointer);
        // displays site type (and arrow)
        cout << st << '\n';
	}
    cout << "********************\n";
}

//! Print Sites in console, with stars pointing to starting sites and arrows pointing to current sites.
void PAHProcess::printSitesMigration() {
    Spointer i;
    std::string st;
    string current_position = "  <--  ";
	string starting_position = "  ***  ";
	string start_current_position = "  ***<--  ";
    cout << "*******************\n";
    cout << "Sites List:\n_____\n";
    // displays total site count
    cout << "Total Site Count: " << m_pah->m_siteList.size() << "\t Total Edge Carbons: " << m_pah->m_carbonList.size() << "\t\t\t Total Internal Carbons: " << m_pah->m_InternalCarbons.size()<<'\n';
	cout << "Total R6 rings: " << m_pah->m_rings << "\t Total R5 rings: " << m_pah->m_rings5_Lone << " lone + "<< m_pah->m_rings5_Embedded << " embedded \t Total R7 rings: " << m_pah->m_rings7_Embedded << '\n';
	cout << "Total walkers identified: " << m_pah->m_R5walker_sites.size() << std::endl;
    for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
        // convert site type into string
        st = kmcSiteName(i->type);
        st.append("\t");
        st.append(kmcSiteName(i->comb));
		for (unsigned int ii = 0;ii != m_pah->m_R5walker_sites.size();ii++){
			Spointer migr_site_start = std::get<0>(m_pah->m_R5walker_sites[ii]);
			int steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
			Spointer site_check = moveIt(migr_site_start, steps);
			Spointer migr_site_start_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
			Spointer site_check_2 = moveIt(migr_site_start_2, steps);
			
			// checks if site i equivalent to stt, if yes put an arrow
			if (steps == 0){
				if (migr_site_start == migr_site_start_2){
					if (migr_site_start == i) {
						st.append(start_current_position);
						st.append(std::to_string(ii));
					}
				}else{
					if (migr_site_start == i) {
						st.append(start_current_position);
						st.append(std::to_string(ii));
						st.append("A");
					}
					if (migr_site_start_2 == i) {
						st.append(start_current_position);
						st.append(std::to_string(ii));
						st.append("B");
					}
				}
			}else{
				if (migr_site_start == migr_site_start_2){
					if (migr_site_start == i) {
						st.append(starting_position);
						st.append(std::to_string(ii));
					}
				}else{
					if (migr_site_start == i) {
						st.append(starting_position);
						st.append(std::to_string(ii));
						st.append("A");
					}
					if (migr_site_start_2 == i) {
						st.append(starting_position);
						st.append(std::to_string(ii));
						st.append("B");
					}
				}
				if (migr_site_start == migr_site_start_2){
					if (site_check == i) {
						st.append(current_position);
						st.append(std::to_string(ii));
					}
				}else{
					if (site_check == i) {
						st.append(current_position);
						st.append(std::to_string(ii));
						st.append("A");
					}
					if (site_check_2 == i) {
						st.append(current_position);
						st.append(std::to_string(ii));
						st.append("B");
					}
				}
			}
		}
		// displays site type (and arrow)
		std::cout << st << '\n';
	} 
    cout << "********************\n";
}

//! Create a copy of the sites type
std::list<std::string> PAHProcess::copySites() const{
	std::list<std::string> BeforeJPSite;
    Spointer i;
    std::string st;
    for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
        // convert site type into string
        st = kmcSiteName(i->type);
        st.append("\t");
        st.append(kmcSiteName(i->comb));
        BeforeJPSite.push_back(st);
    }
	return BeforeJPSite;
}
//! Create a copy of the sites type
std::list<std::string> PAHProcess::copySites(Spointer& stt) const{
    std::list<std::string> BeforeJPSite;
	Spointer i;
    std::string st;
    string pointer = "  <-  ";
    // displays total site count
    for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
        // convert site type into string
        st = kmcSiteName(i->type);
        st.append("\t");
        st.append(kmcSiteName(i->comb));
        // checks if site i equivalent to stt, if yes put an arrow
        if(i == stt) st.append(pointer);
        // displays site type (and arrow)
        BeforeJPSite.push_back(st);
    }
	return BeforeJPSite;
}
//! Print copied sites before JP in console
void PAHProcess::printBeforeSites(std::list<std::string>& before_site_list) const{
    cout << "*******************\n";
    cout << "Sites List Before Jump Process:\n_____\n";
    // displays total site count
    cout << "Total Site Count Before Jump Process: " << m_pah->m_siteList.size() << '\n';
	std::list<std::string>::iterator i;
    for(i=before_site_list.begin(); i!=before_site_list.end(); i++) {
        cout << *i << '\n';
    }
    cout << "********************\n";
}

//! Print sites & site members in console, with arrow pointing at site stt
void PAHProcess::printSitesMemb(Spointer& stt) const{
    Spointer i;
    std::string st;
    string pointer = "  <-  ";
    cout << "*******************\n";
    cout << "Sites List:\n_____\n";
    // displays total site count
    cout << "Total Site Count: " << m_pah->m_siteList.size() << "\t Total Edge Carbons: " << m_pah->m_carbonList.size() << "\t\t\t Total Internal Carbons: " << m_pah->m_InternalCarbons.size()<<'\n';
	cout << "Total R6 rings: " << m_pah->m_rings << "\t Total R5 rings: " << m_pah->m_rings5_Lone << " lone + "<< m_pah->m_rings5_Embedded << " embedded \t Total R7 rings: " << m_pah->m_rings7_Embedded << '\n';
    for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
        // convert site type into string
        cout<< kmcSiteName(i->type);
        cout<<"\t";
        cout<<kmcSiteName(i->comb);
        cout<<'\t'<<i->C1<<'\t'<<i->C2;
        // checks if site i equivalent to stt, if yes put an arrow
        if(i == stt) cout<<pointer;
        // displays site type (and arrow)
        cout <<'\n';
	}
    cout << "********************\n";
}
void PAHProcess::printSitesMemb() const{
    Spointer i;
    std::string st;
    string pointer = "  <-  ";
    cout << "*******************\n";
    cout << "Sites List:\n_____\n";
    // displays total site count
    cout << "Total Site Count: " << m_pah->m_siteList.size() << "\t Total Edge Carbons: " << m_pah->m_carbonList.size() << "\t\t\t Total Internal Carbons: " << m_pah->m_InternalCarbons.size()<<'\n';
	cout << "Total R6 rings: " << m_pah->m_rings << "\t Total R5 rings: " << m_pah->m_rings5_Lone << " lone + "<< m_pah->m_rings5_Embedded << " embedded \t Total R7 rings: " << m_pah->m_rings7_Embedded << '\n';
    for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
        // convert site type into string
        cout<< kmcSiteName(i->type);
        cout<<"\t";
        cout<<kmcSiteName(i->comb);
        cout<<'\t'<<i->C1<<'\t'<<i->C2;
        // checks if site i equivalent to stt, if yes put an arrow
        //if(i == stt) cout<<pointer;
        // displays site type (and arrow)
        cout <<'\n';
    }
    cout << "********************\n";
}
//! Store results in external file, returns success/failure
// Stores structure into a DOT file for Graphviz
bool PAHProcess::saveDOT(const std::string &filename) const {
    bool success = saveDOT(filename, std::string("none"));
    return success;
}
// Stores structure into a DOT file for Graphviz with a title on graph
bool PAHProcess::saveDOT(const std::string &filename, const std::string &title) const {
    ofstream dotfile;
    // creates or opens filename
    //cout << "Attempting to open " << filename << "...\n";//++++
    dotfile.open(filename.c_str());
    if(dotfile.is_open()) { // if successful
        // vector to store indices of C bonded with H atom
        std::vector<int> CH_index;
        // set iterator to starting position
        Cpointer now = m_pah->m_cfirst->C2;
        Cpointer prev = m_pah->m_cfirst;
        //cout << filename << " had been opened, prepare writing..\n";//++++
        // draw undirected graph
        dotfile<< "Graph G {\n";
        // sets max size for drawing structure (in inches)
        dotfile<< "\tgraph [size=15];\n";
        // sets style for C nodes
        dotfile<< "\tnode [shape=circle, style=filled, color=black];\n";
        // sets style for bonds
        dotfile<< "\tedge [weight=2, style=\"setlinewidth(10)\"];\n";
        // sets caption for graph
        if(title!="none") {
            dotfile<< "\tgraph [label=\""<<title<<"\"];\n";
        }
        // get coordinates of current C
        coordtype x = std::get<0>(prev->coords);
		coordtype y = std::get<1>(prev->coords);
        int c=0; // index of current C
        int count=0;
        float CC_CHratio = 10; // length ratio between C-C and C-H bonds
        do {
            //if(now->bridge) {
                //cout<<"Bridge at "<<now<<'\n';
            //}
            c++;
            // write position coordinates of C atom with index c, with ! suffixed to pin node
            dotfile<<"\t"<<c<<" [pos = \""<<x<<','<<y<<"!\"";
            if(c==1) dotfile<<", color=red";
            dotfile<<", height=0.5, width=0.5, fixedsize=shape];\n";
			// stores index of current C
			CH_index.push_back(c);
            if(prev->A=='H') { // if bonded with H
				c++;
                // write position coordinates of H atom with bonded C index, with ! suffixed to pin node
                dotfile<<"\t"<<c<<" [pos = \""<<(x+x_inc(prev->bondAngle2)/CC_CHratio);
                // set style for H node
				dotfile << ',' << (y + y_inc(prev->bondAngle2) / CC_CHratio) << "!\", color=skyblue, fontcolor=skyblue, height=0.25, width=0.25, fixedsize=shape];\n";
            }
            // move iterator to next C
            Cpointer oldnow = now;
            now = moveCPointer(prev, now);
            prev = oldnow;
            // obtain coordinates
			x = std::get<0>(prev->coords);
			y = std::get<1>(prev->coords);

        }
        while
		(!(count != 1 && std::get<0>(prev->coords) == std::get<0>(m_pah->m_cfirst->coords)
		&& std::get<1>(prev->coords) == std::get<1>(m_pah->m_cfirst->coords)));

		// get coordinates of first internal C
		int c_copy = c;
		std::list<cpair>::iterator it;
		for (it = m_pah->m_InternalCarbons.begin(); it != m_pah->m_InternalCarbons.end(); ++it){
			c_copy++;
			coordtype x = std::get<0>(*it);
			coordtype y = std::get<1>(*it);
			dotfile << "\t" << c_copy << " [pos = \"" << x << ',' << y << "!\"";
			dotfile << ", height=0.25, width=0.25, fixedsize=shape];\n";
		}

        // connects each C atom
		for (int i = 1; i != (int)CH_index.size(); i++) {
			dotfile << "\t" << CH_index[i-1] << " -- " << CH_index[i] << ";\n";
		}
		int last = (int)CH_index.size() - 1;
		dotfile << "\t" << CH_index[last] << " -- " << "1;\n"; // connects last C with first C
        // connects each H atom with corresponding C atoms Not needed
        /*for(int i=0; i!=(int)CH_index.size(); i++) {
            dotfile<<"\t"<<CH_index[i]<<" -- "<<CH_index[i]<<";\n";
        }*/
        // close dot file
        dotfile<<"}";
        dotfile.close();
        //cout<<".DOT file writing complete.\n";//++++
        return true;
    }
    // if unable to open/create dot file
    else cout << "unable to open file\n";
    return false;
}

//! Stores structure into a XYZ file.
void PAHProcess::saveXYZ(const std::string &filename, bool optimise) {
	OpenBabel::OBMol mol = passPAH(optimise);
	if (optimise){
		if (HasCurvedMoeities()) mol = optimisePAH(mol, 8000, "mmff94");
		else mol = optimisePAH(mol, 4000, "mmff94");
		//mol = optimisePAH(mol, 4000, "Ghemical");
	}
	ofstream ofs1;
	std::string filename1 = filename;
	filename1.append(".xyz");
	OpenBabel::OBConversion conv;
	OpenBabel::OBFormat *format_out = conv.FindFormat("xyz"); // default output format
	conv.SetInAndOutFormats(format_out, format_out);
	ofs1.open(filename1);
	conv.Write(&mol, &ofs1);
	ofs1.close();
}

//! Stores structure into a XYZ file.
void PAHProcess::save_trajectory_xyz(const std::string &filename, bool optimise) {
	OpenBabel::OBMol mol = passPAH(optimise);
	if (optimise){
		mol = optimisePAH(mol, 4000, "mmff94");
	}
	ofstream ofs1;
	std::string filename1 = filename;
	filename1.append(".xyz");
	OpenBabel::OBConversion conv;
	OpenBabel::OBFormat *format_out = conv.FindFormat("xyz"); // default output format
	conv.SetInAndOutFormats(format_out, format_out);
	ofs1.open(filename1, std::ios::app);
	for (int i=1; i!=5; ++i)	conv.Write(&mol, &ofs1);
	ofs1.close();
}

//! Store results in external file, returns success/failure
// Stores structure into a DOT file for Graphviz
bool PAHProcess::saveDOT3D(const std::string &filename) const {
	bool success = saveDOT3D(filename, std::string("none"));
	return success;
}

// Stores structure into a DOT file for Graphviz with a title on graph for 3D files
bool PAHProcess::saveDOT3D(const std::string &filename, const std::string &title) const {
	ofstream dotfile;
	// creates or opens filename
	//cout << "Attempting to open " << filename << "...\n";//++++
	dotfile.open(filename.c_str());
	if (dotfile.is_open()) { // if successful
		// vector to store indices of C bonded with H atom
		std::vector<int> CH_index;
		// set iterator to starting position
		Cpointer now = m_pah->m_cfirst->C2;
		Cpointer prev = m_pah->m_cfirst;
		//cout << filename << " had been opened, prepare writing..\n";//++++
		// draw undirected graph
		dotfile << "Graph G {\n";
		// sets max size for drawing structure (in inches)
		dotfile << "\tgraph [size=15, dimen=3];\n";
		// sets style for C nodes
		dotfile << "\tnode [shape=point, style=filled, color=black];\n";
		// sets style for bonds
		dotfile << "\tedge [shape=cylinder, color=black];\n";
		// sets caption for graph
		if (title != "none") {
			dotfile << "\tgraph [label=\"" << title << "\"];\n";
		}
		// get coordinates of current C
		coordtype x = std::get<0>(prev->coords);
		coordtype y = std::get<1>(prev->coords);
		coordtype z = std::get<2>(prev->coords);
		int c = 0; // index of current C
		int count = 0;
		float CC_CHratio = 10; // length ratio between C-C and C-H bonds
		do {
			//if(now->bridge) {
			//cout<<"Bridge at "<<now<<'\n';
			//}
			c++;
			// write position coordinates of C atom with index c, with ! suffixed to pin node
			dotfile << "\t" << c << " [pos = \"" << x << ',' << y << ',' << z << "!\"";
			if (c == 1) dotfile << ", color=red";
			dotfile << ", height=0.5, width=0.5, fixedsize=shape];\n";
			// stores index of current C
			CH_index.push_back(c);
			if (prev->A == 'H') { // if bonded with H
				c++;
				// write position coordinates of H atom with bonded C index, with ! suffixed to pin node
				dotfile << "\t" << c << " [pos = \"" << (x + x_inc(prev->bondAngle2) / CC_CHratio);
				// set style for H node
				dotfile << ',' << (y + y_inc(prev->bondAngle2) / CC_CHratio);
				// set style for H node
				dotfile << ',' << (z + 10*sin(30*M_PI/180) / CC_CHratio) << "!\", color=skyblue, fontcolor=skyblue, height=0.25, width=0.25, fixedsize=shape];\n";
			}
			// move iterator to next C
			Cpointer oldnow = now;
			now = moveCPointer(prev, now);
			prev = oldnow;
			// obtain coordinates
			x = std::get<0>(prev->coords);
			y = std::get<1>(prev->coords);
			z = std::get<2>(prev->coords);

		} while
			(!(count != 1 && std::get<0>(prev->coords) == std::get<0>(m_pah->m_cfirst->coords)
			&& std::get<1>(prev->coords) == std::get<1>(m_pah->m_cfirst->coords)
			&& std::get<2>(prev->coords) == std::get<2>(m_pah->m_cfirst->coords)));

		// get coordinates of first internal C
		int c_copy = c;
		std::list<cpair>::iterator it;
		for (it = m_pah->m_InternalCarbons.begin(); it != m_pah->m_InternalCarbons.end(); ++it){
			c_copy++;
			coordtype x = std::get<0>(*it);
			coordtype y = std::get<1>(*it);
			coordtype z = std::get<2>(*it);
			dotfile << "\t" << c_copy << " [pos = \"" << x << ',' << y << ',' << z << "!\"";
			dotfile << ", height=0.25, width=0.25, fixedsize=shape];\n";
		}
		// connects each C atom
		for (int i = 1; i != (int)CH_index.size(); i++) {
			dotfile << "\t" << CH_index[i - 1] << " -- " << CH_index[i] << ";\n";
		}
		int last = (int)CH_index.size() - 1;
		dotfile << "\t" << CH_index[last] << " -- " << "1;\n"; // connects last C with first C
		// connects each H atom with corresponding C atoms
		/*for (int i = 0; i != (int)CH_index.size(); i++) {
			dotfile << "\tH" << CH_index[i] << " -- " << "C" << CH_index[i] << ";\n";
		}*/
		
		// close dot file
		dotfile << "}";
		dotfile.close();
		//cout<<".DOT file writing complete.\n";//++++
		return true;
	}
	// if unable to open/create dot file
	else cout << "unable to open file\n";
	return false;
}

// Protected Read Processes
//! Get other member of the site a particular C atom is a member of
Cpointer PAHProcess::getPair(const Cpointer Carb, bool after) const {/*
    if(after) {
        return Carb->S2->C2;
    }else {
        return Carb->S1->C1;
    }*/
    return NULL;
}
//! Returns the next carbon atom after current, coming from the previous
Cpointer PAHProcess::moveCPointer(Cpointer &previous, Cpointer &current) const {
    //cout << "Moving from " << current << " with " << previous << " preceding..\n";
    if(current == 0x0) {
        cout << "ERROR: Moving from NULL C Pointer\n\n";
        return NULL;
    }
    // check for bridge
    if(!(current->bridge)) {
    // if not bridge continue to next C
        //cout << "current is not bridge, going to " << current->C2 << '\n';
        return current->C2;
        //cout << "moved on to " << current << '\n';
    } else {
        if(previous == current->C1) {
            // if coming from main PAH, move iterator to bridged PAH
            return current->C3;
        }else if(previous == current->C3) {
            // if coming from bridged PAH, move to next on main PAH
            return current->C2;
        }
    }
    return current->C2;
    //cout << "moved on to " << current << '\n';
}
//! Check if process is allowed [returns true for now]
bool PAHProcess::allowed(const Spointer& st, StructureProc proc) const {
    return true;
}
//! Jump to a position coordinate given starting position and angle towards new position
cpair PAHProcess::jumpToPos(const cpair& starting, const angletype& direction, const angletype& direction2, const bondlength& distance) const{
    // calculate new coordinates
	cpair temp = std::make_tuple(std::get<0>(starting) +distance*x_inc(direction), std::get<1>(starting) +distance*y_inc(direction), std::get<2>(starting) +distance*z_inc(direction2));
	
	/*temp.first = starting.first + distance*x_inc(direction);
	temp.second = starting.second + distance*y_inc(direction);*/
    return temp;
}
//! Jump to a position coordinate given starting position and a vector.
cpair PAHProcess::jumpToPos(const cpair& starting, const cpair& direction, const bondlength& distance) const{
    //check if vector is unitary
	cpair temp = scale_vector(direction);
	// calculate new coordinates
	//cpair temp2 = std::make_tuple(std::get<0>(starting) +distance * std::get<0>(direction), std::get<1>(starting) +distance * std::get<1>(direction), std::get<2>(starting) +distance * std::get<2>(direction));
	cpair temp2 = std::make_tuple(std::get<0>(starting) +distance * std::get<0>(temp), std::get<1>(starting) +distance * std::get<1>(temp), std::get<2>(starting) +distance * std::get<2>(temp));
    return temp2;
}
//! Gets a vector between two points.
cpair PAHProcess::get_vector(cpair p1, cpair p2) const{
	cpair temp = std::make_tuple(std::get<0>(p2) - std::get<0>(p1), std::get<1>(p2) - std::get<1>(p1), std::get<2>(p2) - std::get<2>(p1));
	//check if vector is not 0
	double magnitude = sqrt(std::get<0>(temp)*std::get<0>(temp) + std::get<1>(temp)*std::get<1>(temp) + std::get<2>(temp)*std::get<2>(temp));
	if (magnitude < 1e-3){
		//This is an error. Distance between two carbons is too small. Print message and hope that a vector in Z direction makes sense.
		temp = std::make_tuple(0.0, 0.0, 1.0);
		cout << "Error in PAHProcess::get_vector. Passing same coordinate twice. Returning arbitrary vector <0,0,1>.\n";
	}
	//check if vector is unitary
	cpair temp2 = scale_vector(temp);
    return temp2;
}
//! Rescales a vector.
cpair PAHProcess::scale_vector(cpair vec) const{
	//check if vector is unitary
	double tol = 1e-3;
	double magnitude = sqrt(std::get<0>(vec)*std::get<0>(vec) + std::get<1>(vec)*std::get<1>(vec) + std::get<2>(vec)*std::get<2>(vec));
	cpair temp = vec;
	if (magnitude > tol)
	{
		temp = std::make_tuple(std::get<0>(vec)/magnitude, std::get<1>(vec)/magnitude, std::get<2>(vec)/magnitude);
		/*std::get<0>(temp) = std::get<0>(temp)/magnitude;
		std::get<1>(temp) = std::get<1>(temp)/magnitude;
		std::get<2>(temp) = std::get<2>(temp)/magnitude;*/
	}
	else{
		//Vector is really small. Dividing by a very small number risks raising errors.
		temp = std::make_tuple(std::get<0>(vec), std::get<1>(vec), std::get<2>(vec));
	}
	return temp;
}

//! Changes vector direction.
cpair PAHProcess::invert_vector(cpair vec) const{
	//check if vector is unitary
	cpair temp = scale_vector(vec);
	cpair temp2 = std::make_tuple(std::get<0>(temp) * -1.0, std::get<1>(temp) * -1.0, std::get<2>(temp) * -1.0);
	return temp2;
}

//! Returns the resultant of two vectors.
cpair PAHProcess::add_vector(cpair vec1, cpair vec2) const{
    //check if vector is unitary
	cpair vec1_adj = scale_vector(vec1);
	cpair vec2_adj = scale_vector(vec2);
	// calculate resultant
	cpair temp = std::make_tuple(std::get<0>(vec1_adj) + std::get<0>(vec2_adj), std::get<1>(vec1_adj) + std::get<1>(vec2_adj), std::get<2>(vec1_adj) + std::get<2>(vec2_adj));
	cpair temp2 = scale_vector(temp);
    return temp2;
}

//! Returns the normal vector to the face of a ring.
cpair PAHProcess::norm_vector(cpair p1, cpair p2, cpair p3) const{
	cpair vec1 = scale_vector(std::make_tuple(std::get<0>(p2) - std::get<0>(p1), std::get<1>(p2) - std::get<1>(p1), std::get<2>(p2) - std::get<2>(p1)));
	cpair vec2 = scale_vector(std::make_tuple(std::get<0>(p3) - std::get<0>(p2), std::get<1>(p3) - std::get<1>(p2), std::get<2>(p3) - std::get<2>(p2)));
	cpair temp = std::make_tuple(std::get<1>(vec1)*std::get<2>(vec2) - std::get<2>(vec1)*std::get<1>(vec2), std::get<2>(vec1)*std::get<0>(vec2) - std::get<0>(vec1)*std::get<2>(vec2), std::get<0>(vec1)*std::get<1>(vec2) - std::get<1>(vec1)*std::get<0>(vec2));
	double theta = acos(std::get<0>(vec1) * std::get<0>(vec2) + std::get<1>(vec1) * std::get<1>(vec2) + std::get<2>(vec1) * std::get<2>(vec2))*180.0/M_PI;
	if (theta < 180.0){
		temp = std::make_tuple(std::get<0>(temp) * -1.0, std::get<1>(temp) * -1.0, std::get<2>(temp) * -1.0);
	}
	if ( (theta <= 0.1 && theta >= -0.1) || (theta >= 179.9 && theta <= 180.1) ){
		cout << "Error in PAHProcess::norm_vector. Computed the cross product of almost parallel vectors. Returning arbitrary vector. \n";
		//Arbitrarily define the normal vector as 
		temp = std::make_tuple(std::get<0>(temp), std::get<1>(temp), std::get<2>(temp));
	}
	cpair temp2 = scale_vector(temp);
	return temp2;
}

//! Returns the angle in DEGREES between three points.
double PAHProcess::get_angle(cpair p1, cpair p2, cpair p3) const{
	cpair vec1 = scale_vector(std::make_tuple(std::get<0>(p2) - std::get<0>(p1), std::get<1>(p2) - std::get<1>(p1), std::get<2>(p2) - std::get<2>(p1)));
	cpair vec2 = scale_vector(std::make_tuple(std::get<0>(p3) - std::get<0>(p2), std::get<1>(p3) - std::get<1>(p2), std::get<2>(p3) - std::get<2>(p2)));
	cpair temp = std::make_tuple(std::get<1>(vec1)*std::get<2>(vec2) - std::get<2>(vec1)*std::get<1>(vec2), std::get<2>(vec1)*std::get<0>(vec2) - std::get<0>(vec1)*std::get<2>(vec2), std::get<0>(vec1)*std::get<1>(vec2) - std::get<1>(vec1)*std::get<0>(vec2));
	double theta = acos(std::get<0>(vec1) * std::get<0>(vec2) + std::get<1>(vec1) * std::get<1>(vec2) + std::get<2>(vec1) * std::get<2>(vec2))*180.0/M_PI;
	return theta;
}

//! Returns the cross product of two vectors. If v1 is the surface normal vector and v2 goes from C->newC v1xv2 redefines the next growth vector.
cpair PAHProcess::cross_vector (cpair vec1, cpair vec2) const{
	//check if vector is unitary
	cpair vec1_adj = scale_vector(vec1);
	cpair vec2_adj = scale_vector(vec2);
	cpair temp = std::make_tuple(	std::get<1>(vec1_adj) * std::get<2>(vec2_adj) - std::get<2>(vec1_adj) * std::get<1>(vec2_adj),
	 								std::get<2>(vec1_adj) * std::get<0>(vec2_adj) - std::get<0>(vec1_adj) * std::get<2>(vec2_adj), 
									std::get<0>(vec1_adj) * std::get<1>(vec2_adj) - std::get<1>(vec1_adj) * std::get<0>(vec2_adj));
	double magnitude = sqrt(std::get<0>(temp)*std::get<0>(temp) + std::get<1>(temp)*std::get<1>(temp) + std::get<2>(temp)*std::get<2>(temp));
	if (magnitude <= 1e-3){
		//Cross product of two parallel vectors. 
		cout << "Error in PAHProcess::cross_vector. Computed the cross product of almost parallel vectors. Returning arbitrary vector. \n"; //SETBREAKPOINT
		//Arbitrarily define the normal vector as 
		temp = std::make_tuple(std::get<0>(temp), std::get<1>(temp), std::get<2>(temp));
	}
	cpair temp2 = scale_vector(temp);
	return temp2;
}

//! Returns the dot product of two vectors. 
double PAHProcess::dot_vector (cpair vec1, cpair vec2) const{
	//check if vector is unitary
	cpair vec1_adj = scale_vector(vec1);
	cpair vec2_adj = scale_vector(vec2);
	double temp = std::get<0>(vec1_adj) * std::get<0>(vec2_adj) + std::get<1>(vec1_adj) * std::get<1>(vec2_adj) + std::get<2>(vec1_adj) * std::get<2>(vec2_adj);
	return temp;
}

//! Search a particular site (si) from svector associated with stype and erases it.
void PAHProcess::delSiteFromMap(const kmcSiteType& stype, const Spointer& si) {
    // search site vector for site 'si'
    for(int i=0; i<(int)m_pah->m_siteMap[stype].size(); i++) {
        // erase when found
        if(m_pah->m_siteMap[stype][i] == si) {
            m_pah->m_siteMap[stype].erase(m_pah->m_siteMap[stype].begin()+i);
            break;
        }
    }
}
//! Overload: search and erase from svectors associated with all site types in vector v
void PAHProcess::delSiteFromMap(const std::vector<kmcSiteType>& v, const Spointer& si) {
    for(int q=0; q!=(int)v.size(); q++) delSiteFromMap(v[q], si);
}
/*! Check for steric hindrance for a site growth process
  ! O -> positions to check
  ! FE                     ZZ & RFE               AC & RZZ                         
  !   O1--O2                   O1--O2                   O1--O2                 
  !   /     \                 /      \                 /      \                
  ! -C       O3             -C        O3--           -C        C-               
  !   \\    /                 \\     /                 \      /                 
  !    C --O4                  C -- C                   C----C                  
  !   /                       /      \                 /      \                 
  !                                                                              
*/
bool PAHProcess::checkHindrance(const Spointer& st) const {
    // position to check, start at C1 of site
    cpair mpos;
	cpair FEvec = get_vector(st->C1->coords, st->C2->coords);
	cpair ZZvec = get_vector(st->C1->C2->coords, st->C2->coords);
    switch(st->type) {
	case R5:
		mpos = jumpToPos(st->C1->coords, st->C1->bondAngle2, 0, 1.4*pow(3, 0.5)); // position O1
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2 - 60), 0, 1.4*pow(3, 0.5)); // position O2
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2 - 120), 0, 1.4*pow(3, 0.5)); // position O3
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2 - 180), 0, 1.4*pow(3, 0.5)); // position O4
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		break;
    case FE:
		// move position to O1, O2, O3, and O4
		mpos = jumpToPos(st->C1->coords, st->C1->growth_vector, 1.4); // position O1
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, st->C2->growth_vector, 1.4); // position O2
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, FEvec, 1.4); // position O3
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(st->C2->coords, st->C2->growth_vector, 1.4); // position O4
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		/*mpos = jumpToPos(st->C1->coords, st->C1->bondAngle2, 0, 1.4); // position O1
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2 - 60), 0, 1.4); // position O2
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2 - 120), 0, 1.4); // position O3
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2 - 180), 0, 1.4); // position O4
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		break;*/
        break;
    case RFE:
    case ZZ:
		// move position to O1, O2 and O3
		mpos = jumpToPos(st->C1->coords, st->C1->growth_vector, 1.4); // position O1
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, ZZvec, 1.4); // position O2
		if (m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
		mpos = jumpToPos(st->C2->coords, st->C2->growth_vector, 1.4); // position O3
		if (m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
		
		mpos = jumpToPos(st->C1->coords, normAngle(st->C1->bondAngle1 + 120), 0, 1.4*pow(3, 0.5)); // position O1
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle1 + 60), 0, 1.4*pow(3, 0.5)); // position O2
		if (m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
		/*mpos = jumpToPos(st->C1->coords, normAngle(st->C1->bondAngle1 + 120), 0, 1.4); // position O1
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle1 + 60), 0, 1.4); // position O2
		if (m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle1), 0, 1.4); // position O3
		if (m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
		
		mpos = jumpToPos(st->C1->coords, normAngle(st->C1->bondAngle1 + 120), 0, 1.4*pow(3, 0.5)); // position O1
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle1 + 60), 0, 1.4*pow(3, 0.5)); // position O2
		if (m_pah->m_cpositions.count(mpos)) return true; // check if position occupied*/
		break;
    case RZZ:
	case R5R6AC:
    case AC:
		// move position to O1 and O2
		mpos = jumpToPos(st->C1->coords, st->C1->bondAngle2, 0, 1.4); // position O1
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2 - 60), 0, 1.4); // position O2
		if (checkHindrance_C_PAH(mpos)) return true; // check if position occupied
		break;

        /*// move position to O1 and O2
        mpos = jumpToPos(st->C1->coords, st->C1->bondAngle2, 0, 1); // position O1
        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
		mpos = jumpToPos(st->C1->coords, st->C1->bondAngle2, 0, pow(3,0.5)); // position O1
		if (m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
        mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2-60), 0, 1); // position O2
		if (m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
		mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2 - 60), 0, pow(3,0.5)); // position O2
        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
        break;*/
    default:
        return true;
    }
    return false;
}
//! Checks if new C position is already occupied by edge carbons. 
bool PAHProcess::checkHindrance_C_PAH(cpair coords) const {
	double tol = 3e-1;
	for (std::set<cpair>::iterator it = m_pah->m_cpositions.begin(); it != m_pah->m_cpositions.end(); ++it) {
		double dist = getDistance_twoC(coords, *it);
		if (dist <= tol) return true;
	}
	return false;
}
//! Checks if new C position is already occupied by internal carbons. 
cpair PAHProcess::checkHindrance_C_intPAH(cpair coords) const {
	std::list<cpair>::iterator it1, it2;
	double minimal_dist = 1e3;
	for (it1 = m_pah->m_InternalCarbons.begin(); it1 != m_pah->m_InternalCarbons.end(); ++it1) {
		double dist_x = std::get<0>(coords) - std::get<0>(*it1);
		double dist_y = std::get<1>(coords) - std::get<1>(*it1);
		double dist_z = std::get<2>(coords) - std::get<2>(*it1);
		double dist = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
		if (dist < minimal_dist){
			minimal_dist = dist;
			it2 = it1;
		}
		//if (abs(std::get<0>(coords) -std::get<0>(*it1)) < tol && abs(std::get<1>(coords) -std::get<1>(*it1)) < tol && abs(std::get<2>(coords) -std::get<2>(*it1)) < tol){
			//Remove internal carbon from list
			//it2 = it1;
		//}
	}
	cpair temp = *it2;
	m_pah->m_InternalCarbons.erase(it2);
	return temp;
}
//! Checks if new C position is already occupied by edge carbons. 
bool PAHProcess::checkHindrance_newC(Cpointer C_1) const {
	double tol = 5e-1;
	cpair mpos = jumpToPos(C_1->coords, C_1->growth_vector, 1.4);
	for (std::set<cpair>::iterator it = m_pah->m_cpositions.begin(); it != m_pah->m_cpositions.end(); ++it) {
		double dist = getDistance_twoC(mpos, *it);
		if (dist <= tol) return true;
	}
	return false;
}

//! Checks if new C position is already occupied by edge carbons. 
bool PAHProcess::checkHindrance_newCposition(Cpointer C_1) const {
	double tol = 0.25;
	cpair mpos = jumpToPos(C_1->coords, C_1->growth_vector, 1.4);
	if (getDistance_twoCsquared(mpos, C_1->coords)<0.5) return false;
	for (std::set<cpair>::iterator it = m_pah->m_cpositions.begin(); it != m_pah->m_cpositions.end(); ++it) {
		double dist = getDistance_twoCsquared(mpos, *it);
		if (dist <= tol) return true;
	}
	if (HasCurvedMoeities()){
		for (std::list<cpair>::iterator it1 = m_pah->m_InternalCarbons.begin(); it1 != m_pah->m_InternalCarbons.end(); ++it1) {
			double dist = getDistance_twoCsquared(mpos, *it1);
			if (dist <= tol) return true;
		}
	}
	return false;
}

//! Checks if new C position is already occupied by edge carbons. 
bool PAHProcess::checkHindrance_newCposition(cpair mpos) const {
	double tol = 0.25;
	for (std::set<cpair>::iterator it = m_pah->m_cpositions.begin(); it != m_pah->m_cpositions.end(); ++it) {
		double dist = getDistance_twoCsquared(mpos, *it);
		if (dist <= tol) return true;
	}
	for (std::list<cpair>::iterator it1 = m_pah->m_InternalCarbons.begin(); it1 != m_pah->m_InternalCarbons.end(); ++it1) {
		double dist = getDistance_twoCsquared(mpos, *it1);
		if (dist <= tol) return true;
	}
	return false;
}

bool PAHProcess::checkHindrance_twoC(const Cpointer C_1, const Cpointer C_2) const {
	double tol = 5e-1;
	double dist = getDistance_twoC(C_1, C_2);
	if (dist<= tol)	return false;
	else return true;
}

double PAHProcess::getDistance_twoC(const Cpointer C_1, const Cpointer C_2) const {
	double xdist, ydist, zdist;
	xdist = std::get<0>(C_1->coords) - std::get<0>(C_2->coords);
	ydist = std::get<1>(C_1->coords) - std::get<1>(C_2->coords);
	zdist = std::get<2>(C_1->coords) - std::get<2>(C_2->coords);
	return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
}

double PAHProcess::getDistance_twoCsquared(const Cpointer C_1, const Cpointer C_2) const {
	double xdist, ydist, zdist;
	xdist = std::get<0>(C_1->coords) - std::get<0>(C_2->coords);
	ydist = std::get<1>(C_1->coords) - std::get<1>(C_2->coords);
	zdist = std::get<2>(C_1->coords) - std::get<2>(C_2->coords);
	return xdist*xdist + ydist*ydist + zdist*zdist;
}

double PAHProcess::getDistance_twoC(const cpair C_1, const cpair C_2) const {
	double xdist, ydist, zdist;
	xdist = std::get<0>(C_1) - std::get<0>(C_2);
	ydist = std::get<1>(C_1) - std::get<1>(C_2);
	zdist = std::get<2>(C_1) - std::get<2>(C_2);
	return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
}

double PAHProcess::getDistance_twoCsquared(const cpair C_1, const cpair C_2) const {
	double xdist, ydist, zdist;
	xdist = std::get<0>(C_1) - std::get<0>(C_2);
	ydist = std::get<1>(C_1) - std::get<1>(C_2);
	zdist = std::get<2>(C_1) - std::get<2>(C_2);
	return xdist*xdist + ydist*ydist + zdist*zdist;
}

/*! Check steric hindrance for phenyl addition reactions
! O -> positions to check
!
! -C        O2 -- O3
!   \       /      \
!    C -- O1       O4
!   /       \      /
! -C        O6 -- O5
*/
bool PAHProcess::checkHindrancePhenyl(const Cpointer C_1) const {
    // position to check, start at O1
	cpair mpos;
	// This function should compare new atomic coordinates to existing ones. 
	//Vector typespace
	cpair vec1 = get_vector(C_1->C1->coords, C_1->coords);
	cpair vec2 = get_vector(C_1->C2->coords, C_1->coords);
	cpair ivec = invert_vector(C_1->growth_vector);
	cpair ivec2 = invert_vector(vec2);
	mpos = jumpToPos(C_1->coords, C_1->growth_vector, 1.4); // position O1
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, vec2, 1.4); // position O2
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, C_1->growth_vector, 1.4); // position O3
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, vec1, 1.4); //position O4
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, ivec2, 1.4); // position O5
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, ivec, 1.4); // position O6
	if (checkHindrance_C_PAH(mpos)) return true;
	
	//Angle typespace
	/*mpos = jumpToPos(C_1->coords, C_1->bondAngle2, 0, 1.4); // position O1
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2 + 60), 0, 1.4); // position O2
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, C_1->bondAngle2, 0, 1); // position O3
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2 - 60), 0, 1.4); //position O4
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2 - 120), 0, 1.4); // position O5
	if (checkHindrance_C_PAH(mpos)) return true;
	mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2 - 180), 0, 1.4); // position O6
	if (checkHindrance_C_PAH(mpos)) return true;*/
	
    /*cpair mpos = jumpToPos(C_1->coords, C_1->bondAngle2, 0, 1); // position O1
    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
    mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2+60), 0, 1); // position O2
    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
    mpos = jumpToPos(mpos, C_1->bondAngle2, 0, 1); // position O3
    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
    mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2-60), 0, 1); //position O4
    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
    mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2-120), 0, 1); // position O5
    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
    mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2-180), 0, 1); // position O6
    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied*/
    return false;
}
//! Choose random site of a site type st, returns iterator to site
Spointer PAHProcess::chooseRandomSite(kmcSiteType st, rng_type &rng) {
    // to choose any principal site
    if(st == any) {
        // choose any site index from m_pah->m_siteList
        // There does not seem to be any handling of the case that there are no sites. FULLERENES!!

        // Set up an object to generate an integer uniformly distributed on [0, size - 1]
        typedef boost::uniform_smallint<unsigned int> site_index_distrib;
        site_index_distrib siteIndexDistrib(0,  m_pah->m_siteList.size()-1);
        boost::variate_generator<rng_type &, site_index_distrib> siteIndexGenerator(rng, siteIndexDistrib);

        // move iterator to site index and return iterator
        return moveIt(m_pah->m_siteList.begin(), siteIndexGenerator());
    } else if(st == benz) { //to choose sites for phenyl addition
        return chooseRandomSite(PHsites, rng);
	} else if(st == Methyl) {
		std::vector<Spointer> svector;
		for (Spointer site_it = m_pah->m_siteList.begin(); site_it != m_pah->m_siteList.end(); site_it++) {
			if (site_it->C1->A == 'M'){
				svector.push_back(site_it);
			}
		}
		int new_sz = (int)svector.size()-1;
		if (new_sz>0){
			// Set up an object to generate an uniform integer on [0, sz] now that we
			// know that sz > 0 (and can safely be cast to an unsigned type.
			typedef boost::uniform_smallint<unsigned int> site_index_distrib;
			site_index_distrib siteIndexDistrib(0, static_cast<unsigned int>(new_sz));
			boost::variate_generator<rng_type &, site_index_distrib> siteIndexGenerator(rng, siteIndexDistrib);

			const unsigned r = siteIndexGenerator();
			//cout<<"~~Chose "<<r<<"th site..\n";
			return svector[r];
		}
		return svector[0]; // if vector has only one element
	} else {
        //cout << "~~Choosing from " << m_pah->m_siteMap[st].size() << " sites...\n";
        // choose site index from site vector associated with site type st
        int sz = getSiteCount(st) - 1;
		if (sz == ((int)m_pah->m_siteMap[st].size()-1) ){
			// The number of sites coincides with all available sites
			//int sz = ((int) m_pah->m_siteMap[st].size())-1; // size - 1
			if(sz > 0) {
				// Set up an object to generate an uniform integer on [0, sz] now that we
				// know that sz > 0 (and can safely be cast to an unsigned type.
				typedef boost::uniform_smallint<unsigned int> site_index_distrib;
				site_index_distrib siteIndexDistrib(0, static_cast<unsigned int>(sz));
				boost::variate_generator<rng_type &, site_index_distrib> siteIndexGenerator(rng, siteIndexDistrib);

				const unsigned r = siteIndexGenerator();
				//cout<<"~~Chose "<<r<<"th site..\n";
				return m_pah->m_siteMap[st][r];
			}
			return m_pah->m_siteMap[st][0]; // if vector has only one element
		}
		else {
			//There is unavailable sites (Methyl or other type)
			//Create a copy of the site map for this site in a vector 
			std::vector<Spointer> svector = m_pah->m_siteMap[st];
			for(int i=(int)m_pah->m_siteMap[st].size()-1; i>=0; i--) {
				if(m_pah->m_siteMap[st][i]->C1->A == 'M' || m_pah->m_siteMap[st][i]->C2->A == 'M'){
					svector.erase(svector.begin()+i);
				}
			}
			int new_sz = (int)svector.size()-1;
			if (new_sz>0){
				// Set up an object to generate an uniform integer on [0, sz] now that we
				// know that sz > 0 (and can safely be cast to an unsigned type.
				typedef boost::uniform_smallint<unsigned int> site_index_distrib;
				site_index_distrib siteIndexDistrib(0, static_cast<unsigned int>(new_sz));
				boost::variate_generator<rng_type &, site_index_distrib> siteIndexGenerator(rng, siteIndexDistrib);

				const unsigned r = siteIndexGenerator();
				//cout<<"~~Chose "<<r<<"th site..\n";
				return svector[r];
			}
			return svector[0]; // if vector has only one element
		}
    }
}

//! Choose a random site of any site types in vtype
Spointer PAHProcess::chooseRandomSite(std::vector<kmcSiteType> vtype, rng_type &rng) {
    // stores number of sites for each site type
    std::vector<unsigned int> noOfSites;
    // stores total number of sites
    unsigned int sum_s=0;
    for(int i=0; i!= (int) vtype.size(); i++) {
        noOfSites.push_back(getSiteCount(vtype[i]));
        sum_s+=getSiteCount(vtype[i]);
    }
    int sz = sum_s-1; // size - 1
    unsigned r; // random number between 0 to sz
    if(sz > 0) {
        // Set up an object to generate an uniform integer on [0, sz] now that we
        // know that sz > 0 (and can safely be cast to an unsigned type.
        typedef boost::uniform_smallint<unsigned int> site_index_distrib;
        site_index_distrib siteIndexDistrib(0, static_cast<unsigned int>(sz));
        boost::variate_generator<rng_type &, site_index_distrib> siteIndexGenerator(rng, siteIndexDistrib);

        r = siteIndexGenerator();
    }
    else {
        r=0;
    }

    for(int i=0; i!=(int) noOfSites.size(); i++) {
        if(r < noOfSites[i]) return m_pah->m_siteMap[vtype[i]][r];
        else r -= noOfSites[i];
    }
    cout<<"ERROR: chooseRandomSite(vector): cannot choose proper site\n";
    return m_pah->m_siteMap[vtype[0]][0];
}

//! Creates a list of sites where an embedded R5 can migrate to.
std::list<Spointer> PAHProcess::listMigrationSites (Spointer& stt){
	std::list<Spointer> Migr_sites = {stt};
	if (stt->type != ACR5 && stt->type != R5R6 && stt->type != FEACR5) return Migr_sites;
	bool checking_site = true;
	Spointer sFE2 = stt;
	cpair R5coords;
	int s_type = (int)stt->type;
	bool R5R6_b4;
	bool check_left = true;
	bool check_right = true;
	
	//Get R5 internal coordinates
	if (s_type == ACR5) R5coords = findR5internal(stt->C1->C2, stt->C2->C1);
	else if (s_type == R5R6){
		if ( isR5internal(stt->C1->C1, stt->C1)) {
			R5coords = findR5internal(stt->C1->C1, stt->C1);
			R5R6_b4 = true;
		}
		else if ( isR5internal(stt->C2, stt->C2->C2)) {
			R5coords = findR5internal(stt->C2, stt->C2->C2);
			R5R6_b4 = false;
		}
		else {
			//R5 not found
			return Migr_sites;
		}
	}
	else{
		if ( isR5internal(stt->C1->C2, stt->C1->C2->C2)) {
			R5coords = findR5internal(stt->C1->C2, stt->C1->C2->C2);
			R5R6_b4 = true; check_right = false;
		}
		else if ( isR5internal(stt->C2->C1->C1, stt->C2->C1)) {
			R5coords = findR5internal(stt->C2->C1->C1, stt->C2->C1);
			R5R6_b4 = false; check_left = false;
		}
		else {
			//R5 not found
			return Migr_sites;
		}
	}
	
	//Fundamental assumption: R5-R7 pairs cannot move away from each other!
	if (m_pah->m_R7loc.size()>=1){
		for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
			double distR5R7 = getDistance_twoC(*it, R5coords);
			if (distR5R7 < 3.1) {
				m_pah->m_R5loc.push_back(R5coords);
				return Migr_sites;
			}
		}
	}
	
	//Add available sites to the left of stt.
	if (s_type == R5R6 && R5R6_b4) {
		//Check if coupled site can lead to migration
		sFE2 = moveIt(sFE2, -1);
		if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504)  checking_site = false; check_left = false; // This would violate the IPR.
		if ((int)sFE2->type >= 602 && (int)sFE2->type <= 604)  checking_site = false; check_left = false; // This would violate the IPR.
		if ((int)sFE2->type >= 1002 && (int)sFE2->type <= 1004)  checking_site = false; check_left = false; // This would violate the IPR.
		if ((int)sFE2->type >= 2000) checking_site = false; check_left = false;
	}
	if (check_left){
		do{
			sFE2 = moveIt(sFE2, -1);
			Spointer checkR5_1 = moveIt(sFE2,-1);
			Spointer checkR5_2 = moveIt(sFE2,-2);
			//Check for unsupported sites. This section heavily assumes that the Isolated Pentagon Rule is valid.
			if ((int)sFE2->type == 0 && (int)checkR5_1->type == 0 && (int)checkR5_2->type == 0) checking_site = false; // The result would be an indene, not supported. YET!
			if ((int)sFE2->type == 100 || (int)sFE2->type == 101 || (int)sFE2->type == 501 || (int)sFE2->type == 2002) checking_site = false; // This would violate the IPR.
			if ((int)sFE2->type == 9999 || (int)sFE2->type == -1 || sFE2->type == None) checking_site = false;
			if ((int)sFE2->type == 0){
				if ((int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 || (int)checkR5_1->type == 100) checking_site = false;
				if ((int)checkR5_1->type >= 1002 && (int)checkR5_1->type <= 1004) checking_site = false;
				if ((int)checkR5_1->type >= 2002 && (int)checkR5_1->type <= 2204) checking_site = false;
				if ((int)checkR5_1->type >= 2204 && (int)checkR5_1->type <= 2205) checking_site = false;
				if ((int)checkR5_1->type == 0){
					if ((int)checkR5_2->type == 101 || (int)checkR5_2->type == 501) checking_site = false;
					if ((int)checkR5_2->type >= 1002 && (int)checkR5_2->type <= 1004) checking_site = false;
					if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2205) checking_site = false;
				}
			}
			Cpointer CR5_otherside_end = sFE2->C2->C1;
			if (CR5_otherside_end->bridge) checking_site = false;
			if (CR5_otherside_end->C2->bridge) checking_site = false; 
			if (CR5_otherside_end->C1->bridge) checking_site = false;
			if (CR5_otherside_end->C2->C2->bridge) checking_site = false; 
			if (CR5_otherside_end->C1->C1->bridge) checking_site = false;
			
			//Check for other side being valid
			Cpointer thirdC_after = findThirdC(CR5_otherside_end);
			if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
				Spointer opp_site_after = findSite(thirdC_after);
				//I have no clue how to flag an Spointer as error. This can cause seg faults.
				if (opp_site_after != m_pah->m_siteList.end()){
					int os_endtype = opp_site_after->type;
					if (os_endtype >= 200 && os_endtype <= 203) checking_site = false;
					if (os_endtype == 101) checking_site = false;
					if (os_endtype >= 600 && os_endtype <= 603) checking_site = false;
					if (os_endtype >= 1000 && os_endtype <= 1003) checking_site = false;
					if (os_endtype >= 500 && os_endtype <= 504) checking_site = false;
					if (os_endtype >= 2000 && os_endtype <= 2205) checking_site = false;
					if (os_endtype >= 2103 && os_endtype <= 2105) checking_site = false;
					if (os_endtype >= 2204 && os_endtype <= 2205) checking_site = false;
					if (os_endtype == 9999 || os_endtype == -1 || opp_site_after->type == None) checking_site = false;
				}
			}
			
			//check that two pentagons (including internals) will not collide
			cpair R5coords_end = endposR5internal(CR5_otherside_end, CR5_otherside_end->C2);
			if (m_pah->m_R5loc.size()>=1){
				for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
					double distR5s = getDistance_twoC(*it, R5coords_end);
					if (distR5s < 2.8) {
						//This distance is a parameter of this jump process. Might need some more tuning. 
						//2.8 seems appropiate but may reject too many jumps.
						//Two pentagons will be next to each other violating the Isolated Pentagon Rule
						checking_site = false;
						break;
					}
				}
			}
			
			//check that pentagon and heptagon (including internals) will not collide
			if (m_pah->m_R7loc.size()>=1){
				for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
					double distR5R7 = getDistance_twoC(*it, R5coords_end);
					if (distR5R7 < 2.6) {
						//This distance is a parameter of this jump process. Might need some more tuning. 
						//2.8 seems appropiate but may reject too many jumps.
						//Two pentagons will be next to each other violating the Isolated Pentagon Rule
						checking_site = false;
						break;
					}
				}
			}
			if (sFE2==Migr_sites.back() && Migr_sites.size()>1) {
				checking_site = false; check_left=false; check_right=false;
			}
			if ((int)sFE2->type % 10 >=2 && (Migr_sites.front()->type == FE || Migr_sites.front()->type == R5R6)) checking_site = false;
			if (checking_site) Migr_sites.push_front(sFE2);
			if ((int)sFE2->type % 10 >=2) checking_site = false;
			if (sFE2->comb == FE2) checking_site = false;
		} while (checking_site == true);
	}
	//Add available sites to the right of stt
	sFE2 = stt;
	checking_site = true;
	//Add available sites to the left of stt.
	if (s_type == R5R6 && !R5R6_b4) {
		//Check if coupled site can lead to migration
		sFE2 = moveIt(sFE2, +1);
		if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504)  checking_site = false; check_right = false; // This would violate the IPR.
		if ((int)sFE2->type >= 602 && (int)sFE2->type <= 604)  checking_site = false; check_right = false; // This would violate the IPR.
		if ((int)sFE2->type >= 1002 && (int)sFE2->type <= 1004)  checking_site = false; check_right = false; // This would violate the IPR.
		if ((int)sFE2->type >= 2000) checking_site = false; check_right = false;
	}
	if (check_right){
		do{
			sFE2 = moveIt(sFE2, +1);
			Spointer checkR5_1 = moveIt(sFE2,+1);
			Spointer checkR5_2 = moveIt(sFE2,+2);
			//Check for unsupported sites. This section heavily assumes that the Isolated Pentagon Rule is valid.
			if ((int)sFE2->type == 0 && (int)checkR5_1->type == 0 && (int)checkR5_2->type == 0) checking_site = false; // The result would be an indene, not supported. YET!
			if ((int)sFE2->type == 100 || (int)sFE2->type == 101 || (int)sFE2->type == 501 || (int)sFE2->type == 2002) checking_site = false; // This would violate the IPR.
			if ((int)sFE2->type == 9999 || (int)sFE2->type == -1 || sFE2->type == None) checking_site = false;
			if ((int)sFE2->type == 0){
				if ((int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 || (int)checkR5_1->type == 100) checking_site = false;
				if ((int)checkR5_1->type >= 1002 && (int)checkR5_1->type <= 1004) checking_site = false;
				if ((int)checkR5_1->type >= 2002 && (int)checkR5_1->type <= 2204) checking_site = false;
				if ((int)checkR5_1->type >= 2204 && (int)checkR5_1->type <= 2205) checking_site = false;
				if ((int)checkR5_1->type == 0){
					if ((int)checkR5_2->type == 101 || (int)checkR5_2->type == 501) checking_site = false;
					if ((int)checkR5_2->type >= 1002 && (int)checkR5_2->type <= 1004) checking_site = false;
					if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2205) checking_site = false;
				}
			}
			Cpointer CR5_otherside_end = sFE2->C1->C2;
			if (CR5_otherside_end->bridge) checking_site = false;
			if (CR5_otherside_end->C2->bridge) checking_site = false; 
			if (CR5_otherside_end->C1->bridge) checking_site = false;
			if (CR5_otherside_end->C2->C2->bridge) checking_site = false; 
			if (CR5_otherside_end->C1->C1->bridge) checking_site = false;
			
			//Check for other side being valid
			Cpointer thirdC_after = findThirdC(CR5_otherside_end);
			if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
				Spointer opp_site_after = findSite(thirdC_after);
				//I have no clue how to flag an Spointer as error. This can cause seg faults.
				if (opp_site_after != m_pah->m_siteList.end()){
					int os_endtype = opp_site_after->type;
					if (os_endtype >= 200 && os_endtype <= 203) checking_site = false;
					if (os_endtype == 101) checking_site = false;
					if (os_endtype >= 600 && os_endtype <= 603) checking_site = false;
					if (os_endtype >= 1000 && os_endtype <= 1003) checking_site = false;
					if (os_endtype >= 500 && os_endtype <= 504) checking_site = false;
					if (os_endtype >= 2000 && os_endtype <= 2205) checking_site = false;
					if (os_endtype >= 2103 && os_endtype <= 2105) checking_site = false;
					if (os_endtype >= 2204 && os_endtype <= 2205) checking_site = false;
					if (os_endtype == 9999 || os_endtype == -1 || opp_site_after->type == None) checking_site = false;
				}
			}
			
			//check that two pentagons (including internals) will not collide
			cpair R5coords_end = endposR5internal(CR5_otherside_end->C1, CR5_otherside_end,true);
			if (m_pah->m_R5loc.size()>=1){
				for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
					double distR5s = getDistance_twoC(*it, R5coords_end);
					if (distR5s < 2.8) {
						//This distance is a parameter of this jump process. Might need some more tuning. 
						//2.8 seems appropiate but may reject too many jumps.
						//Two pentagons will be next to each other violating the Isolated Pentagon Rule
						checking_site = false;
						break;
					}
				}
			}
			//check that pentagon and heptagon (including internals) will not collide
			if (m_pah->m_R7loc.size()>=1){
				for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
					double distR5R7 = getDistance_twoC(*it, R5coords_end);
					if (distR5R7 < 2.6) {
						//This distance is a parameter of this jump process. Might need some more tuning. 
						//2.8 seems appropiate but may reject too many jumps.
						//Two pentagons will be next to each other violating the Isolated Pentagon Rule
						checking_site = false;
						break;
					}
				}
			}
			if (sFE2 == Migr_sites.front() && Migr_sites.size()>1) {
				checking_site = false; check_left=false; check_right=false;
			}
			if ((int)sFE2->type % 10 >=2 && (Migr_sites.back()->type == FE || Migr_sites.back()->type == R5R6)) checking_site = false;
			if (checking_site) Migr_sites.push_back(sFE2);
			if ((int)sFE2->type % 10 >=2) checking_site = false;
			if (sFE2->comb == FE2) checking_site = false;
		} while (checking_site == true);
	}
	m_pah->m_R5loc.push_back(R5coords);
	return Migr_sites;
}

//! Assigns probabilities to individual sites and randomly selects a finishing site for a migration process.
std::tuple<Spointer,bool,int> PAHProcess::chooseRandomMigrationSite (Spointer& start_site, std::map<std::string,double> rates, std::list<Spointer> avail_end_sites, int N_end_steps, rng_type &rng) {
	int N_end_sites = avail_end_sites.size();
	std::tuple<Spointer,bool,int> Site_direction;
	if (N_end_sites == 1) {
		Site_direction = std::make_tuple(start_site, false, 0);
		return Site_direction;
	}
	std::vector<int> position(N_end_sites,0);
	std::vector<double> left_sites_probs(N_end_sites,0.0);
	std::vector<double> right_sites_probs(N_end_sites,0.0);
	std::vector<double> end_sites_probs(N_end_sites,0.0);
	//Each site has three transition probabilites. To solve the problem we have to calculate them first.
	size_t ii = 0; size_t jj = 0;
	for (std::list<Spointer>::iterator it = avail_end_sites.begin(); it != avail_end_sites.end(); ++it) {
		if ( (*it)->comb == FE2) {
			left_sites_probs[ii] = 0.0;
			right_sites_probs[ii] = 0.0;
			end_sites_probs[ii] = rates["R5"];
		}
		else if ( (*it)->type == FE || (*it)->type == R5R6) {
			left_sites_probs[ii] = rates["Corner"];
			right_sites_probs[ii] = rates["Corner"];
			end_sites_probs[ii] = rates["R5R6"];
		}
		else if ( (*it)->type == ZZ || (*it)->type == ACR5) {
			left_sites_probs[ii] = rates["Middle"];
			right_sites_probs[ii] = rates["Middle"];
			end_sites_probs[ii] = rates["ACR5"];
		}
		else if ( (*it)->type == AC || (*it)->type == FEACR5) {
			std::list<Spointer>::iterator it2 = it;
			it2 ++;
			if (it2 == avail_end_sites.end()) it2--; it2--;
			if ((*it2)->type == FE){
				left_sites_probs[ii] = rates["Corner"];
				right_sites_probs[ii] = rates["Corner"];
				end_sites_probs[ii] = rates["R5R6ZZ"];
			}else{
				left_sites_probs[ii] = rates["Middle"];
				right_sites_probs[ii] = rates["Middle"];
				end_sites_probs[ii] = rates["FEACR5"];
			}
		}
		else if ( (*it)->type == BY5 || (*it)->type == ZZACR5) {
			std::list<Spointer>::iterator it2 = it;
			it2 ++;
			if (it2 == avail_end_sites.end()) it2--; it2--;
			if ((*it2)->type == FE){
				left_sites_probs[ii] = rates["Corner"];
				right_sites_probs[ii] = rates["Corner"];
				end_sites_probs[ii] = rates["R5R6AC"];
			}else{
				left_sites_probs[ii] = rates["Corner"];
				right_sites_probs[ii] = rates["Corner"];
				end_sites_probs[ii] = rates["ZZACR5"];
			}
		}
		else if ( (*it)->type == BY6 || (*it)->type == ACACR5) {
			std::list<Spointer>::iterator it2 = it;
			it2 ++;
			if (it2 == avail_end_sites.end()) it2--; it2--;
			if ((*it2)->type == FE){
				left_sites_probs[ii] = rates["Corner"];
				right_sites_probs[ii] = rates["Corner"];
				end_sites_probs[ii] = rates["R5R6BY5"];
			}else{
				left_sites_probs[ii] = rates["Corner"];
				right_sites_probs[ii] = rates["Corner"];
				end_sites_probs[ii] = rates["ACACR5"];
			}
		}
		else if ( (*it)->type == RZZ) {
			left_sites_probs[ii] = rates["Corner"];
			right_sites_probs[ii] = rates["Corner"];
			end_sites_probs[ii] = rates["R5R6FER"];
		}
		else if ( (*it)->type == RAC) {
			left_sites_probs[ii] = rates["Corner"];
			right_sites_probs[ii] = rates["Corner"];
			end_sites_probs[ii] = rates["R5R6ZZR"];
		}
		else if ( (*it)->type == RBY5) {
			left_sites_probs[ii] = rates["Corner"];
			right_sites_probs[ii] = rates["Corner"];
			end_sites_probs[ii] = rates["R5R6ACR"];
		}
		else if ( (*it)->type == R5R6ZZ) {
			left_sites_probs[ii] = rates["Middle"];
			right_sites_probs[ii] = rates["Middle"];
			end_sites_probs[ii] = rates["R5R6FER5R6"];
		}
		else if ( (*it)->type == R5R6AC) {
			left_sites_probs[ii] = rates["Corner"];
			right_sites_probs[ii] = rates["Corner"];
			end_sites_probs[ii] = rates["R5R6ZZR5R6"];
		}
		else if ( (*it)->type == R5R6BY5) {
			left_sites_probs[ii] = rates["Corner"];
			right_sites_probs[ii] = rates["Corner"];
			end_sites_probs[ii] = rates["R5R6ACR5R6"];
		}
		if ( (*it) == start_site) jj = ii;
		if ( (*it) == avail_end_sites.front()) left_sites_probs[ii] = 0.0;
		if ( (*it) == avail_end_sites.back()) right_sites_probs[ii] = 0.0;
		ii += 1;
	}

	//Calculate the parameter for the random walk length
	double weighted_sum_numerator = 0.0;
	double weighted_sum_denominator = 0.0;
	for (unsigned int kk = 0; kk!=N_end_sites;kk++){
		weighted_sum_numerator += end_sites_probs[kk];
		weighted_sum_denominator += (left_sites_probs[kk] + right_sites_probs[kk])/2.0;
	}
	weighted_sum_numerator = weighted_sum_numerator / (float)N_end_sites;
	weighted_sum_denominator = weighted_sum_denominator / (float)N_end_sites;
	double lambda_param = weighted_sum_numerator / weighted_sum_denominator;

	//Generate exponentially distributed Number of migration steps
	boost::exponential_distribution<double> N_stepsDistrib(lambda_param);
	boost::variate_generator<rng_type &, boost::exponential_distribution<double>> N_stepsGenerator(rng, N_stepsDistrib);
	int N_steps_end_internal = (int)N_stepsGenerator();

	//Perform the random walk.
	int N_steps = 0;
	bool migrating_site = true;
	do{
		std::vector<double> transition_rates{left_sites_probs[jj], right_sites_probs[jj], end_sites_probs[jj]};
		// chooses index from a vector of weights (double number in this case) randomly
		boost::uniform_01<rng_type &, double> uniformGenerator(rng);
		size_t ind = chooseIndex<double>(transition_rates, uniformGenerator);
		if (ind == 0) jj --;
		else if (ind == 1) jj++;
		else migrating_site = false;
		N_steps++;
		if (N_steps >= N_steps_end_internal) migrating_site = false;
		if (N_steps >= N_end_steps) migrating_site = false;
	} while(migrating_site);
    
	//Return pointer to selected site
	Spointer chosenSite;
	ii = 0;
	bool direction = true;
	for (std::list<Spointer>::iterator it = avail_end_sites.begin(); it != avail_end_sites.end(); ++it) {
		if (ii == jj) {
			chosenSite = (*it);
			break;
		}
		if ( (*it) == start_site) direction = false;
		ii += 1;
	}
	Site_direction = std::make_tuple(chosenSite, direction, N_steps);
	return Site_direction;
}

//! Check if site sc has a combined site type cstype
/*bool PAHProcess::checkCombSites(const Spointer& sc, const kmcSiteType& cstype) const {
    for(int i=0; i!= (int) sc->comb.size(); i++) {
        if(sc->comb[i]==cstype) {
            return true;
        }
    }
    return false;
}*/

// Write Processes
//! Creates a lone carbon atom at position 0,0
Cpointer PAHProcess::addC() {
    Cpointer cb;
    // Create new carbon atom at memory pointed by h
    cb = new Carbon;
    if(!m_pah->m_carbonList.insert(cb).second)
        std::cout<<"ERROR: ADDING SAME CARBON POINTER TO SET\n";
    m_pah->m_cpositions.insert(cb->coords); // store coordinates
    addCount(1,0); // add a C count
    return cb;
}
/*! Create a new carbon atom attached next to C1
  ! C1
  !    \ <----------- angle1
  !      C __ (C1->)C2
  !        angle 2
*/
Cpointer PAHProcess::addC(Cpointer C_1, angletype angle1, angletype angle2, bondlength length, bool bulk) {
    Cpointer cb;
    // Create new carbon atom
    cb = new Carbon;
    // Set details of new Carbon
    cb->C1 = C_1;
    cb->C2 = C_1->C2;
    cb->bondAngle1 = normAngle(angle2); // convert angle to +ve/-ve form
    
	// set new coordinates and store
	cpair mpos = jumpToPos(C_1->coords, angle1, 0, length);
	if (bulk){// Carbon is coming from bulk, removing from Internal carbons
		mpos = checkHindrance_C_intPAH(mpos);
	}
	cb->coords = mpos;
    if(!m_pah->m_carbonList.insert(cb).second)
        std::cout<<"ERROR: ADDING SAME CARBON POINTER TO SET\n";
    m_pah->m_cpositions.insert(cb->coords);
    // Edit details of connected carbon(s)
    if(C_1->C2 != NULL) {
        // change member pointer of original neighbour of C_1
        C_1->C2->C1 = cb;
    }
    C_1->bondAngle1 = normAngle(angle1);
    C_1->C2 = cb;
    if(!bulk) addCount(1,0);
    if(C_1 == m_pah->m_clast) m_pah->m_clast = cb;
	
	//Any Carbon addition means that the structure will not have an optimised geometry.
	m_pah->m_optimised = false;
    return cb;
}

/*! Create a new carbon atom attached next to C1 using a vector to tell the position of the next carbon atom.
  ! C1
  !    \ <----------- vector
  !      C __ (C1->)C2
  !        angle 2
*/
Cpointer PAHProcess::addC(Cpointer C_1, cpair direction, bondlength length, bool bulk) {
    Cpointer cb;
    // Create new carbon atom
    cb = new Carbon;
    // Set details of new Carbon
    cb->C1 = C_1;
    cb->C2 = C_1->C2;
    
	// set new coordinates and store
	cpair mpos = jumpToPos(C_1->coords, direction, length);
	if (bulk){// Carbon is coming from bulk, removing from Internal carbons
		mpos = checkHindrance_C_intPAH(mpos);
		double newdist = getDistance_twoC(C_1->coords,mpos);
		if (newdist>1.7){
			std::cout << "WARNING. Returning carbon from internal carbons with distance " << newdist << ".\n";
			std::cout << "Adding from C coords (" << std::get<0>(C_1->coords) << "," << std::get<1>(C_1->coords) << "," << std::get<2>(C_1->coords) << ")\n";
			std::cout << "New carbon coords    (" << std::get<0>(mpos) << "," << std::get<1>(mpos) << "," << std::get<2>(mpos) << ")\n";
		}
	}
	cb->coords = mpos;
    if(!m_pah->m_carbonList.insert(cb).second)
        std::cout<<"ERROR: ADDING SAME CARBON POINTER TO SET\n";
    m_pah->m_cpositions.insert(cb->coords);
	cb->A = 'C';
    // Edit details of connected carbon(s)
    if(C_1->C2 != NULL) {
        // change member pointer of original neighbour of C_1
        C_1->C2->C1 = cb;
    }
    C_1->C2 = cb;
    if(!bulk) addCount(1,0);
    if(C_1 == m_pah->m_clast) m_pah->m_clast = cb;
	//Any Carbon addition means that the structure will not have an optimised geometry.
	m_pah->m_optimised = false;
    return cb;
}

//Overload routine. Passes a OpenBabel molecule so it is modified. 
Cpointer PAHProcess::addC(Cpointer C_1, cpair direction, bondlength length, OpenBabel::OBMol mol, bool bulk) {
    Cpointer cb;
    // Create new carbon atom
    cb = new Carbon;
    // Set details of new Carbon
    cb->C1 = C_1;
    cb->C2 = C_1->C2;
    
	// set new coordinates and store
	cpair mpos = jumpToPos(C_1->coords, direction, length);
	if (bulk){// Carbon is coming from bulk, removing from Internal carbons
		mpos = checkHindrance_C_intPAH(mpos);
	}
	cb->coords = mpos;
    if(!m_pah->m_carbonList.insert(cb).second)
        std::cout<<"ERROR: ADDING SAME CARBON POINTER TO SET\n";
    m_pah->m_cpositions.insert(cb->coords);
    // Edit details of connected carbon(s)
    if(C_1->C2 != NULL) {
        // change member pointer of original neighbour of C_1
        C_1->C2->C1 = cb;
    }
    C_1->C2 = cb;
    if(!bulk) addCount(1,0);
    if(C_1 == m_pah->m_clast) m_pah->m_clast = cb;
	//Add carbon to OpenBabel structure
	OpenBabel::OBAtom *atom  = mol.NewAtom();
	atom->SetAtomicNum(6);
	atom->SetVector(std::get<0>(cb->coords),std::get<1>(cb->coords),std::get<2>(cb->coords));
	atom->SetAromatic();
	atom->SetImplicitValence(3);
	
	double tol = 2e-1;
	for(OpenBabel::OBMolAtomIter     a(mol); a; ++a) {
		double x_pos = a->GetX();
		double y_pos = a->GetY();
		double z_pos = a->GetZ();
		cpair temp = std::make_tuple(x_pos, y_pos, z_pos);
		if (getDistance_twoC(C_1->coords, temp) < tol){
			// a is the starting carbon
			mol.AddBond(a->GetIdx(), atom->GetIdx(),5);
			break;
		}
	}
	//Any Carbon addition means that the structure will not have an optimised geometry.
	m_pah->m_optimised = false;
    return cb;
}

//! Adds bond to OpenBabel structure between two carbons. 
void PAHProcess::addOBbond(Cpointer C_1, Cpointer C_2, OpenBabel::OBMol mol){
	double tol = 2e-1;
	OpenBabel::OBMolAtomIter c1(mol);
	OpenBabel::OBMolAtomIter c2(mol);
	bool r1 = false;
	bool r2 = false;
	for(OpenBabel::OBMolAtomIter     a(mol); a; ++a) {
		double x_pos = a->GetX();
		double y_pos = a->GetY();
		double z_pos = a->GetZ();
		cpair temp = std::make_tuple(x_pos, y_pos, z_pos);
		if (getDistance_twoC(C_1->coords, temp) < tol){
			c1 = a;
			r1 = true;
		}
		if (getDistance_twoC(C_2->coords, temp) < tol){
			c2 = a;
			r2 = true;
		}
		if (r1&&r2) break;
	}
	mol.AddBond(c1->GetIdx(), c2->GetIdx(), 5);
}

void PAHProcess::moveC(Cpointer C_1, cpair newpos) {
	// Remove coordinates of C from m_pah->m_cpositions
	double min_dist = 1e3;
	std::set<cpair>::iterator it2;
	for (std::set<cpair>::iterator it = m_pah->m_cpositions.begin(); it != m_pah->m_cpositions.end(); ++it) {
		double x_dist = std::get<0>(C_1->coords) - std::get<0>(*it);
		double y_dist = std::get<1>(C_1->coords) - std::get<1>(*it);
		double z_dist = std::get<2>(C_1->coords) - std::get<2>(*it);
		double dist = sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);
		if (dist < min_dist) {
			it2 = it;
			min_dist = dist;
		}
		if (dist < 1e-2) break;
	}
	m_pah->m_cpositions.erase(*it2);
	//m_pah->m_cpositions.erase(m_pah->m_cpositions.find(C_1->coords));
	// Set new positions of Carbon C_1
	cpair temp = std::make_tuple(std::get<0>(newpos), std::get<1>(newpos), std::get<2>(newpos));
	/*C_1->coords.first = newpos.first;
	C_1->coords.second = newpos.second;*/
	C_1->coords = temp;
	m_pah->m_cpositions.insert(C_1->coords);
}

void PAHProcess::moveC(Cpointer C_1, Cpointer Cprev, double new_distance) {
	// Remove coordinates of C from m_pah->m_cpositions
	double min_dist = 1e3;
	std::set<cpair>::iterator it2;
	for (std::set<cpair>::iterator it = m_pah->m_cpositions.begin(); it != m_pah->m_cpositions.end(); ++it) {
		double x_dist = std::get<0>(C_1->coords) - std::get<0>(*it);
		double y_dist = std::get<1>(C_1->coords) - std::get<1>(*it);
		double z_dist = std::get<2>(C_1->coords) - std::get<2>(*it);
		double dist = sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);
		if (dist < min_dist) {
			it2 = it;
			min_dist = dist;
		}
		if (dist < 1e-2) break;
	}
	m_pah->m_cpositions.erase(*it2);
	//m_pah->m_cpositions.erase(m_pah->m_cpositions.find(C_1->coords));
	//Calculate distance between two carbons
	double dist = pow(pow(std::get<0>(C_1->coords) - std::get<0>(Cprev->coords), 2) + pow(std::get<1>(C_1->coords) - std::get<1>(Cprev->coords), 2) + pow(std::get<2>(C_1->coords) - std::get<2>(Cprev->coords), 2), 0.5);
	//Calculate new positions
	cpair temp = std::make_tuple(std::get<0>(Cprev->coords) + (std::get<0>(C_1->coords) - std::get<0>(Cprev->coords)) / dist*new_distance, std::get<1>(Cprev->coords) + (std::get<1>(C_1->coords) - std::get<1>(Cprev->coords)) / dist*new_distance, std::get<2>(Cprev->coords) + (std::get<2>(C_1->coords) - std::get<2>(Cprev->coords)) / dist*new_distance);
	/*C_1->coords.first = newpos.first;
	C_1->coords.second = newpos.second;*/
	C_1->coords = temp;
	m_pah->m_cpositions.insert(C_1->coords);
}
//! Move a carbon in z direction.
void PAHProcess::moveC_z(Cpointer C_1, double z_distance) {
	// Remove coordinates of C from m_pah->m_cpositions
	double min_dist = 1e3;
	std::set<cpair>::iterator it2;
	for (std::set<cpair>::iterator it = m_pah->m_cpositions.begin(); it != m_pah->m_cpositions.end(); ++it) {
		double x_dist = std::get<0>(C_1->coords) - std::get<0>(*it);
		double y_dist = std::get<1>(C_1->coords) - std::get<1>(*it);
		double z_dist = std::get<2>(C_1->coords) - std::get<2>(*it);
		double dist = sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);
		if (dist < min_dist) {
			it2 = it;
			min_dist = dist;
		}
		if (dist < 1e-2) break;
	}
	m_pah->m_cpositions.erase(*it2);
	//m_pah->m_cpositions.erase(m_pah->m_cpositions.find(C_1->coords));
	// Set new positions of Carbon C_1
	cpair temp = std::make_tuple(std::get<0>(C_1->coords), std::get<1>(C_1->coords), std::get<2>(C_1->coords) + z_distance);/*C_1->coords.first = newpos.first;
	C_1->coords.second = newpos.second;*/
	C_1->coords = temp;
	m_pah->m_cpositions.insert(C_1->coords);
}
//! Adds an R5 to the list of R5s and R7s
void PAHProcess::addR5internal(Cpointer C_1, Cpointer C_2, bool invert_dir) {
	cpair mpos = endposR5internal(C_1, C_2, invert_dir);
	m_pah->m_R5loc.push_back(mpos);
}
//! Removes an R5 from the list of R5s and R7s
void PAHProcess::removeR5internal(Cpointer C_1, Cpointer C_2) {
	std::list<cpair>::iterator it1, it2;
	double minimal_dist = 1e3;
	for (it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		double dist_x = std::get<0>(C_1->coords) + std::get<0>(C_2->coords) - 2*std::get<0>(*it1);
		double dist_y = std::get<1>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<1>(*it1);
		double dist_z = std::get<2>(C_1->coords) + std::get<2>(C_2->coords) - 2*std::get<2>(*it1);
		double dist = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
		if (dist < minimal_dist){
			minimal_dist = dist;
			it2 = it1;
		}
	}
	m_pah->m_R5loc.erase(it2);
}
//! Removes an R5 from the list of R5s and R7s
void PAHProcess::removeR7internal(Cpointer C_1, Cpointer C_2) {
	std::list<cpair>::iterator it1, it2;
	double minimal_dist = 1e3;
	for (it1 = m_pah->m_R7loc.begin(); it1 != m_pah->m_R7loc.end(); ++it1) {
		double dist_x = std::get<0>(C_1->coords) + std::get<0>(C_2->coords) - 2*std::get<0>(*it1);
		double dist_y = std::get<1>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<1>(*it1);
		double dist_z = std::get<2>(C_1->coords) + std::get<2>(C_2->coords) - 2*std::get<2>(*it1);
		double dist = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
		if (dist < minimal_dist){
			minimal_dist = dist;
			it2 = it1;
		}
	}
	m_pah->m_R7loc.erase(it2);
}
//! Return internal R5 associated to two carbons and deletes it from R5 list.
cpair PAHProcess::findR5internal(Cpointer C_1, Cpointer C_2) {
	std::list<cpair>::iterator it1, it2;
	double minimal_dist = 1e3;
	for (it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		double dist_x = std::get<0>(C_1->coords) + std::get<0>(C_2->coords) - 2*std::get<0>(*it1);
		double dist_y = std::get<1>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<1>(*it1);
		double dist_z = std::get<2>(C_1->coords) + std::get<2>(C_2->coords) - 2*std::get<2>(*it1);
		double dist = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
		if (dist < minimal_dist){
			minimal_dist = dist;
			it2 = it1;
		}
	}
	cpair temp = *it2;
	m_pah->m_R5loc.erase(it2);
	return temp;
}

//! Return internal R5 associated to two carbons and deletes it from R5 list.
cpair PAHProcess::findR7internal(Cpointer C_1, Cpointer C_2) {
	std::list<cpair>::iterator it1, it2;
	double minimal_dist = 1e3;
	for (it1 = m_pah->m_R7loc.begin(); it1 != m_pah->m_R7loc.end(); ++it1) {
		double dist_x = std::get<0>(C_1->coords) + std::get<0>(C_2->coords) - 2*std::get<0>(*it1);
		double dist_y = std::get<1>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<1>(*it1);
		double dist_z = std::get<2>(C_1->coords) + std::get<2>(C_2->coords) - 2*std::get<2>(*it1);
		double dist = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
		if (dist < minimal_dist){
			minimal_dist = dist;
			it2 = it1;
		}
	}
	cpair temp = *it2;
	m_pah->m_R7loc.erase(it2);
	return temp;
}

//! Returns distance from point C_1 to the line formed between points C_2 and C_3
double PAHProcess::getDistance_point_to_line(const cpair C_1, const cpair C_2, const cpair C_3) const {
	cpair vec1 = std::make_tuple(std::get<0>(C_1) - std::get<0>(C_2), std::get<1>(C_1) - std::get<1>(C_2), std::get<2>(C_1) - std::get<2>(C_2));
	cpair vec2 = std::make_tuple(std::get<0>(C_1) - std::get<0>(C_3), std::get<1>(C_1) - std::get<1>(C_3), std::get<2>(C_1) - std::get<2>(C_3));
	cpair vec3 = std::make_tuple(std::get<0>(C_3) - std::get<0>(C_2), std::get<1>(C_3) - std::get<1>(C_2), std::get<2>(C_3) - std::get<2>(C_2));
	cpair temp = std::make_tuple(	std::get<1>(vec1) * std::get<2>(vec2) - std::get<2>(vec1) * std::get<1>(vec2),
	 								std::get<2>(vec1) * std::get<0>(vec2) - std::get<0>(vec1) * std::get<2>(vec2), 
									std::get<0>(vec1) * std::get<1>(vec2) - std::get<1>(vec1) * std::get<0>(vec2));
	double num = std::sqrt(std::get<0>(temp) * std::get<0>(temp) + std::get<1>(temp) * std::get<1>(temp) + std::get<2>(temp) * std::get<2>(temp));
	double den = std::sqrt(std::get<0>(vec3) * std::get<0>(vec3) + std::get<1>(vec3) * std::get<1>(vec3) + std::get<2>(vec3) * std::get<2>(vec3));
	double dist = num / den;
	return dist;
}

//! Are the two carbon atoms members of an R5 with coordinates in R5Internal??
bool PAHProcess::isR5internal(Cpointer C_1, Cpointer C_2) {
	cpair R5vec = get_vector(C_1->coords,C_2->coords);
	//Plane parameters for C_1
	double a_1 = std::get<0>(R5vec);
	double b_1 = std::get<1>(R5vec);
	double c_1 = std::get<2>(R5vec);
	double d_1 = std::get<0>(R5vec)*-1.0*std::get<0>(C_1->coords) + std::get<1>(R5vec)*-1.0*std::get<1>(C_1->coords) + std::get<2>(R5vec)*-1.0*std::get<2>(C_1->coords);
	//Plane parameters for C_2
	double a_2 = -1.0*std::get<0>(R5vec);
	double b_2 = -1.0*std::get<1>(R5vec);
	double c_2 = -1.0*std::get<2>(R5vec);
	double d_2 = -1.0*std::get<0>(R5vec)*-1.0*std::get<0>(C_2->coords) - std::get<1>(R5vec)*-1.0*std::get<1>(C_2->coords) - std::get<2>(R5vec)*-1.0*std::get<2>(C_2->coords);
	//Middle point
	double bond_length = getDistance_twoC(C_1,C_2);

	for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		double dist1 = abs(a_1*std::get<0>(*it1) + b_1*std::get<1>(*it1) + c_1*std::get<2>(*it1) + d_1);
		double dist2 = abs(a_2*std::get<0>(*it1) + b_2*std::get<1>(*it1) + c_2*std::get<2>(*it1) + d_2);
		double dist3 = getDistance_point_to_line(*it1,C_1->coords,C_2->coords);
		if (dist1 <= bond_length && dist2 <= bond_length && dist3 <= bond_length*0.95) return true;
	}
	return false;
	//Previous version - uses endposR5internal
	/*cpair R5_pos_loc = endposR5internal(C_1, C_2, invert_dir);
	std::list<cpair>::iterator it1;
	double minimal_dist = 0.9;
	for (it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		double dist = getDistance_twoC(R5_pos_loc, *it1);
		if (dist <= minimal_dist) return true;
	}
	return false;*/
}

//! Are the two carbon atoms members of an R7 with coordinates in R7Internal??
bool PAHProcess::isR7internal(Cpointer C_1, Cpointer C_2) {
	if(m_pah->m_R7loc.size()==0) return false;
	cpair R7vec = get_vector(C_1->coords,C_2->coords);
	//Plane parameters for C_1
	double a_1 = std::get<0>(R7vec);
	double b_1 = std::get<1>(R7vec);
	double c_1 = std::get<2>(R7vec);
	double d_1 = std::get<0>(R7vec)*-1.0*std::get<0>(C_1->coords) + std::get<1>(R7vec)*-1.0*std::get<1>(C_1->coords) + std::get<2>(R7vec)*-1.0*std::get<2>(C_1->coords);
	//Plane parameters for C_2
	double a_2 = -1.0*std::get<0>(R7vec);
	double b_2 = -1.0*std::get<1>(R7vec);
	double c_2 = -1.0*std::get<2>(R7vec);
	double d_2 = -1.0*std::get<0>(R7vec)*-1.0*std::get<0>(C_2->coords) - std::get<1>(R7vec)*-1.0*std::get<1>(C_2->coords) - std::get<2>(R7vec)*-1.0*std::get<2>(C_2->coords);
	//Middle point
	double bond_length = getDistance_twoC(C_1,C_2);

	for (std::list<cpair>::iterator it1 = m_pah->m_R7loc.begin(); it1 != m_pah->m_R7loc.end(); ++it1) {
		double dist1 = abs(a_1*std::get<0>(*it1) + b_1*std::get<1>(*it1) + c_1*std::get<2>(*it1) + d_1);
		double dist2 = abs(a_2*std::get<0>(*it1) + b_2*std::get<1>(*it1) + c_2*std::get<2>(*it1) + d_2);
		double dist3 = getDistance_point_to_line(*it1,C_1->coords,C_2->coords);
		if (dist1 <= bond_length && dist2 <= bond_length && dist3 <= bond_length*0.95) return true;
	}
	return false;
	//Previous method
	/*cpair R7_pos_loc = endposR7internal(C_1, C_2, invert_dir);
	std::list<cpair>::iterator it1;
	double minimal_dist = 0.9;
	for (it1 = m_pah->m_R7loc.begin(); it1 != m_pah->m_R7loc.end(); ++it1) {
		double dist = getDistance_twoC(R7_pos_loc, *it1);
		if (dist <= minimal_dist) return true;
	}
	return false;*/
}

//! Return coords of final position of an internal R5 based on two carbons
cpair PAHProcess::endposR5internal(Cpointer C_1, Cpointer C_2, bool invert_dir) {
	//C_1 is the carbon that stays in the R5, C_2 the R6 that becomes part of the R5
	double R5_dist = getDistance_twoC(C_1, C_2);
	if (C_1->C2 == C_2 && R5_dist < 2.35){
		double r = R5_dist/(2.0*tan(M_PI/5.0));
		cpair R5dir = get_vector(C_1->coords,C_2->coords);
		cpair normvec;
		if (invert_dir) normvec = invert_vector(norm_vector(C_1->coords, C_2->coords, C_2->C2->coords));
		else normvec = norm_vector(C_1->coords, C_2->coords, C_2->C2->coords);
		cpair crossvec = cross_vector(R5dir, normvec);
		cpair mpos = std::make_tuple(std::get<0>(C_1->coords) + R5_dist/2.0 * std::get<0>(R5dir) + r * std::get<0>(crossvec), 
									 std::get<1>(C_1->coords) + R5_dist/2.0 * std::get<1>(R5dir) + r * std::get<1>(crossvec), 
									 std::get<2>(C_1->coords) + R5_dist/2.0 * std::get<2>(R5dir) + r * std::get<2>(crossvec));
		return mpos;
	}
	else{
		//Assumes that C_2 == C_1->C2->C2
		cpair R5dir = get_vector(C_1->coords,C_2->coords);
		cpair normvec;
		if (invert_dir) normvec = invert_vector(norm_vector(C_1->coords, C_2->coords, C_2->C2->coords));
		else normvec = norm_vector(C_1->coords, C_2->coords, C_2->C2->coords);
		cpair crossvec = cross_vector(R5dir, normvec);
		cpair mpos = std::make_tuple(std::get<0>(C_1->coords) + R5_dist/2.0 * std::get<0>(R5dir) + 0.8 * std::get<0>(crossvec), 
									 std::get<1>(C_1->coords) + R5_dist/2.0 * std::get<1>(R5dir) + 0.8 * std::get<1>(crossvec), 
									 std::get<2>(C_1->coords) + R5_dist/2.0 * std::get<2>(R5dir) + 0.8 * std::get<2>(crossvec));
		return mpos;
	}
}

//! Return coords of final position of an internal R5 based on two carbons
cpair PAHProcess::endposR7internal(Cpointer C_1, Cpointer C_2, bool invert_dir) {
	//C_1 is the carbon that stays in the R5, C_2 the R6 that becomes part of the R5
	double R7_dist = getDistance_twoC(C_1, C_2);
	if (C_1->C2 == C_2 || R7_dist < 2.45){
		cpair R7dir = get_vector(C_1->coords,C_2->coords);
		double r = 1.505;
		cpair normvec;
		if (invert_dir) normvec = invert_vector(norm_vector(C_1->coords, C_2->coords, C_2->C2->coords));
		else normvec = norm_vector(C_1->coords, C_2->coords, C_2->C2->coords);
		cpair crossvec = cross_vector(R7dir, normvec);
		cpair mpos = std::make_tuple(std::get<0>(C_1->coords) + R7_dist/2.0 * std::get<0>(R7dir) + r * std::get<0>(crossvec), 
									 std::get<1>(C_1->coords) + R7_dist/2.0 * std::get<1>(R7dir) + r * std::get<1>(crossvec), 
									 std::get<2>(C_1->coords) + R7_dist/2.0 * std::get<2>(R7dir) + r * std::get<2>(crossvec));
		return mpos;
	}
	else {
		cpair R7dir = get_vector(C_1->coords,C_2->coords);
		double r = 1.505 - 0.629131422 * R7_dist/2.0;
		cpair normvec;
		if (invert_dir) normvec = invert_vector(norm_vector(C_1->coords, C_2->coords, C_2->C2->coords));
		else normvec = norm_vector(C_1->coords, C_2->coords, C_2->C2->coords);
		cpair crossvec = cross_vector(R7dir, normvec);
		cpair mpos = std::make_tuple(std::get<0>(C_1->coords) + R7_dist/2.0 * std::get<0>(R7dir) + r * std::get<0>(crossvec), 
									 std::get<1>(C_1->coords) + R7_dist/2.0 * std::get<1>(R7dir) + r * std::get<1>(crossvec), 
									 std::get<2>(C_1->coords) + R7_dist/2.0 * std::get<2>(R7dir) + r * std::get<2>(crossvec));
		return mpos;
	}
}

int nancounter = 0;
//! Passes a PAH from MOpS to OpenBabel. Returns a mol object.
OpenBabel::OBMol PAHProcess::passPAH(bool detectBonds) {
	std::vector<int> methyl_list_flat;
	//R6 Bay detection
	std::vector<int> R6pairs1, R6pairs2;
	std::vector<int>::iterator R6_iter1, R6_iter2;
	std::vector<Cpointer>CR6_pair1, CR6_pair2;
	std::vector<Cpointer>::iterator resR6, resR62;
	//Second neighbours bond detection
	std::vector<int> C_intlist, first_neighbour, second_neighbour, bridge_neighbour, bridge_neighbour2;
	std::vector<Cpointer>C_list, C_first_neighbour, C_second_neighbour, C_bridge_neighbour, C_bridge_neighbour2;
	std::vector<int>::iterator sn_iter, sn_iter1, sn_iter2, sn_iter3, sn_iter4;
	if (detectBonds){
		for (std::list<Site>::iterator site_it = m_pah->m_siteList.begin(); site_it != m_pah->m_siteList.end(); site_it++) {
			if ( (int)site_it->type % 10 >= 4 ){
				Cpointer CR6_1 = site_it->C1;
				Cpointer CR6_2 = site_it->C2;
				if (getDistance_twoC (CR6_1,CR6_2) <=2.17){
					CR6_pair1.push_back(CR6_1);
					CR6_pair2.push_back(CR6_2);
					R6pairs1.push_back(0);
					R6pairs2.push_back(0);
				}
			}
		}
		Ccontainer::iterator itCsn;
		for (itCsn = m_pah->m_carbonList.begin(); itCsn != m_pah->m_carbonList.end(); ++itCsn) {
			Cpointer C_check = *itCsn;
			C_list.push_back(C_check);
			C_first_neighbour.push_back(C_check->C2);
			C_second_neighbour.push_back(C_check->C2->C2);
			if (C_check->bridge) {
				C_bridge_neighbour.push_back(C_check);
				C_bridge_neighbour2.push_back(C_check->C3);
				bridge_neighbour.push_back(0);
				bridge_neighbour2.push_back(0);
			}
			C_intlist.push_back(0);
			first_neighbour.push_back(0);
			second_neighbour.push_back(0);
		}
	}
	
	//Creates molecule object
	OpenBabel::OBMol mol;
	mol.BeginModify();
	
	//Loads external Carbon Atoms
	int counter = 1;
	Ccontainer::iterator it;
	bool nanflag = false;
	for (it = m_pah->m_carbonList.begin(); it != m_pah->m_carbonList.end(); ++it) {
		bool R6C1 = false; bool R6C2 = false; 
		Cpointer C_change = *it;
		//Edge carbons
		OpenBabel::OBAtom *atom  = mol.NewAtom();
		atom->SetAtomicNum(6);
		atom->SetVector(std::get<0>(C_change->coords),std::get<1>(C_change->coords),std::get<2>(C_change->coords));
		atom->SetAromatic();
		atom->SetImplicitValence(3);
		if ( std::isnan (std::get<0>(C_change->coords))){
			cout << "Error. NaN in coordinates while passing them to OB.\n"; 				//SETBREAKPOINT
			nanflag = true;
		} 
		//atom->SetHyb(2);
		
		//BY6 bond detection.
		if (detectBonds){
			resR6 = std::find(std::begin(CR6_pair1), std::end(CR6_pair1), C_change);
			if (resR6 != std::end(CR6_pair1)){
				int indexR6_1 = std::distance(CR6_pair1.begin(), resR6);
				R6_iter1 = R6pairs1.begin();
				std::advance(R6_iter1, indexR6_1);
				*R6_iter1 = (atom->GetIdx());
				R6C1 = true;
			}
			resR62 = std::find(std::begin(CR6_pair2), std::end(CR6_pair2), C_change);
			if (resR62 != std::end(CR6_pair2)){
				int indexR6_2 = std::distance(CR6_pair2.begin(), resR62);
				R6_iter2 = R6pairs2.begin();
				std::advance(R6_iter2, indexR6_2);
				*R6_iter2 = (atom->GetIdx());
				R6C2 = true;
			}
			
			auto result = std::find(std::begin(C_list), std::end(C_list), C_change);
			int index = std::distance(C_list.begin(), result);
			sn_iter = C_intlist.begin();
			std::advance(sn_iter, index);
			*sn_iter = (atom->GetIdx());
			
			auto result1 = std::find(std::begin(C_first_neighbour), std::end(C_first_neighbour), C_change);
			int index1 = std::distance(C_first_neighbour.begin(), result1);
			sn_iter1 = first_neighbour.begin();
			std::advance(sn_iter1, index1);
			*sn_iter1 = (atom->GetIdx());
			
			auto result2 = std::find(std::begin(C_second_neighbour), std::end(C_second_neighbour), C_change);
			int index2 = std::distance(C_second_neighbour.begin(), result2);
			sn_iter2 = second_neighbour.begin();
			std::advance(sn_iter2, index2);
			*sn_iter2 = (atom->GetIdx());
			
			if(C_change->bridge){
				auto result3 = std::find(std::begin(C_bridge_neighbour), std::end(C_bridge_neighbour), C_change);
				if (result3 != C_bridge_neighbour.end()){
					int index3 = std::distance(C_bridge_neighbour.begin(), result3);
					sn_iter3 = bridge_neighbour.begin();
					std::advance(sn_iter3, index3);
					*sn_iter3 = (atom->GetIdx());
				}
				
				auto result4 = std::find(std::begin(C_bridge_neighbour2), std::end(C_bridge_neighbour2), C_change);
				if (result4 != C_bridge_neighbour2.end()){
					int index4 = std::distance(C_bridge_neighbour2.begin(), result4);
					sn_iter4 = bridge_neighbour2.begin();
					std::advance(sn_iter4, index4);
					*sn_iter4 = (atom->GetIdx());
				}
			}
			
		}
		counter ++;
		
		//Hydrogens
		if (C_change->A == 'H'){
			OpenBabel::OBAtom *atom  = mol.NewAtom();
			atom->SetAtomicNum(1);
			if (R6C1){
				cpair vec1 = get_vector(C_change->C1->coords,C_change->coords);
				cpair vec2 = get_vector(C_change->coords,C_change->C2->coords);
				cpair normvec = cross_vector(vec1, vec2); 
				atom->SetVector(std::get<0>(C_change->coords) + 1.085 * std::get<0>(normvec),std::get<1>(C_change->coords) + 1.085 * std::get<1>(normvec),std::get<2>(C_change->coords) + 1.085 * std::get<2>(normvec));
			}
			else if (R6C2){
				cpair vec1 = get_vector(C_change->C1->coords,C_change->coords);
				cpair vec2 = get_vector(C_change->coords,C_change->C2->coords);
				cpair normvec = cross_vector(vec1, vec2); 
				atom->SetVector(std::get<0>(C_change->coords) - 1.085 * std::get<0>(normvec),std::get<1>(C_change->coords) - 1.085 * std::get<1>(normvec),std::get<2>(C_change->coords) - 1.085 * std::get<2>(normvec));
			}
			else {
				atom->SetVector(std::get<0>(C_change->coords) + 1.085 * std::get<0>(C_change->growth_vector),std::get<1>(C_change->coords) + 1.085 * std::get<1>(C_change->growth_vector),std::get<2>(C_change->coords) + 1.085 * std::get<2>(C_change->growth_vector));
			}
			//atom->SetVector(std::get<0>(C_change->coords) + 1.0*cos(C_change->bondAngle2*M_PI/180),std::get<1>(C_change->coords) + 1.0*sin(C_change->bondAngle2*M_PI/180),std::get<2>(C_change->coords));
			counter ++;
		}
		else if (C_change->A == 'M'){
			OpenBabel::OBAtom *methylc  = mol.NewAtom();
			methylc->SetAtomicNum(6);
			methylc->SetVector(std::get<0>(C_change->coords) + 1.5 * std::get<0>(C_change->growth_vector),std::get<1>(C_change->coords) + 1.5 * std::get<1>(C_change->growth_vector),std::get<2>(C_change->coords) + 1.5 * std::get<2>(C_change->growth_vector));
			cpair Hcoords = std::make_tuple(std::get<0>(C_change->coords) + (1.5 + 1.095*sin(20.8/180.0*M_PI)) * std::get<0>(C_change->growth_vector),std::get<1>(C_change->coords) + (1.5 + 1.095*sin(20.8/180.0*M_PI)) * std::get<1>(C_change->growth_vector),std::get<2>(C_change->coords) + (1.5 + 1.095*sin(20.8/180.0*M_PI)) * std::get<2>(C_change->growth_vector));
			OpenBabel::OBAtom *h1  = mol.NewAtom();
			h1->SetAtomicNum(1);
			OpenBabel::OBAtom *h2  = mol.NewAtom();
			h2->SetAtomicNum(1);
			OpenBabel::OBAtom *h3  = mol.NewAtom();
			h3->SetAtomicNum(1);
			double Radius = 1.095*cos(20.8/180.0*M_PI);
			
			if (abs(std::get<2>(C_change->growth_vector) - 1.0) < 1e-3){
				h1->SetVector(std::get<0>(Hcoords),std::get<1>(Hcoords)+Radius,std::get<2>(Hcoords));
				h2->SetVector(std::get<0>(Hcoords)+Radius*cos(-30.0/180.0*M_PI),std::get<1>(Hcoords)+Radius*sin(-30.0/180.0*M_PI),std::get<2>(Hcoords));
				h3->SetVector(std::get<0>(Hcoords)+Radius*cos(210.0/180.0*M_PI),std::get<1>(Hcoords)+Radius*sin(210.0/180.0*M_PI),std::get<2>(Hcoords));
			}
			else if (abs(std::get<1>(C_change->growth_vector) - 1.0) < 1e-3){
				h1->SetVector(std::get<0>(Hcoords),std::get<1>(Hcoords),std::get<2>(Hcoords)+Radius);
				h2->SetVector(std::get<0>(Hcoords)+Radius*cos(-30.0/180.0*M_PI),std::get<1>(Hcoords),std::get<2>(Hcoords)+Radius*sin(-30.0/180.0*M_PI));
				h3->SetVector(std::get<0>(Hcoords)+Radius*cos(210.0/180.0*M_PI),std::get<1>(Hcoords),std::get<2>(Hcoords)+Radius*sin(210.0/180.0*M_PI));
			}
			else if (abs(std::get<0>(C_change->growth_vector) - 1.0) < 1e-3){
				h1->SetVector(std::get<0>(Hcoords),std::get<1>(Hcoords),std::get<2>(Hcoords)+Radius);
				h2->SetVector(std::get<0>(Hcoords),std::get<1>(Hcoords)+Radius*cos(-30.0/180.0*M_PI),std::get<2>(Hcoords)+Radius*sin(-30.0/180.0*M_PI));
				h3->SetVector(std::get<0>(Hcoords),std::get<1>(Hcoords)+Radius*cos(210.0/180.0*M_PI),std::get<2>(Hcoords)+Radius*sin(210.0/180.0*M_PI));
			}
			else{
				double aa = std::get<0>(C_change->growth_vector);
				double bb = std::get<1>(C_change->growth_vector);
				double cc = std::get<2>(C_change->growth_vector);
				//double rr = sqrt(bb*bb*Radius*Radius*(aa*aa+bb*bb));
				//double h1x = (rr + (aa*aa + bb*bb)*std::get<0>(Hcoords) )/(aa*aa+bb*bb);
				//double h1y = (-aa*rr + (aa*aa*bb + bb*bb*bb)*std::get<1>(Hcoords)) / (aa*aa*bb+bb*bb*bb);
				//double h1z = std::get<2>(Hcoords);
				cpair Angiras_vector = std::make_tuple(cc, 0.0, -1.0*aa);
				Angiras_vector = scale_vector(Angiras_vector);
				cpair Casper_vector = cross_vector(Angiras_vector, C_change->growth_vector);
				
				h1->SetVector(std::get<0>(Hcoords) + Radius*(std::get<0>(Angiras_vector)), std::get<1>(Hcoords) + Radius*(std::get<1>(Angiras_vector)), std::get<2>(Hcoords) + Radius*(std::get<2>(Angiras_vector)));
				h2->SetVector(std::get<0>(Hcoords) + Radius*(-0.5*std::get<0>(Angiras_vector) + sqrt(3.0)/2.0*std::get<0>(Casper_vector)), std::get<1>(Hcoords) + Radius*(-0.5*std::get<1>(Angiras_vector) + sqrt(3.0)/2.0*std::get<1>(Casper_vector)), std::get<2>(Hcoords) + Radius*(-0.5*std::get<2>(Angiras_vector) + sqrt(3.0)/2.0*std::get<2>(Casper_vector)));
				h3->SetVector(std::get<0>(Hcoords) + Radius*(-0.5*std::get<0>(Angiras_vector) - sqrt(3.0)/2.0*std::get<0>(Casper_vector)), std::get<1>(Hcoords) + Radius*(-0.5*std::get<1>(Angiras_vector) - sqrt(3.0)/2.0*std::get<1>(Casper_vector)), std::get<2>(Hcoords) + Radius*(-0.5*std::get<2>(Angiras_vector) - sqrt(3.0)/2.0*std::get<2>(Casper_vector)));
							
				//h1->SetVector( h1x, h1y, h1z);
				//cpair get_vector(Hcoords, std::make_tuple(h1x, h1y, h1z));
				//h2->SetVector( (rr + (aa*aa + bb*bb)*std::get<0>(Hcoords) )/(aa*aa+bb*bb), (-aa*rr + (aa*aa*bb + bb*bb*bb)*std::get<1>(Hcoords)) / (aa*aa*bb+bb*bb*bb), std::get<2>(Hcoords));
			}
			//atom->SetVector(std::get<0>(C_change->coords) + 1.0*cos(C_change->bondAngle2*M_PI/180),std::get<1>(C_change->coords) + 1.0*sin(C_change->bondAngle2*M_PI/180),std::get<2>(C_change->coords));
			counter +=4;
			methyl_list_flat.push_back(methylc->GetIdx()-1);
			methyl_list_flat.push_back(methylc->GetIdx());
			methyl_list_flat.push_back(h1->GetIdx());
			methyl_list_flat.push_back(h2->GetIdx());
			methyl_list_flat.push_back(h3->GetIdx());
		}
	}
	//Internal carbons
	std::list<cpair>::iterator it2;
	for (it2 = m_pah->m_InternalCarbons.begin(); it2 != m_pah->m_InternalCarbons.end(); ++it2){
		OpenBabel::OBAtom *atom  = mol.NewAtom();
		atom->SetAtomicNum(6);
		atom->SetVector(std::get<0>(*it2),std::get<1>(*it2),std::get<2>(*it2));
		atom->SetAromatic();
		atom->SetImplicitValence(3);
		double checkcoord = atom->GetX();
		if (std::isnan(checkcoord)) {
			cout << "Error. NaN in internal coordinates while passing them to OB.\n"; //SETBREAKPOINT
			nanflag = true;
		} 
		//atom->SetHyb(2);
		counter ++;
	}
	
	if (detectBonds){
		//Bond generation
		mol.ConnectTheDots();
		//connectPAH(mol);
		//Deletes bonds longer than 1.85A. 
		vector <unsigned int> bonds_to_delete;
		for (OpenBabel::OBBondIterator bond_iter=mol.BeginBonds(); bond_iter != mol.EndBonds(); bond_iter++){
			OpenBabel::OBBond* bond = *bond_iter;
			double length_bond = bond->GetLength();
			if (length_bond >= 1.85) bonds_to_delete.push_back( (*bond_iter)->GetIdx() );
		}
		if(bonds_to_delete.size() != 0){
			vector <unsigned int>::iterator itb=bonds_to_delete.end();
			itb--;
			for (OpenBabel::OBBondIterator it=mol.EndBonds(); true; ){
				it--;
				if ( (*it)->GetIdx() == (*itb) ){
					mol.DeleteBond((*it), true);
					if (itb == bonds_to_delete.begin()) {break;}
					else {itb--;}
				}
			}
			
		}
		//Adds Bonds within C and its C2. Deletes bonds between C and C->C2->C2. Adds bonds between bridged carbons.
		std::vector<int>::iterator sn_iter;
		std::vector<int>::iterator sn_iter1 = first_neighbour.begin();
		std::vector<int>::iterator sn_iter2 = second_neighbour.begin();
		std::vector<int>::iterator sn_iter3 = bridge_neighbour.begin();
		std::vector<int>::iterator sn_iter4 = bridge_neighbour2.begin();
		for(sn_iter = C_intlist.begin(); sn_iter != C_intlist.end(); ++sn_iter){
			OpenBabel::OBBond* my_bond = mol.GetBond(*sn_iter, *sn_iter1);
			OpenBabel::OBBond* my_bond2 = mol.GetBond(*sn_iter, *sn_iter2);
			if (my_bond == NULL) mol.AddBond(*sn_iter, *sn_iter1,5);
			if (my_bond2 != NULL) mol.DeleteBond(my_bond2);
			my_bond = mol.GetBond(*sn_iter, *sn_iter1);
			if (my_bond->GetLength()>8.0) {
				std::cout << "WARNING. Bond length larger than 8 Armstrongs before optimisation. Deleting bond." << std::endl;
				mol.DeleteBond(my_bond);
			}
			++sn_iter1;
			++sn_iter2;
			auto find_bridge = std::find(std::begin(bridge_neighbour), std::end(bridge_neighbour), *sn_iter);
			if (find_bridge != bridge_neighbour.end()){
				OpenBabel::OBBond* my_bond_bridge = mol.GetBond(*sn_iter3, *sn_iter4);
				if (my_bond_bridge == NULL) mol.AddBond(*sn_iter3, *sn_iter4,5);
				++sn_iter3;
				++sn_iter4;
			}
		}
		
		//Deletes bonds within unclosed BY6.
		std::vector<int>::iterator it_R6pairs1, it_R6pairs2;
		it_R6pairs2 = R6pairs2.begin();
		for(it_R6pairs1 = R6pairs1.begin(); it_R6pairs1 != R6pairs1.end(); ++it_R6pairs1){
			OpenBabel::OBBond* my_bond = mol.GetBond(*it_R6pairs1, *it_R6pairs2);
			if (my_bond != NULL) mol.DeleteBond(my_bond);
			++it_R6pairs2;
		}
		
		//Checks hydrogen bonds //NEEDS DEBUGGING. NOT COMPLETE.
		for(OpenBabel::OBMolAtomIter     a(mol); a; ++a) {
			auto find_methyl = std::find(methyl_list_flat.begin(), methyl_list_flat.end(), a->GetIdx()); 
			if (find_methyl == methyl_list_flat.end()){
				//The atom is not in a methyl group.
				if (a->GetAtomicNum() == 1){
					OpenBabel::OBAtom *prev_atom = mol.GetAtom(a->GetIdx() - 1);
					if (prev_atom -> GetAtomicNum() == 6){
						OpenBabel::OBBond* my_bond = mol.GetBond(a->GetIdx(), a->GetIdx() - 1);
						if (my_bond == NULL){
							vector <OpenBabel::OBBond*> bonds_to_delete;
							for( OpenBabel::OBAtomAtomIter     nbr(*a); nbr; ++nbr ) { 
								OpenBabel::OBBond* my_bond_del = mol.GetBond(a->GetIdx(), nbr->GetIdx());
								bonds_to_delete.push_back(my_bond_del);
							}
							if(bonds_to_delete.size() != 0){
								vector <OpenBabel::OBBond*>::iterator itb;
								for (itb=bonds_to_delete.begin(); itb!=bonds_to_delete.end(); itb++){
									if ( (*itb) != NULL ){
										mol.DeleteBond((*itb), true);
									}
								}
							}
							mol.AddBond(a->GetIdx(), a->GetIdx() - 1,5);
						}
					}
				}
			}
			else {
				//The atom is in a methyl group. We need to find in which methyl group first.
				int index_methyl = std::distance(methyl_list_flat.begin(), find_methyl);
				if (index_methyl%5 == 0){
					//This is the carbon attached to the methyl.
					OpenBabel::OBBond* my_bond = mol.GetBond(a->GetIdx(), a->GetIdx() + 1);
					if (my_bond == NULL) mol.AddBond(a->GetIdx(), a->GetIdx() + 1,5);
				}
				else if (index_methyl%5 == 1){
					//This is the methyl carbon.
					for (int i=1; i<=3; i++){
						OpenBabel::OBBond* my_bond = mol.GetBond(a->GetIdx(), a->GetIdx() + i);
						if (my_bond == NULL) mol.AddBond(a->GetIdx(), a->GetIdx() + i,5);
					}
				}
				else {
					//These are the hydrogens from the methyl group.
					vector <OpenBabel::OBBond*> bonds_to_delete;
					for( OpenBabel::OBAtomAtomIter     nbr(*a); nbr; ++nbr ) {
						auto find_nbr = std::find(std::begin(methyl_list_flat), std::end(methyl_list_flat), nbr->GetIdx());
						if (find_nbr == methyl_list_flat.end()){
							//Hydrogen is attached to other atom. Delete this bond.
							OpenBabel::OBBond* my_bond_del = mol.GetBond(a->GetIdx(), nbr->GetIdx());
							bonds_to_delete.push_back(my_bond_del);
						}
						if(bonds_to_delete.size() != 0){
							vector <OpenBabel::OBBond*>::iterator itb;
							for (itb=bonds_to_delete.begin(); itb!=bonds_to_delete.end(); itb++){
								if ( (*itb) != NULL ){
									mol.DeleteBond((*itb), true);
								}
							}
						}	
					}
				}
			}
		}
		
		vector<int> valence_2_carbons, valence_4_carbons;
		for(OpenBabel::OBMolAtomIter     a(mol); a; ++a) {
			auto find_methyl = std::find(methyl_list_flat.begin(), methyl_list_flat.end(), a->GetIdx()); 
			if (find_methyl == methyl_list_flat.end()){
				//The atom is not in a methyl group.
				if (a->GetAtomicNum() == 6){
					//Check for carbons with valence 2 and add them to list.
					if (a->GetValence() == 2) valence_2_carbons.push_back(a->GetIdx());
					if (a->GetValence() == 4) valence_4_carbons.push_back(a->GetIdx());
				}
			}
		}
		
		//Adds bond between two carbons with valence 2. 
		if(valence_2_carbons.size() != 0){ 
			for (unsigned int ii=0; ii!=valence_2_carbons.size(); ++ii){
				OpenBabel::OBAtom *my_atom  = mol.GetAtom(valence_2_carbons[ii]);
				if (my_atom->GetValence() == 2){
					int kk = ii;
					double min_dist = 1e3;
					for (unsigned int jj =0; jj!=valence_2_carbons.size(); ++jj){
						if (jj!= ii){
							double my_dist = my_atom->GetDistance(valence_2_carbons[jj]);
							if (my_dist < min_dist){
								kk = jj;
								min_dist = my_dist;
							}
						}
					}
					if (kk != ii){
						OpenBabel::OBBond* my_bond = mol.GetBond(valence_2_carbons[ii], valence_2_carbons[kk]);
						if (my_bond == NULL) mol.AddBond(valence_2_carbons[ii], valence_2_carbons[kk],5);
						my_bond = mol.GetBond(valence_2_carbons[ii], valence_2_carbons[kk]);
						if (my_bond->GetLength()>8.0) {
							std::cout << "WARNING. Bond length larger than 8 Armstrongs before optimisation. Deleting bond." << std::endl;
							mol.DeleteBond(my_bond);
						}
					}
				}
			}
		}
		
		//Deletes bond between two carbons with valence 4. 
		if(valence_4_carbons.size() != 0){ 
			for (unsigned int ii=0; ii!=valence_4_carbons.size(); ++ii){
				OpenBabel::OBAtom *my_atom  = mol.GetAtom(valence_4_carbons[ii]);
				if (my_atom->GetValence() == 4){
					int kk = ii;
					double min_dist = 1e3;
					for (unsigned int jj =0; jj!=valence_4_carbons.size(); ++jj){
						if (jj!= ii){
							double my_dist = my_atom->GetDistance(valence_4_carbons[jj]);
							if (my_dist < min_dist){
								kk = jj;
								min_dist = my_dist;
							}
						}
					}
					if (kk != ii){
						OpenBabel::OBBond* my_bond = mol.GetBond(valence_4_carbons[ii], valence_4_carbons[kk]);
						if (my_bond != NULL) mol.DeleteBond(my_bond, true);
					}
				}
			}
		}
		
		mol.SetAromaticPerceived();
		//Comment in to print the details of every bond in the OBmol object
		/*std::cout << "Bond  First Atom     Second Atom   Length    Order   Aromatic\n";
		for (OpenBabel::OBBondIterator bond_iter=mol.BeginBonds(); bond_iter != mol.EndBonds(); bond_iter++){
			OpenBabel::OBBond* bond = *bond_iter;
			int bond_number = bond->GetIdx();
			double length_bond = bond->GetLength();
			int first_idx = bond->GetBeginAtomIdx();
			int second_idx = bond->GetEndAtomIdx();
			int order = bond->GetBondOrder();
			int aromatic, aromatic_a1, aromatic_a2, hyb_a1, hyb_a2;
			if (bond->IsAromatic()) aromatic = 1;
			else aromatic = 0;
			OpenBabel::OBAtom *my_first_atom  = mol.GetAtom(first_idx);
			if (my_first_atom->IsAromatic()) aromatic_a1 = 1;
			else aromatic_a1 = 0;
			hyb_a1 = my_first_atom->GetHyb();
			OpenBabel::OBAtom *my_second_atom  = mol.GetAtom(second_idx);
			if (my_second_atom->IsAromatic()) aromatic_a2 = 1;
			else aromatic_a2 = 0;
			hyb_a2 = my_second_atom->GetHyb();

			std::cout << bond_number << "     " << first_idx << " (" << aromatic_a1 <<") ("<< hyb_a1 << ")       " << second_idx << " (" << aromatic_a2 <<") ("<< hyb_a2 << ")       " << length_bond << "      " << order << "       " << aromatic << "\n";
		}*/

		//PerceiveBondOrders calls several routines that try to identify aromatic and unaromatic parts of a molecule. This is very nice but expensive.
		//If this is not called, OB recognises bonds as single bonds.
		//mol.PerceiveBondOrders();

		//Comment in to print the details of every bond in the OBmol object
		mol.SetAromaticPerceived();
		/*std::cout << "Bond  First Atom     Second Atom     Length    Order   Aromatic\n";
		for (OpenBabel::OBBondIterator bond_iter=mol.BeginBonds(); bond_iter != mol.EndBonds(); bond_iter++){
			OpenBabel::OBBond* bond = *bond_iter;
			int bond_number = bond->GetIdx();
			double length_bond = bond->GetLength();
			int first_idx = bond->GetBeginAtomIdx();
			int second_idx = bond->GetEndAtomIdx();
			int order = bond->GetBondOrder();
			int aromatic, aromatic_a1, aromatic_a2, hyb_a1, hyb_a2;
			if (bond->IsAromatic()) aromatic = 1;
			else aromatic = 0;
			OpenBabel::OBAtom *my_first_atom  = mol.GetAtom(first_idx);
			if (my_first_atom->IsAromatic()) aromatic_a1 = 1;
			else aromatic_a1 = 0;
			hyb_a1 = my_first_atom->GetHyb();
			OpenBabel::OBAtom *my_second_atom  = mol.GetAtom(second_idx);
			if (my_second_atom->IsAromatic()) aromatic_a2 = 1;
			else aromatic_a2 = 0;
			hyb_a2 = my_second_atom->GetHyb();

			std::cout << bond_number << "     " << first_idx << " (" << aromatic_a1 <<") ("<< hyb_a1 << ")       " << second_idx << " (" << aromatic_a2 <<") ("<< hyb_a2 << ")       " << length_bond << "      " << order << "       " << aromatic << "\n";
		}*/
		
		//Checks for sp3 carbon //NEEDS DEBUGGING. NOT COMPLETE.
		/*for(OpenBabel::OBMolAtomIter     a(mol); a; ++a) {
			if (a->GetAtomicNum() == 6 && a->GetValence() >= 4){
				OpenBabel::OBAtomAtomIter wrong_nbr; 
				double maxdist = 0.1;
				for( OpenBabel::OBAtomAtomIter     nbr(*a); nbr; ++nbr ) {
					if (nbr->GetAtomicNum() == 6 && a->GetDistance(nbr->GetIdx()) > maxdist){
						wrong_nbr = nbr;
						maxdist = a->GetDistance(nbr->GetIdx());
					}
				}
				
				OpenBabel::OBBond* my_bond_del = mol.GetBond(a->GetIdx(), wrong_nbr->GetIdx());
				mol.DeleteBond(my_bond_del);
			}
		}*/
		
		//Assign aromatic bond orders to all bonds. //NOT WORKING. DO NOT USE.
		/*for (OpenBabel::OBBondIterator bond_iter=mol.BeginBonds(); bond_iter != mol.EndBonds(); bond_iter++){
			OpenBabel::OBBond* my_bond = *bond_iter;
			int BOrder = my_bond->GetBO();
			int a1 = my_bond->GetBeginAtomIdx(); int a2 = my_bond->GetEndAtomIdx(); 
			cout << a1 << "-" << a2 << " Bond order = " << BOrder << "\n";
		}*/
		
		
		/*for (OpenBabel::OBBondIterator bond_iter=mol.BeginBonds(); bond_iter != mol.EndBonds(); bond_iter++){
			OpenBabel::OBBond* my_bond = *bond_iter;
			OpenBabel::OBAtom *a1 = my_bond->GetBeginAtom(); OpenBabel::OBAtom *a2 = my_bond->GetEndAtom(); 
			//int a1int = my_bond->GetBeginAtomIdx(); int a2int = my_bond->GetEndAtomIdx();
			if ( a1->GetAtomicNum() == 6 && a2->GetAtomicNum() == 6){
				my_bond->SetBO(5);
			}
			int BOrder = my_bond->GetBO();
			//cout << a1int << "-" << a2int << " Bond order = " << BOrder << "\n";
		}
		mol.SetAromaticPerceived();*/
		
	}
	//mol.AddHydrogens();
	mol.EndModify();
	if (nanflag){
		//NaN in coordinates. Output a txt file with all coordinates.
		ofstream ofs1;
		std::string filename = "KMC_DEBUG/NAN_error_";
		filename.append(std::to_string(nancounter));
		filename.append(".xyz");
		OpenBabel::OBConversion conv;
		OpenBabel::OBFormat *format_out = conv.FindFormat("xyz"); // default output format
		conv.SetInAndOutFormats(format_out, format_out);
		ofs1.open(filename);
		conv.Write(&mol, &ofs1);
		ofs1.close();
		cout<<"Saving file " << filename<<"\n";
		printSites();
		nancounter++;
	}
	return mol;
}
//Needed for connect the dots.
bool SortAtomZ(const std::pair<OpenBabel::OBAtom*,double> &a, const std::pair<OpenBabel::OBAtom*,double> &b){
    return (a.second < b.second);
}

//! Connects the atoms in a PAH using OpenBabel routines. Equivalent to OpenBabel::OBMol::ConnectTheDots();
void PAHProcess::connectPAH(OpenBabel::OBMol my_mol) {
    int j,k,max;
    //bool unset = false;
    OpenBabel::OBAtom *atom,*nbr;
    vector<OpenBabel::OBAtom*>::iterator i;
    vector<pair<OpenBabel::OBAtom*,double> > zsortedAtoms;
    vector<double> rad;
    vector<int> zsorted;
    vector<int> bondCount; // existing bonds (e.g., from residues in PDB)

    double *c = new double [my_mol.NumAtoms()*3];
    rad.resize(my_mol.NumAtoms());

    for (j = 0, atom = my_mol.BeginAtom(i) ; atom ; atom = my_mol.NextAtom(i), ++j)
      {
        (atom->GetVector()).Get(&c[j*3]);
        pair<OpenBabel::OBAtom*,double> entry(atom, atom->GetVector().z());
        zsortedAtoms.push_back(entry);
        bondCount.push_back(atom->GetValence());
      }
    sort(zsortedAtoms.begin(), zsortedAtoms.end(), SortAtomZ);

    max = zsortedAtoms.size();

    for ( j = 0 ; j < max ; j++ )
      {
        atom   = zsortedAtoms[j].first;
        rad[j] = OpenBabel::etab.GetCovalentRad(atom->GetAtomicNum());
        zsorted.push_back(atom->GetIdx()-1);
      }

    int idx1, idx2;
    double d2,cutoff,zd;
    for (j = 0 ; j < max ; ++j)
      {
        idx1 = zsorted[j];
        for (k = j + 1 ; k < max ; k++ )
          {
            idx2 = zsorted[k];

            // bonded if closer than elemental Rcov + tolerance
            cutoff = SQUARE(rad[j] + rad[k] + 0.45);

            zd  = SQUARE(c[idx1*3+2] - c[idx2*3+2]);
            if (zd > 25.0 )
              break; // bigger than max cutoff

            d2  = SQUARE(c[idx1*3]   - c[idx2*3]);
            d2 += SQUARE(c[idx1*3+1] - c[idx2*3+1]);
            d2 += zd;

            if (d2 > cutoff)
              continue;
            if (d2 < 0.40)
              continue;

            atom = my_mol.GetAtom(idx1+1);
            nbr  = my_mol.GetAtom(idx2+1);

            if (atom->IsConnected(nbr))
              continue;
            if (atom->IsHydrogen() && nbr->IsHydrogen())
              continue;

            my_mol.AddBond(idx1+1,idx2+1,1);
          }
      }

    // If between BeginModify and EndModify, coord pointers are NULL
    // setup molecule to handle current coordinates
/*
    if (_c == NULL)
      {
        _c = c;
        for (atom = my_mol.BeginAtom(i);atom;atom = my_mol.NextAtom(i))
          atom->SetCoordPtr(&_c);
        _vconf.push_back(c);
        unset = true;
      }
*/
    // Cleanup -- delete long bonds that exceed max valence
    OpenBabel::OBBond *maxbond, *bond;
    double maxlength;
    vector<OpenBabel::OBBond*>::iterator l;
    int valCount;

    for (atom = my_mol.BeginAtom(i);atom;atom = my_mol.NextAtom(i))
      {
        while (atom->BOSum() > static_cast<unsigned int>(OpenBabel::etab.GetMaxBonds(atom->GetAtomicNum()))
               || atom->SmallestBondAngle() < 45.0)
          {
            bond = atom->BeginBond(l);
            maxbond = bond;
            // Fix from Liu Zhiguo 2008-01-26
            // loop past any bonds
            // which existed before ConnectTheDots was called
            // (e.g., from PDB resdata.txt)
            valCount = 0;
            while (valCount < bondCount[atom->GetIdx() - 1]) {
              bond = atom->NextBond(l);
              if (!bond) // so we add an additional check
                break;
              maxbond = bond;
              valCount++;
            }
            if (!bond) // no new bonds added for this atom, just skip it
              break;

            maxlength = maxbond->GetLength();
            for (bond = atom->NextBond(l);bond;bond = atom->NextBond(l))
              {
                if (!bond)
                  break;
                if (bond->GetLength() > maxlength)
                  {
                    maxbond = bond;
                    maxlength = bond->GetLength();
                  }
              }
            my_mol.DeleteBond(maxbond); // delete the new bond with the longest length
          }
      }

    /*if (unset)
      {
        _c = NULL;
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
          atom->ClearCoordPtr();
        _vconf.resize(_vconf.size()-1);
      }*/

    delete [] c;
 }

//! Passes a PAH from OpenBabel to MOpS. Returns a mol object.
void PAHProcess::passbackPAH(OpenBabel::OBMol mol) {
	//Passes new coordinates back to MOPS
	Ccontainer::iterator it = m_pah->m_carbonList.begin();
	std::list<cpair>::iterator it2 = m_pah->m_InternalCarbons.begin();
	for(OpenBabel::OBMolAtomIter     a(mol); a; ++a) {
		if (it != m_pah->m_carbonList.end()){
			double x_pos = a->GetX();
			double y_pos = a->GetY();
			double z_pos = a->GetZ();
			cpair temp = std::make_tuple(x_pos, y_pos, z_pos);
			Cpointer C_back = *it;
			//cout << "BEFORE passing coordinates = " << std::get<0>(C_back->coords) << " " << std::get<1>(C_back->coords) << " " << std::get<2>(C_back->coords) << "\n";
			moveC(C_back, temp);
			//cout << "OPENBABEL coordinates =      " << a->GetX() << " " << a->GetY() << " " << a->GetZ() << "\n";
			//cout << "AFTER passing coordinates =  " << std::get<0>(C_back->coords) << " " << std::get<1>(C_back->coords) << " " << std::get<2>(C_back->coords) << "\n";
			
			if (C_back->A == 'H') {
				cpair p1 = C_back->coords;
				++a;
				double Hx_pos = a->GetX();
				double Hy_pos = a->GetY();
				double Hz_pos = a->GetZ();
				cpair p2 = std::make_tuple(Hx_pos, Hy_pos, Hz_pos);
				cpair vec = get_vector(p1, p2);
				//cout <<"C vector direction = " << std::get<0>(C_back->growth_vector) << " " << std::get<1>(C_back->growth_vector) << " " << std::get<2>(C_back->growth_vector) << "\n";
				C_back->growth_vector = vec;
				//cout <<"C vector direction = " << std::get<0>(C_back->growth_vector) << " " << std::get<1>(C_back->growth_vector) << " " << std::get<2>(C_back->growth_vector) << "\n";
			}
			
			if (C_back->A == 'M') {
				cpair p1 = C_back->coords;
				++a;
				double Hx_pos = a->GetX();
				double Hy_pos = a->GetY();
				double Hz_pos = a->GetZ();
				cpair p2 = std::make_tuple(Hx_pos, Hy_pos, Hz_pos);
				cpair vec = get_vector(p1, p2);
				//cout <<"C vector direction = " << std::get<0>(C_back->growth_vector) << " " << std::get<1>(C_back->growth_vector) << " " << std::get<2>(C_back->growth_vector) << "\n";
				C_back->growth_vector = vec;
				/*for( OpenBabel::OBAtomAtomIter nbr(*a); nbr; ++nbr )
				{
					mol.DeleteAtom(&*nbr);
				}*/
				//cout <<"C vector direction = " << std::get<0>(C_back->growth_vector) << " " << std::get<1>(C_back->growth_vector) << " " << std::get<2>(C_back->growth_vector) << "\n";
				++a; ++a; ++a;
			}
			++it;
		}
		else {
			cpair temp = std::make_tuple(a->GetX(),a->GetY(),a->GetZ());
			*it2 = temp;
			++it2;
		}
	}
	//Modify the R5 and R7 centers and pass them back to MOpS.
	m_pah->m_R5loc.clear(); m_pah->m_R7loc.clear();
	vector<OpenBabel::OBRing*>::iterator iring;
	vector<int>::iterator j;
	vector<OpenBabel::OBRing*> vr;
    vr = mol.GetSSSR();
	//vector<OpenBabel::OBRing*> *rlist = (vector<OpenBabel::OBRing*>*)mol.GetData("RingList");
	for (iring = vr.begin();iring != vr.end();++iring){
		//cout<<(**iring).PathSize()<<"\n";
		if( (**iring).PathSize()==5 || (**iring).PathSize()==7 ){
			OpenBabel::vector3 centre, norm1, norm2;
			bool centre_found = (**iring).findCenterAndNormal(centre, norm1, norm2);
			cpair temp = std::make_tuple(centre.GetX(),centre.GetY(),centre.GetZ());
			if ((**iring).PathSize()==5) m_pah->m_R5loc.push_back(temp);
			else m_pah->m_R7loc.push_back(temp);
		}
	}
	m_pah->m_optimised = true;
}
int total_error_counter = 0;
int forcefield_error_counter = 0;
//! Minimisation of a PAH
OpenBabel::OBMol PAHProcess::optimisePAH(OpenBabel::OBMol mol, int nsteps, std::string forcefield) {
	mol.BeginModify();
	//Defines a forcefield object
	//OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField("Ghemical");
	OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(forcefield);
	if (!pFF) {
    cerr << ": could not find forcefield " << forcefield << "." <<endl;
    exit (-1);
	}
	
	//pFF->SetLogFile(&cerr);
	//pFF->SetLogLevel(OBFF_LOGLVL_MEDIUM);
	/*pFF->SetVDWCutOff(6.0);
	pFF->SetElectrostaticCutOff(10.0);
	pFF->SetUpdateFrequency(10);
	pFF->EnableCutOff(false);*/
	
	//Initialise minimisation
	if (!pFF->Setup(mol)) {
	  /*std::string filename = "KMC_DEBUG/Forcefield_log.txt";
	  ofstream ofs_ff(filename);
	  pFF->SetLogFile(&ofs_ff);
	  pFF->SetLogLevel(OBFF_LOGLVL_MEDIUM);
	  pFF->Setup(mol);*/
      cout << "Error: could not setup force field.\n" << endl;
	  cout << "Sites before calling optimiser:\n"; //SETBREAKPOINT
	  printSites();
	  std::string filename_error = "KMC_DEBUG/Forcefield_error_";
	  filename_error.append(std::to_string(forcefield_error_counter));
	  saveXYZ(filename_error);
	  cout<<"Saving file: "<< filename_error<<".xyz\n";
	  //cerr << ": could not setup force field." << endl;
      //exit (-1);
	  ++forcefield_error_counter;
	  return mol;
    }
	
	/////// Method that uses initialisation. Not recommended by OpenBabel documentation.
	/*bool done = true;
	pFF->SteepestDescentInitialize(nsteps, 1e-7);
	
	//Perform minimisation
	unsigned int totalSteps = 1;
    while (done) {
		done = pFF->SteepestDescentTakeNSteps(1);
		totalSteps++;
		if (pFF->DetectExplosion()) {
			cerr << "explosion has occured!" << endl;
		}
		else
        pFF->GetCoordinates(mol);
	}*/
	//// Method recommended in case minimisation is just used until the nsteps without modifications.
	pFF->SteepestDescent(nsteps, 1e-5);
	
	pFF->GetCoordinates(mol);
	mol.EndModify();
	return mol;
}

/*! Create a carbon atom bridging next to C_1
  !              newC
  !                /  <---- angle (from C1)
  !          R - C1
  !                \
  !                 R
*/
Cpointer PAHProcess::bridgeC(Cpointer C_1) {
    Cpointer cb;
    // Create new carbon atom
    cb = new Carbon;
    // Set details of new Carbon
    cb->C3 = C_1;
	//cb->bondAngle2 = normAngle(C_1->bondAngle2 - 180); // opposite direction of C_1->bondAngle2
    cb->bridge = true;
	cb->growth_vector = C_1->growth_vector;
    // set new coordinates and store
	cb->coords = jumpToPos(C_1->coords, C_1->growth_vector, 1.4);
    //cb->coords = jumpToPos(C_1->coords, C_1->bondAngle2, 0, 1.4);
    if(!m_pah->m_carbonList.insert(cb).second)
        std::cout<<"ERROR: ADDING SAME CARBON POINTER TO SET\n";
    m_pah->m_cpositions.insert(cb->coords);
    // Set details of C_1
    C_1->bridge = true;
    C_1->C3 = cb;
    updateA(C_1, 'C', C_1->growth_vector);
    addCount(1,0);
	//Any Carbon addition means that the structure will not have an optimised geometry.
	m_pah->m_optimised = false;
    return cb;
}/*
//! Creates a bulk carbon atom connected to C_1
//          C2        C2                                                  
//         /           \                                                  
//  C1 - C              C - C1                                            
//         \           /                                                  
//          C3        C3                                                  
Cpointer PAHProcess::addBC(Cpointer C_1) {
    Cpointer cb;
    // Create new carbon atom
    cb = new Carbon;
    // determine if C_1 is C1, C2 or C3
    cb->C1 = NULL; cb->C2 = NULL; cb->C3 = NULL;
    if(C_1->bondAngle2 == 0 || C_1->bondAngle2 == 180 || C_1->bondAngle2 == -180) 
        cb->C1 = C_1;
    else if(C_1->bondAngle2 == 60 || C_1->bondAngle2 == 120)
        cb->C2 = C_1;
    else if(C_1->bondAngle2 == -60 || C_1->bondAngle2 == -120)
        cb->C3 = C_1;
    else std::cout<<"ERROR: addBC(): Edge bondAngle2 invalid!!\n";
    // Set details of new Carbon
    cb->A = 'C';
    // set new coordinates and store
    cb->coords = jumpToPos(C_1->coords, C_1->bondAngle2);
    m_pah->m_cpositions.insert(cb->coords);
    // Set details of C_1
    C_1->C3 = cb;
    C_1->A = 'C';
    addCount(1,0);
    return cb;
}*/
//! Connects a carbon atom to another carbon (to close loop)
void PAHProcess::connectToC(Cpointer C_1, Cpointer C_2) {
    C_1->C2 = C_2;
    C_2->C1 = C_1;
}
//! Removes a carbon atom from perimeter taking into account if the atom becomes bulk carbon within the PAH
void PAHProcess::removeC(Cpointer C_1, bool bulk) {
    if(C_1 == m_pah->m_cfirst) {
        m_pah->m_cfirst = C_1->C2;
    }else if(C_1 == m_pah->m_clast) {
        m_pah->m_clast = C_1->C1;
    }
	double min_dist = 1e3;
	std::set<cpair>::iterator it2;
	for (std::set<cpair>::iterator it = m_pah->m_cpositions.begin(); it != m_pah->m_cpositions.end(); ++it) {
		double x_dist = std::get<0>(C_1->coords) - std::get<0>(*it);
		double y_dist = std::get<1>(C_1->coords) - std::get<1>(*it);
		double z_dist = std::get<2>(C_1->coords) - std::get<2>(*it);
		double dist = sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);
		if (dist < min_dist) {
			it2 = it;
			min_dist = dist;
		}
		if (dist < 1e-2) break;
	}
	if (it2 == m_pah->m_cpositions.end()){
    //if(m_pah->m_cpositions.find(C_1->coords)== m_pah->m_cpositions.end()) {
		cout << "ERROR: removeC: coordinates (" << std::get<0>(C_1->coords) << ',' << std::get<1>(C_1->coords) << ',' << std::get<2>(C_1->coords)<<") not in m_pah->m_cpositions!\n";
        cout<<"Coordinates of nearby 5 C atoms:\n";
        Cpointer now = C_1->C1->C1->C1->C1->C1;
        for(int i=0; i!=11; i++) {
			cout << '(' << std::get<0>(now->coords) << ',' << std::get<1>(now->coords) << ',' << std::get<2>(now->coords) <<")";
            if(now == C_1) cout<<" - ";
            cout<<'\n';
            now = now->C2;
        }
        //saveDOT("KMC_DEBUG/COORD_PROB.dot");
		saveXYZ("KMC_DEBUG/COORD_PROB.dot");
    }
    // Change details of neighbouring Carbons
    if(C_1->C1 != NULL) { // check prev C atom
        C_1->C1->C2 = C_1->C2; // connect prev C atom to next C atom
    }
    if(C_1->C2 != NULL) { // check next C atom
        C_1->C2->C1 = C_1->C1; // connect next C atom to prev C atom
    }
    if(C_1->bridge) { // change details of C atom bridging to it
        //C_1->C3->bondAngle2;// = 0;
        C_1->C3->C3 = NULL;
        C_1->C3->bridge = false;
    }
    // Remove coordinates of C from m_pah->m_cpositions and add them to m_pah->m_InternalCarbons
	if(bulk) m_pah->m_InternalCarbons.push_back(C_1->coords);
	m_pah->m_cpositions.erase(*it2);
    //m_pah->m_cpositions.erase(m_pah->m_cpositions.find(C_1->coords));
    // delete Carbon object
    delete C_1;
    if(m_pah->m_carbonList.erase(C_1) == 0)
        std::cout<<"ERROR: removeC: NOT ERASING ANY POINTERS!\n";
    if(!bulk) { // if desorption, decrease C count accordingly
        addCount(-1, 0); 
    }
	//Any Carbon removal means that the structure will not have an optimised geometry.
	m_pah->m_optimised = false;
}

//! Adds a site with members C_1 & C_2 before site b4site
Spointer PAHProcess::addSite(kmcSiteType stype, Cpointer C_1, Cpointer C_2, Spointer& b4site) {
    // Create new Site
    Site st;
    // Set details of site
    st.type = stype;
    st.C1 = C_1;
    st.C2 = C_2;
    st.comb = None;
    // Add to list
    Spointer newSite = m_pah->m_siteList.insert(b4site, st);
    // Set details of carbon members
    //C_1->S2 = newSite;
    //C_2->S1 = newSite;
    // adds iterator into site vector
    m_pah->m_siteMap[stype].push_back(newSite);
    return newSite;
}
//! Adds a site with members C_1 & C_2 at end of SiteList
Spointer PAHProcess::addSite(kmcSiteType stype, Cpointer C_1, Cpointer C_2) {
    // Create new Site
    Site st;
    // Set details of site
    st.type = stype;
    st.C1 = C_1;
    st.C2 = C_2;
    st.comb = None;
    // Add to list
    m_pah->m_siteList.push_back(st);
    Spointer newSite = --(m_pah->m_siteList.end());
    // Set details of carbon members
    //C_1->S2 = newSite;
    //C_2->S1 = newSite;
    // adds iterator into site vector
    m_pah->m_siteMap[stype].push_back(newSite);
    return newSite;
}
//! Removes a site
//void PAHProcess::removeSite(Spointer st) {}
//! Changes site type into combined site with R5 (e.g. FE -> RFE)
void PAHProcess::addR5toSite(Spointer& st, Cpointer Carb1, Cpointer Carb2) {
    int stype = (int) st->type;
	int index_change = 101;
	updateSites(st, Carb1, Carb2, index_change); 
	/*if ((kmcSiteType)stype == None) {
		std::cout << "ERROR: add5RtoSite: illegal site type to add 5R to.\n"; 
		printSites(st);
		return;
	}
	
	stype += 101;
	
	// removes site from m_pah->m_siteMap (principal site)
	delSiteFromMap(st->type, st);
	// change site type
	st->type = (kmcSiteType)stype;
	// add site to m_pah->m_siteMap
	m_pah->m_siteMap[st->type].push_back(st);
	// update member C
	st->C1 = Carb1;
	st->C2 = Carb2;*/

	/*
    if(stype<5) // i.e. uncombined principal sites (eg FE)
        stype += 6;
	else if (stype < 10) // i.e. principal sites with 5R at one side (eg RFE)
		stype += 4;
	else if (stype >= 20 && stype <= 22) stype += 5;
	else if (stype == 30 || stype == 31 || stype == 35) stype += 5;
	else if (stype >= 56 && stype<= 58) stype += 6;
    else {
        std::cout<<"ERROR: add5RtoSite: illegal site type to add 5R to.\n";
        return;}
    // removes site from m_pah->m_siteMap (principal site)
    delSiteFromMap(st->type, st);
    // change site type
    st->type = (kmcSiteType) stype;
    // add site to m_pah->m_siteMap
    m_pah->m_siteMap[st->type].push_back(st);
    // update member C
    st->C1 = Carb1;
    st->C2 = Carb2;*/
}
//! Changes site type into combined site without R5 (e.g. RFE -> FE)
void PAHProcess::remR5fromSite(Spointer& st, Cpointer Carb1, Cpointer Carb2) {
    //int stype = (int) st->type;
	int index_change = -101;
	updateSites(st, Carb1, Carb2, index_change); 
}
//! Redraws 5 member rings in a zigzag
void PAHProcess::redrawR5(Spointer& st, Cpointer Carb1, Cpointer Carb2) {
	cpair CZZvec = get_vector(Carb2->coords, Carb2->C2->coords);
	cpair normvec = (norm_vector(Carb1->coords, Carb2->coords, Carb2->C2->coords));
	cpair sumvec = get_vector(Carb1->coords, Carb2->coords);
	cpair growvec = cross_vector(normvec, sumvec);
	removeC(Carb1->C2, false); removeC(Carb1->C2, false);
	cpair mpos = jumpToPos(Carb1->coords, CZZvec, 1.4);
	if (!checkHindrance_C_PAH(mpos)) {
		addC(Carb1, CZZvec, 1.4, true);
		convSiteType(st, Carb1, Carb2, ZZ);
		updateA(Carb1, 'H', growvec);
		updateA(Carb2, 'H', growvec);
		proc_G5R_ZZ(st, Carb1, Carb2);
	}
	else {
		cout << "WARNING. PAHProcess:redrawR5 did not find the ZZ carbon. Trying to add it anyways.\n";
		addC(Carb1, CZZvec, 1.4, true);
		convSiteType(st, Carb1, Carb2, ZZ);
		updateA(Carb1, 'H', growvec);
		updateA(Carb2, 'H', growvec);
		proc_G5R_ZZ(st, Carb1, Carb2);
		/*
		//Carb1->bondAngle1 = normAngle(Carb1->C1->bondAngle1-60);
		Cpointer C1_res, C2_res;
		//C1_res = addC(Carb1, normAngle(Carb1->bondAngle1 + 120), 0, 1.4);
		//C2_res = addC(C1_res, normAngle(Carb1->bondAngle1 - 90), normAngle(Carb1->bondAngle1 - 180), 1.4*pow(3, 0.5));
		updateA(Carb1->C1, Carb2->C2, 'H');
		// update sites and neighbours
		convSiteType(st, C1_res, C2_res, R5);
		Spointer S1, S2, S3, S4;
		// neighbours
		S1 = moveIt(st, -1);
		S2 = moveIt(st, 1);
		addR5toSite(S1, S1->C1, C1_res);
		addR5toSite(S2, C2_res, S2->C2);
		// update combined sites
		S3 = moveIt(S1, -1);
		S4 = moveIt(S2, 1);
		updateCombinedSites(st); // update resulting site
		updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
		updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
		// add ring counts
		m_pah->m_rings5_Lone++;
		addR5internal(C1_res,C2_res);*/
	}
}
int convSiteType_error_counter =0;
//! Changes site type into another site type
void PAHProcess::convSiteType(Spointer& st, Cpointer Carb1, Cpointer Carb2, kmcSiteType t) {
	//Checks for sites with similar number of R6s and R5s.
	int stype = (int)t;
	int prev_stype = (int)st->type;
	if (prev_stype == 9999){
		//Site being converted is an SPIRAL. Keep it as an SPIRAL.
		return;
		//stype = 9999;
		//t = (kmcSiteType) stype;
		//st->type = SPIRAL;
	}
	if (stype == 2004 || stype == 2014) {
		////////////////////////////////////////////////////////////
		//		Optimise before deciding which PAH it is.
		// change site type for the optimiser.
		st->type = (kmcSiteType)(stype);
		// update member C so the optimiser does not fail.
		st->C1 = Carb1;
		st->C2 = Carb2;
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		// change site type back to original site.
		st->type = (kmcSiteType)(prev_stype);
		///////////////////////////////////////////////////////////
		//There are two possible sites ZZACR5 and FEACR5FE. Decide which one.
		Cpointer Ccheck = Carb1->C2; Cpointer Ccheck2 = Ccheck->C2; Cpointer Ccheck3 = Carb2->C1->C1;	Cpointer Ccheck4 = Ccheck3->C2;
		bool decide_PAH = false;
		if (isR5internal(Ccheck, Ccheck2) || isR5internal(Ccheck3, Ccheck4)) decide_PAH = true;
		if(decide_PAH == true){
			stype = 2004;
		}
		else {
			stype = 2014;
		}
		t = (kmcSiteType) stype;
	}
	if (stype == 2104 || stype == 2114) {
		////////////////////////////////////////////////////////////
		//		Optimise before deciding which PAH it is.
		// change site type for the optimiser.
		st->type = (kmcSiteType)(stype);
		// update member C so the optimiser does not fail.
		st->C1 = Carb1;
		st->C2 = Carb2;
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		// change site type back to original site.
		st->type = (kmcSiteType)(prev_stype);
		///////////////////////////////////////////////////////////
		//There are two possible sites R5FEACR5 and ACR5R5R6. Decide which one.
		Cpointer Ccheck = Carb1; Cpointer Ccheck2 = Ccheck->C2; Cpointer Ccheck3 = Carb2->C1;	Cpointer Ccheck4 = Ccheck3->C2;
		bool decide_PAH = false;
		if (isR5internal(Ccheck, Ccheck2) || isR5internal(Ccheck3, Ccheck4)) decide_PAH = true;
		if(decide_PAH == true){
			stype = 2104;
		}
		else {
			stype = 2114;
		}
		t = (kmcSiteType) stype;
	}
	if (stype == 2105 || stype == 2115) {
		////////////////////////////////////////////////////////////
		//		Optimise before deciding which PAH it is.
		// change site type for the optimiser.
		st->type = (kmcSiteType)(stype);
		// update member C so the optimiser does not fail.
		st->C1 = Carb1;
		st->C2 = Carb2;
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		// change site type back to original site.
		st->type = (kmcSiteType)(prev_stype);
		///////////////////////////////////////////////////////////
		//There are two possible sites R5ZZACR5 and ACR5R5R6ZZ. Decide which one.
		Cpointer Ccheck = Carb1; Cpointer Ccheck2 = Ccheck->C2; Cpointer Ccheck3 = Carb2->C1;	Cpointer Ccheck4 = Ccheck3->C2;
		bool decide_PAH = false;
		if (isR5internal(Ccheck, Ccheck2) || isR5internal(Ccheck3, Ccheck4)) decide_PAH = true;
		if(decide_PAH == true){
			stype = 2105;
		}
		else {
			stype = 2115;
		}
		t = (kmcSiteType) stype;
	}
	if (stype == 2005 || stype == 2015) {
		////////////////////////////////////////////////////////////
		//		Optimise before deciding which PAH it is.
		// change site type for the optimiser.
		st->type = (kmcSiteType)(stype);
		// update member C so the optimiser does not fail.
		st->C1 = Carb1;
		st->C2 = Carb2;
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		// change site type back to original site.
		st->type = (kmcSiteType)(prev_stype);
		///////////////////////////////////////////////////////////
		//There are two possible sites ZZACR5 and FEACR5FE. Decide which one.
		Cpointer Ccheck = Carb1->C2;
		Cpointer Ccheck2 = Ccheck->C2;
		Cpointer Ccheck3 = Carb2->C1;
		Cpointer Ccheck4 = Ccheck3->C2;
		bool decide_PAH = false;
		if (isR5internal(Ccheck, Ccheck2) || isR5internal(Ccheck3, Ccheck4)) decide_PAH = true;
		if(decide_PAH == true){
			stype = 2005;
		}
		else {
			stype = 2015;
		}
		t = (kmcSiteType) stype;
	}
    // removes site from m_pah->m_siteMap (principal site)
    delSiteFromMap(st->type, st);
    // change site type
    st->type = t;
	if ( !checkSiteValid(st) && ( (int)t%10 >= 5 && (int)t%10 <= 9) ){
		//Assume an spiral has been formed
		st->type = SPIRAL;
		stype = 9999;
		t = (kmcSiteType) stype;
	}
	if (!checkSiteValid(st)) {
		cout << "ERROR. Invalid site convSiteType. This may be alright if the edge is unreactive but it may also be an error. \n";
		cout << "Sweep::PAHProcess::convSiteType. \n"; //SETBREAKPOINT
		printSites(st); 
		if (m_debug_pah) {
			ifstream  src("KMC_DEBUG/BEFORE.xyz");
			std::string filename = "KMC_DEBUG/BEFORE_convSiteType";
			filename.append(std::to_string(convSiteType_error_counter));
			filename.append(".xyz");
			ofstream dst(filename);
			dst << src.rdbuf();
			src.close();
			dst.close();
			cout<<"Saving file: "<< filename<<"\n";
		}
		std::string filename2 = "KMC_DEBUG/KMC_PAH_convSiteType_error_";
		filename2.append(std::to_string(convSiteType_error_counter));
		saveXYZ(filename2);
		cout<<"Saving file: "<< filename2<<".xyz\n";
		++convSiteType_error_counter;
		return;
		//saveXYZ("KMC_DEBUG/convSiteType_invalid");
	}
    // add site to m_pah->m_siteMap
    m_pah->m_siteMap[t].push_back(st);
    // update member C
    st->C1 = Carb1;
    st->C2 = Carb2;
}
//! Remove site
void PAHProcess::removeSite(Spointer& stt) {
    //remove site from siteMap for both elementary and combined site types
    delSiteFromMap(stt->type, stt);
    delSiteFromMap(stt->comb, stt);
    //then remove existence of site altogether
    m_pah->m_siteList.erase(stt);
}
//! Sets the number of counts of C and H
void PAHProcess::setCount(int CCount, int HCount) {
    m_pah->m_counts.first = CCount;
    m_pah->m_counts.second = HCount;
}
//! Adds C and H count
void PAHProcess::addCount(int C_in, int H_in) {
        m_pah->m_counts.first += C_in;
        m_pah->m_counts.second += H_in;
}

//! Update Processes
//! Sets third species bonded to C_1 to a species sp if it is a reactive surface carbon.
void PAHProcess::updateA(Cpointer C, char sp) {
    if(C->bridge) {
        C->A = 'C'; // a bridge is obviously connected to 3 other C atoms
        C->bondAngle2=normAngle(C->bondAngle1+120);
    }else {
        // Check if C is a reactive surface carbon by finding angle it makes with neighbours
        angletype x = normAngle(C->bondAngle1 - C->C1->bondAngle1);
        // Sets A if C is a reactive surface carbon
        if (x<0) {
            C->A = sp;
            C->bondAngle2=normAngle(C->bondAngle1+120);
        } else {
            C->A = 'C';
            C->bondAngle2=normAngle(C->bondAngle1-120);
        }
    }
}
//! Overload function, forces char to C atom.
void PAHProcess::updateA(Cpointer C, char sp, bool flag) {
	if (flag){
		C->A = sp;
		C->bondAngle2=normAngle(C->bondAngle1+120);
	} else updateA(C, sp);
}

//! Overload function, forces char to C atom and defines growth vector.
void PAHProcess::updateA(Cpointer C, char sp, cpair gro_vec) {
	C->A = sp;
	if (sp == 'C') C->growth_vector=std::make_tuple (0.0,0.0,0.0);
	else C->growth_vector=gro_vec;
	if (sp == 'M') addCount(1, 2);
}

//! Overload function, forces char to C atom and defines growth vector, modifies OBmol object.
void PAHProcess::updateA(Cpointer C, char sp, cpair gro_vec, OpenBabel::OBMol mol) {
	C->A = sp;
	if (sp == 'C') {
		C->growth_vector=std::make_tuple (0.0,0.0,0.0);
		OpenBabel::OBMolAtomIter cmol(mol);
		double tol = 2e-1;
		for(OpenBabel::OBMolAtomIter     a(mol); a; ++a) {
			double x_pos = a->GetX();
			double y_pos = a->GetY();
			double z_pos = a->GetZ();
			cpair temp = std::make_tuple(x_pos, y_pos, z_pos);
			if (getDistance_twoC(C->coords, temp) < tol){
				cmol = a;
				break;
			}
		}
		++cmol;
		mol.DeleteAtom(&*cmol);
	}
	else {
		C->growth_vector=gro_vec;
		cpair hloc = jumpToPos(C->coords, gro_vec, 1.085);
		OpenBabel::OBAtom *atom  = mol.NewAtom();
		atom->SetAtomicNum(1);
		atom->SetVector(std::get<0>(hloc),std::get<1>(hloc),std::get<2>(hloc));
	}
}

//! Overload function, updateA for all C from C_1 to C_2 inclusive
void PAHProcess::updateA(Cpointer C_1, Cpointer C_2, char spc) {
    // Iterate through each Carbon atom
    Cpointer now = C_1;
    Cpointer prev;
    // iterator will not cross bridges if C_1 is a bridge
    if(C_1->bridge) prev = C_1->C3;
    else prev = C_1->C1;
    // a count to check if updateA is looping uncontrollably
    unsigned int count = 0;
    do {
        updateA(now, spc); 
        Cpointer oldnow = now;
        now = moveCPointer(prev, now);
        prev = oldnow;
        count++;
		if (count == 1000){
			cout << "WARNING: Loop count has reached 1000 (Sweep::KMC_ARS::PAHProcess::updateA)\n";
			//saveDOT(std::string("KMC_DEBUG/error_count_100.dot"));
		}
        else if(count > 5000){
            cout << "ERROR: Loop count has reached more than 5000. Stopping simulation now.. (Sweep::KMC_ARS::PAHProcess::updateA)\n";
            std::ostringstream msg;
            msg << "ERROR: Possibility of function looping indefinitely."
                << " (Sweep::KMC_ARS::PAHProcess::updateA)";
			//saveDOT(std::string("KMC_DEBUG/error_count.dot"));
            throw std::runtime_error(msg.str());
            assert(false);
        }

    }
    while (prev != C_2);
}
//! Update all principal sites
void PAHProcess::updateSites() {
    //cout<<"Clearing Sites..\n";
    m_pah->m_siteList.clear(); //cout << "m_pah->m_siteList cleared!\n";
    m_pah->m_siteMap.clear();//cout << "m_pah->m_siteMap cleared!\n";
    //check if m_pah->m_cfirst is not bonded to C, if not iterate until reach one with H
    for(int i=0; i!=200;i++) {
        if (m_pah->m_cfirst->A != 'C') {
            break;
        } else {
            m_pah->m_clast = m_pah->m_cfirst;
            m_pah->m_cfirst = m_pah->m_cfirst->C2;
        }
        if (i==199) {
            std::cout << "ERROR: updateSites(): 200th loop reached without finding C-H\n";
            std::ostringstream msg;
            msg << "ERROR: Possibility of function looping indefinitely."
                << " (Sweep::KMC_ARS::PAHProcess::updateSites)";
            throw std::runtime_error(msg.str());
            assert(false);
            return;
        }
    }
    Cpointer prev = m_pah->m_cfirst;
    Cpointer now = m_pah->m_cfirst->C2;
    Cpointer siteC1=prev;
    kmcSiteType sitetype;
    unsigned int bulk = 0;
    do {
        if(now->A != 'C') {
            switch (bulk) {
            case 0:
                sitetype = FE; break;
            case 1:
                sitetype = ZZ; break;
            case 2:
                sitetype = AC; break;
            case 3:
                sitetype = BY5; break;
            case 4:
                sitetype = BY6; break;
            default:
                std::cout << "ERROR: updateSites(): more than 4 bulk carbon atoms detected\n"; break;
            }
			//if (now->A == 'M') sitetype = Methyl;
			//if (siteC1->A == 'M') sitetype = Methyl;
            Spointer end_of_siteList = m_pah->m_siteList.end();
            addSite(sitetype, siteC1, now, end_of_siteList);
            bulk = 0;
            siteC1 = now;
        }else {
            bulk++;
        }
        Cpointer oldnow = now;
        now = moveCPointer(prev, now);
        prev = oldnow;
    } while (prev != m_pah->m_cfirst);
    // Stores iterators into vectors
    //m_pah->m_siteMap.clear();
    //Spointer it;
    //for(it=m_pah->m_siteList.begin(); it!=m_pah->m_siteList.end(); it++) {
        //m_pah->m_siteMap[it->type].push_back(it);
    //}
    setCount(m_pah->m_counts.first, (int) m_pah->m_siteList.size() + 2 * numberOfMethyl() );
    //stericHindrance();
    //cout << "Principal Sites Updated, H count: "<< m_pah->m_counts[1] << "..\n";
}

//! Update all principal sites
void PAHProcess::updateSites(std::string site_list) {
    //cout<<"Clearing Sites..\n";
    m_pah->m_siteList.clear(); //cout << "m_pah->m_siteList cleared!\n";
    m_pah->m_siteMap.clear();//cout << "m_pah->m_siteMap cleared!\n";
	// create a vector from the string
    std::vector<std::string> siteList_strvec;
    Strings::split(site_list, siteList_strvec, std::string(","));
    // convert into vector of siteTypes
    std::vector<kmcSiteType> siteList_vec;
    for(size_t i=0; i<siteList_strvec.size(); i++) {
        kmcSiteType temp = kmcSiteType_str(siteList_strvec[i]);
        if(temp == Inv) {
            std::cout<<"ERROR: Starting Structure site List contains invalid site type"
                <<".. (PAHProcess::initialise)\n\n";
			std::cout<<"Site = " << siteList_strvec[i] << "\n";
            std::ostringstream msg;
            msg << "ERROR: Starting Structure site List contains invalid site type."
                << " (Sweep::KMC_ARS::PAHProcess::initialise)";
                throw std::runtime_error(msg.str());
                assert(false);
        }
        siteList_vec.push_back(temp);
	}

    //check if m_pah->m_cfirst is not bonded to C, if not iterate until reach one with H
    for(int i=0; i!=200;i++) {
        if (m_pah->m_cfirst->A != 'C') {
            break;
        } else {
            m_pah->m_clast = m_pah->m_cfirst;
            m_pah->m_cfirst = m_pah->m_cfirst->C2;
        }
        if (i==199) {
            std::cout << "ERROR: updateSites(): 200th loop reached without finding C-H\n";
            std::ostringstream msg;
            msg << "ERROR: Possibility of function looping indefinitely."
                << " (Sweep::KMC_ARS::PAHProcess::updateSites)";
            throw std::runtime_error(msg.str());
            assert(false);
            return;
        }
    }
    Cpointer prev = m_pah->m_cfirst;
    Cpointer now = m_pah->m_cfirst->C2;
    Cpointer siteC1=prev;
    kmcSiteType sitetype;
	for (unsigned int ii = 0; ii!=siteList_vec.size();++ii){
		int carbons_persite = (int)siteList_vec[ii] % 10;
		sitetype = siteList_vec[ii];
		for (int jj = 0; jj!= carbons_persite; ++jj){
			now = now->C2;
		}
		Spointer end_of_siteList = m_pah->m_siteList.end();
        addSite(sitetype, siteC1, now, end_of_siteList);
		siteC1 = now;
		now = siteC1->C2;
	}
    setCount(m_pah->m_counts.first, (int) m_pah->m_siteList.size() + 2 * numberOfMethyl() );
}

int updatesites_error_counter = 0;
//! Updates particular site
void PAHProcess::updateSites(Spointer& st, // site to be updated
                               Cpointer Carb1, Cpointer Carb2, // new C members
                               int bulkCchange) { // addition to number of bulk C in site
    // check if site type change is valid (as long as site still principal site)
	int stype = (int)st->type;
	int original_stype = stype;
	if (!checkSiteValid(stype)){
		cout << "ERROR: updateSites: Invalid site type before update.\n"; //SETBREAKPOINT
		cout << "Site type: " << kmcSiteName(st->type)<< "\n";
		cout << "Bulk change: " << bulkCchange<< "\n";
		std::ostringstream msg;
		msg << "ERROR: updateSites: Invalid site type before update.\n"
			<< "Site type: " << kmcSiteName(st->type)<< "\n"
			<< "Bulk change: " << bulkCchange<< "\n";
		//saveDOT("KMC_DEBUG/KMC_PAH_X_UPDATE_prev.dot");
		if (m_debug_pah){
			ifstream  src("KMC_DEBUG/BEFORE.xyz");
			std::string filename = "KMC_DEBUG/BEFORE_";
			filename.append(std::to_string(updatesites_error_counter));
			filename.append(".xyz");
			//filename.append(std::to_string(this->m_pah->m_parent->ID()));
			ofstream dst(filename);
			dst << src.rdbuf();
			src.close();
			dst.close();
			cout<<"Saving file: "<< filename<<"\n";
		}
		std::string filename2 = "KMC_DEBUG/KMC_PAH_X_UPDATE_prev_";
		filename2.append(std::to_string(updatesites_error_counter));
		//filename2.append("_");
		//filename2.append(std::to_string(this->m_pah->m_parent->ID()));
		saveXYZ(filename2);
		cout<<"Saving file: "<< filename2<<".xyz\n";
		filename2.append(".txt");
		savePAH_tofile(filename2);
		cout<<"Saving file: "<< filename2<<"\n";
		//throw std::runtime_error(msg.str());
		//assert(false);
		printSites(st);
		++updatesites_error_counter;
		return;
	}
	if ((stype + bulkCchange) < 0){
		cout << "ERROR: updateSites: Bulk C change invalid (Principal)\n"; //SETBREAKPOINT
		cout << "Trying to add " << bulkCchange << " bulk C to a " << kmcSiteName(st->type);
		cout << " (Sweep::KMC_ARS::PAHProcess::updateSites)";
		std::ostringstream msg;
		msg << "ERROR: Bulk C change invalid (Principal). Trying to add "
			<< bulkCchange << " bulk C to a " << kmcSiteName(st->type)
			<< " (Sweep::KMC_ARS::PAHProcess::updateSites)";
		//saveDOT("KMC_DEBUG/KMC_PAH_X_UPDATE.dot");
		if (m_debug_pah){
			ifstream  src("KMC_DEBUG/BEFORE.xyz");
			std::string filename = "KMC_DEBUG/BEFORE_";
			filename.append(std::to_string(updatesites_error_counter));
			filename.append(".xyz");
			//filename.append(std::to_string(this->m_pah->m_parent->ID()));
			ofstream dst(filename);
			dst << src.rdbuf();
			src.close();
			dst.close();
			cout<<"Saving file: "<< filename<<"\n";
		}
		std::string filename2 = "KMC_DEBUG/KMC_PAH_X_UPDATE_";
		filename2.append(std::to_string(updatesites_error_counter));
		saveXYZ(filename2);
		cout<<"Saving file: "<< filename2<<".xyz\n";
		filename2.append(".txt");
		savePAH_tofile(filename2);
		cout<<"Saving file: "<< filename2<<"\n";
		//saveXYZ("KMC_DEBUG/KMC_PAH_X_UPDATE");
		//throw std::runtime_error(msg.str());
		//assert(false);
		printSites(st);
		++updatesites_error_counter;
		return;
	}
	if (stype + bulkCchange == 2001) {
		stype = 501; bulkCchange = 0;
	}
	if (stype + bulkCchange == 500) {
		stype = 100; bulkCchange = 0;
		//removeR5internal(Carb1,Carb2);
	}
	if (stype + bulkCchange == 2102) {
		stype = 1002; bulkCchange = 0;
	}
	//Converts all site types with an extra digit 2014, 2114, ... into 2004, 2104, ... so the next section changes the index and then decides which resulting site is obtained.
	if (stype %100 >= 10 && stype != 9999) {
		stype -= 10;
	}
	if (stype + bulkCchange == 2004 || stype + bulkCchange == 2014) {
		//There are two possible sites ZZACR5 and FEACR5FE. Decide which one.
		Cpointer Ccheck = Carb1->C2; Cpointer Ccheck2 = Ccheck->C2; Cpointer Ccheck3 = Carb2->C1->C1;	Cpointer Ccheck4 = Ccheck3->C2;
		bool decide_PAH = false;
		if (isR5internal(Ccheck, Ccheck2) || isR5internal(Ccheck3, Ccheck4)) decide_PAH = true;
		if(decide_PAH == true){
			stype = 2004; bulkCchange = 0;
		}
		else {
			stype = 2014; bulkCchange = 0;
		}
	}
	if (stype + bulkCchange == 2104 || stype + bulkCchange == 2114) {
		////////////////////////////////////////////////////////////
		//		Optimise before deciding which PAH it is.
		// change site type for the optimiser. 
		st->type = (kmcSiteType)(stype + bulkCchange);
		// update member C so the optimiser does not fail.
		st->C1 = Carb1;
		st->C2 = Carb2;
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		// change site type back to original site.
		st->type = (kmcSiteType)(original_stype);
		///////////////////////////////////////////////////////////
		
		//There are two possible sites R5FEACR5 and ACR5R5R6. Decide which one.
		Cpointer Ccheck = Carb1; Cpointer Ccheck2 = Ccheck->C2; Cpointer Ccheck3 = Carb2->C1;	Cpointer Ccheck4 = Ccheck3->C2;
		bool decide_PAH = false;
		if (isR5internal(Ccheck, Ccheck2) || isR5internal(Ccheck3, Ccheck4)) decide_PAH = true;
		if(decide_PAH == true){
			stype = 2104; bulkCchange = 0;
		}
		else {
			stype = 2114; bulkCchange = 0;
		}
	}
	if (stype + bulkCchange == 2105 || stype + bulkCchange == 2115) {
		////////////////////////////////////////////////////////////
		//		Optimise before deciding which PAH it is.
		// change site type for the optimiser.
		st->type = (kmcSiteType)(stype + bulkCchange);
		// update member C so the optimiser does not fail.
		st->C1 = Carb1;
		st->C2 = Carb2;
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		// change site type back to original site.
		st->type = (kmcSiteType)(original_stype);
		///////////////////////////////////////////////////////////
		//There are two possible sites R5ZZACR5 and ACR5R5R6ZZ. Decide which one.
		Cpointer Ccheck = Carb1; Cpointer Ccheck2 = Ccheck->C2; Cpointer Ccheck3 = Carb2->C1;	Cpointer Ccheck4 = Ccheck3->C2;
		bool decide_PAH = false;
		if (isR5internal(Ccheck, Ccheck2) || isR5internal(Ccheck3, Ccheck4)) decide_PAH = true;
		if(decide_PAH == true){
			stype = 2105; bulkCchange = 0;
		}
		else {
			stype = 2115; bulkCchange = 0;
		}
	}
	if (stype + bulkCchange == 2005 || stype + bulkCchange == 2015) {
		////////////////////////////////////////////////////////////
		//		Optimise before deciding which PAH it is.
		// change site type for the optimiser.
		st->type = (kmcSiteType)(stype + bulkCchange);
		// update member C so the optimiser does not fail.
		st->C1 = Carb1;
		st->C2 = Carb2;
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		// change site type back to original site.
		st->type = (kmcSiteType)(original_stype);
		///////////////////////////////////////////////////////////
		//There are two possible sites ZZACR5 and FEACR5FE. Decide which one.
		Cpointer Ccheck = Carb1->C2;
		Cpointer Ccheck2 = Ccheck->C2;
		Cpointer Ccheck3 = Carb2->C1;
		Cpointer Ccheck4 = Ccheck3->C2;
		bool decide_PAH = false;
		if (isR5internal(Ccheck, Ccheck2) || isR5internal(Ccheck3, Ccheck4)) decide_PAH = true;
		if(decide_PAH == true){
			stype = 2005; bulkCchange = 0;
		}
		else {
			stype = 2015; bulkCchange = 0;
		}
	}
	if (stype + bulkCchange == 2013) {
		stype = 2003; bulkCchange = 0;
	}
	if (st->type == None)
	{
		//This means a None site is modified.
		cout << "Error. Added " << bulkCchange << " carbons to a None site.\n"; //SETBREAKPOINT
		bulkCchange = 0;
		if(m_debug_pah){
			ifstream  src("KMC_DEBUG/BEFORE.xyz");
			std::string filename = "KMC_DEBUG/BEFORE_";
			filename.append(std::to_string(updatesites_error_counter));
			filename.append(".xyz");
			ofstream dst(filename);
			dst << src.rdbuf();
			src.close();
			dst.close();
			cout<<"Saving file: "<< filename<<"\n";
		}
		std::string filename2 = "KMC_DEBUG/KMC_PAH_X_UPDATE_None_site_";
		filename2.append(std::to_string(updatesites_error_counter));
		//filename2.append(std::to_string(this->m_pah->m_parent->ID()));
		saveXYZ(filename2);
		cout<<"Saving file: "<< filename2<<".xyz\n";
		filename2.append(".txt");
		savePAH_tofile(filename2);
		cout<<"Saving file: "<< filename2<<"\n";
		//throw std::runtime_error(msg.str());
		//assert(false);
		printSites(st);
		++updatesites_error_counter;
		return;
	}
	if (st->type == SPIRAL)
	{
		//A spiral site is modified. Currently we cannot modify them back into reality.
		stype = 9999; 
		bulkCchange = 0;
	}
	if (!checkSiteValid(stype + bulkCchange)){
		//There are two options. An invalid transformation or an spiral formation.
		if ( (stype + bulkCchange)%10 >=5 && (stype + bulkCchange)%10 <=7 ){
			//The transformation just formed a 7 to 9 member site that will not react further. Transform it into a SPIRAL site.
			delSiteFromMap(st->type, st);
			stype = 9999; 
			bulkCchange = 0;
			//Optimise to avoid atoms in the spiral to overlap in positions.
			////////////////////////////////////////////////////////////
			//		Optimise before deciding which PAH it is.
			// change site type for the optimiser.
			st->type = (kmcSiteType)(stype + bulkCchange);
			// update member C so the optimiser does not fail.
			st->C1 = Carb1;
			st->C2 = Carb2;
			if (!m_pah->m_optimised){
				OpenBabel::OBMol mol = passPAH();
				mol = optimisePAH(mol);
				passbackPAH(mol);
			}
			// change site type back to original site.
			st->type = (kmcSiteType)(stype);
			///////////////////////////////////////////////////////////
			m_pah->m_siteMap[st->type].push_back(st);
		}
		else {
			//cout << "JP performed: " << pahpjp.getName();
			cout << "ERROR: updateSites: Created undefined site:\n"; //SETBREAKPOINT
			cout << "Type = " << stype << "\n";
			cout << "BulkC change = " << bulkCchange << "\n";
			/*cout << "ERROR: updateSites: Bulk C + Type gave invalid site\n";
			std::ostringstream msg;
			msg << "ERROR: Bulk C + Type gave invalid site. Trying to add "
				<< bulkCchange << " bulk C to a " << kmcSiteName(st->type)
				<< " (Sweep::KMC_ARS::PAHProcess::updateSites)";*/
			//saveDOT("KMC_DEBUG/KMC_PAH_X_UPDATE.dot");
			if (m_debug_pah){
				ifstream  src("KMC_DEBUG/BEFORE.xyz");
				std::string filename = "KMC_DEBUG/BEFORE_";
				filename.append(std::to_string(updatesites_error_counter));
				filename.append(".xyz");
				//filename.append(std::to_string(this->m_pah->m_parent->ID()));
				ofstream dst(filename);
				dst << src.rdbuf();
				src.close();
				dst.close();
				cout<<"Saving file: "<< filename<<"\n";
			}
			std::string filename2 = "KMC_DEBUG/KMC_PAH_X_UPDATE_BulkplusC_";
			filename2.append(std::to_string(updatesites_error_counter));
			//filename2.append("_");
			//filename2.append(std::to_string(this->m_pah->m_parent->ID()));
			saveXYZ(filename2);
			cout<<"Saving file: "<< filename2<<".xyz\n";
			filename2.append(".txt");
			savePAH_tofile(filename2);
			cout<<"Saving file: "<< filename2<<"\n";
			//saveXYZ("KMC_DEBUG/KMC_PAH_X_UPDATE");
			/*throw std::runtime_error(msg.str());
			assert(false);*/
			delSiteFromMap(st->type, st);
			st->type = None;
			m_pah->m_siteMap[st->type].push_back(st);
			st->C1 = Carb1;
			st->C2 = Carb2;
			printSites(st);
			++updatesites_error_counter;
			return;
		}
	}
	// removes site from m_pah->m_siteMap (principal site)
	delSiteFromMap(st->type, st);
	// change site type
	st->type = (kmcSiteType)(stype + bulkCchange);
	// add site to m_pah->m_siteMap
	m_pah->m_siteMap[st->type].push_back(st);
	// update member C
	st->C1 = Carb1;
	st->C2 = Carb2;
}

/*!
 * The two combined sites AC_FE3 and BY5_FE3 represent an armchair site and 5-member bay site, respectively, next to a combined FE3 site
 *
 * We need to update the free edge sites before we can update the rest of the sites; otherwise, we would come up with an incorrect count of AC_FE3 and BY5_FE3
 */
void PAHProcess::updateCombinedSites() {
	std::vector<int> indicesOfFE;
	std::vector<int> indicesOfNonFE;
	int loopIndex = 0;
	
	for(Spointer i=m_pah->m_siteList.begin(); i!= m_pah->m_siteList.end(); i++) {
		if(i->type == FE) {
			indicesOfFE.push_back(loopIndex);
		}
		else {
			indicesOfNonFE.push_back(loopIndex);
		}
		loopIndex++;
    }

	Spointer S1;
    S1 = m_pah->m_siteList.begin();

	for (unsigned int i = 0 ; i < indicesOfFE.size(); ++i) {
        Spointer S2 = moveIt(S1,indicesOfFE[i]);
		updateCombinedSites(S2);
	}

	for (unsigned int i = 0 ; i < indicesOfNonFE.size(); ++i) {
		Spointer S3 = moveIt(S1,indicesOfNonFE[i]);
        updateCombinedSites(S3);
	}
    //cout << "Combined Sites Updated..\n";
}

//! Update Combined Sites
//void PAHProcess::updateCombinedSites() {
//    // updateCombinedSites(st) for all sites
//    for(Spointer i=m_pah->m_siteList.begin(); i!= m_pah->m_siteList.end(); i++) {
//        updateCombinedSites(i);
//    }
//    //cout << "Combined Sites Updated..\n";
//}

//! Combined site for a particular site
void PAHProcess::updateCombinedSites(Spointer& st) {
    kmcSiteType ori_comb = st->comb;
    if(st->type != None) {
        delSiteFromMap(st->comb, st);
    }
	int s_type = (int)st->type;
	if (s_type >= 500) s_type = 500;
	cpair R5coords;
	Spointer S1 = moveIt(st,-1);
	Spointer S2 = moveIt(st,+1);
	bool check_left = true;
	bool check_right = true;
    switch(s_type) {
    case 0:
        // Check for FE3 (if there's FE on each side of the FE)
        if(moveIt(st,1)->type == FE && moveIt(st,1)->C1->A=='H' && moveIt(st,1)->C2->A=='H' && moveIt(st,-1)->type == FE && moveIt(st,-1)->C1->A=='H' && moveIt(st,-1)->C2->A=='H') {
            //if(st->C1->C1->C1->bridge || st->C2->C2->C2->bridge)
            st->comb = FE3;
            m_pah->m_siteMap[FE3].push_back(st);
            Spointer n1, n2;
            n1 = moveIt(st, -1);
            n2 = moveIt(st, 1);
            if(n1->comb != FE3) updateCombinedSites(n1);
            if(n2->comb != FE3) updateCombinedSites(n2);
            break;
        }
        //Check for FE2 (if only one FE is beside this FE)
        else if( (moveIt(st,1)->type == FE && moveIt(st,1)->C1->A=='H' && moveIt(st,1)->C2->A=='H') || (moveIt(st,-1)->type == FE && moveIt(st,-1)->C1->A=='H' && moveIt(st,-1)->C2->A=='H')){
            Spointer S1,S2;
            S1 = moveIt(st,-1); S2 = moveIt(st, 1);
            // Check if that FE is not a FE3
            if(S2->type == FE && moveIt(S2,1)->type != FE) {
                st->comb = FE2;
                m_pah->m_siteMap[FE2].push_back(st);
                // An FE2 site is a combined site where an FE site has an FE site only 
                // For example, for ZZ - FE - FE - ZZ, both of the FE sites has a combi
                //
                // Zig-zag oxidation reactions are based on the number of side-by-side 
                // If these reactions were based on the number of FE2 sites, we would o
                // So we can either calculate the rate based on the number of FE2 sites
                // or - as has been done here - remove half of the FE2 sites.
                //
                //if(S2->comb == FE2) delSiteFromMap(S2->comb, st); //Commented out to allow two FE2 sites next to each other. gl413
                //
                if(S2->comb != FE2) updateCombinedSites(S2);
            } else if(S1->type == FE && moveIt(S1,-1)->type != FE) {
                st->comb = FE2;
                m_pah->m_siteMap[FE2].push_back(st);
                //if(S1->comb == FE2) delSiteFromMap(S1->comb, st); //Commented out to allow two FE2 sites next to each other. gl413
                if(S1->comb != FE2) updateCombinedSites(S1);
            } else
                st->comb = None;
            break;
        }
        // Check for FE_HACA
        else if(moveIt(st,1)->type != FE && moveIt(st,1)->C1->A=='H' && moveIt(st,1)->C2->A=='H' && moveIt(st,-1)->type != FE && moveIt(st,-1)->C1->A=='H' && moveIt(st,-1)->C2->A=='H' && moveIt(st,1)->type != RFE && moveIt(st,-1)->type != RFE) {
            //if(st->C1->C1->bridge || st->C2->C2->bridge)
            st->comb = FE_HACA;
            m_pah->m_siteMap[FE_HACA].push_back(st);
            break;
        }
        else st->comb = None;
        break;
    case 2:
        // Check for AC_FE3
        if((moveIt(st,2)->comb==FE3 || moveIt(st,-2)->comb==FE3) && (moveIt(st,3)->comb!=FE3 || moveIt(st,-3)->comb!=FE3) && !st->C1->C2->bridge) {
			if (moveIt(st,2)->comb==FE3 && st->C1->A=='H' && st->C2->A=='H' && st->C2->C2->A=='H' && st->C2->C2->C2->A=='H' && st->C2->C2->C2->C2->A=='H') {
				st->comb = AC_FE3;
				m_pah->m_siteMap[AC_FE3].push_back(st);
				break;
			}
			else if (moveIt(st,-2)->comb==FE3 && st->C1->A=='H' && st->C2->A=='H' && st->C1->C1->A=='H' && st->C1->C1->C1->A=='H' && st->C1->C1->C1->C1->A=='H') {
				st->comb = AC_FE3;
				m_pah->m_siteMap[AC_FE3].push_back(st);
				break;
			}
			else st->comb = None;
        }else st->comb = None;
        break;
    case 3:
    	// Check for BY5_FE3
        if((moveIt(st,2)->comb==FE3 || moveIt(st,-2)->comb==FE3) && (moveIt(st,3)->comb!=FE3 || moveIt(st,-3)->comb!=FE3) && !st->C1->C2->bridge) {
			if (moveIt(st,2)->comb==FE3 && st->C1->A=='H' && st->C2->A=='H' && st->C2->C2->A=='H' && st->C2->C2->C2->A=='H' && st->C2->C2->C2->C2->A=='H') {
				st->comb = BY5_FE3;
				m_pah->m_siteMap[BY5_FE3].push_back(st);
				break;
			}
			else if (moveIt(st,-2)->comb==FE3 && st->C1->A=='H' && st->C2->A=='H' && st->C1->C1->A=='H' && st->C1->C1->C1->A=='H' && st->C1->C1->C1->C1->A=='H') {
				st->comb = BY5_FE3;
				m_pah->m_siteMap[BY5_FE3].push_back(st);
				break;
			}
			else st->comb = None;
        }else st->comb = None;
        break;
	case 500:
		//Check for R5R7
		if (isR7internal(st->C1,st->C1->C2) || isR7internal(st->C2->C1,st->C2) ){
			cpair R5coords_R7, R7coords;
			if (isR7internal(st->C1,st->C1->C2) ){
				R7coords = findR7internal(st->C1,st->C1->C2);
				R5coords_R7 = findR5internal(st->C1->C2, st->C1->C2->C2);
			} else{
				R7coords = findR7internal(st->C2->C1,st->C2);
				R5coords_R7 = findR5internal(st->C2->C1->C1, st->C2->C1);
			}
			double distR5R7 = getDistance_twoC(R7coords, R5coords_R7);
			if (distR5R7 < 2.6){
				st->comb = R5R7;
				m_pah->m_siteMap[R5R7].push_back(st);
				m_pah->m_R7loc.push_back(R7coords);
				break;
			}
			m_pah->m_R7loc.push_back(R7coords);
		}
		//Get R5 internal coordinates if they have not been cleared!
		if (st->type == ACR5) R5coords = findR5internal(st->C1->C2, st->C1->C2->C2);
		else if (st->type == R5R6){
			if ( isR5internal(st->C1->C1, st->C1)) {
				R5coords = findR5internal(st->C1->C1, st->C1);
				check_left = false;
			}
			else if ( isR5internal(st->C2, st->C2->C2) ) {
				R5coords = findR5internal(st->C2, st->C2->C2);
				check_right = false;
			}
			else {
				//R5 not found
				st->comb = None;
				break;
			}
		}
		else if ((int)st->type>2000 && (int)st->type<=2103){
			if ( isR5internal(st->C1->C2, st->C1->C2->C2) ) {
				R5coords = findR5internal(st->C1->C2, st->C1->C2->C2);
				check_right = false;
			}
			else if ( isR5internal(st->C2->C1->C1, st->C2->C1) ) {
				R5coords = findR5internal(st->C2->C1->C1, st->C2->C1);
				check_left = false;
			}
			else {
				//R5 not found
				st->comb = None;
				break;
			}
		}else{
			st->comb = None;
			break;
		}

		//Fundamental assumption: R5-R7 pairs cannot move away from each other!
		if (m_pah->m_R7loc.size()>=1){
			for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
				double distR5R7 = getDistance_twoC(*it, R5coords);
				if (distR5R7 < 3.1) {
					check_left = false;
					check_right = false;
					st->comb = None;
					break;
				}
			}
		}
		//Check for MIGR, MIGR2 OR R5R6_MIGR
		if (check_left){
			check_left = checkSiteMigration(S1,true);
		}
		if(check_right){
			check_right = checkSiteMigration(S2,false);
		}

		//Add back the R5 location
		m_pah->m_R5loc.push_back(R5coords);
		if (check_left && check_right && st->type==ACR5) {
			st->comb = MIGR2;
			m_pah->m_siteMap[MIGR2].push_back(st);
			break;
		}
		if ( (check_left || check_right) && ( (int)st->type>=2000 && (int)st->type <=2005)) {
			st->comb = MIGR;
			m_pah->m_siteMap[MIGR].push_back(st);
			break;
		}
		if ( (check_left || check_right) && (int)st->type==2103 ) {
			st->comb = MIGR;
			m_pah->m_siteMap[MIGR].push_back(st);
			break;
		}
		if ( (check_left || check_right) && st->type==R5R6) {
			st->comb = R5R6_MIGR;
			m_pah->m_siteMap[R5R6_MIGR].push_back(st);
			break;
		}
		st->comb = None;
		break;
    default:
        st->comb = None;
        break;
    }
    if((ori_comb == FE3) && (moveIt(st,1)->comb!=FE3 || moveIt(st,-1)->comb!=FE3)) {
        Spointer n1, n2;
        n1 = moveIt(st, -2);
        n2 = moveIt(st, 2);
        if((n1->type == AC || n1->type == BY5)) updateCombinedSites(n1);
        if((n2->type == AC || n2->type == BY5)) updateCombinedSites(n2);
    }
}

//! Combined site for a particular site, different rules during migration processes calls.
void PAHProcess::updateCombinedSitesMigration(Spointer& st) {
    kmcSiteType ori_comb = st->comb;
    if(st->type != None) {
        delSiteFromMap(st->comb, st);
    }
	int s_type = (int)st->type;
	if (s_type >= 500) s_type = 500;
	Spointer S1 = moveIt(st,-1);
	Spointer S2 = moveIt(st,+1);
	bool check_left = true;
	bool check_right = true;
	cpair R5coords;
	int steps=-99999;
    switch(s_type) {
    case 0:
        // Check for FE3 (if there's FE on each side of the FE)
        if(moveIt(st,1)->type == FE && moveIt(st,1)->C1->A=='H' && moveIt(st,1)->C2->A=='H' && moveIt(st,-1)->type == FE && moveIt(st,-1)->C1->A=='H' && moveIt(st,-1)->C2->A=='H') {
            //if(st->C1->C1->C1->bridge || st->C2->C2->C2->bridge)
            st->comb = FE3;
            m_pah->m_siteMap[FE3].push_back(st);
            Spointer n1, n2;
            n1 = moveIt(st, -1);
            n2 = moveIt(st, 1);
            if(n1->comb != FE3) updateCombinedSites(n1);
            if(n2->comb != FE3) updateCombinedSites(n2);
            break;
        }
        //Check for FE2 (if only one FE is beside this FE)
        else if( (moveIt(st,1)->type == FE && moveIt(st,1)->C1->A=='H' && moveIt(st,1)->C2->A=='H') || (moveIt(st,-1)->type == FE && moveIt(st,-1)->C1->A=='H' && moveIt(st,-1)->C2->A=='H')){
            Spointer S1,S2;
            S1 = moveIt(st,-1); S2 = moveIt(st, 1);
            // Check if that FE is not a FE3
            if(S2->type == FE && moveIt(S2,1)->type != FE) {
                st->comb = FE2;
                m_pah->m_siteMap[FE2].push_back(st);
                // An FE2 site is a combined site where an FE site has an FE site only 
                // For example, for ZZ - FE - FE - ZZ, both of the FE sites has a combi
                //
                // Zig-zag oxidation reactions are based on the number of side-by-side 
                // If these reactions were based on the number of FE2 sites, we would o
                // So we can either calculate the rate based on the number of FE2 sites
                // or - as has been done here - remove half of the FE2 sites.
                //
                //if(S2->comb == FE2) delSiteFromMap(S2->comb, st); //Commented out to allow two FE2 sites next to each other. gl413
                //
                if(S2->comb != FE2) updateCombinedSites(S2);
            } else if(S1->type == FE && moveIt(S1,-1)->type != FE) {
                st->comb = FE2;
                m_pah->m_siteMap[FE2].push_back(st);
                //if(S1->comb == FE2) delSiteFromMap(S1->comb, st); //Commented out to allow two FE2 sites next to each other. gl413
                if(S1->comb != FE2) updateCombinedSites(S1);
            } else
                st->comb = None;
            break;
        }
        // Check for FE_HACA
        else if(moveIt(st,1)->type != FE && moveIt(st,1)->C1->A=='H' && moveIt(st,1)->C2->A=='H' && moveIt(st,-1)->type != FE && moveIt(st,-1)->C1->A=='H' && moveIt(st,-1)->C2->A=='H' && moveIt(st,1)->type != RFE && moveIt(st,-1)->type != RFE) {
            //if(st->C1->C1->bridge || st->C2->C2->bridge)
            st->comb = FE_HACA;
            m_pah->m_siteMap[FE_HACA].push_back(st);
            break;
        }
        else st->comb = None;
        break;
    case 2:
        // Check for AC_FE3
        if((moveIt(st,2)->comb==FE3 || moveIt(st,-2)->comb==FE3) && (moveIt(st,3)->comb!=FE3 || moveIt(st,-3)->comb!=FE3) && !st->C1->C2->bridge) {
			if (moveIt(st,2)->comb==FE3 && st->C1->A=='H' && st->C2->A=='H' && st->C2->C2->A=='H' && st->C2->C2->C2->A=='H' && st->C2->C2->C2->C2->A=='H') {
				st->comb = AC_FE3;
				m_pah->m_siteMap[AC_FE3].push_back(st);
				break;
			}
			else if (moveIt(st,-2)->comb==FE3 && st->C1->A=='H' && st->C2->A=='H' && st->C1->C1->A=='H' && st->C1->C1->C1->A=='H' && st->C1->C1->C1->C1->A=='H') {
				st->comb = AC_FE3;
				m_pah->m_siteMap[AC_FE3].push_back(st);
				break;
			}
			else st->comb = None;
        }else st->comb = None;
        break;
    case 3:
    	// Check for BY5_FE3
        if((moveIt(st,2)->comb==FE3 || moveIt(st,-2)->comb==FE3) && (moveIt(st,3)->comb!=FE3 || moveIt(st,-3)->comb!=FE3) && !st->C1->C2->bridge) {
			if (moveIt(st,2)->comb==FE3 && st->C1->A=='H' && st->C2->A=='H' && st->C2->C2->A=='H' && st->C2->C2->C2->A=='H' && st->C2->C2->C2->C2->A=='H') {
				st->comb = BY5_FE3;
				m_pah->m_siteMap[BY5_FE3].push_back(st);
				break;
			}
			else if (moveIt(st,-2)->comb==FE3 && st->C1->A=='H' && st->C2->A=='H' && st->C1->C1->A=='H' && st->C1->C1->C1->A=='H' && st->C1->C1->C1->C1->A=='H') {
				st->comb = BY5_FE3;
				m_pah->m_siteMap[BY5_FE3].push_back(st);
				break;
			}
			else st->comb = None;
        }else st->comb = None;
        break;
	case 500:
		//Check for R5R7
		if (isR7internal(st->C1,st->C1->C2) || isR7internal(st->C2->C1,st->C2)){
			cpair R5coords_R7, R7coords;
			if (isR7internal(st->C1,st->C1->C2) ){
				R7coords = findR7internal(st->C1,st->C1->C2);
				Cpointer CR5_otherside_end = st->C1->C2;
				if (CR5_otherside_end->C2->A=='H') R5coords_R7 = endposR5internal(CR5_otherside_end, CR5_otherside_end->C2);
				else R5coords_R7 = endposR5internal(CR5_otherside_end, CR5_otherside_end->C2,true);
			} else{
				R7coords = findR7internal(st->C2->C1,st->C2);
				Cpointer CR5_otherside_end = st->C2->C1;
				if (CR5_otherside_end->A=='H') R5coords_R7 = endposR5internal(CR5_otherside_end->C1, CR5_otherside_end);
				else R5coords_R7 = endposR5internal(CR5_otherside_end->C1, CR5_otherside_end,true);
			}
			double distR5R7 = getDistance_twoC(R7coords, R5coords_R7);
			if (distR5R7 < 3.2){
				st->comb = R5R7;
				m_pah->m_siteMap[R5R7].push_back(st);
				m_pah->m_R7loc.push_back(R7coords);
				break;
			}
			m_pah->m_R7loc.push_back(R7coords);
		}
		//Get R5 internal coordinates if they have not been cleared!
		//For ACR5 both sides need to be checked. Probably it is an overkill since the site we come from is obviously allowed.
		if (st->type==ACR5){
			check_left = true;
			check_right = true;
		}
		else if (st->type == R5R6){
			int dir = coupledSiteDirection(st);
			if (dir==-1) check_left = false; //Coupled site is at the left, check right
			if (dir==1) check_right = false; //Coupled site is at the right, check left
			//dir = 0 check both
		}
		else if ( ((int)st->type>2000 && (int)st->type <= 2005) || ((int)st->type==2103)){ // Possibly site 2103 - R5ACR5 should be allowed to migrate too.
			bool b4;
			check_left = true;
			check_right = true;
			if ((int)st->type!=2002){
				for(unsigned int ii=0;ii!=m_pah->m_R5walker_sites.size();ii++){
					Spointer start_site = std::get<0>(m_pah->m_R5walker_sites[ii]);
					Spointer start_site_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
					int num_steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
					Spointer check_site = moveIt(start_site,num_steps);
					if (check_site==st){
						steps = num_steps;
						break;
					}
				}
				if (steps < 0) check_left = false;
				else if (steps > 0) check_right = false;
				else{
					if ((int)st->type==2103){
						int coupled_site_dir = coupledSiteDirection(st);
						if (coupled_site_dir == -1){
							check_left = false;
							R5coords = findR5internal(st->C2->C1->C1,st->C2->C1);
						}
						else if (coupled_site_dir == 1){
							check_right = false;
							R5coords = findR5internal(st->C1->C2,st->C1->C2->C2);
						}
						else{
							if (S1->type==R5) {
								check_left = false;
								R5coords = findR5internal(st->C2->C1->C1,st->C2->C1);
							}
							else if(S2->type==R5) {
								check_right = false;
								R5coords = findR5internal(st->C1->C2,st->C1->C2->C2);
							} else{
								std::cout << "Could not find coupled site of R5ACR5 site on updateCombinedSitesMigration." << std::endl;
								R5coords = findR5internal(st->C1->C2,st->C1->C2->C2);
							}
						}
					} else{
						if (isR5internal(st->C1->C2,st->C1->C2->C2)) {
							check_right = false;
							R5coords = findR5internal(st->C1->C2,st->C1->C2->C2);
						}
						else if (isR5internal(st->C2->C1->C1,st->C2->C1)) {
							check_left = false;
							R5coords = findR5internal(st->C2->C1->C1,st->C2->C1);
						}
						else{
							if (st->type==ACACR5 || st->type==ZZACR5){
								st->comb = None;
								break;
							}else {
								std::cout << "Could not find R5 on FEACR5 or similar site that has walker with 0 steps." << std::endl;
								st->comb = None;
								break;
							}
						}
					}
				}
			}
		}
		else {
			st->comb = None;
			break;
		}

		//Check for MIGR, MIGR2 OR R5R6_MIGR
		if (check_left) check_left = checkSiteMigration(S1,true);
		if (check_right) check_right = checkSiteMigration(S2,false);

		if (check_left && check_right && st->type==ACR5) {
			st->comb = MIGR2;
			m_pah->m_siteMap[MIGR2].push_back(st);
			break;
		}
		if ( (check_left || check_right) && ( ((int)st->type>=2000 && (int)st->type <=2005) || ((int)st->type==2103))) {
			st->comb = MIGR;
			m_pah->m_siteMap[MIGR].push_back(st);
			if (steps == 0) m_pah->m_R5loc.push_back(R5coords);
			break;
		}
		if ( (check_left || check_right) && st->type==R5R6) {
			st->comb = R5R6_MIGR;
			m_pah->m_siteMap[R5R6_MIGR].push_back(st);
			break;
		}
		st->comb = None;
		if (steps == 0 && (int)st->type>2002 && (int)st->type < 2099) m_pah->m_R5loc.push_back(R5coords);
		else if (steps == 0 && (int)st->type==2103) m_pah->m_R5loc.push_back(R5coords);
		break;
    default:
        st->comb = None;
        break;
    }
    if((ori_comb == FE3) && (moveIt(st,1)->comb!=FE3 || moveIt(st,-1)->comb!=FE3)) {
        Spointer n1, n2;
        n1 = moveIt(st, -2);
        n2 = moveIt(st, 2);
        if((n1->type == AC || n1->type == BY5)) updateCombinedSites(n1);
        if((n2->type == AC || n2->type == BY5)) updateCombinedSites(n2);
    }
}

//! Returns site x steps after i
Spointer PAHProcess::moveIt(Spointer i, int x) {
    Spointer temp = i;
    if (x>=0) {
        for(int a=0; a<x; a++) { // moves forward x times if x>0
            temp++;
            if(temp== m_pah->m_siteList.end()) {
                temp = m_pah->m_siteList.begin(); // if reach end of list go back to beginning
            }
        }
    } else {
        for(int a=0; a>x; a--) {
            if(temp == m_pah->m_siteList.begin()) {
                temp = m_pah->m_siteList.end();// if reach beginning of list go to end
                temp--;
            } else {
                 temp--; // moves backward x times if x<0
            }
        }
    }
    return temp;
}

// Structure change processes
//! Initialisation of structure, specified in derived classes (starting structures)
PAHStructure& PAHProcess::initialise_new(StartingStructure ss){
    // Pair of ring counts (R6 & R5) & chosen site list string
    std::string chosen;
    std::tuple <int, int, int> rings;
    intpair CH;
	int methyls;
	int rings_Embedded;
	std::list <cpair> IntCarbons;
    // Structure for Benzene
	std::string BENZENE_Sites = "FE,FE,FE,FE,FE,FE";
	auto BENZENE_Rings = std::make_tuple (1, 0, 0) ;
	intpair BENZENE_CH(6, 6);
	int BENZENE_methyls = 0;
	int BENZENE_RINGS_EMBEDDED = 0;
	// Structure for Toluene
	std::string TOLUENE_Sites = "FE,FE,FE,FE,FE,FE";
	auto TOLUENE_Rings = std::make_tuple (1, 0, 0) ;
	intpair TOLUENE_CH(7, 8);
	int TOLUENE_methyls = 1;
	int TOLUENE_RINGS_EMBEDDED = 0;
	// Structure for Naphthalene
	std::string NAPHTHALENE_Sites = "FE,FE,FE,ZZ,FE,FE,FE,ZZ";
	auto NAPHTHALENE_Rings = std::make_tuple(2, 0, 0);
	intpair NAPHTHALENE_CH(10, 0);
	int NAPHTHALENE_RINGS_EMBEDDED = 0;
	// Structure for Pyrene
	std::string PYRENE_Sites = "FE,ZZ,FE,FE,ZZ,FE,ZZ,FE,FE,ZZ";
	auto PYRENE_Rings = std::make_tuple (4, 0, 0);
	intpair PYRENE_CH(16, 10);
	int PYRENE_RINGS_EMBEDDED = 0;
	std::list <cpair> PYRENE_intCarbons;
	PYRENE_intCarbons.push_back(std::make_tuple(1.4*cos(-60.0*M_PI/180.0), 1.4*sin(-60.0*M_PI/180.0), 0.0));
	PYRENE_intCarbons.push_back(std::make_tuple(1.4*cos(-60.0*M_PI/180.0) + 1.4, 1.4*sin(-60.0*M_PI/180.0), 0.0));
	// Structure for Methylpyrene
	std::string METHYLPYRENE_Sites = "FE,ZZ,FE,FE,ZZ,FE,ZZ,FE,FE,ZZ";
	auto METHYLPYRENE_Rings = std::make_tuple (4, 0, 0);
	intpair METHYLPYRENE_CH(17, 12);
	int METHYLPYRENE_RINGS_EMBEDDED = 0;
	std::list <cpair> METHYLPYRENE_intCarbons;
	METHYLPYRENE_intCarbons.push_back(std::make_tuple(1.4*cos(-60.0*M_PI/180.0), 1.4*sin(-60.0*M_PI/180.0), 0.0));
	METHYLPYRENE_intCarbons.push_back(std::make_tuple(1.4*cos(-60.0*M_PI/180.0) + 1.4, 1.4*sin(-60.0*M_PI/180.0), 0.0));
	// Structure for Methylene Phenanthrene radical C15H9
	std::string MPHENANTHRENER_Sites = "FE,ZZ,FE,FE,R5R6,R5R6,FE,FE,ZZ";
	auto MPHENANTHRENER_Rings = std::make_tuple (3, 1, 0);
	intpair MPHENANTHRENER_CH(15, 9);
	int MPHENANTHRENER_RINGS_EMBEDDED = 0;
	std::list <cpair> MPHENANTHRENER_intCarbons;
	MPHENANTHRENER_intCarbons.push_back(std::make_tuple(1.4*cos(-60.0*M_PI/180.0), 1.4*sin(-60.0*M_PI/180.0), 0.0));
	MPHENANTHRENER_intCarbons.push_back(std::make_tuple(1.4*cos(-60.0*M_PI/180.0) + 1.4, 1.4*sin(-60.0*M_PI/180.0), 0.0));
	// Structure for BENZOPYRENE
	std::string BENZOPYRENE_Sites = "ZZ,FE,FE,ZZ,FE,ZZ,ZZ,FE,FE,FE,AC,FE";
	auto BENZOPYRENE_Rings = std::make_tuple (5, 0, 0);
	intpair BENZOPYRENE_CH(20, 12);
	int BENZOPYRENE_RINGS_EMBEDDED = 0;
	//Internal Carbons missing
	// Structure for Coronene
	std::string CORONENE_Sites = "FE,ZZ,FE,ZZ,FE,ZZ,FE,ZZ,FE,ZZ,FE,ZZ";
	auto CORONENE_Rings = std::make_tuple(7, 0, 0);
	intpair CORONENE_CH(24, 12);
	int CORONENE_RINGS_EMBEDDED = 0;
	// Test Structure
	std::string TEST_Sites = "FE,AC,FE,ZZ,RFE,R5,RAC,RFE,R5,RFE,FE,AC,FE,FE,FE,FE,BY5,AC,FE,BY5,FE,ZZ,FE,FE";
	auto TEST_Rings = std::make_tuple(8, 2, 0);
	intpair TEST_CH(18, 10);
	int TEST_methyls = 0;
	int TEST_RINGS_EMBEDDED = 0;
	//Internal Carbons missing
    // Choose structure
    switch(ss) {
    case BENZENE_C:
        chosen = BENZENE_Sites;
        rings = BENZENE_Rings;
        CH = BENZENE_CH;
		methyls = BENZENE_methyls;
		rings_Embedded = BENZENE_RINGS_EMBEDDED;
        break;
	case TOLUENE_C:
        chosen = TOLUENE_Sites;
        rings = TOLUENE_Rings;
        CH = TOLUENE_CH;
		methyls = TOLUENE_methyls;
		rings_Embedded = TOLUENE_RINGS_EMBEDDED;
        break;
    case NAPHTHALENE_C:
        chosen = NAPHTHALENE_Sites;
        rings = NAPHTHALENE_Rings;
        CH = NAPHTHALENE_CH;
		methyls = BENZENE_methyls;
		rings_Embedded = NAPHTHALENE_RINGS_EMBEDDED;
        break;
    case PYRENE_C:
        chosen = PYRENE_Sites;
        rings = PYRENE_Rings;
        CH = PYRENE_CH;
		methyls = BENZENE_methyls;
		rings_Embedded = PYRENE_RINGS_EMBEDDED;
		IntCarbons = PYRENE_intCarbons;
        break;
	case METHYLPYRENE_C:
        chosen = METHYLPYRENE_Sites;
        rings = METHYLPYRENE_Rings;
        CH = METHYLPYRENE_CH;
		methyls = TOLUENE_methyls;
		rings_Embedded = METHYLPYRENE_RINGS_EMBEDDED;
		IntCarbons = METHYLPYRENE_intCarbons;
        break;
	case MPHENANTHRENER_C:
        chosen = MPHENANTHRENER_Sites;
        rings = MPHENANTHRENER_Rings;
        CH = MPHENANTHRENER_CH;
		methyls = BENZENE_methyls;
		rings_Embedded = MPHENANTHRENER_RINGS_EMBEDDED;
		IntCarbons = MPHENANTHRENER_intCarbons;
        break;
    case BENZOPYRENE_C:
        chosen = BENZOPYRENE_Sites;
        rings = BENZOPYRENE_Rings;
        CH = BENZOPYRENE_CH;
		methyls = BENZENE_methyls;
		rings_Embedded = BENZOPYRENE_RINGS_EMBEDDED;
        break;
    case CORONENE_C:
        chosen = CORONENE_Sites;
        rings = CORONENE_Rings;
        CH = CORONENE_CH;
		methyls = BENZENE_methyls;
		rings_Embedded = CORONENE_RINGS_EMBEDDED;
        break;
    case TEST_STRUCT:
        chosen = TEST_Sites;
        rings = TEST_Rings;
        CH = TEST_CH;
		methyls = TEST_methyls;
		rings_Embedded = TEST_RINGS_EMBEDDED;
        break;
    default: 
            std::cout<<"ERROR: Starting Structure undefined.. (PAHProcess::initialise)\n\n";
            assert(false);
            abort();
    }
    // Create Structure
	return initialise_sitelist_string(chosen, std::get<0>(rings), std::get<1>(rings), rings_Embedded, std::get<2>(rings), rings_Embedded, std::get<0>(CH), std::get<1>(CH), methyls, IntCarbons);
    //printSites(m_pah->m_siteList.begin());
}

/*!
 * @param[in]    ss    Initialise the PAH structure corresponding to the starting structure, ss
 *
 * @return       Initialized PAH structure
 */
PAHStructure& PAHProcess::initialise(StartingStructure ss){
    if(m_pah == NULL) {
        PAHStructure* pah = new PAHStructure();
        m_pah = pah;
    }else if(m_pah->m_cfirst != NULL)
         m_pah->clear();
	std::list <cpair> P_intCarbons;
	switch(ss) {
        Cpointer newC;
        
    case BENZENE_C:
        //cout << "newC pointer created\n";
        // add first C atom
        m_pah->m_cfirst = addC();
		updateA(m_pah->m_cfirst, 'H', std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0) );
        //cout << "m_cfirst is " << m_cfirst << '\n';
        // adds next C atoms according to structure
        newC = addC(m_pah->m_cfirst, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
        //cout << "newC is " << newC << '\n';
        newC = addC(newC, std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(1.0,0.0,0.0) );
        //cout << "newC is " << newC << '\n';
        newC = addC(newC, std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0) );
        //cout << "newC is " << newC << '\n';
		newC = addC(newC, std::make_tuple(cos(-180.0 *M_PI / 180.0),sin(-180.0 *M_PI / 180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0) );
        //cout << "newC is " << newC << '\n';
        // adds the last C atom, with bond angle towards m_cfirst
		m_pah->m_clast = addC(newC, std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0), 1.4);
		updateA(m_pah->m_clast, 'H', std::make_tuple(-1.0,0.0,0.0) );
        //cout << "m_clast is " << m_clast << '\n';
        // closes structure
        connectToC(m_pah->m_clast, m_pah->m_cfirst);
        //cout << "m_clast and m_cfirst connected, updating A\n";
        // set C & H counts
        setCount(BENZENE_C, BENZENE_H);
        // set ring counts
        m_pah->m_rings = 1;
        m_pah->m_rings5_Lone = 0;
		m_pah->m_rings5_Embedded = 0;
		m_pah->m_rings7_Lone = 0;
		m_pah->m_rings7_Embedded = 0;
		m_pah->m_optimised = false;
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();
        //cout << "Benzene Initialised!\n";
		
		/*//cout << "newC pointer created\n";
        // add first C atom
        m_pah->m_cfirst = addC();
        //cout << "m_cfirst is " << m_cfirst << '\n';
        // adds next C atoms according to structure
        newC = addC(m_pah->m_cfirst, 0, 0, 1.4);
        //cout << "newC is " << newC << '\n';
        newC = addC(newC, -60, 0, 1.4);
        //cout << "newC is " << newC << '\n';
        newC = addC(newC, -120, 0, 1.4);
        //cout << "newC is " << newC << '\n';
        newC = addC(newC, -180, 0, 1.4);
        //cout << "newC is " << newC << '\n';
        // adds the last C atom, with bond angle towards m_cfirst
        m_pah->m_clast = addC(newC, 120, 60, 1.4);
        //cout << "m_clast is " << m_clast << '\n';
        // closes structure
        connectToC(m_pah->m_clast, m_pah->m_cfirst);
        //cout << "m_clast and m_cfirst connected, updating A\n";
        // update H atoms
        updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(BENZENE_C, BENZENE_H);
        // set ring counts
        m_pah->m_rings = 1;
        m_pah->m_rings5_Lone = 0;
		m_pah->m_rings5_Embedded = 0;
		m_pah->m_rings7_Lone = 0;
		m_pah->m_rings7_Embedded = 0;
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();*/
        //cout << "Benzene Initialised!\n";
        break;
	case TOLUENE_C:
        //cout << "newC pointer created\n";
        // add first C atom
        m_pah->m_cfirst = addC();
		updateA(m_pah->m_cfirst, 'H', std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0) );
        //cout << "m_cfirst is " << m_cfirst << '\n';
        // adds next C atoms according to structure
        newC = addC(m_pah->m_cfirst, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
        //cout << "newC is " << newC << '\n';
        newC = addC(newC, std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0), 1.4);
		updateA(newC, 'M', std::make_tuple(1.0,0.0,0.0) );
        //cout << "newC is " << newC << '\n';
        newC = addC(newC, std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0) );
        //cout << "newC is " << newC << '\n';
		newC = addC(newC, std::make_tuple(cos(-180.0 *M_PI / 180.0),sin(-180.0 *M_PI / 180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0) );
        //cout << "newC is " << newC << '\n';
        // adds the last C atom, with bond angle towards m_cfirst
		m_pah->m_clast = addC(newC, std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0), 1.4);
		updateA(m_pah->m_clast, 'H', std::make_tuple(-1.0,0.0,0.0) );
        //cout << "m_clast is " << m_clast << '\n';
        // closes structure
        connectToC(m_pah->m_clast, m_pah->m_cfirst);
        //cout << "m_clast and m_cfirst connected, updating A\n";
        // set C & H counts
        setCount(TOLUENE_C, TOLUENE_H);
        // set ring counts
        m_pah->m_rings = 1;
        m_pah->m_rings5_Lone = 0;
		m_pah->m_rings5_Embedded = 0;
		m_pah->m_rings7_Lone = 0;
		m_pah->m_rings7_Embedded = 0;
		m_pah->m_methyl_counts = 1;
		m_pah->m_optimised = false;
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();
        //cout << "TOLUENE Initialised!\n";
        break;
    case NAPHTHALENE_C:
        // add first C atom
        m_pah->m_cfirst = addC();
        // adds next C atoms according to structure
        newC = addC(m_pah->m_cfirst, 0, 0, 1.4);
        newC = addC(newC, 60, 0, 1.4);
        newC = addC(newC, 0, 0, 1.4);
        newC = addC(newC, -60, 0, 1.4);
        newC = addC(newC, -120, 0, 1.4);
        newC = addC(newC, -180, 0, 1.4);
        newC = addC(newC, -120, 0, 1.4);
        newC = addC(newC, -180, 0, 1.4);
        // adds the last C atom, with bond angle towards m_cfirst
        m_pah->m_clast = addC(newC, 120, 60, 1.4);
        // closes structure
        connectToC(m_pah->m_clast, m_pah->m_cfirst);
        // update H atoms
        updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(10, 8);
        // set ring counts
        m_pah->m_rings = 2;
        m_pah->m_rings5_Lone = 0;
		m_pah->m_rings5_Embedded = 0;
		m_pah->m_rings7_Lone = 0;
		m_pah->m_rings7_Embedded = 0;
		m_pah->m_optimised = false;
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();
        //cout << "Naphthalene Initialised!\n";
        break;
    case PYRENE_C:
        // add first C atom
        m_pah->m_cfirst = addC();
		updateA(m_pah->m_cfirst, 'H', std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0));
        // adds next C atoms according to structure
		newC = addC(m_pah->m_cfirst, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-60.0*M_PI/180.0),sin(-60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-60.0*M_PI/180.0),sin(-60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(1.0,0.0,0.0) );
		newC = addC(newC, std::make_tuple(cos(-120.0*M_PI/180.0),sin(-120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-180.0*M_PI/180.0),sin(-180.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-120.0*M_PI/180.0),sin(-120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-180.0*M_PI/180.0),sin(-180.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(120.0*M_PI/180.0),sin(120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(180.0*M_PI/180.0),sin(180.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(120.0*M_PI/180.0),sin(120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(-1.0,0.0,0.0) );
		newC = addC(newC, std::make_tuple(cos(60.0*M_PI/180.0),sin(60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0) );
        /*newC = addC(m_pah->m_cfirst, 0, 0, 1.4);
		newC = addC(newC, -60, 0, 1.4);
		newC = addC(newC, 0, 0, 1.4);
		newC = addC(newC, -60, 0, 1.4);
		newC = addC(newC, -120, 0, 1.4);
		newC = addC(newC, -180, 0, 1.4);
		newC = addC(newC, -120, 0, 1.4);
		newC = addC(newC, -180, 0, 1.4);
		newC = addC(newC, 120, 0, 1.4);
		newC = addC(newC, 180, 0, 1.4);
		newC = addC(newC, 120, 0, 1.4);
		newC = addC(newC, 60, 0, 1.4);*/
        // adds the last C atom, with bond angle towards m_cfirst
		m_pah->m_clast = addC(newC, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(m_pah->m_clast, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
        // closes structure
        connectToC(m_pah->m_clast, m_pah->m_cfirst);
        // update H atoms
        //updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(PYRENE_C, PYRENE_H);
        // set ring counts
        m_pah->m_rings = 4;
        m_pah->m_rings5_Lone = 0;
		m_pah->m_rings5_Embedded = 0;
		m_pah->m_rings7_Lone = 0;
		m_pah->m_rings7_Embedded = 0;
		m_pah->m_optimised = false;
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();
		//Set Internal Carbons list
		P_intCarbons.push_back(std::make_tuple(0.0, -2*1.4*cos(30.0*M_PI/180.0), 0.0));
		P_intCarbons.push_back(std::make_tuple(1.4, -2*1.4*cos(30.0*M_PI/180.0), 0.0));
		m_pah->m_InternalCarbons = P_intCarbons;
		//cout << "Pyrene Initialised!\n";
        break;
	case METHYLPYRENE_C:
        // add first C atom
        m_pah->m_cfirst = addC();
		updateA(m_pah->m_cfirst, 'H', std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0));
        // adds next C atoms according to structure
		newC = addC(m_pah->m_cfirst, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-60.0*M_PI/180.0),sin(-60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-60.0*M_PI/180.0),sin(-60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'M', std::make_tuple(1.0,0.0,0.0) );
		newC = addC(newC, std::make_tuple(cos(-120.0*M_PI/180.0),sin(-120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-180.0*M_PI/180.0),sin(-180.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-120.0*M_PI/180.0),sin(-120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-180.0*M_PI/180.0),sin(-180.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(120.0*M_PI/180.0),sin(120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(180.0*M_PI/180.0),sin(180.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(120.0*M_PI/180.0),sin(120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(-1.0,0.0,0.0) );
		newC = addC(newC, std::make_tuple(cos(60.0*M_PI/180.0),sin(60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0) );
        /*newC = addC(m_pah->m_cfirst, 0, 0, 1.4);
		newC = addC(newC, -60, 0, 1.4);
		newC = addC(newC, 0, 0, 1.4);
		newC = addC(newC, -60, 0, 1.4);
		newC = addC(newC, -120, 0, 1.4);
		newC = addC(newC, -180, 0, 1.4);
		newC = addC(newC, -120, 0, 1.4);
		newC = addC(newC, -180, 0, 1.4);
		newC = addC(newC, 120, 0, 1.4);
		newC = addC(newC, 180, 0, 1.4);
		newC = addC(newC, 120, 0, 1.4);
		newC = addC(newC, 60, 0, 1.4);*/
        // adds the last C atom, with bond angle towards m_cfirst
		m_pah->m_clast = addC(newC, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(m_pah->m_clast, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
        // closes structure
        connectToC(m_pah->m_clast, m_pah->m_cfirst);
        // update H atoms
        //updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(METHYLPYRENE_C, METHYLPYRENE_H);
        // set ring counts
        m_pah->m_rings = 4;
        m_pah->m_rings5_Lone = 0;
		m_pah->m_rings5_Embedded = 0;
		m_pah->m_rings7_Lone = 0;
		m_pah->m_rings7_Embedded = 0;
		m_pah->m_methyl_counts = 1;
		m_pah->m_optimised = false;
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();
		//Set Internal Carbons list
		P_intCarbons.push_back(std::make_tuple(0.0, -2*1.4*cos(30.0*M_PI/180.0), 0.0));
		P_intCarbons.push_back(std::make_tuple(1.4, -2*1.4*cos(30.0*M_PI/180.0), 0.0));
		m_pah->m_InternalCarbons = P_intCarbons;
		//cout << "Methylpyrene Initialised!\n";
        break;
	case MPHENANTHRENER_C:
        // add first C atom
        m_pah->m_cfirst = addC();
		updateA(m_pah->m_cfirst, 'H', std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0));
        // adds next C atoms according to structure
		newC = addC(m_pah->m_cfirst, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-58.3*M_PI/180.0),sin(-58.3*M_PI/180.0),0.0), 1.426);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-7.83*M_PI/180.0),sin(-7.83*M_PI/180.0),0.0), 1.408);
		updateA(newC, 'H', std::make_tuple(cos(51.95 *M_PI / 180.0),sin(51.95 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-67.63*M_PI/180.0),sin(-67.63*M_PI/180.0),0.0), 1.415);
		updateA(newC, 'H', std::make_tuple(cos(-6.43 *M_PI / 180.0),sin(-6.43 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-125.08*M_PI/180.0),sin(-125.08*M_PI/180.0),0.0), 1.413);
		updateA(newC, 'H', std::make_tuple(cos(-67.20 *M_PI / 180.0),sin(-67.20 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(171.15*M_PI/180.0),sin(171.15*M_PI/180.0),0.0), 1.387);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-142.56*M_PI/180.0),sin(-142.56*M_PI/180.0),0.0), 1.510);
		updateA(newC, 'H', std::make_tuple(cos(-90.0 *M_PI / 180.0),sin(-90.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(142.62*M_PI/180.0),sin(142.62*M_PI/180.0),0.0), 1.510);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		addR5internal(newC->C1, newC);
		newC = addC(newC, std::make_tuple(cos(-171.23*M_PI/180.0),sin(-171.23*M_PI/180.0),0.0), 1.387);
		updateA(newC, 'H', std::make_tuple(cos(-112.84 *M_PI / 180.0),sin(-112.84 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(125.09*M_PI/180.0),sin(125.09*M_PI/180.0),0.0), 1.412);
		updateA(newC, 'H', std::make_tuple(cos(-173.59 *M_PI / 180.0),sin(173.59 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(67.64*M_PI/180.0),sin(67.64*M_PI/180.0),0.0), 1.415);
		updateA(newC, 'H', std::make_tuple(cos(128.02*M_PI / 180.0),sin(128.02*M_PI / 180.0),0.0) );
        // adds the last C atom, with bond angle towards m_cfirst
		m_pah->m_clast = addC(newC, std::make_tuple(cos(7.849*M_PI/180.0),sin(7.849*M_PI/180.0),0.0), 1.408);
		updateA(m_pah->m_clast, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
        // closes structure
        connectToC(m_pah->m_clast, m_pah->m_cfirst);
        // update H atoms
        //updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(MPHENANTHRENER_C, MPHENANTHRENER_H);
        // set ring counts
        m_pah->m_rings = 3;
        m_pah->m_rings5_Lone = 1;
		m_pah->m_rings5_Embedded = 0;
		m_pah->m_rings7_Lone = 0;
		m_pah->m_rings7_Embedded = 0;
		m_pah->m_optimised = false;
        // update all sites and combined sites
        updateSites("FE,ZZ,FE,FE,R5R6,R5R6,FE,FE,ZZ");
        updateCombinedSites();
		//Set Internal Carbons list
		P_intCarbons.push_back(std::make_tuple(0.005, -2.35836, 0.0));
		P_intCarbons.push_back(std::make_tuple(1.39445, -2.35805, 0.0));
		m_pah->m_InternalCarbons = P_intCarbons;
		//cout << "MPHENANTHRENER Initialised!\n";
        break;
    case BENZOPYRENE_C:
        // add first C atom
        m_pah->m_cfirst = addC();
        // adds next C atoms according to structure
		newC = addC(m_pah->m_cfirst, 0, 0, 1.4);
		newC = addC(newC, -60, 0, 1.4);
		newC = addC(newC, 0, 0, 1.4);
		newC = addC(newC, -60, 0, 1.4);
		newC = addC(newC, -120, 0, 1.4);
		newC = addC(newC, -180, 0, 1.4);
		newC = addC(newC, -120, 0, 1.4);
		newC = addC(newC, -180, 0, 1.4);
		newC = addC(newC, 120, 0, 1.4);
		newC = addC(newC, 180, 0, 1.4);
		newC = addC(newC, 120, 0, 1.4);
		newC = addC(newC, 180, 0, 1.4);
		newC = addC(newC, 120, 0, 1.4);
		newC = addC(newC, 60, 0, 1.4);
		newC = addC(newC, 0, 0, 1.4);
		newC = addC(newC, -60, 0, 1.4);
        // adds the last C atom, with bond angle towards m_cfirst
		m_pah->m_clast = addC(newC, 0, 60, 1.4);
        // closes structure
        connectToC(m_pah->m_clast, m_pah->m_cfirst);
        // update H atoms
        updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(BENZOPYRENE_C, BENZOPYRENE_H);
        // set ring counts
        m_pah->m_rings = 5;
        m_pah->m_rings5_Lone = 0;
		m_pah->m_rings5_Embedded = 0;
		m_pah->m_rings7_Lone = 0;
		m_pah->m_rings7_Embedded = 0;
		m_pah->m_optimised = false;
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();
        //cout << "Benzopyrene Initialised!\n";
        break;
	case CORONENE_C:
		// "FE,ZZ,FE,ZZ,FE,ZZ,FE,ZZ,FE,ZZ,FE,ZZ"
		// add first C atom
        m_pah->m_cfirst = addC();
		updateA(m_pah->m_cfirst, 'H', std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0));
        // adds next C atoms according to structure
		newC = addC(m_pah->m_cfirst, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-60.0*M_PI/180.0),sin(-60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-60.0*M_PI/180.0),sin(-60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(1.0,0.0,0.0) );
		newC = addC(newC, std::make_tuple(cos(-120.0*M_PI/180.0),sin(-120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-60.0*M_PI/180.0),sin(-60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(1.0,0.0,0.0) );
		newC = addC(newC, std::make_tuple(cos(-120.0*M_PI/180.0),sin(-120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-180.0*M_PI/180.0),sin(-180.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-120.0*M_PI/180.0),sin(-120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-60.0 *M_PI / 180.0),sin(-60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(-180.0*M_PI/180.0),sin(-180.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(120.0*M_PI/180.0),sin(120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(180.0*M_PI/180.0),sin(180.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(-120.0 *M_PI / 180.0),sin(-120.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(120.0*M_PI/180.0),sin(120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(-1.0,0.0,0.0) );
		newC = addC(newC, std::make_tuple(cos(60.0*M_PI/180.0),sin(60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
		newC = addC(newC, std::make_tuple(cos(120.0*M_PI/180.0),sin(120.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(-1.0,0.0,0.0) );
		newC = addC(newC, std::make_tuple(cos(60.0*M_PI/180.0),sin(60.0*M_PI/180.0),0.0), 1.4);
		updateA(newC, 'H', std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0) );
        // adds the last C atom, with bond angle towards m_cfirst
		m_pah->m_clast = addC(newC, std::make_tuple(1.0,0.0,0.0), 1.4);
		updateA(m_pah->m_clast, 'C', std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0) );
        // closes structure
        connectToC(m_pah->m_clast, m_pah->m_cfirst);
        // update H atoms
        //updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(CORONENE_C, CORONENE_H);
        // set ring counts
        m_pah->m_rings = 7;
        m_pah->m_rings5_Lone = 0;
		m_pah->m_rings5_Embedded = 0;
		m_pah->m_rings7_Lone = 0;
		m_pah->m_rings7_Embedded = 0;
		m_pah->m_optimised = false;
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();
		//Set Internal Carbons list
		P_intCarbons.push_back(std::make_tuple(0.0, -2*1.4*cos(30.0*M_PI/180.0), 0.0));
		P_intCarbons.push_back(std::make_tuple(1.4, -2*1.4*cos(30.0*M_PI/180.0), 0.0));
		P_intCarbons.push_back(std::make_tuple(-0.7, -3*1.4*cos(30.0*M_PI/180.0), 0.0));
		P_intCarbons.push_back(std::make_tuple(2.1, -3*1.4*cos(30.0*M_PI/180.0), 0.0));
		P_intCarbons.push_back(std::make_tuple(0.0, -4*1.4*cos(30.0*M_PI/180.0), 0.0));
		P_intCarbons.push_back(std::make_tuple(1.4, -4*1.4*cos(30.0*M_PI/180.0), 0.0));
		m_pah->m_InternalCarbons = P_intCarbons;
		//cout << "Pyrene Initialised!\n";
        break;
     default: 
            std::cout<<"ERROR: Starting Structure undefined.. (PAHProcess::initialise)\n\n";
            assert(false);
            abort();
    }
    return *m_pah;
}

/*!
 * Initialization of PAH structure from an file. Allows debugging cases for any structure.
 *
 * @param[in]    siteList_str    	Comma-separated string of site types
 * @param[in]    R6_num          	Number of 6-member rings
 * @param[in]    R5_num_Lone     	Number of edge 5-member rings
 * @param[in]    R5_num_Embedded    Number of embedded 5-member rings
 * @param[in]    R7_num_Lone     	Number of edge 7-member rings
 * @param[in]    R7_num_Embedded    Number of embedded 7-member rings
 * @param[in] 	 edgeCarbons		Vector of strings space separated containing edge carbons coordinates, bridge, A, growth vector coordinates, in the way that sitelist runs.
 * @param[in] 	 internalCarbons	Vector of strings space separated containing internal carbon coordinates
 * @param[in] 	 R5_locs			Vector of strings space separated containing R5 positions
 * @param[in] 	 R7_locs			Vector of strings space separated containing R7 positions 
 *
 * @return       Initialized PAH structure
 */
PAHStructure& PAHProcess::initialise_fromfile(){
    if(m_pah == NULL) {
        PAHStructure* pah = new PAHStructure();
        m_pah = pah;
    }else if(m_pah->m_cfirst != NULL)
        m_pah->clear();
	
	ifstream src("InceptedPAH.inx");
	std::string line;
	std::string siteList_str;
	int R6_num, R5_num_Lone, R5_num_Embedded, R7_num_Lone, R7_num_Embedded;
	std::vector<std::string> edgeCarbons;
	std::vector<std::string> internalCarbons;
	std::vector<std::string> R5_locs;
	std::vector<std::string> R7_locs;
	if (src.is_open()){
		std::getline(src,line); //Read first line. Not needed here.
		std::getline(src,line); //Read second line. This contains the site list.
		line.erase(line.end()-1, line.end());
		siteList_str.assign(line);
		std::getline(src,line); //Read third line. This contains the ring information.
		line.erase(line.end()-1, line.end());
		std::stringstream ss(line);
		int ii = 0;
		while(std::getline(ss, line, ' ')){
			if (ii == 0) R6_num = std::stoi(line);
			if (ii == 1) R5_num_Lone = std::stoi(line);
			if (ii == 2) R5_num_Embedded = std::stoi(line);
			if (ii == 3) R7_num_Lone = std::stoi(line);
			if (ii == 4) R7_num_Embedded = std::stoi(line);
			ii++;
		}
		std::getline(src,line); //Read fourth line. This should say Coordinates
		
		
		while(std::getline (src,line) && line != "Internal\r"){
			line.erase(line.end()-1, line.end());
			edgeCarbons.push_back(line);
		}
		
		while(std::getline (src,line) && line != "R5_locs\r"){
			line.erase(line.end()-1, line.end());
			internalCarbons.push_back(line);
		}
		
		while(std::getline (src,line) && line != "R7_locs\r"){
			line.erase(line.end()-1, line.end());
			R5_locs.push_back(line);
		}
		
		while(std::getline (src,line)){
			line.erase(line.end()-1, line.end());
			R7_locs.push_back(line);
		}
	}
	else{
		cout << "Could not open file InceptedPAH.inx.\n";
		std::ostringstream msg;
		msg << "ERROR: Could not open file InceptedPAH.inx."
                << " (Sweep::KMC_ARS::PAHProcess::initialise_fromfile)";
                throw std::runtime_error(msg.str());
                assert(false);
	}
	src.close();
	
    // create a vector from the string
    std::vector<std::string> siteList_strvec;
    Strings::split(siteList_str, siteList_strvec, std::string(","));
    // convert into vector of siteTypes
    std::vector<kmcSiteType> siteList_vec;
	std::vector<int> carb_siteList_vec;
    for(size_t i=0; i<siteList_strvec.size(); i++) {
		bool defined_length_site = false;
		int number_carbs;
		kmcSiteType temp = kmcSiteType_str(siteList_strvec[i]);
        if(temp == Inv) {
			//Check if site has defined length as last char
			int last_char = siteList_strvec[i].back() - '0';
			if(last_char >= 1 && last_char<=9){
				number_carbs = last_char-1;
				std::string my_string = siteList_strvec[i];
				my_string.pop_back();
				temp = kmcSiteType_str(my_string);
				defined_length_site = true;
			}
			if(temp == Inv) {
				std::cout<<"ERROR: Starting Structure site List contains invalid site type"
					<<".. (PAHProcess::initialise)\n\n";
				std::cout<<"Site = " << siteList_strvec[i] << "\n";
				std::ostringstream msg;
				msg << "ERROR: Starting Structure site List contains invalid site type."
					<< " (Sweep::KMC_ARS::PAHProcess::initialise)";
					throw std::runtime_error(msg.str());
					assert(false);
			}
        }
        siteList_vec.push_back(temp);
		if (!defined_length_site){
			if (i==0) number_carbs = (int)temp %10 + 2; // First site accounts for m_cfirst
			else if (i==siteList_strvec.size()-1) number_carbs = (int)temp %10; // Last site closes PAH on m_cfirst.
			else number_carbs = (int)temp %10 + 1;
		} 
		carb_siteList_vec.push_back(number_carbs);
    }
	createPAH_fromfile(siteList_vec, carb_siteList_vec, R6_num, R5_num_Lone, R5_num_Embedded, R7_num_Lone, R7_num_Embedded, edgeCarbons, internalCarbons, R5_locs, R7_locs);
    return *m_pah;
}

// Create Structure from a file 
void PAHProcess::createPAH_fromfile(std::vector<kmcSiteType>& vec, std::vector<int>& carb_vec, int R6, int R5_Lone, int R5_Embedded, int R7_Lone, int R7_Embedded, std::vector<std::string> edCarbons, std::vector<std::string> inCarbs, std::vector<std::string> R5loc, std::vector<std::string> R7loc) {
	Cpointer newC;
	// start drawing..
	// This loops through sites.
	int carb_vec_sum = 0;
	Cpointer site_carb_1;
    for(size_t i=0; i<vec.size(); i++) {
		int vec_carb_number = carb_vec[i];
		
		//This loops through carbons
		for (int j=carb_vec_sum; j<vec_carb_number+carb_vec_sum; j++) {
			cpair carbon_coords;
			std::string carbon_string;
			carbon_string.assign(edCarbons[j]);
			std::vector<std::string> carbon_vector;
			Strings::split(carbon_string, carbon_vector, std::string(" "));
			carbon_coords = std::make_tuple(std::stod(carbon_vector[0]), std::stod(carbon_vector[1]), std::stod(carbon_vector[2]));
			
			if (i == 0 && j == 0){
				//Add first carbon.
				newC = addC();
				m_pah->m_cfirst = newC;
				m_pah->m_clast = NULLC;
				moveC(newC, carbon_coords);
				if (getDistance_twoC(newC, newC->C1)>1.8){
					std::cout << "Warning. Adding carbon from file with bond length > 1.8." << std::endl;
					std::cout << "Carbon numbers: " << j << " and " << j-1 << "." << std::endl;
				}
				if (std::stoi(carbon_vector[3]) == 1) newC->bridge = true;
				else newC->bridge = false;
				char cstr[carbon_vector[4].size()+1];
				std::strcpy(cstr, carbon_vector[4].c_str());
				updateA(newC, *cstr, std::make_tuple(std::stod(carbon_vector[5]), std::stod(carbon_vector[6]), std::stod(carbon_vector[7])));
			}
			else if ( !checkHindrance_C_PAH(carbon_coords)){
				//Only add carbons that are not occupied. Handles bridged atoms.
				newC = addC(newC, std::make_tuple(1.0,0.0,0.0), 7.17);
				moveC(newC, carbon_coords);
				if (std::stoi(carbon_vector[3]) == 1) newC->bridge = true;
				else newC->bridge = false;
				char cstr[carbon_vector[4].size()+1];
				std::strcpy(cstr, carbon_vector[4].c_str());
				updateA(newC, *cstr, std::make_tuple(std::stod(carbon_vector[5]), std::stod(carbon_vector[6]), std::stod(carbon_vector[7])));
				if (newC->bridge && newC->C1->bridge) {
					newC->C3 = newC->C1;
					newC->C1->C3 = newC;
					newC->C1 = NULLC;
					newC->C3->C2 = NULLC;
				}
			}
			else {
				//Adding a carbon for the second time. This means bridge.
				Cpointer bridgedC = findC(carbon_coords);
				if (!newC->bridge){
					newC->C2 = bridgedC;
					bridgedC->C1 = newC;
					newC = bridgedC;
				}
				else {
					newC = bridgedC;
				}
			}
			if (j == 0) site_carb_1 = newC;
		}
		carb_vec_sum += vec_carb_number;
		Cpointer site_carb_2 = newC;
		if (i==vec.size()-1) site_carb_2 = m_pah->m_cfirst;
		addSite(vec[i], site_carb_1, site_carb_2);
		site_carb_1 = site_carb_2;
	}
	connectToC(newC, m_pah->m_cfirst);
	m_pah->m_clast = newC;
	newC = newC->C2;

    m_pah->m_rings = R6;
	m_pah->m_rings5_Lone = R5_Lone;
	m_pah->m_rings5_Embedded = R5_Embedded;
	m_pah->m_rings7_Lone = R7_Lone;
	m_pah->m_rings7_Embedded = R7_Embedded;
	m_pah->m_optimised = false;
	
	for (size_t i = 0; i<inCarbs.size(); i++){
		cpair carbon_coords;
		std::string carbon_string;
		carbon_string.assign(inCarbs[i]);
		std::vector<std::string> carbon_vector;
		Strings::split(carbon_string, carbon_vector, std::string(" "));
		carbon_coords = std::make_tuple(std::stod(carbon_vector[0]), std::stod(carbon_vector[1]), std::stod(carbon_vector[2]));
		m_pah->m_InternalCarbons.push_back(carbon_coords);
	}
	for (size_t i = 0; i<R5loc.size(); i++){
		cpair carbon_coords;
		std::string carbon_string;
		carbon_string.assign(R5loc[i]);
		std::vector<std::string> carbon_vector;
		Strings::split(carbon_string, carbon_vector, std::string(" "));
		carbon_coords = std::make_tuple(std::stod(carbon_vector[0]), std::stod(carbon_vector[1]), std::stod(carbon_vector[2]));
		m_pah->m_R5loc.push_back(carbon_coords);
	}
	for (size_t i = 0; i<R7loc.size(); i++){
		cpair carbon_coords;
		std::string carbon_string;
		carbon_string.assign(R7loc[i]);
		std::vector<std::string> carbon_vector;
		Strings::split(carbon_string, carbon_vector, std::string(" "));
		carbon_coords = std::make_tuple(std::stod(carbon_vector[0]), std::stod(carbon_vector[1]), std::stod(carbon_vector[2]));
		m_pah->m_R7loc.push_back(carbon_coords);
	}
	m_pah->m_methyl_counts = numberOfMethyl();
	//int totalC_num = 2 * m_pah->m_rings + (CarbonListSize() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded) / 2 + numberOfBridges() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded + 1;
	int totalC_num = 2 * m_pah->m_rings + (CarbonListSize() + 3 * m_pah->m_rings5_Lone + 3 * m_pah->m_rings5_Embedded + 5 * m_pah->m_rings7_Lone + 5 * m_pah->m_rings7_Embedded) / 2 + numberOfBridges() + 1 + numberOfMethyl();
	m_pah->setnumofC(totalC_num);
    m_pah->setnumofH((int)vec.size() + 2 * numberOfMethyl());
    updateCombinedSites();
	
	// check if PAH closes correctly
	// This was moved to the end because the PAH was not giving enough information to debug.
	if (m_pah->m_clast == NULLC || !checkHindrance_twoC(newC->C1, m_pah->m_cfirst)) {
        // PAH did not close properly. invalid structure
        cout << "Error: createPAH: PAH did not close properly. Could be problem " //SETBREAKPOINT
            <<"with site list input...\n";
        std::ostringstream msg;
        msg << "ERROR: PAH did not close properly.."
            << " (Sweep::KMC_ARS::PAHProcess::createPAH_fromfile)";
		printSites();
		std:string filename = "KMC_DEBUG/KMC_PAH_X_CLOSE_fromfile";
		cout << "Saving file " << filename << ".xyz\n";
        saveXYZ(filename);
		//saveDOT("KMC_DEBUG/KMC_PAH_X_CLOSE.dot");
        //throw std::runtime_error(msg.str());
        //assert(false);
        //return;
    }
}

//Save a PAH to file with all details from current typespace. Such file can be opened with initialise_fromfile.
void PAHProcess::savePAH_tofile(const std::string &filename){
	ofstream dst(filename);
	if (dst.is_open()){
		dst << std::to_string(getCHCount().first) << " " << std::to_string(getCHCount().second) << "\n";
		std::string site_list_line;
		for(std::list<Site>::iterator i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
			// convert site type into string
			site_list_line = kmcSiteName(i->type);
			dst << site_list_line << ",";
		}
		dst << "\n";
		//Save state
		dst << std::to_string(m_pah->m_rings) << " " << 
			   std::to_string(m_pah->m_rings5_Lone) << " " <<
			   std::to_string(m_pah->m_rings5_Embedded) << " " <<
			   std::to_string(m_pah->m_rings7_Lone) << " " <<
			   std::to_string(m_pah->m_rings7_Embedded) << "\n";
		//Save edge carbon coordinates
		dst << "Coordinates\n";
		Cpointer Cnow = m_pah->m_cfirst;
		Cpointer Cprev = m_pah->m_cfirst;
		do {
			std::string my_line;
			my_line.append(std::to_string(std::get<0>(Cnow->coords)));
			my_line.append(" ");
			my_line.append(std::to_string(std::get<1>(Cnow->coords)));
			my_line.append(" ");
			my_line.append(std::to_string(std::get<2>(Cnow->coords)));
			my_line.append(" ");
			if (Cnow-> bridge) my_line.append("1");
			else my_line.append("0 ");
			if (Cnow->A == 'C'){
				my_line.append("C 0.0 0.0 0.0 \n");
			}
			else{
				my_line.append("H ");
				my_line.append(std::to_string(std::get<0>(Cnow->coords) + 1.085 * std::get<0>(Cnow->growth_vector)));
				my_line.append(" ");
				my_line.append(std::to_string(std::get<1>(Cnow->coords) + 1.085 * std::get<1>(Cnow->growth_vector)));
				my_line.append(" ");
				my_line.append(std::to_string(std::get<2>(Cnow->coords) + 1.085 * std::get<2>(Cnow->growth_vector)));
				my_line.append("\n");
			}
			dst << my_line;
			
			if (Cnow-> bridge && Cprev != Cnow->C3){
				Cprev = Cnow;
				Cnow = Cnow->C3;
			}
			else{
				Cprev = Cnow;
				Cnow = Cnow->C2;
			}
		}while (Cnow != m_pah->m_cfirst);
		//Save internal coordinates
		dst << "Internal\n";
		for(std::list<cpair>::iterator it = m_pah->m_InternalCarbons.begin(); it != m_pah->m_InternalCarbons.end(); ++it){
			std::string my_line;
			my_line.append(std::to_string(std::get<0>(*it)));
			my_line.append(" ");
			my_line.append(std::to_string(std::get<1>(*it)));
			my_line.append(" ");
			my_line.append(std::to_string(std::get<2>(*it)));
			my_line.append("\n");
			dst << my_line;
		}
		//Save R5 locations
		dst << "R5_locs\n";
		for(std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it != m_pah->m_R5loc.end(); ++it){
			std::string my_line;
			my_line.append(std::to_string(std::get<0>(*it)));
			my_line.append(" ");
			my_line.append(std::to_string(std::get<1>(*it)));
			my_line.append(" ");
			my_line.append(std::to_string(std::get<2>(*it)));
			my_line.append("\n");
			dst << my_line;
		}
		//Save R7 locations
		dst << "R7_locs\n";
		for(std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it != m_pah->m_R7loc.end(); ++it){
			std::string my_line;
			my_line.append(std::to_string(std::get<0>(*it)));
			my_line.append(" ");
			my_line.append(std::to_string(std::get<1>(*it)));
			my_line.append(" ");
			my_line.append(std::to_string(std::get<2>(*it)));
			my_line.append("\n");
			dst << my_line;
		}
	}
	else{
		cout << "Could not open destination file PAHProcess::savePAH_tofile.\n";
	}
	dst.close();
}

/*!
 * Initialization of PAH structure from an existing PAH structure (cloning)
 *
 * @param[in]    siteList_str    Comma-separated string of site types
 * @param[in]    R6_num          Number of 6-member rings
 * @param[in]    R5_num          Number of 5-member rings
 *
 * @return       Initialized PAH structure
 */
PAHStructure& PAHProcess::initialise(std::string siteList_str, int R6_num, int R5_num_Lone, int R5_num_Embedded, int R7_num_Lone, int R7_num_Embedded, Ccontainer edgeCarbons, std::list<cpair> internalCarbons, std::list<cpair> R5_locs, std::list<cpair> R7_locs){
    if(m_pah == NULL) {
        PAHStructure* pah = new PAHStructure();
        m_pah = pah;
    }else if(m_pah->m_cfirst != NULL)
        m_pah->clear();
    // create a vector from the string
    std::vector<std::string> siteList_strvec;
    Strings::split(siteList_str, siteList_strvec, std::string(","));
    // convert into vector of siteTypes
    std::vector<kmcSiteType> siteList_vec;
	std::vector<int> carb_siteList_vec;
    for(size_t i=0; i<siteList_strvec.size(); i++) {
        kmcSiteType temp = kmcSiteType_str(siteList_strvec[i]);
        if(temp == Inv) {
            std::cout<<"ERROR: Starting Structure site List contains invalid site type"
                <<".. (PAHProcess::initialise)\n\n";
			std::cout<<"Site = " << siteList_strvec[i] << "\n";
            std::ostringstream msg;
            msg << "ERROR: Starting Structure site List contains invalid site type."
                << " (Sweep::KMC_ARS::PAHProcess::initialise)";
                throw std::runtime_error(msg.str());
                assert(false);
        }
        siteList_vec.push_back(temp);
		int number_carbs = (int)temp %10 + 1;
		carb_siteList_vec.push_back(number_carbs);
    }
	createPAH(siteList_vec, carb_siteList_vec, R6_num, R5_num_Lone, R5_num_Embedded, R7_num_Lone, R7_num_Embedded, edgeCarbons, internalCarbons, R5_locs, R7_locs, m_pah->m_cfirst->coords);
    return *m_pah;
}

int create_pah_counter = 0;
// Create Structure from an existing PAH structure (cloning)
void PAHProcess::createPAH(std::vector<kmcSiteType>& vec, std::vector<int>& carb_vec, int R6, int R5_Lone, int R5_Embedded, int R7_Lone, int R7_Embedded, Ccontainer edCarbons, std::list<cpair> inCarbs, std::list<cpair> R5loc, std::list<cpair> R7loc, cpair first_carbon_coords) {
    //Add first carbon.
	Cpointer newC = addC();
	moveC(newC, first_carbon_coords);
    m_pah->m_cfirst = newC;
	m_pah->m_clast = NULLC;
	
	//Find first carbon in edge carbon container.	
	double mindist = 1e3;
	Ccontainer::iterator edCarbons_it;
    for(Ccontainer::iterator carbiter=edCarbons.begin(); carbiter != edCarbons.end(); ++carbiter) {
		double distance = getDistance_twoC(newC->coords, (*carbiter)->coords);
		if (distance < mindist){
			edCarbons_it = carbiter;
			mindist = distance;
		}	
    }
		
	// current C and coordinates
	moveC(newC, (*edCarbons_it)->coords);
	updateA(newC, (*edCarbons_it)->A, (*edCarbons_it)->growth_vector);
	Cpointer nextC = (*edCarbons_it)->C2;
    
    // start drawing..
    for(size_t i=0; i<vec.size(); i++) {
		int vec_carb_number = carb_vec[i];
		//int vec_carb_number = (int)vec[i] % 10 + 1;
		
		Cpointer site_carb_1 = newC;
		for (int j=0; j<vec_carb_number; j++) {
			if( !checkHindrance_C_PAH(nextC->coords)){
				newC = addC(newC, std::make_tuple(1.0,0.0,0.0), 7.17);
				moveC(newC, nextC->coords);
				updateA(newC, nextC->A, nextC->growth_vector);
				if (!nextC->bridge) nextC = nextC->C2;
				else {
					if (checkHindrance_C_PAH(nextC->C3->coords)) nextC = nextC->C2;
					else nextC = nextC->C3;
				}
			}
			else{ //position already occupied
				Cpointer Cpos = findC(nextC->coords);
				if(Cpos != m_pah->m_cfirst) { // it is a bridged C atom
					Cpointer Cbridge = Cpos->C1;
					Cpos->bridge = true; Cbridge->bridge = true;
					Cpos->C3 = Cbridge; Cbridge->C3 = Cpos;
					connectToC(newC, Cpos);
					Cbridge->C2 = NULLC;
					newC = Cbridge;
					nextC = nextC->C3->C2;
					j++;
				}
				else { // reached end of PAH
					connectToC(newC, m_pah->m_cfirst);
					m_pah->m_clast = newC;
					newC = newC->C2;
				}
			}
		}
		Cpointer site_carb_2 = newC;
        addSite(vec[i], site_carb_1, site_carb_2);
    }
	
    m_pah->m_rings = R6;
	m_pah->m_rings5_Lone = R5_Lone;
	m_pah->m_rings5_Embedded = R5_Embedded;
	m_pah->m_rings7_Lone = R7_Lone;
	m_pah->m_rings7_Embedded = R7_Embedded;
	m_pah->m_InternalCarbons = inCarbs;
	m_pah->m_R5loc = R5loc;
	m_pah->m_R7loc = R7loc;
	m_pah->m_methyl_counts = numberOfMethyl();
	m_pah->m_optimised = false;
	//int totalC_num = 2 * m_pah->m_rings + (CarbonListSize() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded) / 2 + numberOfBridges() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded + 1;
	int totalC_num = 2 * m_pah->m_rings + (CarbonListSize() + 3 * m_pah->m_rings5_Lone + 3 * m_pah->m_rings5_Embedded + 5 * m_pah->m_rings7_Lone + 5 * m_pah->m_rings7_Embedded) / 2 + numberOfBridges() + 1 + numberOfMethyl();
	m_pah->setnumofC(totalC_num);
    m_pah->setnumofH((int)vec.size() + 2 * numberOfMethyl());
    updateCombinedSites();
	
	// check if PAH closes correctly
	// This was moved to the end because the PAH was not giving enough information to debug.
	if (m_pah->m_clast == NULLC || !checkHindrance_twoC(newC->C1, m_pah->m_cfirst)) {
        // PAH did not close properly. invalid structure
        cout << "Error: createPAH: PAH did not close properly. Could be problem " //SETBREAKPOINT
            <<"with site list input...\n";
        std::ostringstream msg;
        msg << "ERROR: PAH did not close properly.."
            << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
		printSites();
		std::string filename = "KMC_DEBUG/KMC_PAH_X_CLOSE_";
		filename.append(std::to_string(create_pah_counter));
		cout << "Saving file " << filename << ".xyz\n";
        saveXYZ(filename);
		//saveDOT("KMC_DEBUG/KMC_PAH_X_CLOSE.dot");
        //throw std::runtime_error(msg.str());
        //assert(false);
        //return;
		create_pah_counter++;
    }
}

/*!
 * Initialization of PAH structure from binary file (deserialize)
 *
 * @return       Initialized PAH structure
 */
PAHStructure& PAHProcess::initialise(std::vector<std::tuple<int, int, cpair, cpair>> siteList_vector,
									 std::map<int, std::vector<std::tuple<int, int, cpair, cpair>>> siteList_map,
									 int R6_num, int R5_num_Lone, int R5_num_Embedded, 
									 int R7_num_Lone, int R7_num_Embedded, 
									 std::vector<std::tuple<cpair, int, int, cpair>> edgeCarbons, 
									 std::list<cpair> internalCarbons, 
									 std::list<cpair> R5_locs, 
									 std::list<cpair> R7_locs){
    if(m_pah == NULL) {
        PAHStructure* pah = new PAHStructure();
        m_pah = pah;
    }else if(m_pah->m_cfirst != NULL)
        m_pah->clear();
	
	//First fill Ccontainer
	if (edgeCarbons.size() <= 0) {
		std::cout<<"ERROR: Structure read from Binary file has no edge carbons."
			<<".. (PAHProcess::initialise)\n\n";
		std::ostringstream msg;
		msg << "ERROR: Structure read from Binary file has no edge carbons."
			<< " (Sweep::KMC_ARS::PAHProcess::initialise)";
			throw std::runtime_error(msg.str());
			assert(false);
	}
	
	Cpointer newC;

	for (std::vector<std::tuple<cpair, int, int, cpair>>::iterator it = edgeCarbons.begin(); it != edgeCarbons.end(); it++){
		if (it == edgeCarbons.begin()){
			//Add first carbon.
			newC = addC();
			m_pah->m_cfirst = newC;
			m_pah->m_clast = NULLC;
			moveC(newC, std::get<0>(*it));
			if (std::get<2>(*it) == 0) updateA(newC, 'C', std::get<3>(*it));
			if (std::get<2>(*it) == 1) updateA(newC, 'H', std::get<3>(*it));
			if (std::get<2>(*it) == 2) updateA(newC, 'M', std::get<3>(*it));
			if (std::get<1>(*it) == 1){
				it++;
				Cpointer Cbridge = addC(newC, std::make_tuple(1.0,0.0,0.0), 7.17);
				moveC(Cbridge, std::get<0>(*it));
				if (std::get<2>(*it) == 0) updateA(Cbridge, 'C', std::get<3>(*it));
				if (std::get<2>(*it) == 1) updateA(Cbridge, 'H', std::get<3>(*it));
				if (std::get<2>(*it) == 2) updateA(Cbridge, 'M', std::get<3>(*it));
				newC->bridge = true; Cbridge->bridge = true;
				newC->C3 = Cbridge; Cbridge->C3 = newC;
				newC->C2 = NULLC; Cbridge->C1 = NULLC;
				newC = Cbridge;
			}
		}
		else {
			if( !checkHindrance_C_PAH(std::get<0>(*it))){
				newC = addC(newC, std::make_tuple(1.0,0.0,0.0), 7.17);
				moveC(newC, std::get<0>(*it));
				if (std::get<2>(*it) == 0) updateA(newC, 'C', std::get<3>(*it));
				if (std::get<2>(*it) == 1) updateA(newC, 'H', std::get<3>(*it));
				if (std::get<2>(*it) == 2) updateA(newC, 'M', std::get<3>(*it));
				//newC->bridge = std::get<1>(*it);

				//if (std::get<1>(*it) == 0) nextC = nextC->C2;
				//else {
				//	if (checkHindrance_C_PAH(nextC->C3->coords)) nextC = nextC->C2;
				//	else nextC = nextC->C3;
				//}
			}
			else{ //position already occupied
				Cpointer Cpos = findC(std::get<0>(*it));
				Cpointer Cbridge;
				if (Cpos == m_pah->m_cfirst->C3){
					connectToC(newC, Cpos);
					newC = m_pah->m_cfirst;
				} else {
					Cbridge = Cpos->C1;
					Cpos->bridge = true; Cbridge->bridge = true;
					Cpos->C3 = Cbridge; Cbridge->C3 = Cpos;
					connectToC(newC, Cpos);
					Cbridge->C2 = NULLC;
					newC = Cbridge;
				}
				it++;
			}
		}
	}
	
	m_pah->m_clast = newC;
	connectToC(m_pah->m_clast, m_pah->m_cfirst);
	//printStruct();

	//Add sites	
	m_pah->m_siteList.clear(); //cout << "m_pah->m_siteList cleared!\n";
    for(std::vector<std::tuple<int, int, cpair, cpair>>::iterator it=siteList_vector.begin(); it!=siteList_vector.end(); it++) {
        kmcSiteType temp = kmcSiteType(std::get<0>(*it));
		kmcSiteType temp_comb = kmcSiteType(std::get<1>(*it));
        if(temp == Inv) {
            std::cout<<"ERROR: Structure read from Binary site list contains invalid site type."
                <<".. (PAHProcess::initialise)\n\n";
			std::cout<<"Site = " << temp << "\n";
            std::ostringstream msg;
            msg << "ERROR: Structure read from Binary site List contains invalid site type."
                << " (Sweep::KMC_ARS::PAHProcess::initialise)";
                throw std::runtime_error(msg.str());
                assert(false);
        }
		Cpointer siteC1 = findC(std::get<2>(*it));
		Cpointer siteC2 = findC(std::get<3>(*it));
		Spointer newSite = addSite(temp, siteC1, siteC2);
		newSite->comb = temp_comb;
    }
	//updateCombinedSites();

	//Create site map
	m_pah->m_siteMap.clear();//cout << "m_pah->m_siteMap cleared!\n";
	std::map<int, std::vector<std::tuple<int, int, cpair, cpair>>>::iterator map_it;
	for(map_it=siteList_map.begin();map_it!=siteList_map.end();map_it++){
		kmcSiteType site_type = (kmcSiteType)map_it->first;
		svector Spointer_vec;
		std::vector<std::tuple<int, int, cpair, cpair>> site_vec_coords = map_it->second;

		for(unsigned int vec_it = 0; vec_it!=site_vec_coords.size();vec_it++){
			Cpointer siteC1 = findC(std::get<2>(site_vec_coords[vec_it]));
			Cpointer siteC2 = findC(std::get<3>(site_vec_coords[vec_it]));
			Spointer Starget;
			bool site_found = false;

			for(Spointer S1 = m_pah->m_siteList.begin();S1!=m_pah->m_siteList.end();S1++){
				//Check if this is the site to be written
				if (S1->type == site_type || S1->comb == site_type){
					if(S1->C1==siteC1 && S1->C2==siteC2){
						//This is the correct site
						Starget = S1;
						Spointer_vec.push_back(Starget);
						site_found = true;
						break;
					}
				}
			}
			if (!site_found){
				std::cout << "Error. Site not found while reading Site Map from binary file.";
				std::ostringstream msg;
            	msg << "ERROR: Site not found while reading Site Map from binary file."
                	<< " (Sweep::KMC_ARS::PAHProcess::initialise)";
                	throw std::runtime_error(msg.str());
                	assert(false);
			}
		}
		m_pah->m_siteMap[site_type] = Spointer_vec;
	}

    m_pah->m_rings = R6_num;
	m_pah->m_rings5_Lone = R5_num_Lone;
	m_pah->m_rings5_Embedded = R5_num_Embedded;
	m_pah->m_rings7_Lone = R7_num_Lone;
	m_pah->m_rings7_Embedded = R7_num_Embedded;
	m_pah->m_InternalCarbons = internalCarbons;
	m_pah->m_R5loc = R5_locs;
	m_pah->m_R7loc = R7_locs;
	m_pah->m_methyl_counts = numberOfMethyl();
	m_pah->m_optimised = false;
	//int totalC_num = 2 * m_pah->m_rings + (CarbonListSize() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded) / 2 + numberOfBridges() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded + 1;
	int totalC_num = 2 * m_pah->m_rings + (CarbonListSize() + 3 * m_pah->m_rings5_Lone + 3 * m_pah->m_rings5_Embedded + 5 * m_pah->m_rings7_Lone + 5 * m_pah->m_rings7_Embedded) / 2 + numberOfBridges() + 1 + numberOfMethyl();
	m_pah->setnumofC(totalC_num);
    m_pah->setnumofH(SiteListSize() + 2 * numberOfMethyl());
    
	return *m_pah;
}



/*!
 * Initialization of PAH structure given a list of site types and number of 6 and 5-member rings
 *
 * @param[in]    siteList_str    Comma-separated string of site types
 * @param[in]    R6_num          Number of 6-member rings
 * @param[in]    R5_num          Number of 5-member rings
 *
 * @return       Initialized PAH structure
 */
PAHStructure& PAHProcess::initialise_sitelist_string(std::string siteList_str, int R6_num, int R5_num_Lone, int R5_num_Embedded, int R7_num_Lone, int R7_num_Embedded, int num_C, int num_H, int num_CH3, std::list<cpair> internalCarbons){
    if(m_pah == NULL) {
        PAHStructure* pah = new PAHStructure();
        m_pah = pah;
    }else if(m_pah->m_cfirst != NULL)
        m_pah->clear();
    // create a vector from the string
    std::vector<std::string> siteList_strvec;
    Strings::split(siteList_str, siteList_strvec, std::string(","));
    // convert into vector of siteTypes
    std::vector<kmcSiteType> siteList_vec;
    for(size_t i=0; i<siteList_strvec.size(); i++) {
        kmcSiteType temp = kmcSiteType_str(siteList_strvec[i]);
        if(temp == Inv) {
            std::cout<<"ERROR: Starting Structure site List contains invalid site type"
                <<".. (PAHProcess::initialise)\n\n";
            std::ostringstream msg;
            msg << "ERROR: Starting Structure site List contains invalid site type."
                << " (Sweep::KMC_ARS::PAHProcess::initialise)";
                //throw std::runtime_error(msg.str());
                //assert(false);
        }
        siteList_vec.push_back(temp);
    }
	createPAH_sitelist_string(siteList_vec, R6_num, R5_num_Lone, R5_num_Embedded, R7_num_Lone, R7_num_Embedded, num_C, num_H, num_CH3);
    return *m_pah;
}

// Create Structure from vector of site types. This function is used to assign PAHs from output files.
void PAHProcess::createPAH_sitelist_string(std::vector<kmcSiteType>& vec, int R6, int R5_Lone, int R5_Embedded, int R7_Lone, int R7_Embedded, int num_C, int num_H, int num_CH3) {
    //This function worked with Serialise and Deserialise to recognise identical PAHs through the simulation and avoid memory requirements.
	//However, the existence of 3D structures makes it extremely hard for the program to recognise repeated structures. Because of this,
	//it was decided that this function would just reproduce the statistics saved from the binary files to be able to compute mass spectra. 
	//If a PAH is desired for postprocessing, defining it in an additional list "tracked_pahs.txt" will produce its output through the OB routines.
	//gl413
	// current C, bondangle and coordinates
    Cpointer newC=addC(); 
    m_pah->m_cfirst = newC;
    m_pah->m_clast = NULLC;
    // number of bulk C to be added
    int bulkC;
	// angle to start drawing PAH.
	double angle = 0.0;
	int counts = 0; int countsR5 = 0;
    // type of site; if type 0, basic site types (FE - BY6); if type 1, R5 and basic sites with
    // a R5 at one side (RFE - RBY5); if type 2, basic sites wit R5 at each side (RFER - RACR)
    unsigned short int site_t;
    // start drawing..
    for(size_t i=0; i<vec.size(); i++) {
        Cpointer S_C1 = newC;
        // get number of bulk C to be added and site type
        if((int)vec[i] <= 4) {
            bulkC = (int) vec[i]; site_t = 0;
        }else if((int)vec[i] >= 100 && (int)vec[i] <= 104) {
            bulkC = (int) vec[i] - 100; site_t = 1;
        }else if((int)vec[i] >= 202 && (int)vec[i] <= 204) {
            bulkC = (int) vec[i] - 200; site_t = 2;
		}else if((int)vec[i] >= 501 && (int)vec[i] <= 504) {
            bulkC = (int) vec[i] - 500; site_t = 1;
        }else if((int)vec[i] >= 602 && (int)vec[i] <= 604) {
            bulkC = (int) vec[i] - 600; site_t = 2;
		}else if((int)vec[i] >= 1002 && (int)vec[i] <= 1004) {
            bulkC = (int) vec[i] - 1000; site_t = 2;
		}else if((int)vec[i] >= 2002 && (int)vec[i] <= 2115) {
            bulkC = (int) vec[i] %10; site_t = 3;
		}else if((int)vec[i] >= 2204 && (int)vec[i] <= 2205) {
            bulkC = (int) vec[i] %10; site_t = 3;
		}
		else if((int)vec[i] == 9999) {
            bulkC = (int) vec[i] %10; site_t = 3;
		}
		else {
            cout << "createPAH: Unrecognised site type. It becomes impossible to track the number of PAHs without defined site types.\n"
                << "The statistics for this PAH will be excluded.\n"; //SETBREAKPOINT
			cout << "Site = " << vec[i] << "\n";
            std::ostringstream msg;
            msg << "ERROR: Unrecognised site types found in list of sites."
                << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
                //throw std::runtime_error(msg.str());
                //assert(false);
            //return;
        }
		//cout << "Angle = " << angle << ". Next site to draw = " <<  (kmcSiteType)vec[i] << "\n";
        kmcSiteType prevType; 
        if(i==0) prevType = vec.back();
        else prevType = vec[i-1];
        switch(site_t) {
        case 0:
            newC = drawType0Site(newC, bulkC, angle); break;
        case 1:
            newC = drawType1Site(newC, bulkC, prevType, angle); break;
        case 2:
            newC = drawType2Site(newC, bulkC, angle); break;
		case 3:
            newC = drawType3Site(newC, bulkC, angle); break;
        default:
            cout << "createPAH: Invalid site_t number...\n";
            std::ostringstream msg;
            msg << "ERROR: invalid site classification."
                << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
            //    throw std::runtime_error(msg.str());
            //    assert(false);
            //return;
        }
        addSite(vec[i], S_C1, newC);
		counts = bulkC;
		if((int)vec[i] == 100 || (int)prevType == 100) countsR5 = 1;
		else countsR5 = 0;
		angle = normAngle( angle -60.0 + counts*60.0 - 1.0*countsR5*30.0) ;
    }
    // check if PAH closes correctly
	if (m_pah->m_clast == NULLC || checkHindrance_twoC(newC, m_pah->m_cfirst)) {
		//In the previous typespace this check was used to see if the PAH made sense. In 2D this is easily achieved. However, in 3D PAHs this is far from trivial.
        //For this reason the PAH is forced to be "closed" so the carbons are not double counted.
		// PAH did not close properly. invalid structure
        //cout << "createPAH: PAH did not close properly. Could be problem "
        //    <<"with site list input...\n";
        //std::ostringstream msg;
        //msg << "ERROR: PAH did not close properly.."
        //   << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
		//cout << msg;
        //saveDOT("KMC_DEBUG/KMC_PAH_X_CLOSE.dot");
        //throw std::runtime_error(msg.str());
        //assert(false);
        //return;
		//This connects the first and last carbons in the sitelist even if they are in unreal positions.
		connectToC(newC, m_pah->m_cfirst);
        m_pah->m_clast = newC;
    }
    m_pah->m_rings = R6; //Read from binary file.
	m_pah->m_rings5_Lone = R5_Lone; //Read from binary file.
	m_pah->m_rings5_Embedded = R5_Embedded; //Read from binary file.
	m_pah->m_rings7_Lone = R7_Lone; //Read from binary file. However, nbot tracked right now.
	m_pah->m_rings7_Embedded = R7_Embedded; //Read from binary file.
	m_pah->m_methyl_counts = num_CH3; //Read from binary file.
	//m_pah->m_InternalCarbons = inCarbs; // Not used.
    /*for(Ccontainer::iterator i = m_pah->m_carbonList.begin(); //Not used.
        i != m_pah->m_carbonList.end(); i++)
        updateA(*i, 'H');*/
	//Euler's characteristic equation for a topology containing 5, 6 and 7 membered rings with an open contiguous face (outside).
	//This was used to recompute the number of carbons, but it is less expensive to just write it in the binary and read it.
	//int totalC_num = 2 * m_pah->m_rings + (CarbonListSize() + 3 * m_pah->m_rings5_Lone + 3 * m_pah->m_rings5_Embedded + 5 * m_pah->m_rings7_Lone + 5 * m_pah->m_rings7_Embedded) / 2 + numberOfBridges() + 1;
	//m_pah->setnumofC(totalC_num);
    //m_pah->setnumofH((int)vec.size());
	m_pah->setnumofC(num_C);
	m_pah->setnumofH(num_H);
    updateCombinedSites();
}

//! For createPAH function: drawing type 0 sites
Cpointer PAHProcess::drawType0Site(Cpointer Cnow, int bulkC, double angle) {
    // draw site
	//cpair prev_direction, starting_direction, Hdirprev, Hdir;
	angle = normAngle(angle - 60);
    //angletype angle = normAngle(Cnow->bondAngle1-60);
    for(int c=0; c<=bulkC; ++c) {
		//Determine directions
		/*if (Cnow->C1 == NULLC) {
			starting_direction = make_tuple( 1.0, 0.0, 0.0);
			Hdirprev = std::make_tuple(cos(120.0 *M_PI / 180.0),sin(120.0 *M_PI / 180.0),0.0);
			Hdir = std::make_tuple(cos(60.0 *M_PI / 180.0),sin(60.0 *M_PI / 180.0),0.0);
		}
		else {
			prev_direction = get_vector(Cnow->C1->coords,Cnow->coords);
			Hdirprev = Cnow->C1->growth_vector;
			starting_direction = add_vector(prev_direction, invert_vector(Hdirprev));
			Hdir = prev_direction;
		}
		*/
        // check if adding on existing C atom (bridge)
        cpair pos = jumpToPos(Cnow->coords, angle, 0, 1.4);
		//cpair pos = jumpToPos(Cnow->coords, starting_direction, 1.4);
        if(checkHindrance_C_PAH(pos)) {
            // this coordinate is filled
            Cpointer Cpos = findC(pos);
            if(Cpos != m_pah->m_cfirst) { // it is a bridged C atom
                Cpointer Cbridge = Cpos->C1;
                //Cnow->bondAngle1 = angle;
                Cpos->bridge = true; Cbridge->bridge = true;
                Cpos->C3 = Cbridge; Cbridge->C3 = Cpos;
                connectToC(Cnow, Cpos);
                Cbridge->C2 = NULLC;
                angle = normAngle(angle + 120);
                //Cbridge->bondAngle1 = angle;
                Cnow = Cbridge;
                --bulkC;
            }else { // reached end of PAH
                //Cnow->bondAngle1 = angle;
                connectToC(Cnow, Cpos);
                m_pah->m_clast = Cnow;
                return Cpos;
            }
        }else {
			//Cnow = addC(Cnow, starting_direction, 1.4, false);
			//if (c==bulkC) updateA(Cnow, 'H', Hdir);
            Cnow = addC(Cnow, angle, angle, 1.4, false);
            angle = normAngle(angle+60);
        }
    }
    return Cnow;
}

//! For createPAH function: drawing type 1 sites
Cpointer PAHProcess::drawType1Site(Cpointer Cnow, int bulkC, kmcSiteType prevType, double angle) {
    //draw R5 site
	if ((int)prevType == 100) angle = normAngle(angle-30);
	//angletype angle = Cnow->bondAngle1-60;
    if(bulkC == 0) {
		angle = normAngle(angle-60);
        angle = normAngle(angle-30);
        cpair pos = jumpToPos(Cnow->coords, angle, 0, 1.4*sqrt(3.0));
        if(checkHindrance_C_PAH(pos)) { // reached end of PAH
            Cpointer Cpos = findC(pos);
            //Cnow->bondAngle1 = angle;
            connectToC(Cnow, Cpos);
            m_pah->m_clast = Cnow;
            return Cpos; // m_cfirst
        }
        else Cnow = addC(Cnow, angle, normAngle(angle-30), 1.4*sqrt(3.0), false);
        return Cnow;
    }else { //draw RXX site
            return drawType0Site(Cnow, bulkC, angle);
    }
}

//! For createPAH function: drawing type 2 sites
Cpointer PAHProcess::drawType2Site(Cpointer Cnow, int bulkC, double angle) {
    //angletype angle = Cnow->bondAngle1;
    //Cnow->bondAngle1 = normAngle(angle-60);
	angle = normAngle(angle-30);
    return drawType0Site(Cnow, bulkC, angle);
}

//! For createPAH function: drawing type 2 sites
Cpointer PAHProcess::drawType3Site(Cpointer Cnow, int bulkC, double angle) {
    //angletype angle = Cnow->bondAngle1;
    //Cnow->bondAngle1 = normAngle(angle-60);
	angle = normAngle(angle-30);
    return drawType0Site(Cnow, bulkC, angle);
}

//! Finds C atom with specific coordinates
Cpointer PAHProcess::findC(cpair coordinates) {
	double tol = 3e-1;
    for(Ccontainer::iterator i=m_pah->m_carbonList.begin(); i != m_pah->m_carbonList.end(); ++i) {
		double distance = getDistance_twoC(coordinates, (*i)->coords);
		if (distance < tol){
			return (*i);
		}
    }
    return NULLC;
}

//! Find Site to the other side of a bridge
Spointer PAHProcess::findSite(Cpointer C_1) {
	Spointer temp;
	Cpointer C_check = C_1;
	for (int k = 0; k != 6; k++){
		for (std::list<Site>::iterator it = m_pah->m_siteList.begin(); it != m_pah->m_siteList.end(); it++) {
			if (C_check == it->C1){
				// C_1 and C_2 are the first and last carbon of site
				temp = it;
				return temp;
			}
		}
		if (C_check->bridge) C_check = C_check->C3->C1;
		else C_check = C_check->C1;
	}
	cout << "Did not find site PAHProcess::findSite. \n";
	saveXYZ("KMC_DEBUG/findSite");
	return m_pah->m_siteList.end();
}

//! Find Site to the other side of a bridge
Cpointer PAHProcess::findThirdC(Cpointer C_1) {
	for(Ccontainer::iterator i=m_pah->m_carbonList.begin(); i != m_pah->m_carbonList.end(); ++i) {
		Cpointer C_check = *i;
		if (C_check != C_1 && C_check != C_1->C1 && C_check != C_1->C1->C1 && C_check != C_1->C2 && C_check != C_1->C2->C2 && C_check != C_1->C3){
			double dist = getDistance_twoC(C_check, C_1);
			if (dist < 1.8) {
				cpair vec1 = get_vector(C_1->C1->coords, C_1->coords);
				cpair vec2 = get_vector(C_1->C2->coords, C_1->coords);
				cpair addvec = add_vector(vec1, vec2);
				addvec = jumpToPos(C_1->coords,addvec,1.4);
				if (getDistance_twoC(C_check->coords, addvec)<1.8) return C_check; 
			}
		}
	}
	return NULLC;
}

// Check to validate if coordinates of C matches bond angles
bool PAHProcess::checkCoordinates() const{
    // start at first C in site list
    Cpointer start = m_pah->m_siteList.begin()->C1;
    Cpointer now = start;
    Cpointer next = now->C2;
    unsigned int count=0;
    do{
        count++;
        Cpointer oldnext = next;
        cpair corr_coords;
        // first check if next is a valid Carbon pointer
        if(next == NULL) {
            cout<<"checkCoordinates() failed at "<<count<<"th C atom..\n"
                <<"Message: This atom is not a valid pointer to Carbon -- "
                <<next<<"\n";
            //printStruct(next);
            return false;
        }
        if(next->bridge) {// if next C is a bridge
            if(now == next->C1) { // check if current C is next->C1
                corr_coords = jumpToPos(now->coords, now->bondAngle1, 0, 1.4);
                next = next->C3; // cross bridge
            }else if(now == next->C3) { // check if current C is from bridge
                corr_coords = jumpToPos(now->coords, now->bondAngle2, 0, 1.4);
                next = next->C2;
            }
        }else { // if next C is not a bridge
            corr_coords = jumpToPos(now->coords, now->bondAngle1, 0, 1.4);
            next = next->C2;
        }
        now = oldnext;
        // check coordinates
        if(now->coords != corr_coords) {
            cout<<"checkCoordinates() failed at "<<count<<"th C atom..\n"
				<< "Coordinates: (" << std::get<0>(now->coords) << ',' << std::get<1>(now->coords) << ',' << std::get<2>(now->coords) <<")\n"
				<< "Correct      : (" << std::get<0>(corr_coords) << ',' << std::get<1>(corr_coords) << ',' << std::get<2>(corr_coords)<<")\n";
            //saveDOT(std::string("KMC_DEBUG/error_struct.dot"));
            return false;
        }
    }while(now != start);
    return true;
}

// Check if site stt is a valid site type.
bool PAHProcess::checkSiteValid(int type) {
	switch (kmcSiteType(type)){
		case FE: return true;
		case ZZ: return true;
		case AC: return true;
		case BY5: return true;
		case BY6: return true;
		case R5: return true;
		case RFE: return true;
		case RZZ: return true;
		case RAC: return true;
		case RBY5: return true;
		case RFER: return true;
		case RZZR: return true;
		case RACR: return true;
		case FE3: return true;
		case AC_FE3: return true;
		case FE_HACA: return true;
		case BY5_FE3: return true;
		case FE2: return true;
		case R5R7: return true;
		case R5R6: return true;
		case R5R6ZZ: return true;
		case R5R6AC: return true;
		case R5R6BY5: return true;
		case R5R6FER: return true;
		case R5R6ZZR: return true;
		case R5R6ACR: return true;
		case R5R6FER5R6: return true;
		case R5R6ZZR5R6: return true;
		case R5R6ACR5R6: return true;
		case ACR5: return true;
		case FEACR5: return true;
		case ZZACR5: return true;
		case FEACR5FE: return true;
		case ACACR5: return true;
		case FEZZACR5: return true;
		case R5ACR5: return true;
		case R5FEACR5: return true;
		case ACR5R5R6: return true;
		case R5ZZACR5: return true;
		case ACR5R5R6ZZ: return true;
		case R5ACR5R5: return true;
		case ACR5RFER: return true;
		case RAC_FE3: return true;
		case MIGR: return true;
		case MIGR2: return true;
		case R5R6_MIGR: return true;
		case SPIRAL: return true;
		case None: return false; //None should return an error
		case Inv: return true;
		case any: return true;
		case benz: return true;
		case Methyl: return true;
	}
	return false;
}

// Check if site stt is a valid site type.
bool PAHProcess::checkSiteValid(Spointer& stt) {
	return checkSiteValid(stt->type);
}

// Check to see if all sites are connected to each other
bool PAHProcess::checkSiteContinuity() const {
    std::list<Site>::const_iterator i=m_pah->m_siteList.begin();
    for(unsigned int k = 0; k!=(unsigned int)m_pah->m_siteList.size();k++) {
        Cpointer lhs = i->C2; i++;
        if(i == m_pah->m_siteList.end()) i = m_pah->m_siteList.begin();
        Cpointer rhs = i->C1;
        if(lhs != rhs) 
            return false;
    }
    return true;
}

// Check to see if site neighbours has a valid combined site type
bool PAHProcess::checkCombinedSiteType(Spointer& stt) {
    int k = 3;
    if(m_pah->m_siteList.size() < 8) k = 1;
    Spointer startSite = moveIt(stt, -k);
    Spointer endSite = moveIt(stt, k+1);
    bool error = false;
    kmcSiteType error_stype_comb = Inv;
    kmcSiteType error_stype = Inv;
    for(Spointer i=startSite; i!=endSite; i = moveIt(i,1)) {
        switch(i->comb) {
            case FE3:
                if(i->type != FE) {
                    error = true;
                    error_stype_comb = FE3;
                    error_stype = i->type;
                }
                break;
            case FE2:
                if(i->type != FE) {
                    error = true;
                    error_stype_comb = FE2;
                    error_stype = i->type;
                }
                break;
            case AC_FE3:
                if(i->type != AC) {
                    error = true;
                    error_stype_comb = AC_FE3;
                    error_stype = i->type;
                }
                break;
            case FE_HACA:
                if(i->type != FE) {
                    error = true;
                    error_stype_comb = FE_HACA;
                    error_stype = i->type;
                }
                break;
            case BY5_FE3:
                if(i->type != BY5) {
                    error = true;
                    error_stype_comb = BY5_FE3;
                    error_stype = i->type;
                }
                break;
			case R5R7:
				if((int)i->type < 500) {
                    error = true;
                    error_stype_comb = R5R7;
                    error_stype = i->type;
                }
            default: break;
            }
        }
    if(error) {
        std::cout<<"ERROR: Invalid combined site type -- Combined site type "
            <<kmcSiteName(error_stype_comb)<<" on a principal site type "
            <<kmcSiteName(error_stype)<<" (PAHProcess::checkCombinedSiteType)\n";
        return false;
    }
    return true;
}
int perform_process_error_counter = 0;
int r5_error_counter = 0;
int r7_error_counter = 0;
//! Structure processes: returns success or failure
bool PAHProcess::performProcess(const JumpProcess& jp, rng_type &rng, int PAH_ID)
{
	total_error_counter = perform_process_error_counter + r5_error_counter + r7_error_counter + 
						forcefield_error_counter + updatesites_error_counter + convSiteType_error_counter;
	if (total_error_counter >=20){
		std::ostringstream msg;
		msg << "ERROR: Over 20 error files have been written. Terminating program."
			<< " (Sweep::KMC_ARS::PAHProcess::performProcess)";
		throw std::runtime_error(msg.str());
		assert(false);
		abort();
	}
    //printStruct();
    //cout << "Start Performing Process..\n";
    kmcSiteType stp = jp.getSiteType();
	Spointer site_perf;
    int id = jp.getID();
    
    // choose random site of type stp to perform process
	site_perf = chooseRandomSite(stp, rng); //cout<<"[random site chosen..]\n";
	if ( site_perf->type != stp && site_perf->comb != stp && stp != benz && stp != Methyl){
		
		cout << "chooseRandomSite ERROR. Choose random site selected incorrect site. Error. Not performing jump process.\n";
		cout << "Jump process performed: " << id << ". " << jp.getName() << "\n";
		printSites(site_perf);
		return false;
	}
    // stores pointers to site Carbon members
    Cpointer site_C1 = site_perf->C1;
    Cpointer site_C2 = site_perf->C2;
    //cout << jp.getName() << '\n';
    //cout<<'\t'<<kmcSiteName(site_perf->type)<<' '<<site_C1<<' '<<site_C2<<'\n';
    // find structure change function
    //std::ostringstream dotname, dotname2;
	//Comment if not debugging
	//////
	//std::string dotname_pre_jp = "KMC_DEBUG/";
	//dotname_pre_jp.append(std::to_string(PAH_ID));
	//cout << PAH_ID << '\n';
	//dotname_pre_jp.append("pre_jp.dot");
	//saveDOT(dotname_pre_jp);
	/*if(PAH_ID==13){
		cout << "R5 coordinate list before JP:\n";
		for(std::list<cpair>::iterator it=m_pah->m_R5loc.begin(); it!=m_pah->m_R5loc.end();++it){
			cout << std::get<0>(*it) <<", "<< std::get<1>(*it) <<", "<< std::get<2>(*it) <<"\n";
		}
		cout << "R7 coordinate list before JP:\n";
		for(std::list<cpair>::iterator it=m_pah->m_R7loc.begin(); it!=m_pah->m_R7loc.end();++it){
			cout << std::get<0>(*it) <<", "<< std::get<1>(*it) <<", "<< std::get<2>(*it) <<"\n";
		}
	}*/
	//Save an XYZ
	if(m_debug_pah) saveXYZ("KMC_DEBUG/BEFORE");
	//Copy site list before performing process
	std::list<std::string> Sitelist_before = copySites(site_perf);
	//Check if the site has the right number of carbons
	if(!SiteRightSize(site_perf) && (id != 24 && id != 34 && id != 66)  ){
		std::cout << "ERROR: Site selected has incorrect number of carbons. Process not performed."
			<< "ID" << id << " Jump process: " << jp.getName() << " on PAH ID: " << PAH_ID << "..."
			<< " (Sweep::KMC_ARS::PAHProcess::performProcess)";
		std::string filename = "KMC_DEBUG/KMC_PAH_performProcess_SiteSize_";
		filename.append(std::to_string(perform_process_error_counter));
		saveXYZ(filename);
		std::ostringstream msg;
		msg << "ERROR: Site selected has incorrect number of carbons. Process not performed."
			<< "ID" << id << " Jump process: " << jp.getName() << " on PAH ID: " << PAH_ID << "..."
			<< " (Sweep::KMC_ARS::PAHProcess::performProcess)";
		//throw std::runtime_error(msg.str());
		//assert(false);
		//abort();
		printBeforeSites(Sitelist_before);
		cout<<"Saving file: "<< filename<<".xyz\n";
		++perform_process_error_counter;
		return false;
	}
	///////
    switch(id) {
        case 1:
            proc_G6R_AC(site_perf, site_C1, site_C2); break;
        case 2:
            proc_G6R_FE(site_perf, site_C1, site_C2); break;
        case 3:
            proc_L6_BY6(site_perf, site_C1, site_C2); break;
        case 4:
            proc_PH_benz(site_perf, site_C1, site_C2, rng); break;
        case 5:
            proc_D6R_FE3(site_perf, site_C1, site_C2); break;
        case 6:
            proc_O6R_FE3_O2(site_perf, site_C1, site_C2); break;
        case 7:
            proc_O6R_FE3_OH(site_perf, site_C1, site_C2); break;
        case 8:
            proc_O6R_FE_HACA(site_perf, site_C1, site_C2); break;
        case 9:
            proc_O6R_FE_HACA_OH(site_perf, site_C1, site_C2); break;
        case 10:
            proc_G5R_ZZ(site_perf, site_C1, site_C2); break;
        case 11:
            proc_D5R_R5(site_perf, site_C1, site_C2); break;
        case 12:
            proc_C6R_AC_FE3(site_perf, site_C1, site_C2, rng); break;
        case 13:
            proc_C5R_RFE(site_perf, site_C1, site_C2); break;
        case 14:
            proc_C5R_RAC(site_perf, site_C1, site_C2); break;
        case 15:
            proc_M5R_RZZ(site_perf, site_C1, site_C2); break;
        case 16:
            proc_C6R_BY5_FE3(site_perf, site_C1, site_C2, rng); break;
        case 17:
            proc_C6R_BY5_FE3violi(site_perf, site_C1, site_C2, rng); break;
        case 18:
            proc_L5R_BY5(site_perf, site_C1, site_C2); break;
        case 19:
            proc_M6R_BY5_FE3(site_perf, site_C1, site_C2, rng); break;
        case 20:
            proc_O6R_FE2(site_perf, site_C1, site_C2); break;
        case 21:
            proc_O6R_FE2(site_perf, site_C1, site_C2); break;
		case 22:
			proc_D6R_FE_AC(site_perf, site_C1, site_C2); break;
		case 23:
			proc_B6R_ACR5(site_perf, site_C1, site_C2); break;
		case 24:
			//proc_M5R_ACR5_multiple_sites(site_perf, site_C1, site_C2, rng); break;
			//proc_M5R_ACR5_ZZ(site_perf, site_C1, site_C2, rng); break;
			//proc_M5R_FEACR5(site_perf, site_C1, site_C2); break;
			proc_M5R_ACR5_ZZ_light(site_perf, site_C1, site_C2); break;
		case 25:
			proc_G6R_RZZ(site_perf, site_C1, site_C2); break;
		case 26:
			proc_G6R_RFER(site_perf, site_C1, site_C2); break;
		case 27:
			proc_G6R_R5(site_perf, site_C1, site_C2); break;
		case 28:
			proc_L6_RBY5(site_perf, site_C1, site_C2); break;
		case 29:
			proc_L6_RACR(site_perf, site_C1, site_C2); break;
		case 30:
			proc_G5R_RFE(site_perf, site_C1, site_C2); break;
		case 31:
			proc_C6R_BY5_FE3(site_perf, site_C1, site_C2, rng); break;
		case 32:
			proc_C6R_BY5_FE3(site_perf, site_C1, site_C2, rng); break;
		case 33:
			proc_C6R_RAC_FE3(site_perf, site_C1, site_C2, rng); break;
		case 34:
			//proc_M5R_R5R6_multiple_sites(site_perf, site_C1, site_C2, rng); break;
			//proc_MR5_R6(site_perf, site_C1, site_C2); break;
			proc_MR5_R6_light(site_perf, site_C1, site_C2); break;
		case 35:
			proc_GR7_R5R6AC(site_perf, site_C1, site_C2); break;
		case 36:
			proc_GR7_FEACR5(site_perf, site_C1, site_C2); break;
		case 37:
			proc_G6R_R5R6ZZ(site_perf, site_C1, site_C2); break;
		case 38:
			proc_L7_ACACR5(site_perf, site_C1, site_C2); break;
		case 39:
			proc_G6R_R5R6FER(site_perf, site_C1, site_C2); break;
		case 40:
			proc_G6R_R5R6FER5R6(site_perf, site_C1, site_C2); break;
		case 41:
			proc_L7_FEZZACR5(site_perf, site_C1, site_C2); break;
		case 42:
			proc_C5R_RZZR(site_perf, site_C1, site_C2, rng); break;
		case 43:
			proc_C5R_R5R6ZZR(site_perf, site_C1, site_C2, rng); break;
		case 44:
			proc_L6_R5R6BY5(site_perf, site_C1, site_C2); break;
		case 45:
			proc_L6_R5R6ACR(site_perf, site_C1, site_C2); break;
		case 46:
			proc_L6_R5R6ACR5R6(site_perf, site_C1, site_C2); break;
		case 47:
			proc_L6_ZZACR5(site_perf, site_C1, site_C2); break;
		case 48:
			proc_L6_R5FEACR5(site_perf, site_C1, site_C2); break;
		case 49:
			proc_L6_FEACR5FE(site_perf, site_C1, site_C2); break;
		case 50:
			proc_L6_R5ACR5R5(site_perf, site_C1, site_C2); break;
		case 51:
			proc_L7_R5ZZACR5(site_perf, site_C1, site_C2); break;
		case 52:
			proc_L6_ACR5R5R6(site_perf, site_C1, site_C2); break;
		case 53:
			proc_L7_ACR5R5R6ZZ(site_perf, site_C1, site_C2); break;
		case 54:
			proc_A_CH3(site_perf, site_C1, site_C2, rng); break;
		case 55:
			proc_D_CH3(site_perf, site_C1, site_C2); break;
		case 56:
			proc_O5R_R5R6(site_perf, site_C1, site_C2, rng); break;
		case 57:
			proc_O5R_R5R6ZZ(site_perf, site_C1, site_C2, rng); break;
		case 58:
			proc_O5R_R5R6AC(site_perf, site_C1, site_C2, rng); break;
		case 59:
			proc_O5R_R5R6BY5(site_perf, site_C1, site_C2, rng); break;
		case 60:
			proc_O5R_R5R6FER(site_perf, site_C1, site_C2, rng); break;
		case 61:
			proc_O5R_R5R6ZZR(site_perf, site_C1, site_C2, rng); break;
		case 62:
			proc_O5R_R5R6ACR(site_perf, site_C1, site_C2, rng); break;
		case 63:
			proc_O5R_R5R6FER5R6(site_perf, site_C1, site_C2, rng); break;
		case 64:
			proc_O5R_R5R6ZZR5R6(site_perf, site_C1, site_C2, rng); break;
		case 65:
			proc_O5R_R5R6ACR5R6(site_perf, site_C1, site_C2, rng); break;
		case 66:
			//proc_M5R_FEACR5_multiple_sites(site_perf, site_C1, site_C2, sFE2, b4, rng); break;
			//proc_M5R_ACR5_ZZ(site_perf, site_C1, site_C2, rng); break;
			proc_M5R_ACR5_ZZ_ZZ_light(site_perf, site_C1, site_C2, rng); break;
		case 67:
			proc_MR5R7_edge(site_perf, site_C1, site_C2, rng); break;
        default:
            cout<<"ERROR: PAHProcess::performProcess: Process not found\n";
            return false;
    }/*
    if(m_pah->m_siteList.begin()->C1 != m_pah->m_cfirst)
        cout<<"WARNING: C1 of Site 1 does not correspond to m_cfirst..\n";
    if(m_pah->m_cfirst->A != 'H')
        cout<<"WARNING: A of m_cfirst is not H..\n";
    //cout<<"----PROCESS PERFORMED!-----\n";*/
	//printSites();
    Spointer S1,S2,S3,S4;
    S1 = moveIt(site_perf, -1); S2 = moveIt(site_perf, 1);
    S3 = moveIt(site_perf, -2); S4 = moveIt(site_perf, 2);
	if (!checkSiteValid(S1) || !checkSiteValid(S2) || !checkSiteValid(S3) || !checkSiteValid(S4)){
		std::cout << "ERROR. Structure produced invalid site type after performing process " 
			<< "ID" << id << " Jump process: " << jp.getName() << " on PAH ID: " << PAH_ID << "...\n"
			<< "*************\nAfter performing process --\n";
		printBeforeSites(Sitelist_before); //SETBREAKPOINT
		printSites(site_perf);
		if (m_debug_pah){
			ifstream  src("KMC_DEBUG/BEFORE.xyz");
			std::string filename = "KMC_DEBUG/BEFORE_performProcess_";
			filename.append(std::to_string(perform_process_error_counter));
			filename.append(".xyz");
			ofstream dst(filename);
			dst << src.rdbuf();
			src.close();
			dst.close();
			cout<<"Saving file: "<< filename<<"\n";
		}
		std::string filename2 = "KMC_DEBUG/KMC_PAH_performProcess_checkSiteValid_";
		filename2.append(std::to_string(perform_process_error_counter));
		saveXYZ(filename2);
		std::ostringstream msg;
		msg << "ERROR: Structure produced invalid combined site type after performing process "
			<< "ID" << id << " Jump process: " << jp.getName() << " on PAH ID: " << PAH_ID << "..."
			<< " (Sweep::KMC_ARS::PAHProcess::performProcess)";
		//throw std::runtime_error(msg.str());
		//assert(false);
		//abort();
		cout<<"Saving file: "<< filename2<<".xyz\n";
		++perform_process_error_counter;
	}
    if(!checkCombinedSiteType(site_perf) || !checkCombinedSiteType(S1)
        || !checkCombinedSiteType(S2) || !checkCombinedSiteType(S3)
        || !checkCombinedSiteType(S4)) {
        std::cout<<"ERROR. Structure produced invalid combined site type after performing process "
            << "ID"<<id << "Jump process: " << jp.getName() << " on PAH ID: "<< PAH_ID <<"...\n"
            <<"*************\nAfter performing process --\n";
        printBeforeSites(Sitelist_before); //SETBREAKPOINT
		printSites(site_perf);
		if(m_debug_pah){
			ifstream  src("KMC_DEBUG/BEFORE.xyz");
			std::string filename = "KMC_DEBUG/BEFORE_performProcess_";
			filename.append(std::to_string(perform_process_error_counter));
			filename.append(".xyz");
			ofstream dst(filename);
			dst << src.rdbuf();
			src.close();
			dst.close();
			cout<<"Saving file: "<< filename<<"\n";
		}
		std::string filename2 = "KMC_DEBUG/KMC_PAH_performProcess_checkCombinedSite_";
		filename2.append(std::to_string(perform_process_error_counter));
		saveXYZ(filename2);
        std::ostringstream msg;
        msg << "ERROR: Structure produced invalid combined site type after performing process "
            << "ID"<<id << "Jump process: " << jp.getName() <<" on PAH ID: "<< PAH_ID <<"..."
            << " (Sweep::KMC_ARS::PAHProcess::performProcess)";
		cout<<"Saving file: "<< filename2<<".xyz\n";
		++perform_process_error_counter;	
        //throw std::runtime_error(msg.str());
        //assert(false);
        //abort();
    }
	if ( (!SiteRightSize(site_perf) || !SiteRightSize(S1)
        || !SiteRightSize(S2) || !SiteRightSize(S3)
        || !SiteRightSize(S4)) && (id != 24 && id != 34 && id != 66) ) {
		
		std::cout <<"ERROR: Site produced has incorrect number of carbons after process was performed."
			<< "ID" << id << " Jump process: " << jp.getName() << " on PAH ID: " << PAH_ID << "..."
			<< " (Sweep::KMC_ARS::PAHProcess::performProcess)";
		std::string filename = "KMC_DEBUG/KMC_PAH_performProcess_SiteSize_";
		filename.append(std::to_string(perform_process_error_counter));
		saveXYZ(filename);
		std::ostringstream msg;
		msg << "ERROR: Site produced has incorrect number of carbons after process was performed."
			<< "ID" << id << " Jump process: " << jp.getName() << " on PAH ID: " << PAH_ID << "..."
			<< " (Sweep::KMC_ARS::PAHProcess::performProcess)";
		//throw std::runtime_error(msg.str());
		//assert(false);
		//abort();
		printBeforeSites(Sitelist_before);
		cout<<"Saving file: "<< filename<<".xyz\n";
		++perform_process_error_counter;
	}
	//int calc_total = 2 * m_pah->m_rings + (CarbonListSize() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded) / 2 + numberOfBridges() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded + 1;
	int calc_total = 2 * m_pah->m_rings + (CarbonListSize() + 3 * m_pah->m_rings5_Lone + 3 * m_pah->m_rings5_Embedded + 5 * m_pah->m_rings7_Lone + 5 * m_pah->m_rings7_Embedded) / 2 + numberOfBridges() + 1 + numberOfMethyl();
	if(calc_total != getCHCount().first) {
        //saveDOT("KMC_DEBUG/KMC_C_Counts_ERROR.dot");
		if(m_debug_pah){
			ifstream  src("KMC_DEBUG/BEFORE.xyz");
			std::string filename = "KMC_DEBUG/BEFORE_performProcess_";
			filename.append(std::to_string(perform_process_error_counter)); //SETBREAKPOINT
			filename.append(".xyz");
			ofstream dst(filename);
			dst << src.rdbuf();
			src.close();
			dst.close();
			cout<<"Saving file: "<< filename<<"\n";
		}
		std::string filename2 = "KMC_DEBUG/KMC_PAH_performProcess_C_Count_";
		filename2.append(std::to_string(perform_process_error_counter));
		saveXYZ(filename2);
        cout<<"ERROR. Calculated total did not tally with double C counts!\n";
        cout<<"Last process performed: "<<jp.getName()<<", ID = "<<jp.getID()<<", PAH ID = "<<PAH_ID<<'\n';
		cout<<"R6s = "<<m_pah->m_rings<<", R5s = "<<m_pah->m_rings5_Lone<<" alone + "<<m_pah->m_rings5_Embedded<<" embedded, R7s = "<<m_pah->m_rings7_Embedded<<"\n";
        std::ostringstream msg;
        msg << "\nCalculated total: "<<calc_total<<'\n'
            << "Real total: " << getCHCount().first << '\n';
        printBeforeSites(Sitelist_before);
		printSites(site_perf);
		//throw std::runtime_error(msg.str());
		cout<<"Saving file: "<< filename2<<".xyz\n";
		++perform_process_error_counter;
    }

	if (id != 24 && id != 34 && id != 66){
		//Not very useful during migration processes
		if (m_pah->m_R5loc.size()+m_pah->m_R5walker_sites.size() - (m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded) != 0){
			//Try to fix this error by an optimisation.
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		if (m_pah->m_R5loc.size()+m_pah->m_R5walker_sites.size() - (m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded) != 0 && SiteListSize() >= 3){
			//Throw error.
			cout << "Error. Number of R5s in m_pah->m_R5loc does not match number of lone and embedded R5s.\n";
			cout << "\t Total R5 rings: " << m_pah->m_rings5_Lone << " lone + "<< m_pah->m_rings5_Embedded << " embedded.\n";
			cout << "R5 coordinates in list = "<< m_pah->m_R5loc.size() << ".\n";
			cout << " Jump process: " << jp.getName() <<" on PAH ID: "<< PAH_ID <<"\n";
			cout << "Printing internal R5 positions:\n"; //SETBREAKPOINT
			for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
				cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
			}
			if(m_debug_pah){
				ifstream  src("KMC_DEBUG/BEFORE.xyz");
				std::string filename = "KMC_DEBUG/BEFORE_R5count_";
				filename.append(std::to_string(r5_error_counter)); //SETBREAKPOINT
				filename.append(".xyz");
				ofstream dst(filename);
				dst << src.rdbuf();
				src.close();
				dst.close();
				cout<<"Saving file: "<< filename<<"\n";
			}
			std::string filename2 = "KMC_DEBUG/AFTER_R5count_";
			filename2.append(std::to_string(r5_error_counter));
			saveXYZ(filename2);
			printBeforeSites(Sitelist_before);
			printSites(site_perf);
			cout<<"Saving file: "<< filename2<<".xyz\n";
			//filename2.append(".inx");
			//savePAH_tofile(filename2);
			//cout<<"Saving file: "<< filename2<<"\n";
			++r5_error_counter;
			//throw std::runtime_error(msg.str());
		}
		
		if (m_pah->m_R7loc.size() - (m_pah->m_rings7_Embedded) != 0){
			//Try to fix this error by an optimisation.
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		if (m_pah->m_R7loc.size() - (m_pah->m_rings7_Embedded) != 0 && SiteListSize() >= 3){
			//Throw error.
			cout << "Error. Number of R7s in m_pah->m_R7loc does not match number of lone and embedded R5s.\n";
			cout << "\t Total R7 rings: " << m_pah->m_rings7_Embedded << " embedded.\n";
			cout << "R7 coordinates in list = "<< m_pah->m_R7loc.size() << ".\n";
			cout << " Jump process: " << jp.getName() <<" on PAH ID: "<< PAH_ID <<"\n";
			cout << "Printing internal R7 positions:\n"; //SETBREAKPOINT
			for (std::list<cpair>::iterator it1 = m_pah->m_R7loc.begin(); it1 != m_pah->m_R7loc.end(); ++it1) {
				cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
			}
			if(m_debug_pah){
				ifstream  src("KMC_DEBUG/BEFORE.xyz");
				std::string filename = "KMC_DEBUG/BEFORE_R7count_";
				filename.append(std::to_string(r7_error_counter)); //SETBREAKPOINT
				filename.append(".xyz");
				ofstream dst(filename);
				dst << src.rdbuf();
				src.close();
				dst.close();
				cout<<"Saving file: "<< filename<<"\n";
			}
			std::string filename2 = "KMC_DEBUG/AFTER_R7count_";
			filename2.append(std::to_string(r7_error_counter));
			saveXYZ(filename2);
			printBeforeSites(Sitelist_before);
			printSites(site_perf);
			//throw std::runtime_error(msg.str());
			cout<<"Saving file: "<< filename2<<".xyz\n";
			++r7_error_counter;
		}
	}
	
	//Save an XYZ
	//saveXYZ("KMC_DEBUG/AFTER");
	/*Cpointer Ccheck = site_perf->C1->C1;
	for (int ii=1; ii!=6; ++ii){
		double dist = getDistance_twoC(Ccheck,Ccheck->C2);*/
	/*if (PAH_ID==5){
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol, 3000);
		passbackPAH(mol);
	}*/
		/*Ccheck = Ccheck->C2;
	}*/
	//printSites(site_perf);
    return true;
}

//--------------------------------------------------------------------
//----------------- STRUCTURE CHANGE PROCESSES -----------------------
//--------------------------------------------------------------------
// For C++ jump process ID X, search for IDX
// For Matlab jump process X, search for ARX
// ************************************************************
// ID1- R6 growth on AC (AR1 on Matlab)
// ************************************************************
void PAHProcess::proc_G6R_AC(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSitesMemb(stt);
    //printStruct();//++++
	//Optimise PAH if needed.
	Spointer S1 = moveIt(stt, -1); 
    Spointer S2 = moveIt(stt, 1);
	Cpointer Ccheck = C_1;
	Cpointer Ccheck2 = Ccheck->C2;
	for (int i=1;i<=3;i++){
		if (getDistance_twoC(Ccheck, Ccheck2) > 1.75 && !(m_pah->m_optimised)){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		if (Ccheck2->bridge && Ccheck2->C3 != Ccheck) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C3;
		}
		else {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
	}	
	if (getDistance_twoC(C_1,C_2) < 2.6 || (int)S1->type % 10 >=4 || (int)S2->type % 10 >=4) {
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
	}
    Cpointer newC1;
    Cpointer newC2;
    if(checkHindrance_newCposition(C_1) || checkHindrance_newCposition(C_2)) {
        /*cout<<"Site hindered, process not performed.\n"*/ return;}
    if(!(C_1->C2->bridge) || !(C_2->C1->bridge)) { // check if bulk C in AC site is a bridge
        // Add and remove C
		cpair FEvector = get_vector(C_1->C2->coords,C_2->C1->coords);
		cpair AGvector = get_vector(C_1->C2->coords,C_1->coords);
		cpair start_direction = add_vector(FEvector, AGvector);
		cpair AJvector = get_vector(C_2->C1->coords,C_2->coords);
		cpair lat_direction = add_vector(AJvector, invert_vector(AGvector));
		double dist = getDistance_twoC(C_1,C_2);
        removeC(C_1->C2, true);
        removeC(C_2->C1, true);
		
		newC1 = addC(C_1, start_direction, dist/2.0);
		updateA(newC1, 'H', AGvector);
		updateA(C_1, 'C', AGvector);
		newC2 = addC(newC1, lat_direction, dist/2.0);
		//newC2 = addC(newC1, FEvector, 1.4);
		updateA(newC2, 'H', AJvector);
		updateA(C_2, 'C', AJvector);
    }else {
		cpair FEvector = get_vector(C_1->C2->coords,C_2->C1->coords);
		cpair AGvector = get_vector(C_1->C2->coords,C_1->coords);
		cpair start_direction = add_vector(FEvector, AGvector);
		cpair AJvector = get_vector(C_2->C1->coords,C_2->coords);
		cpair lat_direction = add_vector(AJvector, invert_vector(AGvector));
		double dist = getDistance_twoC(C_1,C_2);
		// The best assumption is that the PAH over the bridge is orthogonal to the plane.
		// Needs debugging!!
        // update bridges info, both are no longer bridges
        C_1->C2->C1 = C_1->C2->C3;// update neighbour
        C_2->C1->C2 = C_2->C1->C3;
        C_1->C2->bridge = false;
        C_2->C1->bridge = false;
        connectToC(C_1, C_2);
        // Add C
		newC1 = addC(C_1, start_direction, dist/2.0);
		updateA(newC1, 'H', AGvector);
		updateA(C_1, 'C', AGvector);
		newC2 = addC(newC1, lat_direction, dist/2.0);
		//newC2 = addC(newC1, FEvector, 1.4);
		updateA(newC2, 'H', AJvector);
		updateA(C_2, 'C', AJvector);
    }
    //printStruct();
    // neighbouring sites:
	Spointer S3 = moveIt(S1, -1);
	Spointer S4 = moveIt(S2, 1);
    // Update Site and neighbours
	if ( stt->type == RFER ){
		convSiteType(stt, newC1, newC2, FE);
		convSiteType(S1, S1->C1, newC1, R5R6); // neighbours
		convSiteType(S2, newC2, S2->C2, R5R6);
		if ((int)S3->type < 2000) updateSites(S3, S3->C1, S3->C2, +400); // neighbours of neighbours
		if ((int)S4->type < 2000) updateSites(S4, S4->C1, S4->C2, +400);
	}
	else if ( stt->type == R5R6FER ){
		convSiteType(stt, newC1, newC2, FE);
		if (S1->type==R5) {
			convSiteType(S1, S1->C1, newC1, R5R6); // neighbours
			if ((int)S2->type < 610) updateSites(S2, newC2, S2->C2, +1501);
			else if ((int)S2->type < 2000) updateSites(S2, newC2, S2->C2, +1101);
			else updateSites(S2, newC2, S2->C2, +1);
		}
		else if (S2->type==R5){
			convSiteType(S2, newC2, S2->C2, R5R6);
			if ((int)S1->type < 610) updateSites(S1, S1->C1, newC1, +1501);
			else if ((int)S1->type < 2000) updateSites(S1, S1->C1, newC1, +1101);
			else updateSites(S1, S1->C1, newC1, +1);
		}
		else {
			cout<<"R5R6FER not next to R5. Error.\n";
			saveXYZ("KMC_DEBUG/R5R6FER_error");
		}
	}
	else if ( stt->type == R5R6FER5R6 ){
		convSiteType(stt, newC1, newC2, FE);
		if ((int)S2->type < 610) updateSites(S2, newC2, S2->C2, +1501);
		else if ((int)S2->type < 2000) updateSites(S2, newC2, S2->C2, +1101);
		else updateSites(S2, newC2, S2->C2, +1);
		if ((int)S1->type < 610) updateSites(S1, S1->C1, newC1, +1501);
		else if ((int)S1->type < 2000) updateSites(S1, S1->C1, newC1, +1101);
		else updateSites(S1, S1->C1, newC1, +1);
	}
	else{
		updateSites(stt, newC1, newC2, -2);
		updateSites(S1, S1->C1, newC1, 1); // neighbours
		updateSites(S2, newC2, S2->C2, 1);
	}
    // Update combined site for Site and neighbours
    updateCombinedSites(stt);
    updateCombinedSites(S1); updateCombinedSites(S2);
    updateCombinedSites(S3); updateCombinedSites(S4);
    // add ring counts
    m_pah->m_rings++;
    //printSites(stt);
	//Optimise PAH if needed.
	if (getDistance_twoC(newC2,C_2) > 1.7 && !m_pah->m_optimised) {
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
}
// 
// ************************************************************
// ID2- R6 growth on FE (AR2 on Matlab)
// ************************************************************
void PAHProcess::proc_G6R_FE(Spointer& stt, Cpointer C_1, Cpointer C_2) {
//    printSites(stt);
    Cpointer newC1, newC2, newC3, newC4;
	Spointer S1 = moveIt(stt, -1);
	Spointer S2 = moveIt(stt, 1);
    if(checkHindrance_newCposition(C_1) || checkHindrance_newCposition(C_2)) {
        /*cout<<"Site hindered, process not performed.\n"*/ return;}
    // Add C
	cpair start_direction = C_1->growth_vector;
	newC1 = addC(C_1, start_direction, 1.4);
	updateA(newC1, 'H', get_vector(C_2->coords, C_1->coords));
	updateA(C_1, 'C', C_1->growth_vector);
	newC2 = addC(newC1, get_vector(C_1->C1->coords, C_1->coords), 1.4);
	updateA(newC2, 'H', start_direction);
	newC3 = addC(newC2, get_vector(C_1->coords, C_2->coords), 1.4);
	updateA(newC3, 'H', C_2->growth_vector);
	newC4 = addC(newC3, get_vector(C_2->coords, C_2->C2->coords), 1.4);
	updateA(newC4, 'H', get_vector(C_1->coords, C_2->coords));
	updateA(C_2, 'C', start_direction);
	
	/*newC1 = addC(C_1, normAngle(C_1->bondAngle1 + 120), 0, 1.4);
	newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1 - 60), 0, 1.4);
	newC3 = addC(newC2, normAngle(newC2->C1->bondAngle1 - 60), 0, 1.4);
	newC4 = addC(newC3, normAngle(newC3->C1->bondAngle1 - 60), normAngle(newC3->C1->bondAngle1 - 120), 1.4);*/
    //C_1->C3 = C_2; C_2->C3 = C_1;
    // Add and remove H
    //updateA(C_1->C1, C_2->C2, 'H');
    // neighbouring sites:
    //Spointer S1 = moveIt(stt, -1); 
    //Spointer S2 = moveIt(stt, 1);
    // Update Site and neighbours
    updateSites(stt, newC2, newC3, 0);
    updateSites(S1, S1->C1, newC1, 1);
    updateSites(S2, newC4, S2->C2, 1);
    // Add new Sites
    Spointer newS1 = addSite(FE, newC1, newC2, stt);
    Spointer newS2 = addSite(FE, newC3, newC4, S2);
    // Update combined sites for all new sites and original neighbours
    Spointer S3, S4;
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt);
    updateCombinedSites(newS1); updateCombinedSites(newS2); // new sites
    updateCombinedSites(S1); updateCombinedSites(S2); // original neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours of neighbours
    // Add H count
    addCount(0, 2);
    // add ring counts
    m_pah->m_rings++;
	//Optimise if needed
	if (getDistance_twoC(newC4, C_2) > 1.59 && !m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	} 
}
// 
// ************************************************************
// ID3- BY6 closure reaction (AR14 on Matlab)
// ************************************************************
void PAHProcess::proc_L6_BY6(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSites(stt);
    // Remove C
	int ntype_site = (int)stt->type;
    Cpointer now = C_1->C2;
    do{
        Cpointer next;
        if(!now->bridge) {
            next = now->C2;
            //if(now->C3 != NULL) now->C3->C3 = NULL;
			if (ntype_site > 2000) moveC_z(now, 0.1);
            removeC(now, true);
            now = next;
        }
        else {
			m_pah->m_optimised = false;
            next = now->C3->C2; //
            // bridged bulk C will not be removed from edge
            Cpointer b = now->C3; // C bridged to now
            b->C3 = NULL; now->C3 = NULL;
            b->bridge = false; now->bridge = false;
            //b->bondAngle1 = b->bondAngle2;
            //now->bondAngle2 = 0; b->bondAngle2 = 0;
            // connect C_1 to next and the two bulk atoms still remaining
            connectToC(C_1, next);
            connectToC(b, now);
            now = next;
        }
    }while(now!=C_2);
	updateA(C_1, 'C', C_1->growth_vector);
	updateA(C_2, 'C', C_2->growth_vector);
    // Change bond angle between C_1 and C_2
    //C_1->bondAngle1 = normAngle(C_1->bondAngle1+120);
    // Add and remove H
    //updateA(C_1->C1, C_2->C2, 'H');

    //// Remove BY6 site and combine the neighbouring sites. 
    //// First remove all three from site map. Elementary site types first..
    //delSiteFromMap(moveIt(stt, -1)->type, moveIt(stt, -1));
    //delSiteFromMap(moveIt(stt, 1)->type, moveIt(stt, 1));
    //delSiteFromMap(stt->type, stt);
    //// then for combined site types..
    //delSiteFromMap(moveIt(stt, -1)->comb, moveIt(stt, -1));
    //delSiteFromMap(moveIt(stt, 1)->comb, moveIt(stt, 1));
    //delSiteFromMap(stt->comb, stt);
    // Convert the BY6 site into the resulting site after reaction,
    // finding resulting site type:
	int ntype1 = (int) moveIt(stt, -1)->type;
    int ntype2 = (int) moveIt(stt, 1)->type;
	int newType;
    if(ntype1 < 5 && ntype2 < 5) {
		// A bay with two basic neighbour sites will become a basic neighbour type. 
		// The R5 and R7 balances are handled in the individual jump processes.
        newType = (ntype1+ntype2+2);
        // convert site
        if(newType>4) {
            //saveDOT(std::string("BY6ClosureProblem.dot"));
            //std::cerr<<"ERROR: newType is > 4 (PAHProcess::proc_L6_BY6)\n";
			//saveXYZ("KMC_DEBUG/BY6ClosureProblem");
        }
		convSiteType(stt, stt->C1, stt->C2, (kmcSiteType)4);
        updateSites(stt, moveIt(stt,-1)->C1, moveIt(stt,1)->C2,(newType-4));
    }
    else {
		int new_point = 0;
		if (ntype_site == 4) new_point = 2;
		else if (ntype_site == 104) {
			new_point = 502;
			if (ntype1 == 100) {
				ntype1 = 0;
				Spointer S1 = moveIt(stt, -2);
				int stype1 = (int)S1->type;
				if (stype1 <= 1010) stype1 += 400;
				convSiteType(S1, S1->C1, S1->C2, (kmcSiteType)stype1);
				if ( ntype2 < 2000) new_point = 502;
				else new_point = 102;
			}
			else {
				ntype2 = 0;
				Spointer S2 = moveIt(stt, 2);
				int stype2 = (int)S2->type;
				if (stype2 <= 1010) stype2 += 400;
				convSiteType(S2, S2->C1, S2->C2, (kmcSiteType)stype2);
				if ( ntype1 < 2000) new_point = 502;
				else new_point = 102;
			}
		}
		else if (ntype_site == 204) {
			new_point = 1002; ntype1 = 0; ntype2 = 0;
			Spointer S1 = moveIt(stt, -2); Spointer S2 = moveIt(stt, 2);
			int stype1 = (int)S1->type; int stype2 = (int)S2->type;
			if (stype1 <= 1010) stype1 += 400;
			if (stype2 <= 1010) stype2 += 400;
			convSiteType(S1, S1->C1, S1->C2, (kmcSiteType)stype1);
			convSiteType(S2, S2->C1, S2->C2, (kmcSiteType)stype2);
		}
		else if (ntype_site == 504) {
			new_point = 2003;
			if (ntype1 >= 2103 && ntype2 < 501) {
				//The R5R6BY5 site is next to a complex bay site.
				new_point = ntype1 + 1;
				ntype1= 1;
			}
			else if (ntype2 >= 2103 && ntype1 <501) {
				//The R5R6BY5 site is next to a complex bay site.
				new_point = ntype2 + 1;
				ntype2= 1;
			}
			else if (ntype2 >= 501 && ntype1 >=501) {
				//The R5R6BY5 site is next to a complex bay site.
				new_point += 101;
				if (ntype1 >= 501 && ntype1 <= 604) ntype1 -= 501;
				else if (ntype1 >= 1002 && ntype1 <= 1004) ntype1 -= 901;
				if (ntype2 >= 501 && ntype2 <= 604) ntype2 -= 501;
				else if (ntype2 >= 1002 && ntype2 <= 1004) ntype2 -= 901;
			}
			else {
				if (ntype1 >= 501 && ntype1 <= 604) ntype1 -= 501;
				else if (ntype2 >= 501 && ntype2 <= 604) ntype2 -= 501;
				else if (ntype1 >= 1002 && ntype1 <= 1004) ntype1 -= 901;
				else if (ntype2 >= 1002 && ntype2 <= 1004) ntype2 -= 901;
			}
		}
		else if (ntype_site == 604) {
			new_point = 2103;
			if (ntype1 == 100) {
				Spointer S1 = moveIt(stt, -2);
				int stype1 = (int)S1->type;
				if (stype1 <= 1010) stype1 += 400;
				convSiteType(S1, S1->C1, S1->C2, (kmcSiteType)stype1);
				ntype1 = 0;
				if (ntype2 >= 501 && ntype2 <= 604) ntype2 -= 501;
				else if (ntype2 >= 1002 && ntype2 <= 1004) ntype2 -= 901;
			}
			else {
				Spointer S2 = moveIt(stt, 2);
				int stype2 = (int)S2->type;
				if (stype2 <= 1010) stype2 += 400;
				convSiteType(S2, S2->C1, S2->C2, (kmcSiteType)stype2);
				ntype2 = 0;
				if (ntype1 >= 501 && ntype1 <= 604) ntype1 -= 501;
				else if (ntype1 >= 1002 && ntype1 <= 1004) ntype1 -= 901;
			}
		}
		else if (ntype_site == 1004) {
			new_point = 2102;
			ntype1 = ntype1 % 10;
			ntype2 = ntype2 % 10;
			/*if (ntype1 >= 501 && ntype1 <= 504) ntype1 -= 501;
			else if (ntype2 >= 501 && ntype1 <= 504) ntype2 -= 501;*/
		}
		else if (ntype_site == 2014 || ntype_site == 2004) {
			new_point = 2;
			ntype1 = ntype1;
			ntype2 = ntype2;
			/*if (ntype1 >= 501 && ntype1 <= 504) ntype1 -= 501;
			else if (ntype2 >= 501 && ntype1 <= 504) ntype2 -= 501;*/
		}
		else if (ntype_site == 2005 || ntype_site == 2015) {
			new_point = 2;
			ntype1 = ntype1;
			ntype2 = ntype2;
		}
		else if (ntype_site == 2104 || ntype_site == 2105 || ntype_site == 2205) {
			//A site type 2104 has an R5 member to only one side. Detect where is the site.
			new_point = 0;
			Spointer S1 = moveIt(stt, -1);
			Spointer S2 = moveIt(stt, +1);
			if (ntype1 == 100 || ntype2 == 100){
				if (ntype1 == 100 && ntype2 == 100) {
					if (ntype_site == 2104)	{
						cout << "Error. R5FEACR5 site cannot have R5 sites to both sides.\n";
						saveXYZ("KMC_DEBUG/R5FEACR5_error");
						//return;
					}
					else if (ntype_site == 2105)	 {
						cout << "Error. R5ZZACR5 site cannot have R5 sites to both sides.\n";
						saveXYZ("KMC_DEBUG/R5ZZACR5_error");
						//return;
					}
					else {
						cout << "Error. ACR5FER site cannot have R5 sites to both sides.\n";
						saveXYZ("KMC_DEBUG/ACR5FER_error");
						//return;
					}
					
				}
				else if (ntype1 == 100){
					//R5 to the left
					Spointer S3 = moveIt(S1, -1);
					if ((int)S3->type < 2000) updateSites(S3, S3->C1, S3->C2, +400);
					if ( (int)S2->type < 2000) {
						ntype1 = 502;
						ntype2 = ntype2;
					}
					else { 
						ntype1 = 102;
						ntype2 = ntype2;
					}
				}
				else {
					//R5 to the right
					Spointer S4 = moveIt(S2, +1);
					if ((int)S4->type < 2000) updateSites(S4, S4->C1, S4->C2, +400);
					if ( (int)S1->type < 2000) {
						ntype1 = ntype1;
						ntype2 = 502;
					}
					else { 
						ntype1 = ntype1;
						ntype2 = 102;
					}
				}
			}
			else if ( ((int)S2->type >= 501 && (int)S2->type <= 504)  &&  ((int)S1->type >= 501 && (int)S1->type <= 504) ) {
				//Both sites are R5R6 neighbours. Detect which one is has an R5 on the R5ACR5FE site.
				new_point = 2102;
				ntype1 = ntype1 % 10;
				ntype2 = ntype2 % 10;
			}
			else if ( !((int)S2->type >= 501 && (int)S2->type <= 504)  &&  ((int)S1->type >= 501 && (int)S1->type <= 504) ) {
				//R5R6 neighbour to the left. 
				new_point = 2002;
				ntype1 = ntype1 % 10;
				ntype2 = ntype2;
			}
			else if ( ((int)S2->type >= 501 && (int)S2->type <= 504)  &&  !((int)S1->type >= 501 && (int)S1->type <= 504) ) {
				//R5R6 neighbour to the right. 
				new_point = 2002;
				ntype1 = ntype1;
				ntype2 = ntype2 % 10;
			}
			else if ( ((int)S2->type >= 602 && (int)S2->type <= 1004)  &&  ((int)S1->type >= 602 && (int)S1->type <= 1004) ) {
				//Both sites are R5R6 neighbours. Detect which one is has an R5 on the R5ACR5FE site.
				new_point = 2202;
				ntype1 = ntype1 % 10;
				ntype2 = ntype2 % 10;
			}
			else if ( !((int)S2->type >= 602 && (int)S2->type <= 1004)  &&  ((int)S1->type >= 602 && (int)S1->type <= 1004) ) {
				//R5R6 neighbour to the left. 
				new_point = 2102;
				ntype1 = ntype1 % 10;
				ntype2 = ntype2;
			}
			else if ( ((int)S2->type >= 602 && (int)S2->type <= 1004)  &&  !((int)S1->type >= 602 && (int)S1->type <= 1004) ) {
				//R5R6 neighbour to the right. 
				new_point = 2102;
				ntype1 = ntype1;
				ntype2 = ntype2 % 10;
			}
			else if ( ((int)S2->type >2000) ){
				ntype2 = ntype2;
				new_point = 2;
				if ((int)S1->type >= 501 && (int)S1->type <= 604) ntype1 = ntype1 - 501;
				else ntype1 = ntype1; //This is very likely to form a SPIRAL.
			}
			else if ( ((int)S1->type >2000) ){
				ntype1 = ntype1;
				new_point = 2;
				if ((int)S2->type >= 501 && (int)S2->type <= 604) ntype2 = ntype2 - 501;
				else ntype2 = ntype2; //This is very likely to form a SPIRAL.
			}
			else {
				cout << "Error. No R5 or R5R6 neighbour sites to R5FEACR5. \n";
				cout << "Saving file KMC_DEBUG/R5FEACR5_error \n";
				saveXYZ("KMC_DEBUG/R5FEACR5_error");
				//return;
			}
		}
		else if (ntype_site == 2114 || ntype_site == 2115) {
			new_point = 2;
			ntype1 = ntype1;
			ntype2 = ntype2;
			/*if (ntype1 >= 501 && ntype1 <= 504) ntype1 -= 501;
			else if (ntype2 >= 501 && ntype1 <= 504) ntype2 -= 501;*/
		}
		else if (ntype_site == 2204) {
			//A site type 2204 has R5 members to both sides.
			Spointer S1 = moveIt(stt, -1);
			Spointer S2 = moveIt(stt, +1);
			new_point = 0;
			if ((int)S1->type>2000 || (int)S2->type>2000 ){
				new_point = 0;
				if ((int)S1->type>2000 && (int)S2->type>2000 ){
					//Almost certain this will be an SPIRAL
					new_point = 0;
					ntype1=ntype1;
					ntype2=ntype2;
				}
				else if ((int)S1->type>2000){
					new_point = 0;
					ntype1=ntype1;
					if (S2->type == R5){
						Spointer S4 = moveIt(S2, +1);
						if ( (int)S4->type < 2000) updateSites(S4, S4->C1, S4->C2, +400);
						ntype2 = 102;
					}
					else ntype2 = (ntype2 % 10) + 101;
				}
				else if ((int)S2->type>2000){
					new_point = 0;
					ntype2=ntype2;
					if (S1->type == R5){
						Spointer S3 = moveIt(S1, -1);
						if ( (int)S3->type < 2000) updateSites(S3, S3->C1, S3->C2, +400);
						ntype1 = 102;
					}
					else ntype1 = (ntype1 % 10) + 101;
				}
			}
			else if ( ((int)S2->type >= 501 && (int)S2->type <= 1004) && ((int)S1->type >= 501 && (int)S1->type <= 1004) ){
				new_point = 2102;
				if (ntype1 > 600) ntype1 = (ntype1 % 10) + 100;
				else ntype1 = (ntype1 % 10);
				if (ntype2 > 600) ntype2 = (ntype2 % 10) + 100;
				else ntype2 = (ntype2 % 10);
			}
			else if((int)S1->type >= 501 && (int)S1->type <= 1004) {
				if (S2->type != R5){
					new_point = 0;
					if (ntype1 > 600) ntype1 = (ntype1 % 10) + 2102;
					else ntype1 = (ntype1 % 10) + 2002;
				}
				else{
					Spointer S4 = moveIt(S2, +1);
					if ( (int)S4->type < 2000) updateSites(S4, S4->C1, S4->C2, +400);
					ntype2 = 100;
					if (ntype1 > 600) ntype1 = (ntype1 % 10) + 2102;
					else ntype1 = (ntype1 % 10) + 2002;
				}
			}
			else if((int)S2->type >= 501 && (int)S2->type <= 1004) {
				if (S1->type != R5){
					new_point = 0;
					if (ntype2 > 600) ntype2 = (ntype2 % 10) + 2102;
					else ntype2 = (ntype2 % 10) + 2002;
				}
				else{
					Spointer S3 = moveIt(S1, -1);
					if ( (int)S3->type < 2000) updateSites(S3, S3->C1, S3->C2, +400);
					ntype1 = 100;
					if (ntype2 > 600) ntype2 = (ntype2 % 10) + 2102;
					else ntype2 = (ntype2 % 10) + 2002;
				}
			}
			else if (S1->type == R5 || S2->type ==R5){
				if (S1->type == R5){
					Spointer S3 = moveIt(S1, -1);
					if ( (int)S3->type < 2000) updateSites(S3, S3->C1, S3->C2, +400);
					ntype1 = 501;
				}
				else ntype1 = ntype1 % 10;
				if (S2->type == R5){
					Spointer S4 = moveIt(S2, +1);
					if ( (int)S4->type < 2000) updateSites(S4, S4->C1, S4->C2, +400);
					ntype2 = 501;
				}
				else ntype2 = ntype2 % 10;
			}
		}
		
		newType = (new_point + ntype1 + ntype2);
		/*if (ntype1 >= 501 && ntype1 <= 504) ntype1 -= 501;
		else if (ntype2 >= 501 && ntype1 <= 504) ntype2 -= 501;*/
		if (newType == 9999 || ntype1 == 9999 || ntype2 == 9999){
			//Neighbour is a SPIRAL, convert product site in SPIRAL
			newType = 9999;
		}
		if ((kmcSiteType)newType == None) {
			//saveDOT(std::string("BY6ClosureProblem.dot"));
			std::cerr << "ERROR: newType is None (PAHProcess::proc_L6_BY6)\n";
		}
		
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);

		/*
        int newType = 0;
		if ((ntype1 - 56) >= 0) { // R5R6 neighbour to one side.
			ntype1 -= 56;
			newType = 56;
		}
        else if((ntype1-6)>=0) { // R5 neighbour to one side.
            ntype1 -= 6;
            newType = 6;
        }
		if ((ntype2 - 56) >= 0) { // R5R6 neighbour to the other side.
			ntype2 -= 56;
			if (newType > 0 && newType < 100) {
				newType = 66; // R5R6 neighbour to one side and R5 neighbour to the other side. 
				//saveDOT(std::string("BY6ClosureProblem.dot"));
				//std::cerr << "Problem with Bay Closure\n";
			}
			else if (newType > 100) {
				newType = 60; // R5R6 neighbours to both sides. Resulting site not defined.
				saveDOT(std::string("BY6ClosureProblem.dot"));
				std::cerr << "Problem with Bay Closure\n";
			}
			else newType = 56; // R5R6 neighbour to the second side only.
		}
		else if ((ntype2 - 6) >= 0) { // R5 neighbour to the other side
			ntype2 -= 6;
			if (newType>0 && newType<100) newType = 10; // R5 neighbours to both sides.
			else if (newType > 100) {
				newType = 66; // R5R6 neighbour to one side and R5 neighbour to the other side. 
				//saveDOT(std::string("BY6ClosureProblem.dot"));
				//std::cerr << "Problem with Bay Closure\n";
			}
			else newType = 6; // R5 neighbour to the other side only.
		}
        newType += ntype1+ntype2+2;
        convSiteType(stt, moveIt(stt,-1)->C1, moveIt(stt,1)->C2, (kmcSiteType) newType);*/
    }    
    // erase the existence of the neighbouring sites
    Spointer Srem1 = moveIt(stt,-1);
    Spointer Srem2 = moveIt(stt, 1);
    removeSite(Srem1);
    removeSite(Srem2);
    // update combined sites and neighbours
    Spointer S1 = moveIt(stt,-1); Spointer S2 = moveIt(stt,1);
    //Spointer S3 = moveIt(S1,-1); Spointer S4 = moveIt(S2,1);
    updateCombinedSites(stt);
    updateCombinedSites(S1); updateCombinedSites(S2);
    //updateCombinedSites(S3); updateCombinedSites(S4);

    //printSites(stt);
    // update H count
    addCount(0,-2);
    // add ring counts
    m_pah->m_rings++;
	//Optimise PAH if needed.
	if (getDistance_twoC(C_1,C_2) > 1.6 || ntype_site > 2000) {
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol, 8000);
			passbackPAH(mol);
		}
	}
}
// 
// ************************************************************
// ID4- phenyl addition (AR15 in Matlab)
// ************************************************************
void PAHProcess::proc_PH_benz(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
    Cpointer chosen;
    bool before; // true if C1 of site is chosen, false if C2
    // choose one of the C atoms if site type is FE/AC/ZZ
    if(stt->type == FE || stt->type == AC || stt->type == ZZ) {
        // Define a distribution that has two equally probably outcomes
        boost::bernoulli_distribution<> choiceDistrib;
        // Now build an object that will generate a sample using rng
        boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
        if(choiceGenerator()) {
            chosen = C_1;
			Spointer s_adjacent = moveIt(stt, -1);
			if ( (int)s_adjacent->type %10 >=3) return;
            before=true;
        }
        else {
            chosen = C_2;
			Spointer s_adjacent = moveIt(stt, +1);
			if ( (int)s_adjacent->type %10 >=3) return;
            before=false;
        }
    }
	else if (stt->type == BY5 || stt->type == BY6) return; // Jump process not allowed for BY5 or BY6. YET!!!!
    else { //if not, choose the C which is not part of a 5-membered ring
        if(moveIt(stt,-1)->type == R5) {
            chosen = C_2;
            before=false;
        }
        else {
            chosen = C_1;
            before=true;
        }
    }
    // neighbouring site to be updated:
    Spointer neighbour;
    if(before) {
        neighbour = moveIt(stt,-1);
    }
    else {
        neighbour = moveIt(stt,1);
    }
    // check hindrance
    if(checkHindrancePhenyl(chosen)) return;
    // add C atoms
    Cpointer newC;
	cpair vec1 = get_vector(chosen->C1->coords, chosen->coords);
	cpair vec2 = get_vector(chosen->C2->coords, chosen->coords);
	cpair vec3 = chosen->growth_vector;
	//cpair ivec1 = invert_vector(vec1);
	cpair ivec2 = invert_vector(vec2);
	cpair ivec3 = invert_vector(vec3);
    newC = bridgeC(chosen);
	newC = addC(newC, vec2, 1.4); 
	updateA(newC, 'H', invert_vector(vec1));
	updateA(newC->C1,'C',newC->C1->growth_vector);
	newC = addC(newC, vec3, 1.4); 
	updateA(newC, 'H', vec2);
	newC = addC(newC, vec1, 1.4); 
	updateA(newC, 'H', vec3);
	newC = addC(newC, ivec2, 1.4);
	updateA(newC, 'H', vec1);
	newC = addC(newC, ivec3, 1.4);
	updateA(newC, 'H', ivec2);
	
	/*newC = addC(newC, normAngle(chosen->bondAngle2 + 60), 0, 1.4);
	newC = addC(newC, normAngle(newC->C1->bondAngle1 - 60), 0, 1.4);
	newC = addC(newC, normAngle(newC->C1->bondAngle1 - 60), 0, 1.4);
	newC = addC(newC, normAngle(newC->C1->bondAngle1 - 60), 0, 1.4);
	newC = addC(newC, normAngle(newC->C1->bondAngle1 - 60), normAngle(newC->C1->bondAngle1 - 120), 1.4);*/
    connectToC(newC, chosen->C3);
    // which site do we add the new sites before?
    Spointer addBefore;
    if(before) addBefore = stt;
    else addBefore = neighbour;
    // Add new sites (4 new FE sites)
    Cpointer Cnow = chosen->C3->C2;
    Spointer nS1 = addSite(FE, Cnow, Cnow->C2, addBefore); Cnow = Cnow->C2;
    Spointer nS2 = addSite(FE, Cnow, Cnow->C2, addBefore); Cnow = Cnow->C2;
    Spointer nS3 = addSite(FE, Cnow, Cnow->C2, addBefore); Cnow = Cnow->C2;
    Spointer nS4 = addSite(FE, Cnow, Cnow->C2, addBefore); Cnow = Cnow->C2;
    // add and remove H
    //updateA(chosen->C1, chosen->C2, 'H');
    // update sites and neighbour
    // new member C for stt and neighbour
    Cpointer s_C1, s_C2, n_C1, n_C2;
    Spointer n1, n2; // new neighbours for updated sites
    if(before) {
        s_C1 = moveIt(stt, -1)->C2;
        s_C2 = stt->C2;
        n_C1 = neighbour->C1;
        n_C2 = moveIt(neighbour, 1)->C1;
        n1 = moveIt(stt, 1);
        n2 = moveIt(neighbour, -1);
    }else {
        s_C1 = stt->C1;
        s_C2 = moveIt(stt, 1)->C1;
        n_C1 = moveIt(neighbour, -1)->C2;
        n_C2 = neighbour->C2;
        n1 = moveIt(stt, -1);
        n2 = moveIt(neighbour, 1);
    }
    updateSites(stt, s_C1, s_C2, 2);
    updateSites(neighbour, n_C1, n_C2, 2);
    // update combined sites for all new sites and neighbours (and their neighbours)
    updateCombinedSites(stt); updateCombinedSites(neighbour);
    updateCombinedSites(nS1); updateCombinedSites(nS2); updateCombinedSites(nS3); updateCombinedSites(nS4);
    updateCombinedSites(n1); updateCombinedSites(n2);
    // update H count
    addCount(0,4);
    // add ring counts
    m_pah->m_rings++;
}
// 
// ************************************************************
// ID5- R6 desorption at FE (AR8 in Matlab)
// ************************************************************
void PAHProcess::proc_D6R_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSites(stt);
    // cannot happen if R6 is next to bridge
    if(stt->C1->C1->C1->bridge || stt->C2->C2->C2->bridge) return;
	// cannot happen if neighbour is a pentagon.
	Spointer checkR5_1 = moveIt(stt,-2);
	Spointer checkR5_2 = moveIt(stt,+2);
	if ( (int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 || (int)checkR5_1->type == 2002) return;
	if ( (int)checkR5_2->type == 101 || (int)checkR5_2->type == 501 || (int)checkR5_2->type == 2002) return;
    // member C atoms of resulting FE site
    Cpointer C1_new, C2_new;
    C1_new = C_1->C1->C1;
    C2_new = C_2->C2->C2;
    // delete neighbouring FE sites from site map. Identify them first:
    Spointer S1 = moveIt(stt,-1);
    Spointer S2 = moveIt(stt,1);
    // then remove them
    removeSite(S1);
    removeSite(S2);
	// Resulting H vectors
	cpair Hdir1 = get_vector(C1_new->coords, C1_new->C2->coords);
	cpair Hdir2 = get_vector(C2_new->coords, C2_new->C1->coords);
	// remove C
	for (int i = 0; i<4; i++) {
		removeC(C1_new->C2, false);
	}
	updateA(C1_new, 'H', Hdir1);
	updateA(C2_new, 'H', Hdir2);
    // update site stt and new neighbours
    // new neighbours:
    S1 = moveIt(stt,-1); S2 = moveIt(stt,1);
	if (isR5internal(C1_new, C2_new) ){
		//FE3 desorption allows pentagon going to the edge
		if ( (int)S1->type >= 2000 && (int)S1->type <= 2100 ) updateSites(S1, S1->C1, C1_new, -1901);
		else if ( (int)S1->type >= 2100 ) {
			Spointer S3 = moveIt(S1, -1);
			if ( (int)S3->type == 100) updateSites(S1, S1->C1, C1_new, -1901);
			else updateSites(S1, S1->C1, C1_new, -1501);
		}
		else {
			cout << "FE3 desorption next to an R5 with wrong neightbour site. S1 site type = " << S1->type << "\n";
			if (m_debug_pah){
				ifstream src("KMC_DEBUG/BEFORE.xyz");
				std::string filename = "KMC_DEBUG/BEFORE_FE3-R5-S1_error";
				filename.append(".xyz");
				ofstream dst(filename);
				dst << src.rdbuf();
				src.close();
				dst.close();
				cout << "Saving file KMC_DEBUG/BEFORE_FE3-R5-S1_error.\n";
			}
			saveXYZ("KMC_DEBUG/FE3-R5-S1 desorption_error");
			cout << "Saving file KMC_DEBUG/FE3-R5-S1 desorption_error.\n"; 
			printSites(S1);
		}
		if ( (int)S2->type >= 2000 && (int)S2->type <= 2100 ) updateSites(S2, C2_new, S2->C2, -1901);
		else if ( (int)S2->type >= 2100 ) {
			Spointer S4 = moveIt(S2, +1);
			if ( (int)S4->type == 100) updateSites(S2, C2_new, S2->C2, -1901);
			else updateSites(S2, C2_new, S2->C2, -1501);
		}
		else {
			cout << "FE3 desorption next to an R5 with wrong neightbour site. S2 site type = " << S2->type << "\n";
			if(m_debug_pah){
				ifstream src("KMC_DEBUG/BEFORE.xyz");
				std::string filename = "KMC_DEBUG/BEFORE_FE3-R5-S2_error";
				filename.append(".xyz");
				ofstream dst(filename);
				dst << src.rdbuf();
				src.close();
				dst.close();
				cout << "Saving file KMC_DEBUG/BEFORE_FE3-R5-S2_error.\n";
			}
			saveXYZ("KMC_DEBUG/FE3-R5-S2 desorption_error");
			cout << "Saving file KMC_DEBUG/FE3-R5-S2 desorption_error.\n";
			printSites(S2);
		}
		updateSites(stt, C1_new, C2_new, 100);
	}
	else {
		//Flat PAH desorption
		updateSites(stt, C1_new, C2_new, 0); // only change C1 and C2, site still FE
		updateSites(S1, S1->C1, C1_new, -1);
		updateSites(S2, C2_new, S2->C2, -1);
		// add H atoms
		//updateA(C1_new->C1, C2_new->C2, 'H');
	}
    // update combined sites
    Spointer S3, S4;
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt); // update resulting site
    updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
    
    // update H count
    addCount(0,-2);
    // add ring counts
    m_pah->m_rings--;
}
// ************************************************************
// ID6- R6 oxidation at FE by O2 (AR10 in Matlab)
// ************************************************************
void PAHProcess::proc_O6R_FE3_O2(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    proc_D6R_FE3(stt, C_1, C_2);
}
// ************************************************************
// ID7- R6 oxidation at FE by OH (AR11 in Matlab)
// ************************************************************
void PAHProcess::proc_O6R_FE3_OH(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    proc_D6R_FE3(stt, C_1, C_2);
}

// ************************************************************
// ID8- R6 oxidation at FE-HACA (AR12 in Matlab)
// ************************************************************
void PAHProcess::proc_O6R_FE_HACA(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//This jump process used to remove a full ring. However, the work by Singh 2016 showed that it probably forms an R5 which is more interesting.
	//Converts an FE_HACA into an R5.
    //printSites(stt);
    Spointer S1, S2;
	Cpointer CRem = C_2;
	Cpointer CRem_before = C_1;
	Cpointer CRem_next = C_2->C2;
    S1 = moveIt(stt, -1);
    S2 = moveIt(stt, 1);
	
	//check that two pentagons (including internals) will not collide
	if (m_pah->m_R5loc.size()>=1){
		cpair R5coords_end = endposR5internal(CRem_before, CRem);
		for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
			double distR5s = getDistance_twoC(*it, R5coords_end);
			if (distR5s < 2.8) {
				//This assumes that the R6 is adjacent to an R5. The correct route according to both Raj2013 and Singh2016 is the removal of two carbon atoms.
				//Structural effects on the oxidation of soot particles by O2: Experimentaland theoretical study. Raj2013
				// Oxidation of Graphene-Edge Six- and Five-Member Rings by Molecular Oxygen. Singh2016
				//This distance is a parameter of this jump process. Might need some more tuning. 
				//2.8 seems appropiate but may reject too many jumps.
				//Two pentagons will be next to each other violating the Isolated Pentagon Rule
				//return;
				//proc_O6R_FE_HACA_double(stt, C_1, C_2);
				return;
			}
		}
	}
	
    Cpointer C_bulk;
    for(C_bulk = S1->C1; C_bulk != S1->C2; C_bulk=C_bulk->C2) {
        if(C_bulk->bridge) return;
    }
    for(C_bulk = S2->C1; C_bulk != S2->C2; C_bulk=C_bulk->C2) {
        if(C_bulk->bridge) return;
    }
	//Get location of remaining carbon
	cpair newpos;
	cpair start_dir = get_vector(CRem_before->C1->coords, CRem->C2->coords);
	cpair normvec = invert_vector(norm_vector(CRem_before->coords, CRem->coords, CRem->C2->coords));
	cpair crossvec = cross_vector(start_dir, normvec);
	cpair Hdir = crossvec;
	bool optimise_flag = false;
	double dist = getDistance_twoC(CRem_before->C1,CRem->C2);
	if (dist < 3.4){
		dist = dist /2.0;
		double dist2 = sqrt(2.89 - dist * dist);
		newpos = std::make_tuple(std::get<0>(CRem_before->C1->coords) + dist * std::get<0>(start_dir) + dist2 * std::get<0>(crossvec), std::get<1>(CRem_before->C1->coords) + dist * std::get<1>(start_dir) + dist2 * std::get<1>(crossvec), std::get<2>(CRem_before->C1->coords) + dist * std::get<2>(start_dir) + dist2 * std::get<2>(crossvec));
	} else{
		dist = dist /2.0;
		newpos = std::make_tuple(std::get<0>(CRem_before->C1->coords) + dist * std::get<0>(start_dir) + dist/2.0 * std::get<0>(crossvec), std::get<1>(CRem_before->C1->coords) + dist * std::get<1>(start_dir)  + dist/2.0 * std::get<1>(crossvec), std::get<2>(CRem_before->C1->coords) + dist * std::get<2>(start_dir)  + dist/2.0 * std::get<2>(crossvec));
		optimise_flag = true;
	}
	//Find if ring is exposed on other side
	Cpointer thirdC = findThirdC(CRem_before);
	Cpointer thirdC2 = findThirdC(CRem_next);
	bool other_side_bool = false;
	Spointer other_side;
	if (thirdC != NULLC || thirdC2 != NULLC) {
		other_side_bool = true;
		if (thirdC == NULLC) other_side = findSite(thirdC2);
		else other_side = findSite(thirdC);
	}
	addR5internal(CRem_before,CRem_next, true);
	//Remove carbon
	removeC(CRem, false);
	//Move remaining carbon to new position
	moveC(CRem_before, newpos);
	updateA(CRem_before, 'H', Hdir);
	
	if (optimise_flag && !m_pah->m_optimised){
		//saveXYZ("KMC_DEBUG/Beforeoptimoxid");
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
    // Update Sites and neighbouring sites
	removeSite(stt);
	if ( (int)S1->type < 2000) updateSites(S1, S1->C1, CRem_before, +500);
	else updateSites(S1, S1->C1, CRem_before, +100);
	if ( (int)S2->type < 2000) updateSites(S2, CRem_before, S2->C2, +500);
	else updateSites(S2, CRem_next, S2->C2, +100);
	stt = S1;
	if (other_side_bool) {
		if (other_side != m_pah->m_siteList.end()) {
			updateSites(other_side, other_side->C1, other_side->C2, +2000);
			updateCombinedSites(other_side);
		}
	}
    Spointer S3, S4;
    // update combined sites for all sites and their neighbours
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(S1); updateCombinedSites(S2);
    updateCombinedSites(S3); updateCombinedSites(S4);
    addCount(0,-1);
    m_pah->m_rings--;
	m_pah->m_rings5_Lone++;
}


// ************************************************************
// ID8- R6 oxidation at FE-HACA (AR12 in Matlab)
// ************************************************************
void PAHProcess::proc_O6R_FE_HACA_double(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//This jump process removes a FE and assumes you form an AC or another site from the carbon inside.
	//This causes several issues because recognising what site will be formed is not trivial.
    //printSites(stt);
    // member C atoms of resulting AC site
    Cpointer C1_res, C2_res;
	//printStruct(C_1);
    C1_res = C_1->C1;
    C2_res = C_2->C2;
	Cpointer Ccheck = C1_res->C1;
	for (int ii = 0; ii<4; ii++){
		double Rdist = getDistance_twoC(Ccheck, Ccheck->C2);
		if (Rdist < 1.2 || Rdist > 1.6){
			if (!m_pah->m_optimised){
				OpenBabel::OBMol mol = passPAH();
				mol = optimisePAH(mol);
				passbackPAH(mol);
			}
		}
		Ccheck = Ccheck->C2;
	}
		
	cpair Cdir = get_vector(C_2->coords,C_2->C2->coords);
	cpair FEdir = get_vector(C_1->coords,C_2->coords);
	cpair Hdir1 = C_2->growth_vector;
	cpair Hdir2 = C_1->growth_vector;
	// Check for not allowed JPs
	Spointer S1 = moveIt(stt, -1);
	Spointer S2 = moveIt(stt, +1);
	if (S1->type== RFE || S2->type == RFE) return; //This would expose three edges of a pentagon. Not supported YET!
    // check if process will result in a bridge
	cpair pos = jumpToPos(C1_res->coords, Cdir, 1.4);
    //cpair pos = jumpToPos(C1_res->coords, normAngle(C1_res->bondAngle1-120), 0, 1.4);
	bool bridge = checkHindrance_C_PAH((pos));
	if (bridge) {
		Cpointer thirdC = findThirdC(C1_res);
		Spointer other = findSite(thirdC);
		if (other != m_pah->m_siteList.end()){
			if ((int)S1->type > 2000 || (int)S2->type > 2000 || S1->type== R5 || S2->type == R5 || S1->type== RFE || S2->type == RFE || (int)other->type > 2000) return; //This would expose three edges of a pentagon. Not supported YET!
		}
	}
	bool hept_bool = false;
	if  (isR7internal(C_1, C_2) ) hept_bool = true;
    // check if site is next to a bridge
    if (C1_res->bridge || C2_res->bridge) return;
    //if(bridge) return;
    // remove C
    removeC(C_1, false);
    removeC(C_2, false);
	if (bridge){
		Cpointer C_k = findC(pos);
		if (C_k == NULLC) {
			std::ostringstream msg;
            msg << "Error finding Carbon in Oxidation to bridge. \n";
			saveXYZ("KMC_DEBUG/Oxidation_to_bridge_error_find_carbon");
            //throw std::runtime_error(msg.str());
		}			
		Cpointer C_kk = C_k->C1;
		C1_res->C2 = C_k;
		C_k->C1 = C1_res;
		updateA(C1_res, 'H', Hdir1);
		//C1_res->bondAngle1 = normAngle(normAngle(C1_res->C1->bondAngle1) - 60);
		C_k->bridge = true;
		//C_k->bondAngle1 = normAngle(C1_res->bondAngle1 - 60);
		//C_k->bondAngle2 = normAngle(C1_res->bondAngle1 + 60);
		C_k->C3 = C_kk;
		C_kk->C2 = C2_res;
		C_kk->bridge = true;
		C_kk->C3 = C_k;
		//C_kk->bondAngle1 = normAngle(C_kk->C1->bondAngle1 - 60);
		//C_kk->bondAngle2 = normAngle(C_k->bondAngle2 - 180);
		C2_res->C1 = C_kk;
		updateA(C2_res,'H', Hdir2);
		//C2_res->bondAngle1 = normAngle(C2_res->C1->bondAngle1 - 60);
		//updateA(C1_res, 'H');
		//updateA(C2_res, 'H');
	}
	else {
		Cpointer Cnew1 = addC(C1_res, Cdir, 1.4, true);
		updateA(C1_res, 'H', Hdir1);
		Cpointer Cnew2 = addC(Cnew1, FEdir, 1.4, true);
		updateA(C2_res,'H', Hdir2);
		//saveXYZ("KMC_DEBUG/oxidation_issue");
		
		if (hept_bool){
			//Assume that this site was in an heptagon.
			removeR7internal(C1_res, C2_res);
			addC(Cnew2, FEdir, 1.4, true);
		}
    }
    // update sites and neighbours
    Spointer S3, S4;
    S1 = moveIt(stt,-1); S2 = moveIt(stt,1); S3 = moveIt(S1, -1); S4 = moveIt(S2, +1);
	if (!hept_bool){
		bool leftR5 = (isR5internal(C1_res->C1, C1_res));
		bool rightR5 = (isR5internal(C2_res, C2_res->C2));
		int leftsite = (int)S1->type;
		int rightsite = (int)S2->type;
		int zerosite = 2;
		
		if (leftR5){
			//R5 to the left
			if (leftsite > 500 && leftsite < 2000) {
				zerosite += 100;
				updateSites(S1, S1->C1, C1_res, -401);
				if ((int)S3->type < 2000) updateSites(S3, S3->C1, S3->C2, -400);
			}
			else if (leftsite > 2000) {
				zerosite += 500;
				if (leftsite == 2103) {
					if((int)S3->type == 100) updateSites(S1, S1->C1, C1_res, -1501);
					else updateSites(S1, S1->C1, C1_res, -1101);
				}
				else if (leftsite == 2114 || leftsite == 2115) {
					updateSites(S1, S1->C1, C1_res, -11);
				}
				else if (leftsite == 2205) updateSites(S1, S1->C1, C1_res, -1);
				else updateSites(S1, S1->C1, C1_res, -1501);
			}
			else updateSites(S1, S1->C1, C1_res, -1);
		} else{
			updateSites(S1, S1->C1, C1_res, -1);
		}
		if (rightR5){
			//R5 to the right
			if (rightsite > 500 && rightsite < 2000) {
				zerosite += 100;
				updateSites(S2, C2_res, S2->C2, -401);
				if ((int)S4->type < 2000) updateSites(S4, S4->C1, S4->C2, -400);
			}
			else if (rightsite > 2000) {
				zerosite += 500;
				if (rightsite == 2103) {
					if((int)S4->type == 100) updateSites(S2, C2_res, S2->C2, -1501);
					else updateSites(S2, C2_res, S2->C2, -1101);
				}
				else if (rightsite == 2114 || rightsite == 2115) {
					updateSites(S2, C2_res, S2->C2, -11);
				}
				else if (rightsite == 2205) updateSites(S2, C2_res, S2->C2, -1);
				else updateSites(S2, C2_res, S2->C2, -1501);
			}
			else updateSites(S2, C2_res, S2->C2, -1);
		} else{
			updateSites(S2, C2_res, S2->C2, -1);
		}
		
		if ( (isR5internal(C1_res->C2, C2_res->C1)) && !leftR5 && !rightR5 ) {
			zerosite = 2002;
			Cpointer thirdC = findThirdC(C1_res->C2);
			Cpointer thirdC2 = findThirdC(C2_res->C1);
			if (thirdC == NULLC && thirdC2 == NULLC){
				//Both carbons connected to the R5 that got exposed were embedded.
				m_pah->m_rings5_Embedded--;
				m_pah->m_rings5_Lone++;
			}
		}
		
		updateSites(stt, C1_res, C2_res, zerosite);
		
		
		/*if (S1->type == R5R6 && (S2->type == R5R6)) {
			S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
			updateSites(stt, C1_res, C2_res, 202);
			updateSites(S1, S1->C1, C1_res, -1);
			updateSites(S2, C2_res, S2->C2, -1);
			if ((int)S3->type > 500 && (int)S3->type < 2000) updateSites(S3, S3->C1, S3->C2, -400);
			if ((int)S4->type > 500 && (int)S4->type < 2000) updateSites(S4, S4->C1, S4->C2, -400);
			if (!m_pah->m_optimised){
				OpenBabel::OBMol mol = passPAH();
				mol = optimisePAH(mol);
				passbackPAH(mol);
			}
		}	
		else if (S1->type == R5R6 || (S2->type == R5R6)) {
			if (S1->type == R5R6){
				S3 = moveIt(S1, -1);
				if (isR5internal(C2_res, C2_res->C2, true) || isR5internal(C2_res, C2_res->C2, false)) {
					if ((int)S2->type>2000) {
						updateSites(stt, C1_res, C2_res, 602);
					else updateSites(stt, C1_res, C2_res, 602);
				}
				else updateSites(stt, C1_res, C2_res, 102);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
				if ((int)S3->type <2000)  updateSites(S3, S3->C1, S3->C2, -400);
			}
			else {
				S4 = moveIt(S2, 1);
				if (isR5internal(C1_res->C1, C1_res, true) || isR5internal(C1_res->C1, C1_res, false)) updateSites(stt, C1_res, C2_res, 502);
				else updateSites(stt, C1_res, C2_res, 102);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
				if ((int)S4->type <2000) updateSites(S4, S4->C1, S4->C2, -400);
			}
			if (!m_pah->m_optimised){
				OpenBabel::OBMol mol = passPAH();
				mol = optimisePAH(mol);
				passbackPAH(mol);
			}
		}
		else if (S1->type == R5ACR5 || (S2->type == R5ACR5)) {
			if (S1->type == R5ACR5 && S2->type == R5ACR5){
				updateSites(stt, C1_res, C2_res, 1002);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
				S3 = moveIt(S1, -1);
				if (S3->type==R5) updateSites(S1, S1->C1, S1->C2, -400);
				S4 = moveIt(S2, 1);
				if (S4->type==R5) updateSites(S2, S2->C1, S2->C2, -400);
			}
			else if (S1->type == R5ACR5){
				S3 = moveIt(S1, -1);
				if (isR5internal(C2_res, C2_res->C2, true) || isR5internal(C2_res, C2_res->C2, false)) updateSites(stt, C1_res, C2_res, 1002);
				else updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
				if (S3->type==R5) updateSites(S1, S1->C1, S1->C2, -400);
			}
			else {
				S4 = moveIt(S2, 1);
				if (isR5internal(C1_res->C1, C1_res, true) || isR5internal(C1_res->C1, C1_res, false)) updateSites(stt, C1_res, C2_res, 1002);
				else updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
				if (S4->type==R5) updateSites(S2, S2->C1, S2->C2, -400);
			}
			if (!m_pah->m_optimised){
				OpenBabel::OBMol mol = passPAH();
				mol = optimisePAH(mol);
				passbackPAH(mol);
			}
		}
		else if (S1->type == ACR5 && S2->type == ACR5) {
			updateSites(stt, C1_res, C2_res, 1002);
			updateSites(S1, S1->C1, C1_res, -1501);
			updateSites(S2, C2_res, S2->C2, -1501);
		}
		else if (S1->type == ACR5 || S2->type == ACR5) {
			if (S1->type == ACR5){
				if (isR5internal(C2_res, C2_res->C2, true) || isR5internal(C2_res, C2_res->C2, false)) updateSites(stt, C1_res, C2_res, 1002);
				else updateSites(stt, C1_res, C2_res, 502);
			}
			else {
				if (isR5internal(C1_res->C1, C1_res, true) || isR5internal(C1_res->C1, C1_res, false)) updateSites(stt, C1_res, C2_res, 1002);
				else updateSites(stt, C1_res, C2_res, 502);
			}
			updateSites(S1, S1->C1, C1_res, -1);
			updateSites(S2, C2_res, S2->C2, -1);
		}
		else if (((int)S1->type >= 2003 && (int)S1->type <= 2005) && ((int)S2->type >= 2003 && (int)S2->type <= 2005)){
			cpair R5check1 = endposR5internal(C1_res->C1, C1_res);
			cpair R5check2 = endposR5internal(C2_res, C2_res->C2, true);
			bool check_1 = false; bool check_2 = false;
			for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
				double distR51 = getDistance_twoC(*it, R5check1);
				double distR52 = getDistance_twoC(*it, R5check1);
				if (distR51 < 0.5) check_1 = true;
				if (distR52 < 0.5) check_2 = true;
			}
			if ( check_1 && check_2){
				updateSites(stt, C1_res, C2_res, 2);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else if (check_1){
				updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1501);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else if (check_2){
				updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1501);
			}
			else {
				updateSites(stt, C1_res, C2_res, 1003);
				updateSites(S1, S1->C1, C1_res, -1501);
				updateSites(S2, C2_res, S2->C2, -1501);
			}
		}
		else if ((int)S1->type >= 2003 && (int)S1->type <= 2005) {
			cpair R5check1 = endposR5internal(C1_res->C1, C1_res);
			bool check_1 = false; 
			for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
				double distR51 = getDistance_twoC(*it, R5check1);
				if (distR51 < 0.5) check_1 = true;
			}
			if (check_1){
				updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1501);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else{
				updateSites(stt, C1_res, C2_res, 2);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
			}
		}
		else if ((int)S2->type >= 2003 && (int)S2->type <= 2005) {
			cpair R5check1 = endposR5internal(C1_res, C1_res->C1, true);
			bool check_1 = false; 
			for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
				double distR51 = getDistance_twoC(*it, R5check1);
				if (distR51 < 0.5) check_1 = true;
			}
			if (check_1){
				updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1501);
			}
			else{
				updateSites(stt, C1_res, C2_res, 2);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
			}
		}
		/*else if (S1->type == R5R6 || (S2->type == R5R6)) {
			if (S1->type == R5R6 && (S2->type == R5R6)) {
				convSiteType(stt, C1_res, C2_res, FE);
				S3 = moveIt(S1, -1);
				S4 = moveIt(S2, +1);
				convSiteType(S3, S3->C1, S3->C2, FE);
				convSiteType(S4, S4->C1, S4->C2, FE);
				convSiteType(S1, S1->C1, S1->C1->C2, R5);
				convSiteType(S2, S2->C1, S2->C1->C2, R5);
				removeR5internal(S1->C1, S1->C1->C2); redrawR5(S1, S1->C1->C1, S1->C2->C2);
				removeR5internal(S2->C1, S2->C1->C2); redrawR5(S2, S2->C1->C1, S2->C2->C2);
				m_pah->m_rings5_Embedded--; m_pah->m_rings5_Embedded--;
			}
			else if (S1->type == R5R6){
				convSiteType(stt, C1_res, C2_res, ZZ);
				S3 = moveIt(S1, -1);
				convSiteType(S3, S3->C1, S3->C2, FE);
				convSiteType(S1, S1->C1, S1->C1->C2, R5);
				removeR5internal(S1->C1, S1->C1->C2); redrawR5(S1, S1->C1->C1, S1->C2->C2);
				updateSites(S2, C2_res, S2->C2, -1);
				m_pah->m_rings5_Embedded--;
			}
			else {
				convSiteType(stt, C1_res, C2_res, ZZ);
				S4 = moveIt(S2, +1);
				convSiteType(S4, S4->C1, S4->C2, FE);
				convSiteType(S2, S2->C1, S2->C1->C2, R5);
				removeR5internal(S2->C1, S2->C1->C2); redrawR5(S2, S2->C1, S2->C1->C2);
				updateSites(S1, S1->C1, C1_res, -1);
				m_pah->m_rings5_Embedded--;
				/*convSiteType(stt, C1_res, C2_res, ZZ);
				S4 = moveIt(S2, +1);
				convSiteType(S4, S4->C1, S4->C2, FE);
				convSiteType(S2, S2->C1, S2->C1->C2, R5);
				redrawR5(S2, S2->C1->C1, S2->C2->C2);
				updateSites(S1, S1->C1, C1_res, -1);
				m_pah->m_rings5_Embedded--;
			}
		}
		else {
			if (isR5internal(C1_res->C2, C2_res->C1) || isR5internal(C1_res->C2, C2_res->C1, true)) {
				updateSites(stt, C1_res, C2_res, 2002);
				Cpointer thirdC = findThirdC(C1_res->C2);
				Cpointer thirdC2 = findThirdC(C2_res->C1);
				if (thirdC == NULLC && thirdC2 == NULLC){
					//Both carbons connected to the R5 that got exposed were embedded.
					m_pah->m_rings5_Embedded--;
					m_pah->m_rings5_Lone++;
				}
			}
			else updateSites(stt, C1_res, C2_res, 2);
			updateSites(S1, S1->C1, C1_res, -1);
			updateSites(S2, C2_res, S2->C2, -1);
		}*/
		// update combined sites
		S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
		updateCombinedSites(stt); // update resulting site
		updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
		updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
		// add ring counts
		m_pah->m_rings--;
	}
	else {
		//Assuming that the site removed was an heptagon
		bool R5_bool0, R5_bool1, R5_bool2;
		if (isR5internal(C1_res, C1_res->C2)) R5_bool0 = true;
		else R5_bool0 = false;
		if ( isR5internal(C1_res->C2, C1_res->C2->C2) || isR5internal(C2_res->C1->C1, C2_res->C1) ) R5_bool1 = true;
		else R5_bool1 = false;
		if (isR5internal(C2_res->C1, C2_res) ) R5_bool2 = true;
		else R5_bool2 = false;
		
		if (R5_bool0 && R5_bool2){
			if ((int)S1->type >= 2002 && (int)S2->type >= 2002) updateSites(stt, C1_res, C2_res, 1003);
			else if ((int)S1->type >= 2002 || (int)S2->type >= 2002) updateSites(stt, C1_res, C2_res, 603);
			else updateSites(stt, C1_res, C2_res, 203);
		}
		if ( R5_bool0 && R5_bool1) updateSites(stt, C1_res, C2_res, 2103);
		if ( R5_bool2 && R5_bool1) updateSites(stt, C1_res, C2_res, 2103);
		if ( !R5_bool0 && !R5_bool2 && R5_bool1) updateSites(stt, C1_res, C2_res, 2003);
		if (R5_bool0 && ! R5_bool2 && !R5_bool1){
			if ((int)S1->type >= 2002) updateSites(stt, C1_res, C2_res, 503);
			else updateSites(stt, C1_res, C2_res, 103);
		}
		if (!R5_bool0 && R5_bool2 && !R5_bool1){
			if ((int)S2->type >= 2002) updateSites(stt, C1_res, C2_res, 503);
			else updateSites(stt, C1_res, C2_res, 103);
		}
		
		if (S1->type == ACR5 && S2->type == ACR5) {
			//updateSites(stt, C1_res, C2_res, 2204);
			updateSites(S1, S1->C1, C1_res, -1501);
			updateSites(S2, C2_res, S2->C2, -1501);
		}
		else if (S1->type == ACR5 || S2->type == ACR5) {
			//updateSites(stt, C1_res, C2_res, 2103);
			updateSites(S1, S1->C1, C1_res, -1);
			updateSites(S2, C2_res, S2->C2, -1);
		}
		else if (((int)S1->type >= 2003 && (int)S1->type <= 2005) && ((int)S2->type >= 2003 && (int)S2->type <= 2005)){
			double dist1 = getDistance_twoC(C1_res, C1_res->C1);
			double dist2 = getDistance_twoC(C2_res, C2_res->C2);
			if (dist1 > 1.1 && dist2 > 1.1){
				//updateSites(stt, C1_res, C2_res, 2002);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else if (dist1 > 1.1){
				//updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1501);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else if (dist2 > 1.1){
				//updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1501);
			}
			else {
				//updateSites(stt, C1_res, C2_res, 1003);
				updateSites(S1, S1->C1, C1_res, -1501);
				updateSites(S2, C2_res, S2->C2, -1501);
			}
		}
		else if ((int)S1->type >= 2003 && (int)S1->type <= 2005) {
			double dist2 = getDistance_twoC(C1_res, C1_res->C1);
			if (dist2 > 1.1){
				//updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1501);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else{
				//updateSites(stt, C1_res, C2_res, 2);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
			}
		}
		else if ((int)S2->type >= 2003 && (int)S2->type <= 2005) {
			double dist2 = getDistance_twoC(C2_res, C2_res->C2);
			if (dist2 > 1.1){
				//updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1501);
			}
			else{
				//updateSites(stt, C1_res, C2_res, 2);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
			}
		}
		else if (S1->type == R5R6 || (S2->type == R5R6)) {
			if (S1->type == R5R6 && (S2->type == R5R6)) {
				convSiteType(stt, C1_res, C2_res, FE);
				S3 = moveIt(S1, -1);
				S4 = moveIt(S2, +1);
				convSiteType(S3, S3->C1, S3->C2, FE);
				convSiteType(S4, S4->C1, S4->C2, FE);
				convSiteType(S1, S1->C1, S1->C1->C2, R5);
				convSiteType(S2, S2->C1, S2->C1->C2, R5);
				removeR5internal(S1->C1, S1->C1->C2); redrawR5(S1, S1->C1->C1, S1->C2->C2);
				removeR5internal(S2->C1, S2->C1->C2); redrawR5(S2, S2->C1->C1, S2->C2->C2);
				m_pah->m_rings5_Embedded--; m_pah->m_rings5_Embedded--;
			}
			else if (S1->type == R5R6){
				convSiteType(stt, C1_res, C2_res, ZZ);
				S3 = moveIt(S1, -1);
				convSiteType(S3, S3->C1, S3->C2, FE);
				convSiteType(S1, S1->C1, S1->C1->C2, R5);
				removeR5internal(S1->C1, S1->C1->C2); redrawR5(S1, S1->C1->C1, S1->C2->C2);
				updateSites(S2, C2_res, S2->C2, -1);
				m_pah->m_rings5_Embedded--;
			}
			else {
				convSiteType(stt, C1_res, C2_res, ZZ);
				S4 = moveIt(S2, +1);
				convSiteType(S4, S4->C1, S4->C2, FE);
				convSiteType(S2, stt->C2, stt->C2->C2, R5);
				removeR5internal(S2->C1, S2->C1->C2); redrawR5(S2, S2->C1->C1, S2->C2->C2);
				updateSites(S1, S1->C1, C1_res, -1);
				m_pah->m_rings5_Embedded--;
				/*convSiteType(stt, C1_res, C2_res, ZZ);
				S4 = moveIt(S2, +1);
				convSiteType(S4, S4->C1, S4->C2, FE);
				convSiteType(S2, S2->C1, S2->C1->C2, R5);
				redrawR5(S2, S2->C1->C1, S2->C2->C2);
				updateSites(S1, S1->C1, C1_res, -1);
				m_pah->m_rings5_Embedded--;*/
			}
		}
		else {
			//updateSites(stt, C1_res, C2_res, 2003);
			updateSites(S1, S1->C1, C1_res, -1);
			updateSites(S2, C2_res, S2->C2, -1);
		}
		// update combined sites
		S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
		updateCombinedSites(stt); // update resulting site
		updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
		updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
		// add ring counts
		m_pah->m_rings7_Embedded--;
	}
	//if (bridge) saveXYZ("KMC_DEBUG/Oxidation_to_bridge_after");
}
// ************************************************************
// ID9- R6 oxidation at AC by OH (AR13 in Matlab)
// ************************************************************
void PAHProcess::proc_O6R_FE_HACA_OH(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    proc_O6R_FE_HACA(stt, C_1, C_2);
}
// ************************************************************
// ID10- R5 growth on ZZ (AR3 in Matlab)
// ************************************************************
void PAHProcess::proc_G5R_ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSites(stt);
	if (getDistance_twoC(C_1,C_2) < 2.1 && !m_pah->m_optimised) {
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
    if(checkHindrance_newCposition(C_1) || checkHindrance_newCposition(C_2)) {
       /* cout<<"Site hindered, process not performed.\n";*/ return;}
    
    // member C atoms of resulting R5 site
	cpair Cdir = get_vector(C_1->coords,C_2->coords);
	cpair normvec = (norm_vector(C_1->coords, C_1->C2->coords, C_2->coords));
	cpair crossvec = cross_vector(Cdir, normvec);
	cpair Hdir1 = cross_vector( get_vector(C_1->C2->coords, C_2->coords),normvec );
	cpair Hdir2 = cross_vector( get_vector(C_1->coords, C_1->C2->coords),normvec );
	double dist = getDistance_twoC(C_1,C_2);
	double dist2 = dist / 2.35 * 1.47;
	double theta = acos( (dist - dist2)/2.0 / dist2 );
	cpair resultantvec = std::make_tuple((dist - dist2)/2.0 * std::get<0>(Cdir) + 1.5*sin(theta) * std::get<0>(crossvec), (dist - dist2)/2.0 * std::get<1>(Cdir) + 1.5*sin(theta) * std::get<1>(crossvec), (dist - dist2)/2.0 * std::get<2>(Cdir)+ 1.5*sin(theta) * std::get<2>(crossvec));
	cpair intdir = scale_vector(resultantvec);
    removeC(C_2->C1, true);
    Cpointer C1_res, C2_res;
	C1_res = addC(C_1, intdir, dist/2.35*1.5);
	updateA(C1_res, 'H', Hdir1);
	updateA(C_1, 'C', Hdir1);
	C2_res = addC(C1_res, Cdir, dist2); //1.4*sqt(3.0)
	updateA(C2_res, 'H', Hdir2);
	updateA(C_2, 'C', Hdir2);
    /*C1_res = addC(C_1, normAngle(C_1->bondAngle1+120), 0, 1.4);
	C2_res = addC(C1_res, normAngle(C_1->bondAngle1 - 90), normAngle(C_1->bondAngle1 - 180), 1.4*pow(3, 0.5));
    updateA(C_1->C1, C_2->C2, 'H');*/
    // update sites and neighbours
    convSiteType(stt, C1_res, C2_res, R5);
	addR5internal(C1_res,C2_res);
    Spointer S1, S2, S3, S4;
    // neighbours
    S1 = moveIt(stt,-1); 
    S2 = moveIt(stt,1);
    addR5toSite(S1, S1->C1, C1_res);
    addR5toSite(S2, C2_res, S2->C2);
    // update combined sites
    S3 = moveIt(S1, -1); 
    S4 = moveIt(S2, 1);
    updateCombinedSites(stt); // update resulting site
    updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
    // add ring counts
    m_pah->m_rings5_Lone++;
	if ( (dist > 2.7 && !m_pah->m_optimised) || (getDistance_twoC(C2_res, C2_res->C2->C2)<2.4&& !m_pah->m_optimised) ) {
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
}
// ************************************************************
// ID11- R5 desorption (AR7 in Matlab)
// ************************************************************
void PAHProcess::proc_D5R_R5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSites(stt);
    // member C atoms of resulting R5 site
	removeR5internal(C_1,C_2);
    Cpointer C1_res, C2_res;
    C1_res = C_1->C1;
    C2_res = C_2->C2;
	cpair Cdir = add_vector( invert_vector(C_1->growth_vector),C_2->growth_vector);
	double Cdist = getDistance_twoC(C1_res, C2_res);
	cpair Cnewpos = jumpToPos(C1_res->coords,Cdir,Cdist/2.35*1.4);
	Cnewpos = checkHindrance_C_intPAH(Cnewpos);
	double newdist = getDistance_twoC(C1_res->coords,Cnewpos);
	if (newdist >= 1.7){
		m_pah->m_InternalCarbons.push_back(Cnewpos);
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
		Cdir = add_vector( invert_vector(C_1->growth_vector),C_2->growth_vector);
	} else m_pah->m_InternalCarbons.push_back(Cnewpos);
    removeC(C_1, false);
    removeC(C_2, false);
	Cdist = getDistance_twoC(C1_res, C2_res);
	Cpointer newC = addC(C1_res, Cdir, Cdist/2.35*1.4, true);
	updateA(newC, 'C', newC->coords);
	cpair Hdir = add_vector( get_vector(newC->coords, C1_res->coords),get_vector(newC->coords, C2_res->coords) );
	updateA(C1_res, 'H', Hdir);
	updateA(C2_res, 'H', Hdir);
    // update sites and neighbours
    Spointer S1, S2, S3, S4;
    // neighbours
    S1 = moveIt(stt,-1); S2 = moveIt(stt,1);
    convSiteType(stt, C1_res, C2_res, ZZ);
    remR5fromSite(S1, S1->C1, C1_res);
    remR5fromSite(S2, C2_res, S2->C2);
    // update combined sites
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt); // update resulting site
    updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
    // add ring counts
    m_pah->m_rings5_Lone--;
}
// ************************************************************
// ID12- R6 conversion to R5 (AR9 in Matlab)
// ************************************************************
void PAHProcess::proc_C6R_AC_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
    //printSites(stt);
	//printStruct();
	if (C_1->bridge || C_1->C2->bridge || C_2->C1->bridge || C_2->bridge){
		// Process is not possible  since AC is a bridge
		return;}
    if(checkHindrance_newCposition(C_1) || checkHindrance_newCposition(C_2)) {
        /*cout<<"Site hindered, process not performed.\n"*/ return;}
    // check if FE3 is before or after AC
    bool b4 = false;
    if(moveIt(stt,-2)->comb == FE3 && moveIt(stt,2)->comb == FE3) {
        // Define a distribution that has two equally probably outcomes
        boost::bernoulli_distribution<> choiceDistrib;
        // Now build an object that will generate a sample using rng
        boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);

        if(choiceGenerator())
            b4 = true; // if FE3 on both sides, choose a random one
    }
    else {
        if(moveIt(stt,-2)->comb == FE3)
            b4 = true;
    }
	if (b4 && (int)moveIt(stt, -4)->type == 100) return;
	if (b4 && (int)moveIt(stt, -4)->type >= 500 && (int)moveIt(stt, -4)->type <= 504) return; // Two pentagons will collide if this is the case
	if (b4 && (int)moveIt(stt, -4)->type >= 2000 && (int)moveIt(stt, -4)->type <= 2205) return; // Two pentagons will collide if this is the case
	if (!b4 && (int)moveIt(stt, +4)->type == 100) return;
	if (!b4 && (int)moveIt(stt, +4)->type >= 500 && (int)moveIt(stt, +4)->type <= 504) return; // Two pentagons will collide if this is the case
	if (!b4 && (int)moveIt(stt, +4)->type >= 2000 && (int)moveIt(stt, +4)->type <= 2205) return; // Two pentagons will collide if this is the case
    Cpointer C1_res, C2_res, C1_R5, C2_R5, C_xR5;
    Spointer FE_res;
	cpair Hdir1, Hdir2, normvec, perpdir, crossvec, intdir, resultantvec, Cdir;
	double dist, dist2, theta;
    // removing R6 first. save the two member C for resulting FE
    if(b4) {
        FE_res = moveIt(stt, -2); // resulting FE (the FE3)
        C1_res = FE_res->C1->C1->C1;
        C2_res = C_1->C2;
		C_xR5 = C_2;
		perpdir = get_vector(C1_res->coords, C1_res->C2->coords);
		Cdir = get_vector(C_1->C2->coords,C_2->coords);
		normvec = (norm_vector(C_1->C2->coords, C_2->C1->coords, C_2->coords));
		crossvec = cross_vector(Cdir, normvec);
		Hdir1 = cross_vector( get_vector(C_2->C1->coords, C_2->coords),normvec );
		Hdir2 = cross_vector( get_vector(C_1->C2->coords, C_1->C2->C2->coords),normvec );
		dist = getDistance_twoC(C_1->C2,C_2);
		dist2 = dist / 2.35 * 1.47;
		theta = acos( (dist - dist2)/2.0 / 1.5 );
		resultantvec = std::make_tuple((dist - dist2)/2.0 * std::get<0>(Cdir) + dist/2.35*1.5*sin(theta) * std::get<0>(crossvec), (dist - dist2)/2.0 * std::get<1>(Cdir) + dist/2.35*1.5*sin(theta) * std::get<1>(crossvec), (dist - dist2)/2.0 * std::get<2>(Cdir)+ dist/2.35*1.5*sin(theta) * std::get<2>(crossvec));
		intdir = scale_vector(resultantvec);
    }else {
        FE_res = moveIt(stt, 2);
        C1_res = C_2->C1;
        C2_res = FE_res->C2->C2->C2;
		C_xR5 = C_1;
		perpdir = get_vector(C2_res->coords, C2_res->C1->coords);
		Cdir = get_vector(C_1->coords,C_2->C1->coords);
		normvec = (norm_vector(C_1->coords, C_1->C2->coords, C_2->C1->coords));
		crossvec = cross_vector(Cdir, normvec);
		Hdir1 = cross_vector( get_vector(C_1->C2->coords, C_2->C1->coords),normvec );
		Hdir2 = cross_vector( get_vector(C_1->coords, C_1->C2->coords),normvec );
		dist = getDistance_twoC(C_1->C2,C_2);
		dist2 = dist / 2.35 * 1.47;
		theta = acos( (dist - dist2)/2.0 / 1.5 );
		cpair resultantvec = std::make_tuple((dist - dist2)/2.0 * std::get<0>(Cdir) + 1.5*sin(theta) * std::get<0>(crossvec), (dist - dist2)/2.0 * std::get<1>(Cdir) + 1.5*sin(theta) * std::get<1>(crossvec), (dist - dist2)/2.0 * std::get<2>(Cdir)+ 1.5*sin(theta) * std::get<2>(crossvec));
		intdir = scale_vector(resultantvec);
    }
    // remove C atoms
    for(int i=0; i!=4; i++) removeC(C1_res->C2, false);
    //C1_res->bondAngle1 = normAngle(C1_res->bondAngle1-120);
    // now add R5 on resulting ZZ site (used to be AC)
    Cpointer Cstart;
    if(b4) Cstart = C2_res;
    else Cstart = C_1;
    removeC(Cstart->C2, true);
	C1_R5 = addC(Cstart, intdir, dist / 2.35*1.5 );
	updateA(Cstart, 'C', Hdir1);
	updateA(C1_R5, 'H', Hdir1);
	C2_R5 = addC(C1_R5, Cdir, dist2);
	updateA(C2_R5, 'H', Hdir2);
	updateA(C_xR5, 'C', Hdir2);
	if (b4) updateA(C1_res, 'H', perpdir);
	else updateA(C2_res, 'H', perpdir);
	// update H atoms
    //if(b4) updateA(C1_res->C1, C_2->C2, 'H');
    //else updateA(C_1->C1, C2_res->C2, 'H');
    // updating the sites, first remove the two neighbouring FE in the FE3
    Spointer FE1, FE2;
    FE1 = moveIt(FE_res, -1);
    FE2 = moveIt(FE_res, 1);
    // remove the sites
    removeSite(FE1);
    removeSite(FE2);
    // convert the AC site, FE site, and update the neighbours
    Spointer S1, S2, S3, S4;
    if(b4){
        S1 = moveIt(FE_res, -1);
        S2 = moveIt(stt, 1);
        convSiteType(stt, C1_R5, C2_R5, R5); // convert the AC site to R5
        addR5toSite(FE_res, C1_res, C1_R5); // convert the FE site to RFE
        addR5toSite(S2, C2_R5, S2->C2); // update neighbour of resulting R5
        updateSites(S1, S1->C1, C1_res, -1); // update neighbour of resulting RFE
    }else {
        S1 = moveIt(stt, -1);
        S2 = moveIt(FE_res, 1);
        convSiteType(stt, C1_R5, C2_R5, R5); // convert the AC site to R5
        addR5toSite(FE_res, C2_R5, C2_res); // convert the FE site to RFE
        addR5toSite(S1, S1->C1, C1_R5); // update neighbour of resulting R5
        updateSites(S2, C2_res, S2->C2, -1); // update neighbour of resulting RFE
    }
    // update combined sites for all the sites involved and their neighbours
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt); updateCombinedSites(FE_res);// update resulting sites
    updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
    // update H count
    addCount(0, -2);
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
    // add ring counts
    m_pah->m_rings--;
    m_pah->m_rings5_Lone++;
	//addR5internal(C1_R5, C2_R5); // Not needed since the optimisation is being called.
}
// ************************************************************
// ID13- R5 conversion to R6 on FE (AR5 in Matlab)
// ************************************************************
void PAHProcess::proc_C5R_RFE(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSites(stt);
    if(checkHindrance_newCposition(C_1) || checkHindrance_newCposition(C_2)) {
        /*cout<<"Site hindered, process not performed.\n"*/ return;}
    
    // check if R5 is before or after the site
    bool b4=false;
    if(moveIt(stt,-1)->type == R5) b4 = true;
    // R5 will be removed and replaced by ZZ. first identify the C member of the
    // resulting AC site, C_AC which is from bulk. Identify where the R5 site is too.
    Spointer sR5;
    Cpointer C_AC, Cstart;
	cpair starting_direction, Cdir, ZZdir, Hdir;
    if(b4) {
        sR5 = moveIt(stt, -1);
        C_AC = sR5->C1->C1;
		starting_direction = invert_vector( get_vector(C_2->coords, C_2->C2->coords) );
		Cdir = get_vector(C_1->C2->coords, C_2->coords);
		ZZdir = Cdir;
		Hdir = starting_direction;
		Cstart = C_AC;
    }else {
        sR5 = moveIt(stt, 1);
        C_AC = sR5->C2->C2;
		starting_direction = add_vector( get_vector(C_1->C1->coords, C_1->coords), get_vector(C_1->C2->coords,C_1->coords));
		Cdir = get_vector(C_1->coords, C_1->C2->coords);
		ZZdir = invert_vector(starting_direction);
		Hdir = get_vector(C_1->C1->coords, C_1->coords);
		Cstart = C_2->C1;
    }
    // remove R5 first, leaving a ZZ site (adding another C atom after removing C)
	removeR5internal(sR5->C1, sR5->C2);
	for(int i=0; i!=2; i++) removeC(Cstart->C2, false);
    addC(Cstart, ZZdir, 1.4, true);
	updateA(C_AC, 'H', Hdir);
	
	//addC(Cstart, normAngle(Cstart->bondAngle1-120), normAngle(Cstart->bondAngle1-60), 1.4, true);
    // this new C atom is irrelevant. Next add a R6 on the resulting FE (from RFE)
    Cpointer C1_new, C2_new, C3_new, C4_new; // save all new C atoms
	double dist;
    if(b4) {
        Cstart = C_2->C1;
		cpair opp_vec = get_vector(Cstart->C1->coords, Cstart->coords);
		dist = getDistance_twoC(Cstart, Cstart->C2);
		C1_new = addC(Cstart, starting_direction, dist);
		updateA(C1_new, 'H', invert_vector(Cdir));
		C2_new = addC(C1_new, opp_vec, dist);
		updateA(C2_new, 'H', starting_direction);
		C3_new = addC(C2_new, Cdir, dist);
		updateA(C3_new, 'H', opp_vec);
		C4_new = addC(C3_new, invert_vector(starting_direction), dist);
		updateA(C4_new, 'H', Cdir);
		updateA(C_2, 'C', Cdir);
	}
    else {
        Cstart = C_1;
		dist = getDistance_twoC(Cstart, Cstart->C2);
		C1_new = addC(Cstart, starting_direction, dist);
		updateA(C1_new, 'H', invert_vector(Cdir));
		C2_new = addC(C1_new, Hdir, dist);
		updateA(C2_new, 'H', starting_direction);
		C3_new = addC(C2_new, Cdir, dist);
		updateA(C3_new, 'H', Hdir);
		C4_new = addC(C3_new, invert_vector(starting_direction), dist);
		updateA(C4_new, 'H', Cdir);
		updateA(C_1, 'C', Cdir);
	}
    /*C1_new = addC(Cstart, normAngle(Cstart->bondAngle1+120), 0, 1.4);
    C2_new = addC(C1_new, normAngle(C1_new->C1->bondAngle1-60), 0, 1.4);
    C3_new = addC(C2_new, normAngle(C2_new->C1->bondAngle1-60), 0, 1.4);
    C4_new = addC(C3_new, normAngle(C3_new->C1->bondAngle1-60), normAngle(C3_new->C1->bondAngle1-120), 1.4);*/
    // edit sites. first identify the neighbouring sites of resulting AC & FE3
    Spointer S1, S2, S3, S4;
    if(b4) {
        S1 = moveIt(sR5, -1); // neighbour of R5
        S2 = moveIt(stt, 1); // neighbour of RFE (stt)
        convSiteType(sR5, C_AC, C1_new, AC); // convert R5 to AC
        remR5fromSite(S1, S1->C1, C_AC); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        remR5fromSite(stt, C2_new, C3_new); // convert RFE to FE
        addSite(FE, C1_new, C2_new, stt); // add new FE site before the original RFE
        addSite(FE, C3_new, C4_new, S2); // add new FE site after the original RFE
        updateSites(S2, C4_new, S2->C2, +1); // update resulting FE3 neighbour
    }else {
        S1 = moveIt(stt, -1); // neighbour of RFE (stt)
        S2 = moveIt(sR5, 1); // neighbour of R5
        convSiteType(sR5, C4_new, C_AC, AC); // convert R5 to AC
        remR5fromSite(S2, C_AC, S2->C2); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        remR5fromSite(stt, C2_new, C3_new); // convert RFE to FE
        addSite(FE, C1_new, C2_new, stt); // add new FE site before the original RFE
        addSite(FE, C3_new, C4_new, sR5); // add new FE site after the original RFE
        updateSites(S1, S1->C1, C1_new, +1); // update resulting FE3 neighbour
    }
    // update H atoms
    //if(b4) updateA(C_AC->C1, C_2->C2, 'H');
    //else updateA(C_1, C_AC->C2, 'H');
    // update combined sites for all sites involved and their neighbours
    // (excluding new FE sites, since their combined site type will still be None)
    S3 = moveIt(S1, -1);
    S4 = moveIt(S2, 1);
    updateCombinedSites(stt); updateCombinedSites(sR5); // new FE3 and AC
    updateCombinedSites(S1); updateCombinedSites(S2); 
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours
    // update H count
    addCount(0, 2);
    // add ring counts
    m_pah->m_rings++;
    m_pah->m_rings5_Lone--;
}
// ************************************************************
// ID14- R5 conversion to R6 on AC (AR4 in Matlab)
// ************************************************************
void PAHProcess::proc_C5R_RAC(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSites(stt);
	//printStruct();
	if(checkHindrance_newCposition(C_1) || checkHindrance_newCposition(C_2)) {
        /*cout<<"Site hindered, process not performed.\n"*/ return;}
	bool b4;
	//check if R5 in RAC is reactive or R5 is already partially embedded.
	if (moveIt(stt, -1)->type == R5 || moveIt(stt, +1)->type == R5){
		// check if R5 is before or after RAC
		b4 = false;
	}
	else return; //R5 is shouldered and cannot move.
	// check if R5 is before or after RAC
	b4 = false;
    if(moveIt(stt,-1)->type == R5) b4 = true;
    // R5 will be removed and replaced by AC. first identify the C member of the
    // resulting AC site, C_AC. Identify where the R5 site is too.
    Spointer sR5;
    Cpointer C_AC;
	cpair starting_direction, Hdir1, Hdir2, FEdir, Hdir, CZZdir, R5normvec, R6normvec;
    if(b4) {
        sR5 = moveIt(stt, -1);
		R5normvec = norm_vector(sR5->C1->coords, sR5->C2->coords, sR5->C2->C2->coords);
		R6normvec = norm_vector(C_2->C1->coords, C_2->coords, C_2->C2->coords);
        C_AC = sR5->C1->C1;
		cpair FEvector = get_vector(C_1->C2->C2->coords,C_2->C1->coords);
		cpair AGvector = get_vector(C_1->C2->C2->coords,C_1->C2->coords);
		starting_direction = add_vector(FEvector, AGvector);
		cpair AJvector = get_vector(C_2->C1->coords,C_2->coords);
		FEdir = add_vector(AJvector, invert_vector(AGvector));
		//starting_direction = get_vector(C_2->C1->coords, C_2->coords);
		Hdir1 = get_vector(C_1->C2->C2->coords, C_1->C2->coords);
		//FEdir = get_vector(C_1->C2->coords, C_2->coords);
		Hdir2 = starting_direction;
		Hdir = starting_direction;
		CZZdir = get_vector(C_1->C2->coords, C_1->C2->C2->coords);
    }else {
        sR5 = moveIt(stt, 1);
		R5normvec = norm_vector(sR5->C1->coords, sR5->C2->coords, sR5->C2->C2->coords);
		R6normvec = norm_vector(C_1->C1->coords, C_1->coords, C_1->C2->coords);
        C_AC = sR5->C2->C2;
		cpair FEvector = get_vector(C_1->C2->coords,C_2->C1->C1->coords);
		cpair AGvector = get_vector(C_1->C2->coords,C_1->coords);
		starting_direction = add_vector(FEvector, AGvector);
		cpair AJvector = get_vector(C_2->C1->C1->coords,C_2->C2->coords);
		FEdir = add_vector(AJvector, invert_vector(AGvector));
		//starting_direction = get_vector(C_2->C1->C1->coords, C_2->C1->coords);
		Hdir1 = get_vector(C_1->C2->coords, C_1->coords);
		//FEdir = get_vector(C_1->coords, C_2->C1->coords);
		Hdir2 = starting_direction;
		Hdir = Hdir1;
		CZZdir = FEdir;
    }
	double angle_norm_vectors = dot_vector(R5normvec, R6normvec);
	//saveXYZ("KMC_DEBUG/RAC");
    // check if there's a bridge in the BY5
    bool bridge = false;
    if(b4) bridge = C_2->C1->bridge;
    else bridge = C_1->C2->bridge;
    // remove R5 first, leaving a ZZ site (adding another C atom after removing C)
    Cpointer Cstart;
    if(b4) Cstart = C_AC; 
    else Cstart = C_2->C1;
	removeR5internal(sR5->C1, sR5->C2);
    for(int i=0; i!=2; i++) removeC(Cstart->C2, false);
	addC(Cstart, CZZdir, 1.4, true);
	updateA(C_AC, 'H', Hdir);
    // this new C atom is irrelevant. Next add a R6 on the resulting AC (from RAC)
    Cpointer C1_new, C2_new; // save all new C atoms
	if(b4) Cstart = C_AC->C2->C2;
	else Cstart = C_1;
	if(!bridge){
        for(int i=0; i!=2; i++) removeC(Cstart->C2, true);
    } else {//else if there's a bridge, convert the bridge atoms to normal edge atoms
        Cpointer b1 = Cstart->C2;
        Cpointer b2 = b1->C3;
        connectToC(Cstart, b2->C2);
        b1->C3 = NULL;
        b2->C3 = NULL;
        b1->bridge = false;
        b2->bridge = false;
        b1->C1 = b2;
        b2->C2 = b1;
        //angletype a = b2->bondAngle1;
        //b2->bondAngle1 = b2->bondAngle2;
        //b2->bondAngle2 = a;
    }
    C1_new = addC(Cstart, starting_direction, 1.4);
	updateA(C1_new, 'H', Hdir1);
	updateA(Cstart, 'C', Hdir1);
	C2_new = addC(C1_new, FEdir, 1.4);
	updateA(C2_new, 'H', Hdir2);
	updateA(C_2, 'C', FEdir);
    //C1_new = addC(Cstart, normAngle(Cstart->bondAngle1+120), 0, 1.4);
    //C2_new = addC(C1_new, normAngle(C1_new->C1->bondAngle1-60), normAngle(C1_new->C1->bondAngle1-120), 1.4);
    // edit sites. first identify the neighbouring sites of resulting AC & FE3
    Spointer S1, S2, S3, S4;
    if(b4) {
        S1 = moveIt(sR5, -1); // neighbour of R5
        S2 = moveIt(stt, 1); // neighbour of RAC (stt)
        convSiteType(sR5, C_AC, C1_new, AC); // convert R5 to AC
        remR5fromSite(S1, S1->C1, C_AC); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, C1_new, C2_new, FE); // convert RAC to FE
        updateSites(S2, C2_new, S2->C2, +1); // update resulting FE neighbour
    } else {
        S1 = moveIt(stt, -1); // neighbour of RAC (stt)
        S2 = moveIt(sR5, 1); // neighbour of R5
        convSiteType(sR5, C2_new, C_AC, AC); // convert R5 to AC
        remR5fromSite(S2, C_AC, S2->C2); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, C1_new, C2_new, FE); // convert RAC to FE
        updateSites(S1, S1->C1, C1_new, +1); // update resulting FE neighbour
    }
    // update H atoms
    //if(b4) updateA(C_AC->C1, C_2->C2, 'H');
    //else updateA(C_1->C1, C_AC->C2, 'H');
    // update combined sites for all sites involved and their neighbours
    // (excluding new FE sites, since their combined site type will still be None)
	//saveXYZ("KMC_DEBUG/RAC_beforeoptim");
	if (angle_norm_vectors < 0.30){
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
		//saveXYZ("KMC_DEBUG/RAC_afteroptim");
	}
    S3 = moveIt(S1, -1);
    S4 = moveIt(S2, 1);
    updateCombinedSites(stt); updateCombinedSites(sR5); // new FE and AC
    updateCombinedSites(S1); updateCombinedSites(S2); 
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours
    // add ring counts
    m_pah->m_rings++;
    m_pah->m_rings5_Lone--;
}
// ************************************************************
// ID15- R5 migration to neighbouring ZZ (AR6 in Matlab)
// ************************************************************
void PAHProcess::proc_M5R_RZZ(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSites(stt);
    if(checkHindrance_newCposition(C_1) || checkHindrance_newCposition(C_2)) {
        /*cout<<"Site hindered, process not performed.\n"*/ return;}
    
    // check if R5 is before or after RZZ
    bool b4 = false;
    if(moveIt(stt,-1)->type == R5) b4 = true;
    // R5 will be removed and replaced by ZZ. first identify the C member of the
    // resulting RZZ site, C_RZZ. Identify where the R5 site is too.
    Spointer sR5;
	Cpointer C_RZZ, Cstart, Ccheck;
	cpair ZZCdir, Hdir1, Hdir2, Hdir, R5vec;
	double ZZdist;
    if(b4) {
        sR5 = moveIt(stt, -1);
		Hdir1 = add_vector(get_vector(sR5->C1->C1->coords, sR5->C1->coords), get_vector(sR5->C1->C2->coords, sR5->C1->coords) ); //sR5->C1->growth_vector;
		Hdir2 = add_vector(get_vector(sR5->C2->C1->coords, sR5->C2->coords), get_vector(sR5->C2->C2->coords, sR5->C2->coords) ); //sR5->C2->growth_vector;
        C_RZZ = sR5->C1->C1;
		ZZCdir = get_vector(C_1->C2->coords, C_1->C2->C2->coords );
		Hdir = add_vector(get_vector(sR5->C1->C1->C1->coords, sR5->C1->C1->coords), get_vector(sR5->C2->C2->C2->coords, sR5->C2->C2->coords) ); //Hdir2;
		ZZdist = getDistance_twoC(C_1->C2, C_2);
		R5vec = get_vector(C_1->C2->coords, C_2->coords);
    }else {
        sR5 = moveIt(stt, 1);
		Hdir1 = add_vector(get_vector(sR5->C1->C1->coords, sR5->C1->coords), get_vector(sR5->C1->C2->coords, sR5->C1->coords) ); //sR5->C1->growth_vector;
		Hdir2 = add_vector(get_vector(sR5->C2->C1->coords, sR5->C2->coords), get_vector(sR5->C2->C2->coords, sR5->C2->coords) ); //sR5->C2->growth_vector;
        C_RZZ = sR5->C2->C2;
		ZZCdir = get_vector(C_1->coords, C_1->C2->coords );
		Hdir = add_vector(get_vector(sR5->C1->C1->C1->coords, sR5->C1->C1->coords), get_vector(sR5->C2->C2->C2->coords, sR5->C2->C2->coords) ); //Hdir1;
		ZZdist = getDistance_twoC(C_1, C_2->C1);
		R5vec = get_vector(C_1->coords, C_2->C1->coords);
    }
	cpair starting_direction = get_vector(sR5->C1->C1->coords, sR5->C1->coords);
	//Check for ZZ site
	Ccheck = sR5->C1->C1;
	cpair mpos = jumpToPos(Ccheck->coords, ZZCdir, 1.4);
	//cpair mpos = jumpToPos(Ccheck->coords, normAngle(Ccheck->C1->bondAngle1 - 60), 0, 1.4);
    //Remove the R5 from internal coordinates list
	findR5internal(Ccheck->C2, Ccheck->C2->C2);
	// remove R5 first, leaving a ZZ site (adding another C atom after removing C)
    for(int i=0; i!=2; i++) removeC(Ccheck->C2, false);
	Cpointer C1_R5, C2_R5; // save all new C atoms
	if (!checkHindrance_C_PAH(mpos)) {
		//Create ZZ carbon
		addC(Ccheck, ZZCdir, 1.4, true);
		updateA(C_RZZ, 'H', Hdir);
		//addC(Cstart, normAngle(Cstart->bondAngle1 - 120), normAngle(Cstart->bondAngle1 - 60), 1.4, true);
		// Next add a R5 on the neighbouring ZZ
		if (b4) Cstart = C_RZZ->C2->C2;
		else Cstart = C_1;
		//Remove neighbouring ZZ carbon
		removeC(Cstart->C2, true);
		//Add new R5
		C1_R5 = addC(Cstart, starting_direction, ZZdist/2.35*1.5);
		updateA(C1_R5, 'H', Hdir1);
		updateA(Cstart, 'C', Hdir1);
		C2_R5 = addC(C1_R5, R5vec, ZZdist/2.35*1.47);
		updateA(C2_R5, 'H', Hdir2);
		updateA(C2_R5->C2, 'C', Hdir1);
		
		if (getDistance_twoC(C2_R5, C2_R5->C2) > 1.8){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
	}
	else { // migration over bridge.
		if (b4){
			//Adjust angle of previous carbon
			//Cstart->bondAngle1 = normAngle(Cstart->C1->bondAngle1 - 60);
			//Connect carbon to existing ZZ carbon
			saveXYZ("KMC_DEBUG/Migration_over_bridge");
			Cpointer C_newbridge, C_sharedbridge, C_oldbridge;
			if (b4) Cstart = C_RZZ->C2->C2;
			else Cstart = C_1;
			Cstart->C2 = findC(mpos);
			C_newbridge = Cstart->C2;
			C_sharedbridge = Cstart->C2->C1;
			C_oldbridge = Cstart->C2->C1->C3;
			//Adjust bridges
			C_newbridge->bridge = true;
			C_newbridge->C3 = C_sharedbridge;
			C_newbridge->C1 = Cstart;
			//C_newbridge->bondAngle2 = normAngle(Cstart->bondAngle1 + 60);
			C_sharedbridge->C1 = C_oldbridge;
			C_sharedbridge->C2 = C_2;
			C_sharedbridge->C3 = C_newbridge;
			//C_sharedbridge->bondAngle1 = C_newbridge->bondAngle2 + 60;
			//C_sharedbridge->bondAngle2 = C_newbridge->bondAngle2 - 180;
			C_oldbridge->bridge = false;
			C_oldbridge->C2 = C_sharedbridge;
			C_oldbridge->C3 = NULL;
			//C_oldbridge->bondAngle1 = C_oldbridge->C1->bondAngle1 + 60;
			//C_oldbridge->bondAngle2 = 0;
			updateA(Cstart, 'H', Hdir);
			// Next add a R5 on the neighbouring ZZ
			Cstart = C_sharedbridge;
			//Add new R5
			C1_R5 = addC(Cstart, starting_direction, ZZdist/2.35*1.5);
			updateA(C1_R5, 'H', Hdir1);
			updateA(C_2, 'C', Hdir1);
			C2_R5 = addC(C1_R5, R5vec, ZZdist/2.35*1.47 );
			updateA(C2_R5, 'H', Hdir2);
			//C2_R5->C2 = C_2;
			//Update 'H'
			//updateA(C1_R5, 'H'); updateA(C2_R5, 'H'); updateA(C2_R5->C2, 'C');
		}
		else{
			//Adjust angle of previous carbon
			//Cstart->bondAngle1 = normAngle(Cstart->C1->bondAngle1 - 60);
			//Connect carbon to existing ZZ carbon
			saveXYZ("KMC_DEBUG/Migration_over_bridge");
			Cpointer C_newbridge, C_sharedbridge, C_oldbridge;
			Ccheck->C1 = findC(mpos);
			C_newbridge = Ccheck->C1;
			C_sharedbridge = C_newbridge->C2;
			C_oldbridge = C_sharedbridge->C3;
			//Adjust bridges
			C_newbridge->bridge = true;
			C_newbridge->C3 = C_sharedbridge;
			C_newbridge->C2 = Ccheck;
			//C_newbridge->bondAngle1 = normAngle(C_newbridge->C1->bondAngle1 - 60);
			//C_newbridge->bondAngle2 = normAngle(C_newbridge->C1->bondAngle1 + 60);
			C_sharedbridge->C2 = C_oldbridge;
			C_sharedbridge->C1 = C_1;
			C_1->C2 = C_sharedbridge;
			C_sharedbridge->C3 = C_newbridge;
			//C_sharedbridge->bondAngle1 = C_newbridge->bondAngle2 + 60;
			//C_sharedbridge->bondAngle2 = C_newbridge->bondAngle2 - 180;
			C_oldbridge->bridge = false;
			C_oldbridge->C1 = C_sharedbridge;
			C_oldbridge->C3 = NULL;
			//C_oldbridge->bondAngle2 = 0;
			updateA(Ccheck, 'H', Hdir);
			// Next add a R5 on the neighbouring ZZ
			Cstart = C_1;
			//Add new R5
			//Cstart->bondAngle1 = normAngle(Cstart->C1->bondAngle1 + 60);
			C1_R5 = addC(Cstart, starting_direction, ZZdist/2.35*1.5);
			updateA(C1_R5, 'H', Hdir1);
			updateA(C_1, 'C', Hdir1);
			C2_R5 = addC(C1_R5, R5vec, ZZdist/2.35*1.47 );
			updateA(C2_R5, 'H', Hdir2);
			//C1_R5 = addC(Cstart, normAngle(Cstart->bondAngle1), 0, 1.4);
			//C2_R5 = addC(C1_R5, normAngle(C1_R5->C1->bondAngle1 - 90), normAngle(C1_R5->C1->bondAngle1 - 180), 1.4*pow(3, 0.5));
			//C2_R5->C2 = C_2;
			//Update 'H'
			//updateA(C1_R5, 'H'); updateA(C2_R5, 'H'); updateA(C2_R5->C2, 'C');
		}
		if (getDistance_twoC(C2_R5, C2_R5->C2) > 1.8 || getDistance_twoC(C2_R5, C2_R5->C2->C2) < 2.4){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
	}
	//Add R5 to internal coordinates after the migration.
    if (!m_pah->m_optimised) addR5internal(C1_R5, C2_R5);
    // edit sites. first identify the neighbouring sites of resulting RZZ & R5
    Spointer S1, S2, S3, S4, S5, S6;
    if(b4) {
        S1 = moveIt(sR5, -1); // neighbour of R5
        S2 = moveIt(stt, 1); // neighbour of RAC (stt)
        convSiteType(sR5, C_RZZ, C1_R5, RZZ); // convert R5 to RZZ
        remR5fromSite(S1, S1->C1, C_RZZ); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, C1_R5, C2_R5, R5); // convert RZZ to R5
        addR5toSite(S2, C2_R5, S2->C2); // update S2
    } else {
        S1 = moveIt(stt, -1); // neighbour of RAC (stt)
        S2 = moveIt(sR5, 1); // neighbour of R5
        convSiteType(sR5, C2_R5, C_RZZ, RZZ); // convert R5 to RZZ
        addR5toSite(S1, S1->C1, C1_R5); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, C1_R5, C2_R5, R5); // convert RZZ to R5
        remR5fromSite(S2, C_RZZ, S2->C2); // update S2
    }
    // update combined sites for all sites involved and their neighbours
    // (excluding new FE sites, since their combined site type will still be None)
    S3 = moveIt(S1, -1);
    S4 = moveIt(S2, 1);
	S5 = moveIt(S1, -2);
    S6 = moveIt(S2, 2);
    updateCombinedSites(stt); updateCombinedSites(sR5); // new FE and AC
    updateCombinedSites(S1); updateCombinedSites(S2); 
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours
	updateCombinedSites(S5); updateCombinedSites(S6); // neighbours
    //printSites(stt);
    //cout<<sp.None;
}

// ************************************************************
// ID16- R6 migration & conversion to R5 at BY5 (pyrene+R5; pathway 1; AR22 in Matlab)
// ************************************************************
void PAHProcess::proc_C6R_BY5_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
    //printSites(stt);
    // check if there are any bridges in the BY5, cancel process is yes
    Cpointer now=C_1->C2;
    for(int i=0; i!=3; i++) {
        if(now->bridge) return;
        else now = now->C2;
    }
    // check if FE3 is before or after BY5
    bool b4 = false;
    if(moveIt(stt,-2)->comb == FE3 && moveIt(stt,2)->comb == FE3) {
        // Define a distribution that has two equally probably outcomes
        boost::bernoulli_distribution<> choiceDistrib;
        // Now build an object that will generate a sample using rng
        boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);

        b4= choiceGenerator(); // if FE3 on both sides, choose a random one
    }
    else {
        if(moveIt(stt,-2)->comb == FE3)
            b4 = true;
    }
	//check is site neighbouring FE3 has an embedded R5: Two pentagons together are not allowed.
	Spointer neighbour;
	if (b4) {
		neighbour = moveIt(stt, -4);
		if ((int)neighbour->type >= 2002 && (int)neighbour->type <= 2005) return;
	}
	else {
		neighbour = moveIt(stt, 4);
		if ((int)neighbour->type >= 2002 && (int)neighbour->type <= 2005) return;
	}

    // first the R6 will be changed to R5 by removing a C. locate where the FE3 is and which
    // C is to be removed.
    Spointer sFE3;
    Cpointer CRem;// C to be removed from original FE3 site
	cpair Cdir = get_vector(C_1->C2->coords, C_1->coords);
    if(b4) {
        sFE3 = moveIt(stt, -2); 
        CRem = sFE3->C1;
    }else {
        sFE3 = moveIt(stt, 2);
        CRem = sFE3->C2;
    }
	// close BY5 to form R6. First remove all bulk C
	for (int i = 0; i != 3; i++) removeC(C_1->C2, true);
	// add a C atom
	Cpointer Cnew = addC(C_1, C_1->growth_vector, 1.4);
	updateA(C_1, 'C', C_1->growth_vector);
	updateA(Cnew, 'H', Cdir);
	updateA(C_2, 'C', C_2->growth_vector);
	//addC(C_1, normAngle(C_1->bondAngle1 + 120), normAngle(C_1->bondAngle1 + 60), 1.4);
	Cpointer CR5_1 = CRem->C1;
	removeC(CRem, false);
	//CR5_1->bondAngle1 = CR5_1->bondAngle1 - 30;

    // delete neighbouring (to sFE3) FE sites from site map. Identify them first:
    Spointer Srem1 = moveIt(sFE3,-1);
    Spointer Srem2 = moveIt(sFE3,1);
    // then remove them
    removeSite(Srem1);
    removeSite(Srem2);
    // edit sites. first identify the neighbouring sites of resulting RFE & R5
    Spointer S1, S2, S3, S4;
    if(b4) {
        S1 = moveIt(sFE3, -1); // neighbour of FE3
        S2 = moveIt(stt, 1); // neighbour of BY5 (stt)
		convSiteType(sFE3, CR5_1, CR5_1->C2, R5); // convert FE3 to R5
		convSiteType(stt, CR5_1->C2, C_1->C2, RFE); // convert BY5 to RFE
        updateSites(S2, C_1->C2, S2->C2, +1); // add a bulk C to S2
		if ((int)S1->type <= 2000) updateSites(S1, S1->C1, CR5_1, -1);
		addR5internal(CR5_1, CR5_1->C2);
        addR5toSite(S1, S1->C1, CR5_1); // remove a bulk C from S1 and add R5
    } else {
        S2 = moveIt(sFE3, 1); // neighbour of FE3
        S1 = moveIt(stt, -1); // neighbour of BY5 (stt)
		convSiteType(sFE3, CR5_1, CR5_1->C2, R5); // convert FE3 to R5
		convSiteType(stt, C_1->C2, CR5_1, RFE); // convert BY5 to RFE
        updateSites(S1, S1->C1, C_1->C2, +1); // add a bulk C to S1
		if ((int)S2->type <= 2000) updateSites(S2, CR5_1->C2, S2->C2, -1);
		addR5internal(CR5_1, CR5_1->C2);
        addR5toSite(S2, CR5_1->C2, S2->C2); // remove a bulk C from S2 and add R5
    }
    
    // update H atoms
    //updateA(S1->C2->C1, S2->C1->C2, 'H');
    // update combined sites for all sites involved and their neighbours
    // (excluding new FE sites, since their combined site type will still be None)
    S3 = moveIt(S1, -1);
    S4 = moveIt(S2, 1);
    updateCombinedSites(stt); updateCombinedSites(sFE3); // new R5 and RFE
    updateCombinedSites(S1); updateCombinedSites(S2); 
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours

    addCount(0, -2);
    // add ring counts
    m_pah->m_rings5_Lone++;
    //printSites(stt);
    // cout<<sp.None;
    //Currently the pentagons formed by this jump process are too large. Use optimiser to fix it. Not ideal but quick solution.
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
}
// ************************************************************
// ID17- R6 migration & conversion to R5 at BY5 (pyrene+R5; pathway 2-violi; AR24 in Matlab)
// ************************************************************
void PAHProcess::proc_C6R_BY5_FE3violi(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
    proc_C6R_BY5_FE3(stt, C_1, C_2, rng);
}
// ************************************************************
// ID18- BY5 closure (AR16 in Matlab)
// ************************************************************
void PAHProcess::proc_L5R_BY5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//printSites(stt);
	// Remove C
	Cpointer now = C_1->C2;
	Spointer other_side;
	bool bridged_before = false;
	int ostype;
	//Detect bridges
	if (C_1->C2->bridge || C_1->C2->C2->bridge || C_2->C1->bridge){
		bridged_before = true;
		Cpointer bridgecheck = now;
		do {
			if (!bridgecheck->bridge) {
				bridgecheck = bridgecheck->C2;
			}
			else {
				other_side = findSite(bridgecheck->C3->C1);
				ostype = (int)other_side->type;
				if (ostype % 10 == 2 && ostype!= 2) {
					//this violates the isolated pentagon rule.
					return;
				}
				if (ostype >= 2002) {
					//This would create two five member rings next to each other.
					return;
				}
				break;
			}
		}while (bridgecheck != C_2);
	}
	
	do{
		Cpointer next;
		if (!now->bridge) {
			next = now->C2;
			//if(now->C3 != NULL) now->C3->C3 = NULL;
			removeC(now, true);
			now = next;
		}
		else {
			next = now->C3->C2; //
			// bridged bulk C will not be removed from edge
			Cpointer b = now->C3; // C bridged to now
			b->C3 = NULL; now->C3 = NULL;
			b->bridge = false; now->bridge = false;
			//b->bondAngle1 = b->bondAngle2;
			//now->bondAngle2 = 0; b->bondAngle2 = 0;
			// connect C_1 to next and the two bulk atoms still remaining
			connectToC(C_1, next);
			connectToC(b, now);
			now = next;
		}
	} while (now != C_2);
	// Change bond angle between C_1 and C_2
	//Connect C_1 and C_2
	connectToC(C_1, C_2);
	// Add and remove H
	updateA(C_1, 'C', C_1->growth_vector);
	updateA(C_2, 'C', C_2->growth_vector);
	//updateA(C_1->C1, C_2->C2, 'H');
	// Convert the BY6 site into the resulting site after reaction,
	// finding resulting site type:
	//int ntype_site = (int)stt->type;
	int ntype1 = (int)moveIt(stt, -1)->type;
	int ntype2 = (int)moveIt(stt, 1)->type;
	int newType;
	if (ntype1 < 5 && ntype2 < 5) {
		newType = (2002 + ntype1 + ntype2);
		// convert site
		if (newType>2005) {
			//saveDOT(std::string("BY5ClosureProblem.dot"));
			//saveXYZ("KMC_DEBUG/BY5ClosureProblem");
			//std::cerr << "ERROR: newType is > 65 (PAHProcess::proc_L5R_BY5)\n";
		}
		//Spointer Srem1 = moveIt(stt, -1);
		//Spointer Srem2 = moveIt(stt, 1);
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}
	else if (ntype1 == 501 && ntype2 == 501) convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)2204);
	else if (ntype1 == 501 || ntype2 == 501){
		if (ntype1 == 501) ntype1 = 0;
		if (ntype2 == 501) ntype2 = 0;
		newType = 2103 + ntype1 + ntype2;
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}
	else if ( (ntype1 >= 501 && ntype1 <= 504) || (ntype2 >= 501 && ntype2 <= 504)){
		newType = 2002;
		if ( (ntype1 >= 501 && ntype1 <= 504) ) {
			ntype1 = ntype1 % 500;
			newType = newType + 100;
		}
		if ( (ntype2 >= 501 && ntype2 <= 504) ){
			ntype2 = ntype2 % 500;
			newType = newType + 100;
		}
		newType = newType + ntype1 + ntype2;
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}
	else if(ntype1 > 2000 && ntype2 > 2000){
		newType = 9999;
		ntype1 = 0;
		ntype2 = 0;
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}
	else if ( (ntype1 >= 2002 && ntype1 <= 2015) || (ntype2 >= 2002 && ntype2 <= 2015) ){
		if (ntype1 % 2000 >= 10) ntype1 -= 10;
		if (ntype2 % 2000 >= 10) ntype2 -= 10;
		ntype1 = ntype1 % 2002;
		ntype2 = ntype2 % 2002;
		newType = 2104 + ntype1 + ntype2;
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}
	else if ( (ntype1 >= 2103 && ntype1 <= 2115) || (ntype2 >= 2103 && ntype2 <= 2115) ){
		if (ntype1 % 2100 >= 10) ntype1 -= 10;
		if (ntype2 % 2100 >= 10) ntype2 -= 10;
		ntype1 = ntype1 % 2103;
		ntype2 = ntype2 % 2103;
		newType = 2205 + ntype1 + ntype2;
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}
	else if (ntype1 == 9999 || ntype2 == 9999){
		newType = 9999;
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}
	else {
		newType = (ntype1 + ntype2 + 2002);
		if ((kmcSiteType)newType == None) {
			//saveDOT(std::string("BY5ClosureProblem.dot"));
			std::cerr << "ERROR: newType is None (PAHProcess::proc_L5R_BY5)\n";
		}
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}

	//Optimise if needed
	if (getDistance_twoC(C_1, C_2) > 2.6){
		if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		else addR5internal(C_1,C_2, true);
	}
	else addR5internal(C_1,C_2, true);
	
	if (bridged_before){
		if (ostype >= 100) {
			ostype = ostype % 10;
			ostype += 2100;
		}
		else ostype += 2000;
		if (other_side != m_pah->m_siteList.end()){
			convSiteType(other_side, other_side->C1, other_side->C2, (kmcSiteType)ostype);
			updateCombinedSites(other_side);
		}
	}

	
	
	// erase the existence of the neighbouring sites
	Spointer Srem1 = moveIt(stt, -1);
	Spointer Srem2 = moveIt(stt, 1);
	removeSite(Srem1);
	removeSite(Srem2);
	// update combined sites and neighbours
	Spointer S1 = moveIt(stt, -1); Spointer S2 = moveIt(stt, 1);
	//Spointer S3 = moveIt(S1,-1); Spointer S4 = moveIt(S2,1);
	updateCombinedSites(stt);
	updateCombinedSites(S1); updateCombinedSites(S2); 
	//updateCombinedSites(S3); updateCombinedSites(S4);

	//printSites(stt);
	// update H count
	addCount(0, -2);
	// add ring counts
	m_pah->m_rings5_Lone++;
}
// ************************************************************
// ID19- R6 desorption at bay -> pyrene (AR21 in Matlab)
// ************************************************************
void PAHProcess::proc_M6R_BY5_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
    //printSites(stt);
    // check if there are any bridges in the BY5, cancel process if yes
    Cpointer now=C_1->C2;
    for(int i=0; i!=3; i++) {
        if(now->bridge) return;
        else now = now->C2;
    }
    // check if FE3 is before or after BY5
    bool b4 = false;
    if(moveIt(stt,-2)->comb == FE3 && moveIt(stt,2)->comb == FE3) {
        // Define a distribution that has two equally probably outcomes
        boost::bernoulli_distribution<> choiceDistrib;
        // Now build an object that will generate a sample using rng
        boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);

        b4 = choiceGenerator(); // if FE3 on both sides, choose a random one
    }
    else if(moveIt(stt,-2)->comb == FE3) b4 = true;
    // BY5 will be closed to form R6. firstly locate where the FE3 is
    Spointer sFE3;
	cpair Hdir;
    if(b4) {
        sFE3 = moveIt(stt, -2);
		Hdir = sFE3->C1->growth_vector;
    }else {
        sFE3 = moveIt(stt, 2);
		Hdir = sFE3->C2->growth_vector;
    }
	cpair Cdir = get_vector(C_1->C2->coords, C_1->coords);
    // then remove all bulk C in the BY5
    for(int i=0; i!=3; i++) removeC(C_1->C2, true);
    // then add a C to convert it to R6
	Cpointer Cnew = addC(C_1, C_1->growth_vector, 1.4);
	updateA(C_1, 'C', C_1->growth_vector);
	updateA(Cnew, 'H', Cdir);
	updateA(C_2, 'C', C_1->growth_vector);
    //addC(C_1, normAngle(C_1->bondAngle1+120), normAngle(C_1->bondAngle1+60), 1.4);
    // next remove R6 neighbouring with the original R5. removing bulk C:
    Cpointer Cstart, Cend;
	cpair CZZdir;
    if(b4) {
		Cstart = C_1->C1->C1->C1->C1;
		Cend = C_1;
		CZZdir = get_vector(C_1->coords,C_1->C2->coords);
		}
    else {
		Cstart = C_2;
		Cend = C_2->C2->C2->C2->C2;
		CZZdir = invert_vector(Cdir);
	}
    for(int i=0; i!=3; i++) removeC(Cstart->C2, false);
    // then add a C
	addC(Cstart, CZZdir, 1.4, true);
	updateA(Cstart, 'H', Hdir);
	updateA(Cend, 'H', Hdir);
    // delete neighbouring (to sFE3) FE sites from site map. Identify them first:
    Spointer Srem1 = moveIt(sFE3,-1);
    Spointer Srem2 = moveIt(sFE3,1);
    // then remove them
    removeSite(Srem1);
    removeSite(Srem2);
    // edit sites. first identify the neighbouring sites of original BY5 and R6
    Spointer S1, S2, S3, S4;
    if(b4) {
        S1 = moveIt(sFE3, -1); // neighbour of FE3
        S2 = moveIt(stt, 1); // neighbour of BY5 (stt)
        updateSites(sFE3, Cstart, C_1, +1); // convert FE3 (FE) to ZZ
        convSiteType(stt, C_1, C_2->C1, FE); // convert BY5 to FE
        updateSites(S1, S1->C1, Cstart, -1); // remove a bulk C from S1
        updateSites(S2, C_2->C1, S2->C2, +1); // add a bulk C to S2
		//if( (int)S1->type >= 501 && (int) S1->type <= 504) updateSites(sFE3, Cstart, C_1, +500); // convert FE3 to R5R6
    } else {
        S1 = moveIt(stt, -1); // neighbour of BY5 (stt)
        S2 = moveIt(sFE3, 1); // neighbour of FE3
        updateSites(sFE3, Cstart, C_2->C2->C2, +1); // convert FE3 (FE) to ZZ
        convSiteType(stt, C_1->C2, C_2, FE); // convert BY5 to FE
        updateSites(S1, S1->C1, C_1->C2, +1); // add a bulk C to S1
        updateSites(S2, sFE3->C2, S2->C2, -1); // remove a bulk C from S2
		//if((int)S2->type >= 501 && (int)S2->type <= 504) updateSites(sFE3, Cstart, C_2->C2->C2, +500); // convert FE3 to R5R6
    }
    
    // update H atoms
    //updateA(S1->C2->C1, S2->C1->C2, 'H');
    // update combined sites for all sites involved and their neighbours
    // (excluding new FE sites, since their combined site type will still be None)
    S3 = moveIt(S1, -1);
    S4 = moveIt(S2, 1);
    updateCombinedSites(stt); updateCombinedSites(sFE3); // new FE and ZZ
    updateCombinedSites(S1); updateCombinedSites(S2); 
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours

    addCount(0, -2);
    //printSites(stt);
   // cout<<sp.None;
   double distrr = getDistance_twoC(Cnew, Cnew->C2);
   if (distrr < 1.15 || distrr > 1.8){
	   if (!m_pah->m_optimised){
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
	   }
   }
}
// ************************************************************
// ID20 & ID21- R6 oxidation at FE2 site
// ************************************************************
void PAHProcess::proc_O6R_FE2(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//This jump process used to remove a full ring. However, the work by Singh 2016 showed that it probably forms an R5 which is more interesting.
	//Converts an FE2 into an R5.
    //printSites(stt);
    Spointer other = moveIt(stt, -1);
    Spointer S1, S2;
	Cpointer CRem, CRem_before, CRem_next;
    if(other->type == FE) {
        S1 = moveIt(other, -1);
        S2 = moveIt(stt, 1);
		CRem = C_1;
    }
    else {
        other = moveIt(stt, 1);
        S2 = moveIt(other, 1);
        S1 = moveIt(stt, -1);
		CRem = C_2;
    }
	CRem_before = CRem->C1;
	CRem_next = CRem->C2;
	if (S1->type==RFE || S2->type==RFE || S1->type==R5R6 || S2->type==R5R6) return;
	
	//This section assumes that the oxidation of an FE2 site on an heptagon forms a six-member ring. Needs confirmation.
	bool heptagon_flag = false;
	//check that two pentagons (including internals) will not collide
	if (m_pah->m_R7loc.size()>=1){
		cpair R7coords_end = endposR7internal(CRem_before, CRem_next);
		for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
			double distR7s = getDistance_twoC(*it, R7coords_end);
			if (distR7s < 1.7) {
				//FE2 oxidation on an R7
				heptagon_flag = true;
				break;
			}
		}
	}
	
	//check that two pentagons (including internals) will not collide
	if (m_pah->m_R5loc.size()>=1 && !(heptagon_flag)){
		cpair R5coords_end = endposR5internal(CRem_before, CRem);
		for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
			double distR5s = getDistance_twoC(*it, R5coords_end);
			if (distR5s < 3.0) {
				//This distance is a parameter of this jump process. Might need some more tuning. 
				//2.8 seems appropiate but may reject too many jumps.
				//Two pentagons will be next to each other violating the Isolated Pentagon Rule
				return;
			}
		}
	}
	
    Cpointer C_bulk;
    for(C_bulk = S1->C1; C_bulk != S1->C2; C_bulk=C_bulk->C2) {
        if(C_bulk->bridge) return;
    }
    for(C_bulk = S2->C1; C_bulk != S2->C2; C_bulk=C_bulk->C2) {
        if(C_bulk->bridge) return;
    }
	//Remove carbon
	removeC(CRem, false);
    // Remove one of the FE in FE2
    removeSite(other);
    // Update Sites and neighbouring sites
	if (!heptagon_flag){
		convSiteType(stt, CRem_before, CRem_next, R5);
		updateSites(S1, S1->C1, CRem_before, +100);
		updateSites(S2, CRem_next, S2->C2, +100);
		m_pah->m_rings--;
		m_pah->m_rings5_Lone++;
		addR5internal(CRem_before,CRem_next);
	}
	else {
		convSiteType(stt, CRem_before, CRem_next, FE);
		m_pah->m_rings7_Embedded--;
		m_pah->m_rings++;
		removeR7internal(CRem_before, CRem_next);
	}
    Spointer S3, S4;
    // update combined sites for all sites and their neighbours
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt);
    updateCombinedSites(S1); updateCombinedSites(S2);
    updateCombinedSites(S3); updateCombinedSites(S4);
    addCount(0,-1);
}

// ************************************************************
// ID22- R6 desorption from FE to form AC
// ************************************************************
void PAHProcess::proc_D6R_FE_AC(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_O6R_FE_HACA_double(stt, C_1, C_2);
}

//
// ************************************************************
// ID22- Bay-capping
// ************************************************************
void PAHProcess::proc_B6R_ACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//printSitesMemb(stt);
	//printStruct();//++++
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol,2000);
		passbackPAH(mol);
	}
	//saveXYZ("After_1stmin");
	Cpointer newC1;
	Cpointer newC2;
	cpair CR5_1, CR5_2;
	cpair vec1, vec2, Hdir1, Hdir2;
	
	// check if ACR5 has an opposite site.
	Spointer opp_site, opp_site_second;
	bool opp_site_bool = false;
	Cpointer thirdC = findThirdC(C_1->C2);
	Cpointer thirdC2 = findThirdC(C_1->C2->C2);
	if (thirdC != NULLC || thirdC2 != NULLC) {
		opp_site_bool = true;
	}
	if(checkHindrance_newCposition(C_1) || checkHindrance_newCposition(C_2)) {
        /*cout<<"Site hindered, process not performed.\n"*/ return;}
	if (!(C_1->C2->bridge) || !(C_2->C1->bridge)) { // check if bulk C in AC site is a bridge
		//Getting vectors
		cpair FEvector = get_vector(C_1->C2->coords,C_2->C1->coords);
		cpair AGvector = get_vector(C_1->C2->coords,C_1->coords);
		cpair start_direction = add_vector(FEvector, AGvector);
		cpair AJvector = get_vector(C_2->C1->coords,C_2->coords);
		cpair lat_direction = add_vector(AJvector, invert_vector(AGvector));
		double dist = getDistance_twoC(C_1,C_2);
		moveC_z(C_1->C2, 0.2);
		moveC_z(C_2->C1, 0.2);
        removeC(C_1->C2, true);
        removeC(C_2->C1, true);
		
		newC1 = addC(C_1, start_direction, dist/2.0);
		updateA(newC1, 'H', AGvector);
		updateA(C_1, 'C', AGvector);
		newC2 = addC(newC1, lat_direction, dist/2.0);
		//newC2 = addC(newC1, FEvector, 1.4);
		updateA(newC2, 'H', AJvector);
		updateA(C_2, 'C', AJvector);
		//newC2->bondAngle2 = normAngle(newC2->bondAngle1+120);
		
		/*// Add and remove C
		removeC(C_1->C2, true);
		newC1 = addC(C_1, normAngle(C_1->bondAngle1 + 90), 0, 1.4);
		removeC(newC1->C2, true);
		newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1 - 60), normAngle(newC1->C1->bondAngle1 -120), 3.45);
		newC2->bondAngle2 = normAngle(newC2->bondAngle1+120);
		CR5_3 = newC1->coords;
		// Add and remove H
		updateA(C_1, 'C');
		updateA(C_1->C2, 'H', true);
		updateA(C_1->C2->C2, 'H', true);
		CR5_4 = newC2->coords;*/
	}
	else {
		// update bridges info, both are no longer bridges
		C_1->C2->C1 = C_1->C2->C3;// update neighbour
		C_2->C1->C2 = C_2->C1->C3;
		C_1->C2->bridge = false;
		C_2->C1->bridge = false;
		//angletype a = C_2->C1->bondAngle1;
		//C_2->C1->bondAngle1 = C_2->C1->bondAngle2;
		//C_2->C1->bondAngle2 = a;
		//C_1->C2->C3 = NULL;
		//C_2->C1->C3 = NULL;
		// connect C_1 and C_2
		connectToC(C_1, C_2);
		// Add C
		// Add and remove C
		newC1 = addC(C_1, vec1, 1.4);
		updateA(newC1, 'H', Hdir1);
		updateA(C_1, 'C', Hdir1);
		newC2 = addC(newC1, vec2, 1.4);
		updateA(newC2, 'H', Hdir2);
		updateA(C_2, 'C', Hdir2);
		//newC2->bondAngle2 = normAngle(newC2->bondAngle1+120);
		//This makes sure that embedded pentagons does not account this process as embedding
		m_pah->m_rings5_Lone++;
		m_pah->m_rings5_Embedded--;
	}
	CR5_1 = newC1->coords;
	CR5_2 = newC2->coords;
	//printStruct();
	// neighbouring sites:
	Spointer S1 = moveIt(stt, -1);
	Spointer S2 = moveIt(stt, 1);
	// Update Site and neighbours
	convSiteType(stt, newC1, newC2, (kmcSiteType)0);
	updateSites(S1, S1->C1, newC1, 1); // neighbours
	updateSites(S2, newC2, S2->C2, 1);
	// Update combined site for Site and neighbours
	Spointer S3, S4;
	S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
	updateCombinedSites(stt);
	updateCombinedSites(S1); updateCombinedSites(S2);
	updateCombinedSites(S3); updateCombinedSites(S4);
	// add ring counts
	if (opp_site_bool) {
		m_pah->m_rings5_Lone++;
		m_pah->m_rings5_Embedded--;
	}
	m_pah->m_rings++;
	m_pah->m_rings5_Lone--;
	m_pah->m_rings5_Embedded++;
	//printSites(stt);
	//saveXYZ("Before_second_min");
	if (!m_pah->m_optimised){
		OpenBabel::OBMol newmol = passPAH();
		newmol = optimisePAH(newmol, 4000);
		passbackPAH(newmol);
	}
	//optimisePAH(true, true, "curved_guy");
	//includeCurvature(CR5_1, CR5_2, CR5_3, CR5_4, true, "curved_guy");
}

// ************************************************************
// ID23- Embedded 5-member ring migration to ZZ
// ************************************************************
void PAHProcess::proc_M5R_ACR5_ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	bool b4 = false;
	//printStruct();
	//printSites();
	//Spointer S1 = moveIt(stt, -1);
	if (moveIt(stt, -1)->comb == FE2 || moveIt(stt, 1)->comb == FE2) {
		if (moveIt(stt, -1)->comb == FE2 && moveIt(stt, 1)->comb == FE2) {
			// Define a distribution that has two equally probably outcomes
			boost::bernoulli_distribution<> choiceDistrib;
			// Now build an object that will generate a sample using rng
			boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
			b4 = choiceGenerator(); // if FE3 on both sides, choose a random one
		}
		if (moveIt(stt, -1)->comb == FE2){
			// Define a distribution that has two thirds to one third probably outcome
			boost::bernoulli_distribution<> choiceDistrib(1.0/3.0);
			// Now build an object that will generate a sample using rng
			boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
			b4 = choiceGenerator(); 
		}
	}
	else{
		// Define a distribution that has two equally probably outcomes
		boost::bernoulli_distribution<> choiceDistrib;
		// Now build an object that will generate a sample using rng
		boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
		b4 = choiceGenerator(); // if FE3 on both sides, choose a random one
	}

	Spointer sFE2, checkR5_1, checkR5_2;
	Cpointer CRem, CFE, CRem_before, CRem_next, CR5_otherside_1, CR5_otherside_2;
	
	if (b4) {
		sFE2 = moveIt(stt, -1);
		checkR5_1 = moveIt(stt, -2);
		checkR5_2 = moveIt(stt, -3);
		CFE = C_2->C1->C1;
		CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
		CR5_otherside_1 = C_2->C1;
		CR5_otherside_2 = C_2->C1->C1;
	}
	else {
		sFE2 = moveIt(stt, 1);
		checkR5_1 = moveIt(stt, 2);
		checkR5_2 = moveIt(stt, 3);
		CFE = C_1->C2;
		CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
		CR5_otherside_1 = C_1->C2;
		CR5_otherside_2 = C_1->C2->C2;
	}
	//Check for unsupported sites. This section heavily assumes that the Isolated Pentagon Rule is valid.
	if ((int)sFE2->type == 0 && (int)checkR5_1->type == 0 && (int)checkR5_2->type == 0) return; // The result would be an indene, not supported. YET!
	if ((int)sFE2->type == 101 || (int)sFE2->type == 501 || (int)sFE2->type == 2002) return; // This would violate the IPR.
	if ((int)sFE2->type == 0){
		if ((int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 || (int)checkR5_1->type == 100) return;
		//if ((int)checkR5_1->type >= 501 && (int)checkR5_1->type <= 65) return;
		if ((int)checkR5_1->type >= 1002 && (int)checkR5_1->type <= 1004) return;
		if ((int)checkR5_1->type >= 2002 && (int)checkR5_1->type <= 2204) return;
		if ((int)checkR5_1->type >= 2204 && (int)checkR5_1->type <= 2205) return;
		if ((int)checkR5_1->type == 0){
			//if ((int)checkR5_2->type >= 501 && (int)checkR5_2->type <= 65) return;
			if ((int)checkR5_2->type == 101 || (int)checkR5_2->type == 501) return;
			if ((int)checkR5_2->type >= 1002 && (int)checkR5_2->type <= 1004) return;
			if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2205) return;
		}
	}
	if (CRem_next->bridge) return;
	if (CRem_next->C2->bridge) return; if (CRem_next->C1->bridge) return;
	if (CRem_next->C2->C2->bridge) return; if (CRem_next->C1->C1->bridge) return;
	
	//There are two main cases. The R5 of the ACR5 is shown on a short or a long side.
	double bond_distance = 1.4;
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	bool optimised = false;
	if (R5_dist < 1.6){
		//The ACR5 site is on the "short" side of an R5. 
		optimised = true;
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol, 250);
		passbackPAH(mol);
	}
	
	// check if ACR5 has an opposite site.
	Spointer opp_site, opp_site_second, opp_site_after;
	bool opp_site_bool = false; bool opp_site_bool_second = false; bool opp_site_bool_after = false;
	Cpointer thirdC = findThirdC(CR5_otherside_1);
	Cpointer thirdC2 = findThirdC(CR5_otherside_2);
	Cpointer thirdC_after = findThirdC(CRem_next);
	// Seven cases:
	// 1. One pentagon has one exposed edge  and migrates to a location where it will have one exposed edge. Normal migration.
	// 2. One pentagon has one exposed edge  and migrates to a location where it will have two exposed edges. 
	// 3. One pentagon has two exposed edges and migrates to a location where it will have two exposed edges. 
	// 4. One pentagon has two exposed edges and migrates to a location where it will have one exposed edge.
	// 5. One pentagon has two exposed edges and migrates to a location where it will have three exposed edges.
	// 6. One pentagon has three exposed edges and migrates to a location where it will have three exposed edges.
	// 7. One pentagon has three exposed edges and migrates to a location where it will have two exposed edges.
	
	if (thirdC != NULLC) {
		opp_site = findSite(thirdC);
		if (opp_site != m_pah->m_siteList.end()) opp_site_bool = true;
	}
	if (thirdC2 != NULLC) {
		opp_site_second = findSite(thirdC2);
		if (opp_site_second != m_pah->m_siteList.end() && opp_site_second!=opp_site) opp_site_bool_second = true;
	}
	if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
		opp_site_after = findSite(thirdC_after);
		if (opp_site_after != m_pah->m_siteList.end()) {
			opp_site_bool_after = true;
			int os_endtype = opp_site_after->type;
			if (os_endtype >= 200 && os_endtype <= 203) return;
			if (os_endtype == 101) return;
			if (os_endtype >= 600 && os_endtype <= 603) return;
			if (os_endtype >= 1000 && os_endtype <= 1003) return;
			if (os_endtype >= 500 && os_endtype <= 504) return;
			if (os_endtype >= 2000 && os_endtype <= 2205) return;
			if (os_endtype >= 2103 && os_endtype <= 2105) return;
			if (os_endtype >= 2204 && os_endtype <= 2205) return;
		}
	}
	
	//Fundamental assumption: R5-R7 pairs cannot move away from each other!
	cpair R5coords = findR5internal(C_1->C2, C_2->C1);
	if (m_pah->m_R7loc.size()>=1){
		for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
			double distR5R7 = getDistance_twoC(*it, R5coords);
			if (distR5R7 < 3.1) {
				m_pah->m_R5loc.push_back(R5coords);
				return;
			}
		}
	}
	//check that two pentagons (including internals) will not collide
	cpair R5coords_end = endposR5internal(CRem->C1, CRem);
	if (m_pah->m_R5loc.size()>=1){
		for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
			double distR5s = getDistance_twoC(*it, R5coords_end);
			if (distR5s < 2.8) {
				//This distance is a parameter of this jump process. Might need some more tuning. 
				//2.8 seems appropiate but may reject too many jumps.
				//Two pentagons will be next to each other violating the Isolated Pentagon Rule
				m_pah->m_R5loc.push_back(R5coords);
				return;
			}
		}
	}
	R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2;
	if (R5_dist < 2.5) dist2 = 1.4;
	else dist2 = R5_dist / 2.7 * 1.5;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir = get_vector(C_1->C2->coords,C_2->C1->coords);
	cpair normvec = (norm_vector(C_1->C2->coords, C_1->C2->C2->coords, C_1->C2->C2->C2->coords));
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, bond_distance);
	updateA(Cnew, 'H', crossvec);
	removeC(CRem, false);
	if (optimised){
		//The ACR5 site is on the "short" side of an R5. 
		OpenBabel::OBMol newmol = passPAH();
		newmol = optimisePAH(newmol, 250);
		passbackPAH(newmol);
	}
	
	if ((int)checkR5_1->type == 0 && (int)sFE2->type == 0){
		Cpointer C1_new, C2_new;
		if (b4) C1_new = checkR5_1->C1->C1;
		else C1_new = C_1->C2->C2->C2;
		C2_new = C1_new->C2->C2->C2;
		//R5 has moved to the edge and will now be free.
		//removeC(C1_new->C2, false); removeC(C1_new->C2, false);
		//addC(C1_new, normAngle(C1_new->C1->bondAngle1 - 60), normAngle(C1_new->C1->bondAngle1), 1, true);
		if (b4) {
			sFE2->C1 = C_2->C1->C1->C1;
			sFE2->C2 = C_2->C1->C1;
			convSiteType(stt, sFE2->C2, C_2, ZZ);
		}
		else {
			sFE2->C1 = C_1->C2->C2;
			sFE2->C2 = C_1->C2->C2->C2;
			convSiteType(stt, C_1, sFE2->C1, ZZ);
		}
		if (!optimised){
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FE);
			//convSiteType(checkR5_1, C1_new, C1_new->C2->C2, ZZ);
			updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, -1);
			m_pah->m_rings5_Lone--;
			redrawR5(checkR5_1, C1_new, C2_new);
			//proc_G5R_ZZ(checkR5_1, checkR5_1->C1, checkR5_1->C2);
		}
		else {
			if (b4) {
				convSiteType(sFE2, sFE2->C1->C1, sFE2->C2, RFE);
				convSiteType(checkR5_1, C1_new->C2, C2_new->C1, R5);
				updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
			}
			else {
				convSiteType(sFE2, sFE2->C1, sFE2->C2->C2, RFE);
				convSiteType(checkR5_1, C1_new->C2, C2_new->C1, R5);
				updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
			}
		}
	}
	else {
		//Reassign connectivity at next site
		if (b4) sFE2->C2 = C_2->C1->C1;
		else sFE2->C1 = C_1->C2->C2;
	}

	// edit sites. first identify the neighbouring sites of resulting RFE & R5
	Spointer S1, S2, S3, S4;
	if (b4) {
		S1 = moveIt(sFE2, -1);
		S2 = moveIt(stt, +1);
		if ((int)sFE2->type == 0 && (int)S1->type == 0){ //sFE2 is a FE
			//This means that the pentagon has migrated to the edge. This is handled above.
			//convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			//updateSites(S1, S1->C1, sFE2->C1, 500);
		}
		else if ((int)sFE2->type == 0 ){ 
			//This means that the pentagon has migrated to the edge but will have a carbon out of the structure.
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			updateSites(S1, S1->C1, sFE2->C1, 500);
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C2->C1);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
			if (!optimised) addR5internal(sFE2->C2->C1, sFE2->C2);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
			if (!optimised) addR5internal(sFE2->C2->C1, sFE2->C2);
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
			if (!optimised) addR5internal(sFE2->C2->C1, sFE2->C2);
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
			if (!optimised) addR5internal(sFE2->C2->C1, sFE2->C2);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
			if (!optimised) addR5internal(sFE2->C2->C1->C1, sFE2->C2->C1);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
			if (!optimised) addR5internal(sFE2->C2->C1->C1, sFE2->C2->C1);
		}
		convSiteType(stt, sFE2->C2, stt->C2, ZZ);
	}
	else {
		S1 = moveIt(stt, -1); 
		S2 = moveIt(sFE2, 1); 
		if ((int)sFE2->type == 0 && (int)S2->type == 0){ //sFE2 is a FE
			//This means that the pentagon has migrated to the edge. This is handled above.
			//convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			//updateSites(S2, sFE2->C2, S2->C2, 500);
		}
		else if ((int)sFE2->type == 0){
			//This means that the pentagon has migrated to the edge but will have a carbon out of the structure.
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			updateSites(S2, sFE2->C2, S2->C2, 500);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C2);
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C2->C1);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C1->C2->C2);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C1->C2->C2);
		}
		convSiteType(stt, stt->C1, sFE2->C1, ZZ);
	}
	// update H atoms
	/*if (b4){
		updateA(S1->C2->C1, C_2, 'H');
	}
	else{
		updateA(C_1, S2->C1->C2, 'H');
	}*/
	if (opp_site_bool && !opp_site_bool_second && !opp_site_bool_after) {
		if ( (int)opp_site->type >= 2100) {
			Spointer S1_opp_site = moveIt(opp_site, -1);
			Spointer S2_opp_site = moveIt(opp_site, +1);
			if (S1_opp_site->type==R5 || S2_opp_site->type==R5){
				updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
			}
			else updateSites(opp_site, opp_site->C1, opp_site->C2, -100);
		}
		else updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		Spointer S1_opp_site = moveIt(opp_site, -1);
		Spointer S2_opp_site = moveIt(opp_site, +1);
		updateCombinedSites(opp_site); updateCombinedSites(S1_opp_site);  updateCombinedSites(S2_opp_site); 
	}
	else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		if (opp_site_after != S1 && opp_site_after != S2) {
			if ((int)opp_site_after->type >= 500 && (int)opp_site_after->type <= 700) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -400);
			else if ((int)opp_site_after->type >= 1000 && (int)opp_site_after->type <= 2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -800);
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			Spointer S1_opp_site = moveIt(opp_site_after, -1);
			Spointer S2_opp_site = moveIt(opp_site_after, +1);
			updateCombinedSites(opp_site_after); updateCombinedSites(S1_opp_site);  updateCombinedSites(S2_opp_site); 
		}
		updateCombinedSites(opp_site);
	}
	else if (opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
		//convSiteType(opp_site_second, opp_site_second->C1, opp_site_second->C2, (kmcSiteType)new_stype);
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, +1500);
		updateCombinedSites(opp_site);
		updateCombinedSites(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateCombinedSites(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, -1500);
		updateCombinedSites(opp_site_second);
		if (opp_site_after != S1 && opp_site_after != S2) {
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
			updateCombinedSites(opp_site_after);
		}
	}
	else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		if (opp_site_after != S1 && opp_site_after != S2) {
			if ( (int)opp_site_after->type >=2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +100);
			else updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			updateCombinedSites(opp_site_after);
		}
	}
	else if (opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
		updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
		updateCombinedSites(opp_site);
		updateCombinedSites(opp_site_second);
		updateCombinedSites(opp_site_after);
	}
	//printStruct();
	// update combined sites for all sites involved and their neighbours
	// (excluding new FE sites, since their combined site type will still be None)
	if (b4){
		S3 = moveIt(S1, -1);
	}
	else{
		S4 = moveIt(S2, 1);
	}
	updateCombinedSites(stt); updateCombinedSites(sFE2); 
	if (b4){
		updateCombinedSites(S1);
		updateCombinedSites(S3);
	}
	else{
		updateCombinedSites(S2);
		updateCombinedSites(S4); // neighbours
	}
}
int G6RRZZ_error_counter = 0;
// ************************************************************
// ID24 - R6 growth on RZZ 
// ************************************************************
void PAHProcess::proc_G6R_RZZ(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	Cpointer newC1;
	Cpointer newC2;
	if (checkHindrance_newCposition(C_1)|| checkHindrance_newCposition(C_2)) {
		/*cout<<"Site hindered, process not performed.\n"*/ return;
	}
	Spointer S1 = moveIt(stt, -1);
	Spointer S2 = moveIt(stt, 1);
	Spointer S1check = moveIt(stt, -1);
	Spointer S2check = moveIt(stt, 1);
	Spointer SR5;
	bool b4 = false;
	if (stt->type == R5R6ZZ) { //G6R_R5R6ZZ
		if ( (((int) S1->type >= 501 && (int) S1->type <= 504) || ((int) S1->type >= 1002 && (int) S1->type <= 1004) || ((int) S1->type >= 2103 && (int) S1->type <= 2205 ) || S1->type == SPIRAL) && (((int) S2->type >= 501 && (int) S2->type <= 504) || ((int) S2->type >= 1002 && (int) S2->type <= 1004) || ((int) S2->type >= 2103 && (int) S2->type <= 2205 ) || S2->type == SPIRAL) ) {
			if ( (isR5internal(C_1, C_1->C2) ) && !(isR5internal(C_2->C1, C_2) )){
				SR5 = S1;
				S1 = moveIt(SR5, -1);
				b4 = true;
			}
			else if ( !(isR5internal(C_1, C_1->C2)) && (isR5internal(C_2->C1, C_2) ) ){
				SR5 = S2;
				S2 = moveIt(SR5, +1);
				b4 = false;
			}
			else {
				//Use site loop to detect right site as second try.
				for (int i = 1; i<=3; i++){
					//Both sides of the R5R6ZZ seem to be possible. Proceed to identify them.
					int S1_before = (int)moveIt(S1check, -1)->type;
					int S2_after = (int)moveIt(S2check, +1)->type;
					if ( ((S1_before >= 501 && S1_before <= 504) || (S1_before >= 1002 && S1_before <= 1004) || (S1_before >= 2103 && S1_before <= 2205 ) ) && ( (S2_after >= 501 && S2_after <= 504) || (S2_after >= 1002 && S2_after <= 1004) || (S2_after >= 2103 && S2_after <= 2205 ) ) ){
						//Both second neighbour sites could be the coupled site to the R5R6ZZ. Move to next neighbours.
					}
					else if ( ((S1_before >= 501 && S1_before <= 504) || (S1_before >= 1002 && S1_before <= 1004) || (S1_before >= 2103 && S1_before <= 2205 ) ) && !( (S2_after >= 501 && S2_after <= 504) || (S2_after >= 1002 && S2_after <= 1004) || (S2_after >= 2103 && S2_after <= 2205 ) ) ){
						SR5 = S2;
						S2 = moveIt(SR5, +1);
						b4 = false;
						break;
					}
					else if ( !((S1_before >= 501 && S1_before <= 504) || (S1_before >= 1002 && S1_before <= 1004) || (S1_before >= 2103 && S1_before <= 2205 ) ) && ( (S2_after >= 501 && S2_after <= 504) || (S2_after >= 1002 && S2_after <= 1004) || (S2_after >= 2103 && S2_after <= 2205 ) ) ){
						SR5 = S1;
						S1 = moveIt(SR5, -1);
						b4 = true;
						break;
					}
					else {
						//No associated site. Odd. Throw error.
						cout << "G6R_R5R6ZZ has a site with no correct neighbours. Error.\n";
						std::list<std::string> Sitelist_before = copySites(stt);
						printBeforeSites(Sitelist_before);
						printSites(stt);
						std::string filename = "KMC_DEBUG/RZZ_neighbour_error_";
						filename.append(std::to_string(G6RRZZ_error_counter));
						cout << "Saving file " << filename << ".xyz\n";
						saveXYZ(filename);
						G6RRZZ_error_counter++;
						return;
					}
					S1check = moveIt(S1check, -1);
					S2check = moveIt(S2check, +1);
				}
			}
		}
		else {
			if ( ((int) S1->type >= 501 && (int) S1->type <= 504 ) || ((int) S1->type >= 602 && (int) S1->type <= 604 ) || ((int) S1->type >= 1002 && (int) S1->type <= 1004) || ((int) S1->type >= 2103 && (int) S1->type <= 2205 )) {
				SR5 = S1;
				S1 = moveIt(SR5, -1);
				b4 = true;
			}
			else if ( ((int) S2->type >= 501 && (int) S2->type <= 504) || ((int) S2->type >= 602 && (int) S2->type <= 604) || ((int) S2->type >= 1002 && (int) S2->type <= 1004 ) || ((int) S2->type >= 2103 && (int) S2->type <= 2205 )) {
				SR5 = S2;
				S2 = moveIt(SR5, +1);
				b4 = false;
			}
			else if ((int)S1->type == 9999 || (int)S2->type == 9999){
				//SPIRAL site, reject.
				return;
			}
			else {
				cout << "R5R6ZZ not neighbouring an R5R6 site. Error.\n";
				std::list<std::string> Sitelist_before = copySites(stt);
				printBeforeSites(Sitelist_before);
				printSites(stt);
				std::string filename = "KMC_DEBUG/R5R6ZZ_neighbour_error_";
				filename.append(std::to_string(G6RRZZ_error_counter));
				cout << "Saving file " << filename << ".xyz\n";
				saveXYZ(filename);
				G6RRZZ_error_counter++;
				return;
			}
		}
	}
	else { // G6R_RZZ
		if ( S1->type == R5 && S2->type != R5) {
			SR5 = S1;
			S1 = moveIt(SR5, -1);
			b4 = true;
		}
		else if ( S2->type == R5 && S1->type != R5) {
			SR5 = S2;
			S2 = moveIt(SR5, +1);
			b4 = false;
		}
		else {
			cout << "RZZ not neighbouring an R5 site. Error.\n";
			std::list<std::string> Sitelist_before = copySites(stt);
			printBeforeSites(Sitelist_before);
			printSites(stt);
			std::string filename = "KMC_DEBUG/RZZ_neighbour_error_";
			filename.append(std::to_string(G6RRZZ_error_counter));
			cout << "Saving file " << filename << ".xyz\n";
			saveXYZ(filename);
			G6RRZZ_error_counter++;
			return;
		}
	}
	cpair FEvector = get_vector(C_1->C2->coords,C_2->C1->coords);
	cpair AGvector = get_vector(C_1->C2->coords,C_1->coords);
	cpair start_direction = add_vector(FEvector, AGvector);
	cpair AJvector = get_vector(C_2->C1->coords,C_2->coords);
	cpair lat_direction = add_vector(AJvector, invert_vector(AGvector));
	double dist = getDistance_twoC(C_1,C_2);
	if (!(C_1->C2->bridge) || !(C_2->C1->bridge)) { // check if bulk C in AC site is a bridge
		// Add and remove C
		//if(C_1->C2->C3 != NULL) C_1->C2->C3->C3 = NULL;
		//if(C_2->C1->C3 != NULL) C_2->C1->C3->C3 = NULL;
		removeC(C_1->C2, true);
		removeC(C_2->C1, true);
		//newC1 = addC(C_1, normAngle(C_1->bondAngle1 + 120), 0, 1.4);
		//newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1 - 60), normAngle(newC1->C1->bondAngle1 - 120), 1.4);
	}
	else {
		// update bridges info, both are no longer bridges
		C_1->C2->C1 = C_1->C2->C3;// update neighbour
		C_2->C1->C2 = C_2->C1->C3;
		C_1->C2->bridge = false;
		C_2->C1->bridge = false;
		//angletype a = C_2->C1->bondAngle1;
		//C_2->C1->bondAngle1 = C_2->C1->bondAngle2;
		//C_2->C1->bondAngle2 = a;
		// connect C_1 and C_2
		connectToC(C_1, C_2);
		//newC1 = addC(C_1, normAngle(C_1->bondAngle1 + 120), 0, 1.4);
		//newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1 - 60), normAngle(newC1->C1->bondAngle1 - 120), 1.4);
	}
	// Add C & H
	newC1 = addC(C_1, start_direction, dist/2.0);
	updateA(newC1, 'H', AGvector);
	updateA(C_1, 'C', C_1->growth_vector);
	newC2 = addC(newC1, lat_direction, dist/2.0);
	updateA(newC2, 'H', AJvector);
	updateA(C_2, 'C', C_1->growth_vector);
	
	//printStruct();
	// Add and remove HR
	//updateA(C_1->C1, C_2->C2, 'H');
	// neighbouring sites:
	
	Spointer S3, S4;
	S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
	// Update Site and neighbours
	convSiteType(stt, newC1, newC2, FE); // Converts RZZ to FE
	if (b4){
		if (SR5->type == R5) {
			updateSites(SR5, SR5->C1, stt->C1, +401);
			if ( (int)S1->type < 2000) updateSites(S1, S1->C1, S1->C2, +400);
			//else updateSites(S1, S1->C1, stt->C1, +1); THIS WAS WRONG.
		}
		else if ((int) SR5->type >= 501 && (int) SR5->type <= 504 ){
			updateSites(SR5, SR5->C1, stt->C1, +1501);
		}
		else if ((int) SR5->type >= 602 && (int) SR5->type <= 604 ){
			updateSites(SR5, SR5->C1, stt->C1, +1501);
		}
		else if ((int) SR5->type >= 1002 && (int) SR5->type <= 1004) {
			updateSites(SR5, SR5->C1, stt->C1, +1101);
		}
		else updateSites(SR5, SR5->C1, stt->C1, +1);
		updateSites(S2, stt->C2, S2->C2, +1);
	}
	else{
		if (SR5->type == R5) {
			updateSites(SR5, stt->C2, SR5->C2, +401);
			if ( (int)S2->type < 2000) updateSites(S2, S2->C1, S2->C2, +400);
			//else updateSites(S2, stt->C2, S2->C2, +1);  THIS WAS WRONG.
		}
		else if ((int) SR5->type >= 501 && (int) SR5->type <= 504 ){
			updateSites(SR5, stt->C2, SR5->C2, +1501);
		}
		else if ((int) SR5->type >= 602 && (int) SR5->type <= 604 ){
			updateSites(SR5, stt->C2, SR5->C2, +1501);
		}
		else if ((int) SR5->type >= 1002 && (int) SR5->type <= 1004) {
			updateSites(SR5, stt->C2, SR5->C2, +1101);
		}
		else updateSites(SR5, stt->C2, SR5->C2, +1);
		updateSites(S1, S1->C1, stt->C1, +1);
	}

	// Update combined site for Site and neighbours
	updateCombinedSites(stt); updateCombinedSites(SR5);
	updateCombinedSites(S1); updateCombinedSites(S2);
	updateCombinedSites(S3); updateCombinedSites(S4);
	// add ring counts
	m_pah->m_rings++;
	//printSites(stt);
	//printStruct();
	if (getDistance_twoC(stt->C1->C1, stt->C2->C2) >= 3.0 && !m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
}

// ************************************************************
// ID25 - R6 growth on RFER 
// ************************************************************
void PAHProcess::proc_G6R_RFER(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_G6R_AC(stt, C_1, C_2);
}

// ************************************************************
// ID26 - R6 growth on R5
// ************************************************************
void PAHProcess::proc_G6R_R5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	
	//    printSites(stt);
	Cpointer newC1, newC2, newC3, newC4;
	if (checkHindrance(stt)) {
		//cout<<"Site hindered, process not performed.\n"
		return;
	}
	//First we need to resize the R5 bond length or our R6s will be huge.
	//moveC(C_1, C_1->bondAngle1, C_1->bondAngle2, 0.5*(pow(3,0.5)-1));
	//moveC(C_2, C_1->bondAngle1, C_1->bondAngle2, -0.5*(pow(3, 0.5) - 1));

	// Add C
	newC1 = addC(C_1, normAngle(C_1->bondAngle1 + 120), 0, 1.4);
	newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1 - 60), 0, 1.4);
	newC3 = addC(newC2, normAngle(newC2->C1->bondAngle1 - 60), 0, 1.4);
	newC4 = addC(newC3, normAngle(newC3->C1->bondAngle1 - 60), normAngle(newC3->C1->bondAngle1 - 120), 1.4);
	//C_1->C3 = C_2; C_2->C3 = C_1;
	// Add and remove H
	updateA(C_1->C1, C_2->C2, 'H');
	// neighbouring sites:
	Spointer S1 = moveIt(stt, -1);
	Spointer S2 = moveIt(stt, 1);
	// Update Site and neighbours
	updateSites(stt, newC2, newC3, 0);
	updateSites(S1, S1->C1, newC1, 1);
	updateSites(S2, newC4, S2->C2, 1);
	// Add new Sites
	Spointer newS1 = addSite(FE, newC3, newC4, S2);
	Spointer newS2 = addSite(FE, newC2, newC3, newS1);
	Spointer newS3 = addSite(FE, newC1, newC2, newS2);
	
	// Update combined sites for all new sites and original neighbours
	Spointer S3, S4;
	S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
	updateCombinedSites(stt);
	updateCombinedSites(newS1); updateCombinedSites(newS2); updateCombinedSites(newS3); // new sites
	updateCombinedSites(S1); updateCombinedSites(S2); // original neighbours
	updateCombinedSites(S3); updateCombinedSites(S4); // neighbours of neighbours
	// Add H count
	addCount(0, 2);
	// add ring counts
	m_pah->m_rings++;
}

// ************************************************************
// ID27 - RBY5 closure reaction
// ************************************************************
void PAHProcess::proc_L6_RBY5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_L6_BY6(stt, C_1, C_2);
}

// ************************************************************
// ID28 - RACR closure reaction
// ************************************************************
void PAHProcess::proc_L6_RACR(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_L6_BY6(stt, C_1, C_2);
}

// ************************************************************
// ID29 - R5 growth on RFE 
// ************************************************************
void PAHProcess::proc_G5R_RFE(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_G5R_ZZ(stt, C_1, C_2);
}

// ************************************************************
// ID30 - R6 migration & conversion to R5 at RAC
// ************************************************************
void PAHProcess::proc_C6R_RAC_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_C6R_BY5_FE3(stt, C_1, C_2, rng);
}

// ************************************************************
// ID31 - R6 migration & conversion to R5 at RAC
// ************************************************************
void PAHProcess::proc_C6R_RAC_FE3violi(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_C6R_BY5_FE3(stt, C_1, C_2, rng);
}

// ************************************************************
// ID32 - R6 desorption at RAC -> pyrene
// ************************************************************
void PAHProcess::proc_M6R_RAC_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_M6R_BY5_FE3(stt, C_1, C_2, rng);
}

// ************************************************************
// ID34- R5 exchange with R6 
// ************************************************************
void PAHProcess::proc_MR5_R6(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//printStruct();
	/*OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol, 1000);
	passbackPAH(mol);*/
	for (Cpointer Ccheck = C_1->C1; Ccheck != C_2; Ccheck = Ccheck->C2){
		if (getDistance_twoC(Ccheck, Ccheck->C2)>1.6 || getDistance_twoC(Ccheck, Ccheck->C2)<1.2 ){
			if (!m_pah->m_optimised){
				OpenBabel::OBMol mol = passPAH();
				mol = optimisePAH(mol);
				passbackPAH(mol);
			}
		}
	}
	//First check if R6 is to the left or the right of R5
	bool b4 = false;
	Spointer sFE2, checkR5_1, checkR5_2;
	Cpointer CRem, CRem_before, CRem_next, CFE, CR5_otherside_1, CR5_otherside_2;
	if (( isR5internal(C_1->C1, C_1) ) && ( isR5internal(C_2,C_2->C2) )){
		//Pentagons to both sides, JP not allowed
		return;
	}
	else if ( isR5internal(C_1->C1, C_1) ) b4 = false;
	else if ( isR5internal(C_2, C_2->C2) ) b4 = true;
	else return;
	if (b4) {
		sFE2 = moveIt(stt, -1);
		checkR5_1 = moveIt (stt, -2);
		checkR5_2 = moveIt (stt, -3);
		CFE = C_1->C2;
		CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
		CR5_otherside_1 = C_2->C2;
		CR5_otherside_2 = C_2->C1;
	}
	else {
		sFE2 = moveIt(stt, 1);
		checkR5_1 = moveIt (stt, 2);
		checkR5_2 = moveIt (stt, 3);
		CFE = C_1;
		CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
		CR5_otherside_1 = C_1->C1;
		CR5_otherside_2 = C_1->C2;
	}
	//Check for unsupported sites. This section heavily assumes that the Isolated Pentagon Rule is valid.
	if ((int)sFE2->type == 0 && (int)checkR5_1->type == 0 && (int)checkR5_2->type == 0) return; // The result would be an indene, not supported. YET!
	if ((int)sFE2->type == 101 || (int)sFE2->type == 501 || (int)sFE2->type == 2002 || (int)sFE2->type == 1002) return; // This would violate the IPR.
	if ((int)sFE2->type == 0){
		if ((int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 || (int)checkR5_1->type == 100) return;
		//if ((int)checkR5_1->type >= 501 && (int)checkR5_1->type <= 65) return;
		if ((int)checkR5_1->type >= 1002 && (int)checkR5_1->type <= 1004) return;
		if ((int)checkR5_1->type >= 2002 && (int)checkR5_1->type <= 2204) return;
		if ((int)checkR5_1->type >= 2204 && (int)checkR5_1->type <= 2205) return;
		if ((int)checkR5_1->type == 0){
			//if ((int)checkR5_2->type >= 501 && (int)checkR5_2->type <= 65) return;
			if ((int)checkR5_2->type == 101 || (int)checkR5_2->type == 501) return;
			if ((int)checkR5_2->type >= 1002 && (int)checkR5_2->type <= 1004) return;
			if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2205) return;
		}
	}
	if (CRem_next->bridge) return;
	if (CRem_next->C2->bridge) return; if (CRem_next->C1->bridge) return;
	if (CRem_next->C2->C2->bridge) return; if (CRem_next->C1->C1->bridge) return;

	// check if R5R6 has an opposite site.
	Spointer opp_site, opp_site_second, opp_site_after;
	bool opp_site_bool = false; bool opp_site_bool_second = false; bool opp_site_bool_after = false;
	Cpointer thirdC = findThirdC(CR5_otherside_1);
	Cpointer thirdC2 = findThirdC(CR5_otherside_2);
	Cpointer thirdC_after = findThirdC(CRem_next);
	// Seven cases:
	// 1. One pentagon has one exposed edge  and migrates to a location where it will have one exposed edge. Normal migration.
	// 2. One pentagon has one exposed edge  and migrates to a location where it will have two exposed edges. 
	// 3. One pentagon has two exposed edges and migrates to a location where it will have two exposed edges. 
	// 4. One pentagon has two exposed edges and migrates to a location where it will have one exposed edge.
	// 5. One pentagon has two exposed edges and migrates to a location where it will have three exposed edges.
	// 6. One pentagon has three exposed edges and migrates to a location where it will have three exposed edges.
	// 7. One pentagon has three exposed edges and migrates to a location where it will have two exposed edges.
	
	if (thirdC != NULLC) {
		opp_site = findSite(thirdC);
		if (opp_site != m_pah->m_siteList.end()) opp_site_bool = true;
	}
	if (thirdC2 != NULLC) {
		opp_site_second = findSite(thirdC2);
		if (opp_site_second != m_pah->m_siteList.end() && opp_site_second!=opp_site) opp_site_bool_second = true;
	}
	if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
		opp_site_after = findSite(thirdC_after);
		if (opp_site_after != m_pah->m_siteList.end()){
			opp_site_bool_after = true;
			int os_endtype = opp_site_after->type;
			if (os_endtype >= 200 && os_endtype <= 203) return;
			if (os_endtype == 101) return;
			if (os_endtype >= 600 && os_endtype <= 603) return;
			if (os_endtype >= 1000 && os_endtype <= 1003) return;
			if (os_endtype >= 500 && os_endtype <= 504) return;
			if (os_endtype >= 2000 && os_endtype <= 2205) return;
			if (os_endtype >= 2103 && os_endtype <= 2105) return;
			if (os_endtype >= 2204 && os_endtype <= 2205) return;
		}
	}
	
	//Fundamental assumption: R5-R7 pairs cannot move away from each other!
	cpair R5coords = findR5internal(CFE, CFE->C2);
	if (m_pah->m_R7loc.size()>=1){
		for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
			double distR5R7 = getDistance_twoC(*it, R5coords);
			if (distR5R7 < 3.1) {
				m_pah->m_R5loc.push_back(R5coords);
				return;
			}
		}
	}
	//check that two pentagons (including internals) will not collide
	cpair R5coords_end = endposR5internal(CRem->C1, CRem);
	if (m_pah->m_R5loc.size()>=1){
		for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
			double distR5s = getDistance_twoC(*it, R5coords_end);
			if (distR5s < 2.8) {
				//This distance is a parameter of this jump process. Might need some more tuning. 
				//2.8 seems appropiate but may reject too many jumps.
				//Two pentagons will be next to each other violating the Isolated Pentagon Rule
				m_pah->m_R5loc.push_back(R5coords);
				return;
			}
		}
	}
	
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2 = 1.4;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir = get_vector(CFE->coords, CFE->C2->coords);
	cpair normvec;
	if (b4) normvec = invert_vector((norm_vector(CFE->coords, CFE->C2->coords, CFE->C2->C2->coords)));
	else normvec = (norm_vector(CFE->coords, CFE->C2->coords, CFE->C2->C2->coords));
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, 1.4);
	updateA(Cnew, 'H', crossvec);
	removeC(CRem, false);
	if (!m_pah->m_optimised){
		OpenBabel::OBMol newmol = passPAH();
		newmol = optimisePAH(newmol,250);
		passbackPAH(newmol);
	}
	
	//addC(CFE, normAngle(CFE->bondAngle1 + 30), normAngle(CFE->bondAngle1 - 30), 1.4);
	//CRem->C1->bondAngle1 = normAngle(CRem->C1->bondAngle1 - 30);
	//printStruct(CRem);
	if (b4) sFE2->C2 = C_2->C1;
	else sFE2->C1 = C_1->C2;

	// edit sites. first identify the neighbouring sites of resulting RFE & R5
	Spointer S1, S2, S3, S4;
	if (b4) {
		S1 = moveIt(sFE2, -1);
		if ((int)sFE2->type == 0){ //sFE2 is a FE
			if ((int)S1->type == 0){ //S1 is a FE
				Cpointer C1_new, C2_new;
				C1_new = checkR5_1->C1->C1;
				C2_new = C1_new->C2->C2->C2;
				sFE2->C1 = C_2->C1->C1;
				sFE2->C2 = C_2->C1;
				convSiteType(stt, sFE2->C2, C_2, FE);
				convSiteType(sFE2, sFE2->C1->C1, sFE2->C2, RFE);
				convSiteType(checkR5_1, C1_new->C2, C2_new->C1, R5);
				updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
			}
			else{
				convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
				int stype = S1->type;
				if (stype >= 2002 && stype <= 2204) {
					stype = stype + 100;
					convSiteType(S1, S1->C1, sFE2->C1, (kmcSiteType)stype);
				}
				else {
					stype = stype + 500;
					convSiteType(S1, S1->C1, sFE2->C1, (kmcSiteType)stype);
					//updateSites(S1, S1->C1, sFE2->C1, 5);
				}
			}
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
		}
		convSiteType(stt, stt->C2->C1, stt->C2, FE); //stt is originally the R5R6 site that will become the new FE site
		S2 = moveIt(stt, 1); // neighbour of stt
		int stype;
		if ( (int)S2->type >= 2103 && (int)S2->type <= 2205 ) stype = S2->type - 100;
		else stype = S2->type - 500;
		convSiteType(S2, stt->C2, S2->C2, (kmcSiteType)stype);
		/*if ((int)S2->type >= 10 && (int)S2->type <= 12){ //The site in the middle had 2R5s to each site, so lucky!
			updateSites(S2, stt->C2, S2->C2, -3); //convert the neighbour to its same version but next to an R5
		}
		else{
			int stype = S2->type - 55;
			convSiteType(S2, stt->C2, S2->C2, (kmcSiteType)stype);
			//updateSites(S2, stt->C2, S2->C2, -5); //convert the neighbour to its same version but next to an R5
		}
		*/
	}
	else {
		S2 = moveIt(sFE2, 1);
		if ((int)sFE2->type == 0){ //sFE2 is a FE
			if ((int)S2->type == 0){ //S2 is a FE
				Cpointer C1_new, C2_new;
				C1_new = C_1->C2->C2;
				C2_new = C1_new->C2->C2->C2;
				sFE2->C1 = C_1->C2;
				sFE2->C2 = C_1->C2->C2;
				convSiteType(stt, C_1, sFE2->C1, ZZ);
				convSiteType(sFE2, sFE2->C1, sFE2->C2->C2, RFE);
				convSiteType(checkR5_1, C1_new->C2, C2_new->C1, R5);
				updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
			}
			else {
				convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
				int stype = S2->type;
				if (stype >= 2002 && stype <= 2204) {
					stype = stype + 100;
					convSiteType(S2, sFE2->C2, S2->C2, (kmcSiteType)stype);
				}
				else {
					stype = stype + 500;
					convSiteType(S2, sFE2->C2, S2->C2, (kmcSiteType)stype);
					//updateSites(S2, sFE2->C2, S2->C2, 5);
				}
			}
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
		}
		convSiteType(stt, stt->C1, stt->C1->C2, FE); //stt is originally the R5R6 site that will become the new FE site	
		S1 = moveIt(stt, -1); // neighbour of stt
		int stype;
		if ( (int)S1->type >= 2103 && (int)S1->type <= 2205 ) stype = S1->type - 100;
		else stype = S1->type - 500;
		convSiteType(S1, S1->C1, stt->C1, (kmcSiteType)stype);
		/*if ((int)S1->type >= 10 && (int)S1->type <= 12){ //The site in the middle had 2R5s to each site, so lucky!
			updateSites(S1, S1->C1, stt->C1, -3); //convert the neighbour to its same version but next to ONLY ONE R5
		}
		else{
			int stype = S1->type - 55;
			convSiteType(S1, S1->C1, stt->C1, (kmcSiteType)stype);
			//updateSites(S1, S1->C1, stt->C1, -5); //convert the neighbour to its same version but NOT next to an R5
		}*/
	}
	if (opp_site_bool && !opp_site_bool_second && !opp_site_bool_after) {
		if ( (int)opp_site->type >= 2100) {
			Spointer S1_opp_site = moveIt(opp_site, -1);
			Spointer S2_opp_site = moveIt(opp_site, +1);
			if (S1_opp_site->type==R5 || S2_opp_site->type==R5){
				updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
			}
			else updateSites(opp_site, opp_site->C1, opp_site->C2, -100);
		}
		else updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		Spointer S1_opp_site = moveIt(opp_site, -1);
		Spointer S2_opp_site = moveIt(opp_site, +1);
		updateCombinedSites(opp_site); updateCombinedSites(S1_opp_site); updateCombinedSites(S2_opp_site);
	}
	else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		if (opp_site_after != S1 && opp_site_after != S2) {
			if ((int)opp_site_after->type >= 500 && (int)opp_site_after->type <= 700) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -400);
			else if ((int)opp_site_after->type >= 1000 && (int)opp_site_after->type <= 2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -800);
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			Spointer S1_opp_site = moveIt(opp_site_after, -1);
			Spointer S2_opp_site = moveIt(opp_site_after, +1);
			updateCombinedSites(opp_site_after); updateCombinedSites(S1_opp_site); updateCombinedSites(S2_opp_site);
		}
		updateCombinedSites(opp_site);
	}
	else if (opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		if (opp_site != opp_site_second){
			updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
			//convSiteType(opp_site_second, opp_site_second->C1, opp_site_second->C2, (kmcSiteType)new_stype);
			updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, +1500);
			updateCombinedSites(opp_site);
			updateCombinedSites(opp_site_second);
		}
	}
	else if (!opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateCombinedSites(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, -1500);
		updateCombinedSites(opp_site_second);
		if (opp_site_after != S1 && opp_site_after != S2) {
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
			updateCombinedSites(opp_site_after);
		}
	}
	else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		if (opp_site_after != S1 && opp_site_after != S2) {
			if ( (int)opp_site_after->type >=2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +100);
			else updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			updateCombinedSites(opp_site_after);
		}
	}
	else if (opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		if ( (opp_site != opp_site_after) && (opp_site != opp_site_second) ){
			updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
			updateCombinedSites(opp_site);
			updateCombinedSites(opp_site_second);
			updateCombinedSites(opp_site_after);
		}
	}
	
	// update H atoms
	/*if (b4){
		updateA(S1->C2->C1, C_2, 'H');
	}
	else{
		updateA(C_1, S2->C1->C2, 'H');
	}*/
	//printStruct();
	// update combined sites for all sites involved and their neighbours
	// (excluding new FE sites, since their combined site type will still be None)
	if (b4){
		S3 = moveIt(S1, -1);
	}
	else{
		S4 = moveIt(S2, 1);
	}
	updateCombinedSites(stt); updateCombinedSites(sFE2); updateCombinedSites(S1); updateCombinedSites(S2);
	if (b4){
		updateCombinedSites(S3);
	}
	else{
		updateCombinedSites(S4); // neighbours
	}
}

int GR7_R5R6AC_error_counter = 0;
// ************************************************************
// ID35- R7 growth on embedded-obstructed R5
// ************************************************************
void PAHProcess::proc_GR7_R5R6AC(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//printSitesMemb(stt);
	//printStruct();//++++
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH(); 
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	Cpointer newC1;
	Cpointer newC2;
	cpair Hdir1 = get_vector(C_1->C2->coords, C_1->coords);
	cpair Hdir2 = get_vector(C_2->C1->coords, C_2->coords);
	cpair starting_direction = get_vector(C_1->C2->C2->coords,C_2->coords);
	cpair FEdir = get_vector(C_1->coords,C_2->coords);
	if(checkHindrance_newCposition(C_1) || checkHindrance_newCposition(C_2)) {
		/*cout<<"Site hindered, process not performed.\n"*/ return;
	}
	if (!(C_1->C2->bridge) && !(C_2->C1->bridge)) { // check if bulk C in AC site is a bridge
		// Add and remove C
		//if(C_1->C2->C3 != NULL) C_1->C2->C3->C3 = NULL;
		//if(C_2->C1->C3 != NULL) C_2->C1->C3->C3 = NULL;
		moveC_z(C_1->C2, 0.2);
		removeC(C_1->C2, true);
		moveC_z(C_1->C2, 0.2);
		removeC(C_1->C2, true);
		moveC_z(C_2->C1, 0.2);
		removeC(C_2->C1, true);
		newC1 = addC(C_1, starting_direction, 1.55);
		updateA(C_1,'C', C_1->growth_vector);
		updateA(newC1, 'H', Hdir1);
		newC2 = addC(newC1, FEdir, 1.4);
		updateA(C_2,'C', C_2->growth_vector);
		updateA(newC2, 'H', Hdir2);
		//addOBbond(newC2, C_2, mol);
		/*removeC(C_1->C2, true);
		removeC(C_1->C2, true);
		removeC(C_2->C1, true);
		//moveC(C_1, C_1->C1, 1.5);
		newC1 = addC(C_1, normAngle(C_1->bondAngle1 + 120 - atan(0.5/1.5) * 180.0/M_PI), 0, pow(2.5,0.5));
		newC2 = addC(newC1, normAngle(newC1->C1->C1->bondAngle1 - 30), normAngle(newC1->C1->C1->bondAngle1 - 90), 1.4*2.0);*/
	}
	else {
		// update bridges info, both are no longer bridges
		Cpointer Cbridge, Cbridge2, Cnotbridge;
		if (C_1->C2->bridge) {
			Cbridge = C_1->C2;
			Cbridge2 = C_1->C2->C3;
			Cnotbridge = C_2->C1;
		} else{
			Cbridge2 = C_2->C1;
			Cbridge = C_2->C1->C3;
			Cnotbridge = C_1->C2;
		}
		Cbridge->C1 = Cbridge2;// update neighbour
		Cbridge->bridge = false;
		Cbridge->C3 = NULL;
		Cbridge2->C2 = Cbridge;
		Cbridge2->bridge = false;
		Cbridge2->C3 = NULL;
		Cnotbridge->C1 = C_1;
		Cnotbridge->C2 = C_2;
		removeC(Cnotbridge, true);
		// connect C_1 and C_2
		connectToC(C_1, C_2);
		// Add C
		//newC1 = addC(C_1, normAngle(C_1->bondAngle1 + 90), 0, 1.4*1.5);
		//double distR7 = getDistance_twoC(C_1, C_2);
		newC1 = addC(C_1, starting_direction, 1.4);
		updateA(C_1,'C', C_1->growth_vector);
		updateA(newC1, 'H', Hdir1);
		newC2 = addC(newC1, FEdir, 1.4);
		updateA(C_2,'C', C_2->growth_vector);
		updateA(newC2, 'H', Hdir2);
		//addOBbond(newC2, C_2);
	}
	if (!m_pah->m_optimised){
		OpenBabel::OBMol newmol = passPAH();
		newmol = optimisePAH(newmol);
		passbackPAH(newmol);
	}
	//printStruct();
	// Add and remove H
	//updateA(C_1->C1, C_2->C2, 'H');
	// neighbouring sites:
	Spointer S1 = moveIt(stt, -1);
	Spointer S2 = moveIt(stt, 1);
	Spointer SR5;
	int stype = (int)stt->type;
	// Update Site and neighbours
	convSiteType (stt, newC1, newC2, FE); //FE of an heptagon
	if (stype == 503){
		if ( isR5internal(C_1->C1, C_1) ){
			//R5 to the left.
			if ( (int)S1->type >= 501 && (int)S1->type <= 504){ 
			//Side of an R5R6AC
			int stype = (int)S1->type;
			stype = (stype % 10) + 2001;
			convSiteType(S1, S1->C1, newC1, (kmcSiteType)stype); // neighbours
			}
			else if ( (int)S1->type >= 602 && (int)S1->type <= 1004){ 
				int stype = (int)S1->type;
				stype = (stype % 10) + 2101;
				convSiteType(S1, S1->C1, newC1, (kmcSiteType)stype); // neighbours
			}		
			else updateSites(S1, S1->C1, newC1, 1); // neighbours
			updateSites(S2, newC2, S2->C2, 1);
		}
	
		else if ( isR5internal(C_2,C_2->C2) ){
			//R5 to the right.
			if ( (int)S2->type >= 501 && (int)S2->type <= 504){ 
				//Side of an R5R6AC
				int stype = (int)S2->type;
				stype = (stype % 10) + 2001;
				convSiteType(S2, newC2, S2->C2, (kmcSiteType)stype); // neighbours
			}
			else if ( (int)S2->type >= 602 && (int)S2->type <= 1004){ 
				int stype = (int)S2->type;
				stype = (stype % 10) + 2101;
				convSiteType(S2, newC2, S2->C2, (kmcSiteType)stype); // neighbours
			}	
			else updateSites(S2, newC2, S2->C2, 1); // neighbours
			updateSites(S1, S1->C1, newC1, 1);
		}
		else {
			if ( (int)S1->type >= 501 && (int)S2->type <= 501){
				//R5 to the left.
				if ( (int)S1->type >= 501 && (int)S1->type <= 504){ 
				//Side of an R5R6AC
				int stype = (int)S1->type;
				stype = (stype % 10) + 2001;
				convSiteType(S1, S1->C1, newC1, (kmcSiteType)stype); // neighbours
				}
				else if ( (int)S1->type >= 602 && (int)S1->type <= 1004){ 
					int stype = (int)S1->type;
					stype = (stype % 10) + 2101;
					convSiteType(S1, S1->C1, newC1, (kmcSiteType)stype); // neighbours
				}		
				else updateSites(S1, S1->C1, newC1, 1); // neighbours
				updateSites(S2, newC2, S2->C2, 1);
			}
			else if ( (int)S2->type >= 501 && (int)S1->type <= 501){
				//R5 to the right.
				if ( (int)S2->type >= 501 && (int)S2->type <= 504){ 
					//Side of an R5R6AC
					int stype = (int)S2->type;
					stype = (stype % 10) + 2001;
					convSiteType(S2, newC2, S2->C2, (kmcSiteType)stype); // neighbours
				}
				else if ( (int)S2->type >= 602 && (int)S2->type <= 1004){ 
					int stype = (int)S2->type;
					stype = (stype % 10) + 2101;
					convSiteType(S2, newC2, S2->C2, (kmcSiteType)stype); // neighbours
				}	
				else updateSites(S2, newC2, S2->C2, 1); // neighbours
				updateSites(S1, S1->C1, newC1, 1);
			}
			else {
				cout << "R5 not found on site R5R6AC for GR7_R5R6AC.\n";
				if(m_debug_pah){
					ifstream  src("KMC_DEBUG/BEFORE.xyz");
					std::string filename = "KMC_DEBUG/BEFORE_GR7_R5R6AC_error_";
					filename.append(std::to_string(GR7_R5R6AC_error_counter));
					filename.append(".xyz");
					ofstream dst(filename);
					dst << src.rdbuf();
					src.close();
					dst.close();
					cout << "Saving file " << filename <<"\n";
				}
				std::string fileout = "KMC_DEBUG/GR7_R5R6AC_error_";
				fileout.append(std::to_string(GR7_R5R6AC_error_counter));
				cout << "Saving file " << fileout << ".xyz\n";
				saveXYZ(fileout); //SETBREAKPOINT
				GR7_R5R6AC_error_counter++;
				cout << "Printing internal R5 positions:.\n";
				for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
					cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
				}
				cpair endR5_1 = endposR5internal(C_1->C1, C_1, true);
				cpair endR5_2 = endposR5internal(C_1->C1, C_1, false);
				cpair endR5_3 = endposR5internal(C_2, C_2->C2, true);
				cpair endR5_4 = endposR5internal(C_2, C_2->C2, false);
				cout << "Looking in positions:\n";
				cout <<	std::get<0>(endR5_1) << ", " << std::get<1>(endR5_1) << ", " << std::get<2>(endR5_1) <<"\n";
				cout <<	std::get<0>(endR5_2) << ", " << std::get<1>(endR5_2) << ", " << std::get<2>(endR5_2) <<"\n";
				cout <<	std::get<0>(endR5_3) << ", " << std::get<1>(endR5_3) << ", " << std::get<2>(endR5_3) <<"\n";
				cout <<	std::get<0>(endR5_4) << ", " << std::get<1>(endR5_4) << ", " << std::get<2>(endR5_4) <<"\n";
			}
		}
	}
	else {
		// GR7 on FEACR5
		updateSites(S1, S1->C1, newC1, 1);
		updateSites(S2, newC2, S2->C2, 1);
	}
	// Update combined site for Site and neighbours
	Spointer S3, S4;
	S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
	updateCombinedSites(stt);
	updateCombinedSites(S1); updateCombinedSites(S2);
	updateCombinedSites(S3); updateCombinedSites(S4);
	// add ring counts
	m_pah->m_rings7_Embedded++;
	//printSites(stt);
}

// ************************************************************
// ID35- R7 growth on embedded-obstructed R5
// ************************************************************
void PAHProcess::proc_GR7_FEACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_GR7_R5R6AC(stt, C_1, C_2);
	m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
}

// ************************************************************
// ID37- R6 growth on R5R6ZZ
// ************************************************************
void PAHProcess::proc_G6R_R5R6ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_G6R_RZZ(stt, C_1, C_2);
}

// 
// ************************************************************
// ID38- R7 bay closure on ACACR5
// ************************************************************
void PAHProcess::proc_L7_ACACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {

	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	proc_L6_BY6(stt, C_1, C_2);
	m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	m_pah->m_rings--; m_pah->m_rings7_Embedded++;
}

// ************************************************************
// ID39 - R6 growth on R5R6FER 
// ************************************************************
void PAHProcess::proc_G6R_R5R6FER(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_G6R_AC(stt, C_1, C_2);
}

// ************************************************************
// ID40 - R6 growth on R5R6FER5R6 
// ************************************************************
void PAHProcess::proc_G6R_R5R6FER5R6(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_G6R_AC(stt, C_1, C_2);
}

// ************************************************************
// ID41- R7 bay closure on FEZZACR5
// ************************************************************
void PAHProcess::proc_L7_FEZZACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	proc_L6_BY6(stt, C_1, C_2);
	m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	m_pah->m_rings--; m_pah->m_rings7_Embedded++;
}

// ************************************************************
// ID42 - R5 conversion to R6 on RZZR (2R5 collision)
// ************************************************************
void PAHProcess::proc_C5R_RZZR(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
    //printSites(stt);
	//printStruct();
	bool b4;
	if (stt->type == RZZR){
		// Define a distribution that has two equally probably outcomes
		boost::bernoulli_distribution<> choiceDistrib;
		// Now build an object that will generate a sample using rng
		boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
		b4 = choiceGenerator();
	}
	else{
		//Site must be R5R6ZZR
		Spointer S1 = moveIt(stt, -1);
		Spointer S2 = moveIt(stt, +1);
		if (S1->type == R5 && (int)S2->type >= 501) b4 = true;
		else if ( (int)S1->type >= 501 && S2->type == R5) b4 = false;
		else {
			cout<<"proc_C5R_RZZR has inconsistent sites. Process not performed.\n";
			printSites(stt);
			return;
		}
	}
    Spointer sR5;
    Cpointer C_AC;
	cpair starting_direction, Hdir1, Hdir2, FEdir, Hdir, CZZdir;
    if(b4) {
        sR5 = moveIt(stt, -1);
        C_AC = sR5->C1->C1;
		starting_direction = get_vector(C_2->C1->coords, C_2->coords);
		Hdir1 = get_vector(C_1->C2->C2->coords, C_1->C2->coords);
		FEdir = get_vector(C_1->C2->coords, C_2->coords);
		Hdir2 = starting_direction;
		Hdir = starting_direction;
		CZZdir = get_vector(C_1->C2->coords, C_1->C2->C2->coords);
    }else {
        sR5 = moveIt(stt, 1);
        C_AC = sR5->C2->C2;
		starting_direction = get_vector(C_2->C1->C1->coords, C_2->C1->coords);
		Hdir1 = get_vector(C_1->C2->coords, C_1->coords);
		FEdir = get_vector(C_1->coords, C_2->C1->coords);
		Hdir2 = starting_direction;
		Hdir = Hdir1;
		CZZdir = FEdir;
    }
    // check if there's a bridge in the BY5
    bool bridge = false;
    if(b4) bridge = C_2->C1->bridge;
    else bridge = C_1->C2->bridge;
    // remove R5 first, leaving a ZZ site (adding another C atom after removing C)
    Cpointer Cstart;
    if(b4) Cstart = C_AC; 
    else Cstart = C_2->C1;
	removeR5internal(sR5->C1, sR5->C2);
    for(int i=0; i!=2; i++) removeC(Cstart->C2, false);
	addC(Cstart, CZZdir, 1.4, true);
	updateA(C_AC, 'H', Hdir);
    // this new C atom is irrelevant. Next add a R6 on the resulting AC (from RAC)
    Cpointer C1_new, C2_new; // save all new C atoms
	if(b4) Cstart = C_AC->C2->C2;
	else Cstart = C_1;
	if(!bridge){
        for(int i=0; i!=2; i++) removeC(Cstart->C2, true);
    } else {//else if there's a bridge, convert the bridge atoms to normal edge atoms
        Cpointer b1 = Cstart->C2;
        Cpointer b2 = b1->C3;
        connectToC(Cstart, b2->C2);
        b1->C3 = NULL;
        b2->C3 = NULL;
        b1->bridge = false;
        b2->bridge = false;
        b1->C1 = b2;
        b2->C2 = b1;
    }
    C1_new = addC(Cstart, starting_direction, 1.4);
	updateA(C1_new, 'H', Hdir1);
	updateA(Cstart, 'C', Hdir1);
	C2_new = addC(C1_new, FEdir, 1.4);
	updateA(C2_new, 'H', Hdir2);
	updateA(C_2, 'C', FEdir);
    // edit sites. first identify the neighbouring sites of resulting AC & FE3
    Spointer S1, S2, S3, S4;
    if(b4) {
        S1 = moveIt(sR5, -1); // neighbour of R5
        S2 = moveIt(stt, 1); // neighbour of RZZR (stt)
		S4 = moveIt(stt, 2); // neighbour of R5 (stt)
        convSiteType(sR5, C_AC, C1_new, AC); // convert R5 to AC
        remR5fromSite(S1, S1->C1, C_AC); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, C1_new, C2_new, FE); // convert RZZR to FE
		if (S2->type != R5) {
			if ((int)S2->type >=501 && (int)S2->type <=605) updateSites(S2, C2_new, S2->C2, +1501);
			else if ((int)S2->type >=602 && (int)S2->type <=604) updateSites(S2, C2_new, S2->C2, +1501);
			else if ((int)S2->type >=1002 && (int)S2->type <=1004) updateSites(S2, C2_new, S2->C2, +1101);
			else updateSites(S2, C2_new, S2->C2, +1);
		}
        else {
			updateSites(S2, C2_new, S2->C2, +401); // update resulting FE neighbour
			if ((int)S4->type < 2000) updateSites(S4, S4->C1, S4->C2, +400); // update resulting FE neighbour
		}
    } else {
        S1 = moveIt(stt, -1); // neighbour of RZZR (stt)
        S2 = moveIt(sR5, 1); // neighbour of R5
		S3 = moveIt(stt, -2); // neighbour of R5 (stt)
        convSiteType(sR5, C2_new, C_AC, AC); // convert R5 to AC
        remR5fromSite(S2, C_AC, S2->C2); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, C1_new, C2_new, FE); // convert RAC to FE
		if (S1->type != R5) {
			if ((int)S1->type >=501 && (int)S1->type <=605) updateSites(S1, S1->C1, C1_new, +1501);
			else if ((int)S1->type >=602 && (int)S1->type <=604) updateSites(S1, S1->C1, C1_new, +1501);
			else if ((int)S1->type >=1002 && (int)S1->type <=1004) updateSites(S1, S1->C1, C1_new, +1101);
			else updateSites(S1, S1->C1, C1_new, +1);
		}
        else {
			updateSites(S1, S1->C1, C1_new, +401); // update resulting FE neighbour
			if ((int)S3->type < 2000) updateSites(S3, S3->C1, S3->C2, +400); // update resulting FE neighbour
		}
    }
    S3 = moveIt(S1, -1);
    S4 = moveIt(S2, 1);
    updateCombinedSites(stt); updateCombinedSites(sR5); // new FE and AC
    updateCombinedSites(S1); updateCombinedSites(S2); 
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours
    // add ring counts
    m_pah->m_rings++;
    m_pah->m_rings5_Lone--;
}

// ************************************************************
// ID43 - R5 conversion to R6 on R5R6ZZR (2R5 collision)
// ************************************************************
void PAHProcess::proc_C5R_R5R6ZZR(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_C5R_RZZR(stt, C_1, C_2, rng);
}

// ************************************************************
// ID44 - R5R6BY5 closure reaction
// ************************************************************
void PAHProcess::proc_L6_R5R6BY5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	proc_L6_BY6(stt, C_1, C_2);
}

// ************************************************************
// ID45 - R5R6ACR closure reaction
// ************************************************************
void PAHProcess::proc_L6_R5R6ACR(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	proc_L6_BY6(stt, C_1, C_2);
}

// ************************************************************
// ID46 - R5R6ACR5R6 closure reaction
// ************************************************************
void PAHProcess::proc_L6_R5R6ACR5R6(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	proc_L6_BY6(stt, C_1, C_2);
}

// ************************************************************
// ID47 - ZZACR5 closure reaction
// ************************************************************
void PAHProcess::proc_L6_ZZACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	// check if ZZACR5 has an opposite site of the R5.
	Cpointer Ccheck = C_1->C2;
	Cpointer Ccheck2;
	if (Ccheck->bridge) Ccheck2 = Ccheck->C3;
	else Ccheck2 = Ccheck->C2;
	bool opp_site_bool = false;
	do {
		if( (isR5internal(Ccheck, Ccheck2)) && opp_site_bool==false){
			Cpointer thirdC = findThirdC(Ccheck);
			Cpointer thirdC2 = findThirdC(Ccheck2);
			if (thirdC != NULLC || thirdC2 != NULLC) opp_site_bool = true;
		}
		if (Ccheck->bridge && Ccheck2->bridge) {
			if (Ccheck->C3 == Ccheck2){
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C2;
			}
			else {
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C3;
			}
		}
		else if (Ccheck->bridge && !Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
		else if (!Ccheck->bridge && Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C3;
		}
		else{
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
	}while (Ccheck != C_2->C1);
	
	//cout << "Printing internal R5 positions:.\n";
	/*for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
	}*/
	proc_L6_BY6(stt, C_1, C_2); 
	if (opp_site_bool == false){
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	}
}

// ************************************************************
// ID48 - R5FEACR5 closure reaction
// ************************************************************
void PAHProcess::proc_L6_R5FEACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	// check if ZZACR5 has an opposite site of the R5.
	Cpointer Ccheck = C_1->C2;
	Cpointer Ccheck2;
	if (Ccheck->bridge) Ccheck2 = Ccheck->C3;
	else Ccheck2 = Ccheck->C2;
	bool opp_site_bool = false;
	do {
		if( (isR5internal(Ccheck, Ccheck2) ) && opp_site_bool==false){
			Cpointer thirdC = findThirdC(Ccheck);
			Cpointer thirdC2 = findThirdC(Ccheck2);
			if (thirdC != NULLC || thirdC2 != NULLC) opp_site_bool = true;
		}
		if (Ccheck->bridge && Ccheck2->bridge) {
			if (Ccheck->C3 == Ccheck2){
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C2;
			}
			else {
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C3;
			}
		}
		else if (Ccheck->bridge && !Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
		else if (!Ccheck->bridge && Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C3;
		}
		else{
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
	}while (Ccheck != C_2->C1);
	
	//cout << "Printing internal R5 positions:.\n";
	/*for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
	}*/
	proc_L6_BY6(stt, C_1, C_2); 
	if (opp_site_bool == false){
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	}
}

// ************************************************************
// ID49 - FEACR5FE closure reaction
// ************************************************************
void PAHProcess::proc_L6_FEACR5FE(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	// check if ZZACR5 has an opposite site of the R5.
	Cpointer Ccheck = C_1->C2;
	Cpointer Ccheck2;
	if (Ccheck->bridge) Ccheck2 = Ccheck->C3;
	else Ccheck2 = Ccheck->C2;
	bool opp_site_bool = false;
	do {
		if( (isR5internal(Ccheck, Ccheck2) ) && opp_site_bool==false){
			Cpointer thirdC = findThirdC(Ccheck);
			Cpointer thirdC2 = findThirdC(Ccheck2);
			if (thirdC != NULLC || thirdC2 != NULLC) opp_site_bool = true;
		}
		if (Ccheck->bridge && Ccheck2->bridge) {
			if (Ccheck->C3 == Ccheck2){
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C2;
			}
			else {
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C3;
			}
		}
		else if (Ccheck->bridge && !Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
		else if (!Ccheck->bridge && Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C3;
		}
		else{
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
	}while (Ccheck != C_2->C1);
	
	//cout << "Printing internal R5 positions:.\n";
	/*for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
	}*/
	proc_L6_BY6(stt, C_1, C_2); 
	if (opp_site_bool == false){
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	}
}

// ************************************************************
// ID50 - R5ACR5R5 closure reaction
// ************************************************************
void PAHProcess::proc_L6_R5ACR5R5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	// check if ZZACR5 has an opposite site of the R5.
	Cpointer Ccheck = C_1->C2;
	Cpointer Ccheck2;
	if (Ccheck->bridge) Ccheck2 = Ccheck->C3;
	else Ccheck2 = Ccheck->C2;
	bool opp_site_bool = false;
	do {
		if( (isR5internal(Ccheck, Ccheck2) ) && opp_site_bool==false){
			Cpointer thirdC = findThirdC(Ccheck);
			Cpointer thirdC2 = findThirdC(Ccheck2);
			if (thirdC != NULLC || thirdC2 != NULLC) opp_site_bool = true;
		}
		if (Ccheck->bridge && Ccheck2->bridge) {
			if (Ccheck->C3 == Ccheck2){
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C2;
			}
			else {
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C3;
			}
		}
		else if (Ccheck->bridge && !Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
		else if (!Ccheck->bridge && Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C3;
		}
		else{
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
	}while (Ccheck != C_2->C1);
	
	//cout << "Printing internal R5 positions:.\n";
	/*for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
	}*/
	proc_L6_BY6(stt, C_1, C_2); 
	if (opp_site_bool == false){
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	}
}

// ************************************************************
// ID51 - R7 bay closure on R5ZZACR5
// ************************************************************
void PAHProcess::proc_L7_R5ZZACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	// check if ZZACR5 has an opposite site of the R5.
	Cpointer Ccheck = C_1->C2;
	Cpointer Ccheck2;
	if (Ccheck->bridge) Ccheck2 = Ccheck->C3;
	else Ccheck2 = Ccheck->C2;
	bool opp_site_bool = false;
	do {
		if( (isR5internal(Ccheck, Ccheck2) ) && opp_site_bool==false){
			Cpointer thirdC = findThirdC(Ccheck);
			Cpointer thirdC2 = findThirdC(Ccheck2);
			if (thirdC != NULLC || thirdC2 != NULLC) opp_site_bool = true;
		}
		if (Ccheck->bridge && Ccheck2->bridge) {
			if (Ccheck->C3 == Ccheck2){
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C2;
			}
			else {
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C3;
			}
		}
		else if (Ccheck->bridge && !Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
		else if (!Ccheck->bridge && Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C3;
		}
		else{
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
	}while (Ccheck != C_2->C1);
	
	//cout << "Printing internal R5 positions:.\n";
	/*for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
	}*/
	proc_L6_BY6(stt, C_1, C_2); 
	if (opp_site_bool == false){
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	}
	m_pah->m_rings--; m_pah->m_rings7_Embedded++;
	
}

// ************************************************************
// ID52 - ACR5R5R6 closure reaction
// ************************************************************
void PAHProcess::proc_L6_ACR5R5R6(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	// check if ZZACR5 has an opposite site of the R5.
	Cpointer Ccheck = C_1->C2;
	Cpointer Ccheck2;
	if (Ccheck->bridge) Ccheck2 = Ccheck->C3;
	else Ccheck2 = Ccheck->C2;
	int opp_site_bool = 0;
	do {
		if( (isR5internal(Ccheck, Ccheck2) ) && opp_site_bool<=1){
			Cpointer thirdC = findThirdC(Ccheck);
			Cpointer thirdC2 = findThirdC(Ccheck2);
			if (thirdC != NULLC || thirdC2 != NULLC) opp_site_bool += 1;
		}
		if (Ccheck->bridge && Ccheck2->bridge) {
			if (Ccheck->C3 == Ccheck2){
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C2;
			}
			else {
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C3;
			}
		}
		else if (Ccheck->bridge && !Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
		else if (!Ccheck->bridge && Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C3;
		}
		else{
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
	}while (Ccheck != C_2->C1);
	
	//cout << "Printing internal R5 positions:.\n";
	/*for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
	}*/
	proc_L6_BY6(stt, C_1, C_2); 
	if (opp_site_bool == 0){
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	}
	else if (opp_site_bool == 1) {
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	}
}

// ************************************************************
// ID53 - R7 bay closure on ACR5R5R6ZZ
// ************************************************************
void PAHProcess::proc_L7_ACR5R5R6ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if (!m_pah->m_optimised){
		OpenBabel::OBMol mol = passPAH();
		mol = mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	// check if ZZACR5 has an opposite site of the R5.
	Cpointer Ccheck = C_1->C2;
	Cpointer Ccheck2;
	if (Ccheck->bridge) Ccheck2 = Ccheck->C3;
	else Ccheck2 = Ccheck->C2;
	int opp_site_bool = 0;
	do {
		if( (isR5internal(Ccheck, Ccheck2) ) && opp_site_bool<=1){
			Cpointer thirdC = findThirdC(Ccheck);
			Cpointer thirdC2 = findThirdC(Ccheck2);
			if (thirdC != NULLC || thirdC2 != NULLC) opp_site_bool += 1;
		}
		if (Ccheck->bridge && Ccheck2->bridge) {
			if (Ccheck->C3 == Ccheck2){
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C2;
			}
			else {
				Ccheck = Ccheck2;
				Ccheck2 = Ccheck2->C3;
			}
		}
		else if (Ccheck->bridge && !Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
		else if (!Ccheck->bridge && Ccheck2->bridge) {
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C3;
		}
		else{
			Ccheck = Ccheck2;
			Ccheck2 = Ccheck2->C2;
		}
	}while (Ccheck != C_2->C1);
	
	//cout << "Printing internal R5 positions:.\n";
	/*for (std::list<cpair>::iterator it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		cout << std::get<0>(*it1) << ", " << std::get<1>(*it1) << ", " << std::get<2>(*it1) <<"\n";
	}*/
	proc_L6_BY6(stt, C_1, C_2); 
	if (opp_site_bool == 0){
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	}
	else if (opp_site_bool == 1) {
		m_pah->m_rings5_Lone--; m_pah->m_rings5_Embedded++;
	}
	m_pah->m_rings--; m_pah->m_rings7_Embedded++;
}

// ************************************************************
// ID54 - CH3 addition
// ************************************************************
void PAHProcess::proc_A_CH3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	Cpointer chosen;
    bool before; // true if C1 of site is chosen, false if C2
	if(stt->C1->A != 'H' && stt->C2->A != 'H') return;
	if(stt->C1->A == 'H' && stt->C2->A != 'H') {
		before = true;
		chosen = C_1;
	}
	else if (stt->C1->A != 'H' && stt->C2->A == 'H') {
		before = false;
		chosen = C_2;
	}
	else if ((int)stt->type>=100){
		Spointer S1 = moveIt(stt, -1);
		if ((int)S1->type>=100){
			before = false;
			chosen = C_2;
		}
		else {
			before = true;
			chosen = C_1;
		}
	}
	else {
		// choose one of the C atoms if site type is FE/AC/ZZ
		// Define a distribution that has two equally probably outcomes
		boost::bernoulli_distribution<> choiceDistrib;
		// Now build an object that will generate a sample using rng
		boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
		if(choiceGenerator()) {
			chosen = C_1;
			Spointer s_adjacent = moveIt(stt, -1);
			if ( (int)s_adjacent->type %10 >=3 ) return;
			before=true;
		}
		else {
			chosen = C_2;
			Spointer s_adjacent = moveIt(stt, +1);
			if ( (int)s_adjacent->type %10 >=3 ) return;
			before=false;
		}
    }
    // check hindrance
    if(checkHindrance_newCposition(chosen)) return;
    // add C atoms
	updateA(chosen, 'M', chosen->growth_vector);
	// neighbouring site to be updated:
    Spointer neighbour, neighbour2, prev;
    if(before) {
		prev = moveIt(stt,1);
        neighbour = moveIt(stt,-1);
		neighbour2 = moveIt(stt,-2);
    }
    else {
		prev = moveIt(stt,-1);
        neighbour = moveIt(stt,1);
		neighbour2 = moveIt(stt,2);
    }
	m_pah->m_methyl_counts += 1;

    // update combined sites for all new sites and neighbours (and their neighbours)
    updateCombinedSites(stt); updateCombinedSites(prev); updateCombinedSites(neighbour); updateCombinedSites(neighbour2);
}

// ************************************************************
// ID55 - CH3 desorption
// ************************************************************
void PAHProcess::proc_D_CH3(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	if(stt->C1->A != 'M') {
		cout << "Error. No methyl group in site for CH3 desorption.\n";
		return;
	}
	updateA(stt->C1, 'H', stt->C1->growth_vector);
	updateSites(stt,C_1,C_2,0);
    Spointer neighbour, neighbour2, prev, prev2;
	prev2 = moveIt(stt,-2);
	prev = moveIt(stt,-1);
	neighbour = moveIt(stt,1);
	neighbour2 = moveIt(stt,2);
	addCount(-1, -2);
	m_pah->m_methyl_counts -= 1;

    // update combined sites for all new sites and neighbours (and their neighbours)
    updateCombinedSites(stt); updateCombinedSites(prev); updateCombinedSites(prev2); updateCombinedSites(neighbour); updateCombinedSites(neighbour2);
}

// ************************************************************
// ID56 - Oxidation of R5R6 site
// ************************************************************
void PAHProcess::proc_O5R_R5R6(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	//printStruct();
	for (Cpointer Ccheck = C_1->C1; Ccheck != C_2; Ccheck = Ccheck->C2){
		if (getDistance_twoC(Ccheck, Ccheck->C2)>1.58 || getDistance_twoC(Ccheck, Ccheck->C2)<1.2 ){
			if (!m_pah->m_optimised){
				OpenBabel::OBMol mol = passPAH();
				mol = mol = optimisePAH(mol);
				passbackPAH(mol);
			}
		}
	}
	//saveXYZ("KMC_DEBUG/Oxidation_PAH");
	//First check if R6 is to the left or the right of R5
	bool b4 = false;
	Spointer other;
	Cpointer CRem, CRem_before, CRem_next;
	//The pentagons are at one or both ends of the PAH. First identify which atom to remove.
	if ((int)stt->type > 1000 && (int)stt->type < 1005 ){
		//Site with embedded R5 to both sides.
		// Define a distribution that has two equally probably outcomes
		boost::bernoulli_distribution<> choiceDistrib;
		// Now build an object that will generate a sample using rng
		boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
		if(choiceGenerator()) {
			b4 = true;
			other = moveIt(stt, +1);
			CRem = stt->C2;
		}
		else {
			other = moveIt(stt, -1);
			CRem = stt->C1;
		}
	}
	else if ((int)stt->type > 600 && (int)stt->type < 605){
		//Partially embedded R5 to one side and R5 to the other side.
		Spointer next_site = moveIt(stt, -1);
		if(next_site->type == R5){
			b4 = true;
			other = moveIt(stt, +1);
			CRem = stt->C2;
		}
		else{
			other = moveIt(stt, -1);
			CRem = stt->C1;
		}
	}
	else {
		//Partially embedded R5 at one end or the other.
		if ( isR5internal(C_1->C1, C_1) ) {
			other = moveIt(stt, -1);
			CRem = stt->C1;
		}
		else if ( isR5internal(C_2,C_2->C2) ) {
			b4 = true;
			other = moveIt(stt, +1);
			CRem = stt->C2;
		}
		else if((int)moveIt(stt,-1)->type == 501){
			other = moveIt(stt, -1);
			CRem = stt->C1;
		}
		else if((int)moveIt(stt,+1)->type == 501){
			b4 = true;
			other = moveIt(stt, +1);
			CRem = stt->C2;
		}
		else{
			//Embedded R5R6 not found. Flag error.
			std::cout << "Embedded R5 not found on R5R6 oxidation.\n";
			return;
		}
	}

	CRem_next = CRem->C2;
	CRem_before = CRem->C1;
	Cpointer thirdC = findThirdC(CRem_before);
	Cpointer thirdC2 = findThirdC(CRem_next);
	Spointer other_side;
	bool bridged = false;
	cpair Cdir1, Cdir2, normdir, Cdir0;
	if (thirdC != NULLC && thirdC2 != NULLC) {
		bridged = true;
		other_side = findSite(thirdC);
	}
	//Get normal vector
	normdir = norm_vector(CRem_before->coords, CRem->coords, CRem_next->coords);
	//Remove carbon first
	removeC(CRem, false);
	removeR5internal(CRem_before,CRem_next);
	if (bridged == false) {
		Cdir1 = cross_vector(get_vector(CRem_before->coords, CRem_next->coords), normdir);
		Cdir2 = get_vector(CRem_before->coords, CRem_next->coords);
		double Cdist = getDistance_twoC(CRem_before, CRem_next);
		//Add carbons from internal positions
		Cpointer newC = addC(CRem_before, Cdir1, Cdist/2.40*1.47, true);
		newC = addC(newC, Cdir2, 1.45, true);
	}
	else{
		thirdC->C3 = thirdC2;
		thirdC2->C3 = thirdC;
		thirdC2->C2 = CRem_next;
		thirdC->C1 = CRem_before;
		CRem_before->C2 = thirdC;
		CRem_next->C1 = thirdC2;
		thirdC2->bridge = true;
		thirdC->bridge = true;
	}
	//Add Hydrogen
	Cdir1 = get_vector(CRem_before->coords, CRem_before->C2->coords);
	Cdir0  = get_vector(CRem_before->C1->coords, CRem_before->coords);
	cpair Hdir1 = cross_vector(add_vector(Cdir0, Cdir1),invert_vector(normdir));
	updateA(CRem_before, 'H', Hdir1);
	Cdir1 = get_vector(CRem_next->coords, CRem_next->C2->coords);
	Cdir0  = get_vector(CRem_next->C1->coords, CRem_next->coords);
	Hdir1 = cross_vector(add_vector(Cdir0, Cdir1), invert_vector(normdir));
	updateA(CRem_next, 'H', Hdir1);
	//New site appearing
	Spointer newSite;
	if (b4){
		if ((int)stt->type < 2000) updateSites(stt, stt->C1, CRem_before, -501);
		else updateSites(stt, stt->C1, CRem_before, -101);
		if ((int)other->type < 2000) updateSites(other, CRem_next, other->C2, -501);
		else updateSites(other, CRem_next, other->C2, -101);
		newSite = addSite(AC, CRem_before, CRem_next, other);
	}
	else{
		if ((int)stt->type < 2000) updateSites(stt, CRem_next, stt->C2, -501);
		else updateSites(stt, CRem_next, stt->C2, -101);
		if ((int)other->type < 2000) updateSites(other, other->C1, CRem_before, -501);
		else updateSites(other, other->C1, CRem_before, -101);
		newSite = addSite(AC, CRem_before, CRem_next, stt);
	}
	
	if (bridged){
		if (other_side != m_pah->m_siteList.end()){
			updateSites(other_side, other_side->C1, other_side->C2, -2000);
			Spointer S1_os = moveIt(other_side, -1); Spointer S2_os = moveIt(other_side, +1);
			updateCombinedSites(S1_os); updateCombinedSites(S2_os);
		}
	}
	
	Spointer S1 = moveIt(newSite, -1); Spointer S3 = moveIt(newSite, -2); Spointer S5 = moveIt(newSite, -3);
	Spointer S2 = moveIt(newSite, +1); Spointer S4 = moveIt(newSite, +2); Spointer S6 = moveIt(newSite, +3);
    // update combined sites for all new sites and neighbours (and their neighbours)
    updateCombinedSites(stt); updateCombinedSites(other); updateCombinedSites(newSite); 
	updateCombinedSites(S1); updateCombinedSites(S2); updateCombinedSites(S3); updateCombinedSites(S4); updateCombinedSites(S5); updateCombinedSites(S6);
	m_pah->m_rings5_Lone--;
	addCount(0,+1);
}

// ************************************************************
// ID57 - Oxidation of R5R6ZZ site
// ************************************************************
void PAHProcess::proc_O5R_R5R6ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_O5R_R5R6(stt, C_1, C_2, rng);
}

// ************************************************************
// ID58 - Oxidation of R5R6AC site
// ************************************************************
void PAHProcess::proc_O5R_R5R6AC(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_O5R_R5R6(stt, C_1, C_2, rng);
}

// ************************************************************
// ID59 - Oxidation of R5R6BY5 site
// ************************************************************
void PAHProcess::proc_O5R_R5R6BY5(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_O5R_R5R6(stt, C_1, C_2, rng);
}

// ************************************************************
// ID60 - Oxidation of R5R6FER site
// ************************************************************
void PAHProcess::proc_O5R_R5R6FER(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_O5R_R5R6(stt, C_1, C_2, rng);
}

// ************************************************************
// ID61 - Oxidation of R5R6ZZR site
// ************************************************************
void PAHProcess::proc_O5R_R5R6ZZR(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_O5R_R5R6(stt, C_1, C_2, rng);
}

// ************************************************************
// ID62 - Oxidation of R5R6ACR site
// ************************************************************
void PAHProcess::proc_O5R_R5R6ACR(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_O5R_R5R6(stt, C_1, C_2, rng);
}

// ************************************************************
// ID63 - Oxidation of R5R6FER5R6 site
// ************************************************************
void PAHProcess::proc_O5R_R5R6FER5R6(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_O5R_R5R6(stt, C_1, C_2, rng);
}

// ************************************************************
// ID64 - Oxidation of R5R6ZZR5R6 site
// ************************************************************
void PAHProcess::proc_O5R_R5R6ZZR5R6(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_O5R_R5R6(stt, C_1, C_2, rng);
}

// ************************************************************
// ID65 - Oxidation of R5R6BY5 site
// ************************************************************
void PAHProcess::proc_O5R_R5R6ACR5R6(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	proc_O5R_R5R6(stt, C_1, C_2, rng);
}

// ************************************************************
// ID66- Five member ring migration around the corner
// ************************************************************
void PAHProcess::proc_M5R_ACR5_around_corner(Spointer& stt, Cpointer C_1, Cpointer C_2, Spointer& sFE2, bool b4, int migr_index) {
	//This process is called when the R5 walker moves to another edge. 
	//It heals the vacancy present on the current edge and creates a vacancy on the next edge.
	//It also resets the walker position to the new edge.
	//It assumes the starting site is an ACR5 site or similar.
	// The pentagon migrated N times and ended at the same position.
	int ii = migr_index;
	if(stt==sFE2){
		//Walker came back to starting position from other side without resetting.
		if ( ((int)sFE2->type >= 1 && (int)sFE2->type <= 4) || ((int)sFE2->type >= 102 && (int)sFE2->type <= 104) ){
			//This means that the pentagon has migrated to a basic site.
			if (b4) {
				updateSites(sFE2, sFE2->C1, sFE2->C2, 2000 + (int)sFE2->type);
				if ((int)sFE2->type<2000) {
					Spointer S4 = moveIt(sFE2, +1);
					updateSites(S4, S4->C1, S4->C2,+500);
				}
			}
			else {
				updateSites(sFE2, sFE2->C1, sFE2->C2, 2000 + (int)sFE2->type);
				if ((int)sFE2->type<2000) {
					Spointer S3 = moveIt(sFE2, -1);
					updateSites(S3, S3->C1, S3->C2,+500);
				}
			}
		} else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){
			//This means that the pentagon has migrated to neighbour an edge R5 edge.
			Spointer S1_left = moveIt(sFE2, -1);
			Spointer S2_right = moveIt(sFE2, +1);
			if (b4) {
				updateSites(sFE2, sFE2->C1, sFE2->C2, 1600 + (int)sFE2->type);
				if ((int)S2_right->type<2000) {
					Spointer S4 = moveIt(S2_right, +1);
					updateSites(S4, S4->C1, S4->C2,+500);
				}
			}
			else {
				updateSites(sFE2, sFE2->C1, sFE2->C2, 1600 + (int)sFE2->type);
				if ((int)sFE2->type<2000) {
					Spointer S3 = moveIt(S1_left, -1);
					updateSites(S3, S3->C1, S3->C2,+500);
				}
			}
		} else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){
			//This means that the pentagon has migrated to a bay containing an R5.
			if (b4) {
				updateSites(sFE2, sFE2->C1, sFE2->C2, 100 + (int)sFE2->type);
				if ((int)sFE2->type<2000) {
					Spointer S4 = moveIt(sFE2, +1);
					updateSites(S4, S4->C1, S4->C2,+500);
				}
			}
			else {
				updateSites(sFE2, sFE2->C1, sFE2->C2, 100 + (int)sFE2->type);
				if ((int)sFE2->type<2000) {
					Spointer S3 = moveIt(sFE2, -1);
					updateSites(S3, S3->C1, S3->C2,+500);
				}
			}
		}
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0; // Walker came back to pos 0
	
		//Update combined sites
		Spointer S1, S2, S3, S4, S5, S6;
		S1 = moveIt(sFE2, -1);
		S2 = moveIt(sFE2, +1);
		S3 = moveIt(S1, -1);
		S4 = moveIt(S2, +1);
		S5 = moveIt(S1, -2);
		S6 = moveIt(S2, +2);
		updateCombinedSitesMigration(S1); updateCombinedSitesMigration(S2);
		updateCombinedSitesMigration(S3); updateCombinedSitesMigration(S4); updateCombinedSitesMigration(S5); updateCombinedSitesMigration(S6);
	}

	int steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
	bool dir;
	if (steps < 0) dir = true;
	else dir = false;

	//Remove R5coords from m_pah->m_R5loc. This is done by starter function/
	//findR5internal(C_1->C2, C_2->C1);
	// First select carbons and sites affected.
	Cpointer CFE, CRem, CRem_next, CRem_before;
	Spointer checkR5_1, checkR5_2;
	if (b4) {
		CFE = C_1->C2;
		checkR5_1 = moveIt(sFE2, -1);
		checkR5_2 = moveIt(sFE2, -2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C1;
		else CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
	}
	else {
		CFE = C_2->C1->C1;
		checkR5_1 = moveIt(sFE2, +1);
		checkR5_2 = moveIt(sFE2, +2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C2;
		else CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
	}
	//int end_site_type = (int)sFE2->type;
	
	//Add a new carbon between current R5 carbons of ACR5
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2;
	if (R5_dist < 2.5) dist2 = 1.4;
	else dist2 = R5_dist / 2.7 * 1.5;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir = get_vector(C_1->C2->coords,C_2->C1->coords);
	cpair normvec = (norm_vector(C_1->C2->coords, C_1->C2->C2->coords, C_1->C2->C2->C2->coords));
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, 1.4);
	updateA(Cnew, 'H', crossvec);
	
	//Remove carbon from end site 
	removeC(CRem, false);
	/*OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol, 500);
	passbackPAH(mol);*/
	
	//First adjust starting site and add new site if needed.
	Spointer stt_coupled, newSite;
	if (b4) stt_coupled = moveIt(stt,-1);
	else stt_coupled = moveIt(stt,+1);
	if (b4) {
		//convSiteType(stt, Cnew, stt->C2, ZZ);
		updateSites(stt, Cnew, stt->C2, 0);
		newSite = addSite(ZZ, Cnew->C1->C1, Cnew, stt);
	} else{
		//convSiteType(stt, stt->C1, Cnew, ZZ);
		updateSites(stt, stt->C1, Cnew, 0);
		newSite = addSite(ZZ, Cnew, Cnew->C2->C2, stt_coupled);
	}
	
	Spointer S1 = moveIt(sFE2, -1);
	Spointer S2 = moveIt(sFE2, +1);
	
	//Adjust sites for ending position.
	//Migration outside.
	if ((int)checkR5_1->type == 0 && (int)sFE2->type == 0){
		if (b4) {
			Spointer sFE2_right = moveIt(sFE2, +1);
			updateSites(sFE2_right, CRem_before, sFE2_right->C2, +100);
			convSiteType(checkR5_1, CRem_next, CRem_before, R5);
			updateSites(checkR5_2, checkR5_2->C1, CRem_next, +100);
		}
		else {
			Spointer sFE2_left = moveIt(sFE2, -1);
			updateSites(sFE2_left, sFE2_left->C1, CRem_before, +100);
			convSiteType(checkR5_1, CRem_before, CRem_next, R5);
			updateSites(checkR5_2, CRem_next, checkR5_2->C2, +100);
		}
		removeSite(sFE2);
		std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type == 0){
		//This means that the pentagon has migrated to the edge but will have a single carbon out of the structure.
		if (b4) {
			updateSites(S1, S1->C1, sFE2->C1, 500);
			updateSites(S2, sFE2->C1, S2->C2, 500);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 500);
			updateSites(S2, sFE2->C2, S2->C2, 500);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2);
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type >= 1 && (int)sFE2->type <= 4){
		//This means that the pentagon has migrated to a basic site.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 2000 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
			if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 2000 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S1;
			if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2); 
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 2000 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
			if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 2000 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S1;
			if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2);
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		Spointer S1 = moveIt(sFE2, -1);
		Spointer S2 = moveIt(sFE2, +1);
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 1600 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
			if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 1600 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S1;
			if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2);
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){
		//This means that the pentagon has migrated to a bay containing an R5.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 100 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
			if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 100 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S1;
			if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2);
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}

	Spointer finsite = std::get<0>(m_pah->m_R5walker_sites[ii]);
	if((int)finsite->type>=2003) {
		if(dir) addR5internal(finsite->C2->C1->C1,finsite->C2->C1,true);
		else addR5internal(finsite->C1->C2,finsite->C1->C2->C2,true);
	}
	
	//Update combined sites
	Spointer stt_1, stt_2, end_site, origin_site, S3, S4, S5, S6;
	if (b4){
		stt_1 = moveIt(newSite, +1);
		stt_2 = moveIt(stt, +1);
		end_site = S2;
		origin_site = moveIt(end_site,+1);
		S1 = moveIt(end_site,-1);
		S2 = moveIt(origin_site,1);
	}else{
		stt_1 = moveIt(stt, -1);
		stt_2 = moveIt(newSite, -1);
		end_site = S1;
		origin_site = moveIt(end_site,-1);
		S1 = moveIt(origin_site,-1);
		S2 = moveIt(end_site,1);
	}
	S3 = moveIt(S1, -1);
	S4 = moveIt(S2, +1);
	S5 = moveIt(S1, -2);
	S6 = moveIt(S2, +2);
	updateCombinedSitesMigration(stt_1); updateCombinedSitesMigration(stt); updateCombinedSitesMigration(newSite); updateCombinedSitesMigration(stt_2);
	updateCombinedSitesMigration(origin_site); updateCombinedSitesMigration(end_site);
	updateCombinedSitesMigration(S1); updateCombinedSitesMigration(S2);
	updateCombinedSitesMigration(S3); updateCombinedSitesMigration(S4); updateCombinedSitesMigration(S5); updateCombinedSitesMigration(S6);
}

// ************************************************************
// ID67- Recursive partially embedded 5-member ring migration
// ************************************************************
void PAHProcess::proc_M5R_R5R6_out_of_corner(Spointer& stt, Cpointer C_1, Cpointer C_2, Spointer& sFE2, bool b4, int migr_index) {
	//This process is called when the R5 walker moves out of the corner into an edge. 
	//It heals the vacancy present on the current corner and creates a vacancy on the next edge.
	//It also resets the walker position to the new edge.
	//It assumes the starting site is an R5R6 site or similar.
	// First select carbons and sites affected.
	int ii = migr_index;
	Cpointer CFE, CRem, CRem_next, CRem_before;
	Spointer checkR5_1, checkR5_2;
	if (b4) {
		CFE = C_2->C1;
		checkR5_1 = moveIt(sFE2, -1);
		checkR5_2 = moveIt(sFE2, -2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C1;
		else CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
	}
	else {
		CFE = C_1;
		checkR5_1 = moveIt(sFE2, +1);
		checkR5_2 = moveIt(sFE2, +2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C2;
		else CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
	}
	//int end_site_type = (int)sFE2->type;
	//Add a new carbon between current R5 carbons of ACR5
	Cpointer Cnew;
	if (b4) {
		cpair R5dirFE = get_vector(CFE->coords,CFE->C2->C2->coords);
		double R5dirFEdist = getDistance_twoC(CFE->coords,CFE->C2->C2->coords);
		double bond_dist = R5dirFEdist / 2.0;
		double magnFE = bond_dist*sin(M_PI / 3.0);
		cpair normvec = invert_vector(norm_vector(CFE->coords,CFE->C2->coords,CFE->C2->C2->coords));
		cpair crossvecFE = cross_vector(R5dirFE, normvec);
		cpair resultantvecFE = std::make_tuple(R5dirFEdist/4.0 * std::get<0>(R5dirFE) + magnFE * std::get<0>(crossvecFE),
											R5dirFEdist/4.0 * std::get<1>(R5dirFE) + magnFE * std::get<1>(crossvecFE), 
											R5dirFEdist/4.0 * std::get<2>(R5dirFE)+ magnFE * std::get<2>(crossvecFE));
		cpair CnewdirFE = scale_vector(resultantvecFE);
		cpair mposFE = jumpToPos(CFE->coords,CnewdirFE,bond_dist);
		mposFE = jumpToPos(mposFE,R5dirFE,bond_dist);
		moveC(CFE->C2,mposFE);
		Cnew = addC(CFE,CnewdirFE,bond_dist);

		double dotprod = std::get<0>(CnewdirFE)*std::get<0>(R5dirFE) + std::get<1>(CnewdirFE)*std::get<1>(R5dirFE) + std::get<2>(CnewdirFE)*std::get<2>(R5dirFE);
		cpair Hdir = std::make_tuple(std::get<0>(CnewdirFE) - 2.0*dotprod*std::get<0>(R5dirFE),
									std::get<1>(CnewdirFE) - 2.0*dotprod*std::get<1>(R5dirFE),
									std::get<2>(CnewdirFE) - 2.0*dotprod*std::get<2>(R5dirFE));
		updateA(Cnew,'H',Hdir);
		updateA(Cnew->C2,'H',CnewdirFE);
	}
	else {
		cpair R5dirFE = get_vector(CFE->C1->coords,CFE->C2->coords);
		double R5dirFEdist = getDistance_twoC(CFE->C1->coords,CFE->C2->coords);
		double bond_dist = R5dirFEdist / 2.0;
		double magnFE = bond_dist*sin(M_PI / 3.0);
		cpair normvec = invert_vector(norm_vector(CFE->C1->coords,CFE->coords,CFE->C2->coords));
		cpair crossvecFE = cross_vector(R5dirFE, normvec);
		cpair resultantvecFE = std::make_tuple(R5dirFEdist/4.0 * std::get<0>(R5dirFE) + magnFE * std::get<0>(crossvecFE),
											R5dirFEdist/4.0 * std::get<1>(R5dirFE) + magnFE * std::get<1>(crossvecFE), 
											R5dirFEdist/4.0 * std::get<2>(R5dirFE)+ magnFE * std::get<2>(crossvecFE));
		cpair CnewdirFE = scale_vector(resultantvecFE);
		cpair mposFE = jumpToPos(CFE->C1->coords,CnewdirFE,bond_dist);
		moveC(CFE,mposFE);
		double dotprod = std::get<0>(CnewdirFE)*std::get<0>(R5dirFE) + std::get<1>(CnewdirFE)*std::get<1>(R5dirFE) + std::get<2>(CnewdirFE)*std::get<2>(R5dirFE);
		cpair Hdir = std::make_tuple(std::get<0>(CnewdirFE) - 2.0*dotprod*std::get<0>(R5dirFE),
									std::get<1>(CnewdirFE) - 2.0*dotprod*std::get<1>(R5dirFE),
									std::get<2>(CnewdirFE) - 2.0*dotprod*std::get<2>(R5dirFE));
		updateA(CFE,'H',Hdir);
		Cnew = addC(CFE,R5dirFE,bond_dist);
		updateA(Cnew,'H',CnewdirFE);
	}
	
	//Remove carbon from end site 
	removeC(CRem, false);
	
	//First adjust starting site and add new site if needed.
	Spointer stt_coupled, newSite;
	if (b4) stt_coupled = moveIt(stt,+1);
	else stt_coupled = moveIt(stt,-1);
	if (b4) {
		convSiteType(stt, Cnew->C1->C1, Cnew, ZZ);
		newSite = addSite(FE, Cnew, Cnew->C2, stt_coupled);
	} else{
		convSiteType(stt, Cnew, Cnew->C2->C2, ZZ);
		newSite = addSite(FE, Cnew->C1, Cnew, stt);
	}
	/*if (b4) {
		//convSiteType(stt, Cnew->C1->C1, Cnew, ZZ);
		sFE2->C1 = sFE2->C1;
		sFE2->C2 = Cnew;
		//newSite = addSite(R5R6, Cnew, stt_coupled->C1, stt_coupled);
		stt->C1 = Cnew;
		stt->C2 = stt_coupled->C1;
	} else{
		//convSiteType(stt, Cnew, Cnew->C2->C2, ZZ);
		sFE2->C1 = Cnew;
		sFE2->C2 = sFE2->C2;
		//newSite = addSite(R5R6, stt_coupled->C2, Cnew, stt);
		stt->C1 = stt_coupled->C2;
		stt->C2 = Cnew;

	}*/
	//Adjust sites for ending position.
	Spointer S1 = moveIt(sFE2, -1);
	Spointer S2 = moveIt(sFE2, +1);
	//Migration outside.
	if ((int)checkR5_1->type == 0 && (int)sFE2->type == 0){
		if (b4) {
			Spointer sFE2_right = moveIt(sFE2, +1);
			updateSites(sFE2_right, CRem_before, sFE2_right->C2, +100);
			convSiteType(checkR5_1, CRem_next, CRem_before, R5);
			updateSites(checkR5_2, checkR5_2->C1, CRem_next, +100);
		}
		else {
			Spointer sFE2_left = moveIt(sFE2, -1);
			updateSites(sFE2_left, sFE2_left->C1, CRem_before, +100);
			convSiteType(checkR5_1, CRem_before, CRem_next, R5);
			updateSites(checkR5_2, CRem_next, checkR5_2->C2, +100);
		}
		removeSite(sFE2);
		std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
		std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type == 0){
		//This means that the pentagon has migrated to the edge but will have a single carbon out of the structure.
		if (b4) {
			updateSites(S1, S1->C1, sFE2->C1, 500);
			updateSites(S2, sFE2->C1, S2->C2, 500);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 500);
			updateSites(S2, sFE2->C2, S2->C2, 500);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2);
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type >= 1 && (int)sFE2->type <= 4){
		//This means that the pentagon has migrated to a basic site.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 2000 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
			if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 2000 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S1;
			if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2); 
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 2000 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
			if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 2000 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S1;
			if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2);
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		Spointer S1 = moveIt(sFE2, -1);
		Spointer S2 = moveIt(sFE2, +1);
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 1600 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
			if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 1600 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S1;
			if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2);
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}
	else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){
		//This means that the pentagon has migrated to a bay containing an R5.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 100 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S2;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S2;
			if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 100 + (int)sFE2->type);
			std::get<0>(m_pah->m_R5walker_sites[ii]) = S1;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = S1;
			if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}
		}
		checkRemR5Walkers(ii, b4, sFE2);
		removeSite(sFE2);
		std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
	}

	Spointer finsite = std::get<0>(m_pah->m_R5walker_sites[ii]);
	if((int)finsite->type>=2003) {
		if(b4) addR5internal(finsite->C2->C1->C1,finsite->C2->C1,true);
		else addR5internal(finsite->C1->C2,finsite->C1->C2->C2,true);
	}

	//Update combined sites
	Spointer S3, S4, S5, S6;
	if (b4){
		S1 = moveIt(stt,-1);
		S2 = moveIt(newSite,+1);
		S3 = moveIt(S1, -1);
		S4 = moveIt(S2, 1);
		S5 = moveIt(S1, -2);
		S6 = moveIt(S2, 2);

		/*stt_1 = moveIt(newSite, +1);
		stt_2 = moveIt(stt, +1);
		stt_3 = moveIt(newSite, +2);
		stt_4 = moveIt(newSite, +3);*/
	}else{
		S1 = moveIt(newSite,-1);
		S2 = moveIt(stt,+1);
		S3 = moveIt(S1, -1);
		S4 = moveIt(S2, 1);
		S5 = moveIt(S1, -2);
		S6 = moveIt(S2, 2);
		/*stt_1 = moveIt(stt, -1);
		stt_2 = moveIt(newSite, -1);
		stt_3 = moveIt(newSite, -2);
		stt_4 = moveIt(newSite, -3);*/
	}
	//S3 = moveIt(S1, -1);
	//S4 = moveIt(S2, +1);
	//updateCombinedSitesMigration(stt);
	if (b4){
		updateCombinedSitesMigration(newSite);
		updateCombinedSitesMigration(stt);
		updateCombinedSitesMigration(S1);
		updateCombinedSitesMigration(S2);
		updateCombinedSitesMigration(S3);
		updateCombinedSitesMigration(S4);
		updateCombinedSitesMigration(S5);
		updateCombinedSitesMigration(S6);
	}
	else{
		updateCombinedSitesMigration(newSite);
		updateCombinedSitesMigration(stt);
		updateCombinedSitesMigration(S1);
		updateCombinedSitesMigration(S2);
		updateCombinedSitesMigration(S3);
		updateCombinedSitesMigration(S4);
		updateCombinedSitesMigration(S5);
		updateCombinedSitesMigration(S6);
	}
}

// ************************************************************
// ID66- Termination of ACR5 migration
// ************************************************************
void PAHProcess::proc_M5R_ACR5_termination(Spointer& stt, Cpointer C_1, Cpointer C_2, Spointer& sFE2, bool b4) {
	// The pentagon migrated N times and ended at the same position.
	if (sFE2 == stt) return;
	//Remove R5coords from m_pah->m_R5loc. This is done by starter function/
	//findR5internal(C_1->C2, C_2->C1);
	// First select carbons and sites affected.
	Cpointer CFE, CRem, CRem_next, CRem_before, CR5_otherside_1, CR5_otherside_2;
	Spointer checkR5_1, checkR5_2;
	if (b4) {
		//CFE = C_2->C1->C1;
		CFE = C_1->C2;
		//CR5_otherside_1 = C_2->C1;
		CR5_otherside_1 = C_1->C2->C2;
		//CR5_otherside_2 = C_2->C1->C1;
		CR5_otherside_2 = C_1->C2;
		checkR5_1 = moveIt(sFE2, -1);
		checkR5_2 = moveIt(sFE2, -2);
		CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
	}
	else {
		//CFE = C_1->C2;
		CFE = C_2->C1->C1;
		//CR5_otherside_1 = C_1->C2;
		CR5_otherside_1 = C_2->C1->C1;
		//CR5_otherside_2 = C_1->C2->C2;
		CR5_otherside_2 = C_2->C1;
		checkR5_1 = moveIt(sFE2, +1);
		checkR5_2 = moveIt(sFE2, +2);
		CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
	}
	//int end_site_type = (int)sFE2->type;

	//Add a new carbon between current R5 carbons of ACR5
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2;
	if (R5_dist < 2.5) dist2 = 1.4;
	else dist2 = R5_dist / 2.7 * 1.5;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir = get_vector(C_1->C2->coords,C_2->C1->coords);
	cpair normvec = (norm_vector(C_1->C2->coords, C_1->C2->C2->coords, C_1->C2->C2->C2->coords));
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, 1.4);
	updateA(Cnew, 'H', crossvec);
	
	//Remove carbon from end site 
	removeC(CRem, false);
	
	//First adjust starting site and add new site if needed.
	Spointer stt_coupled, newSite;
	if (b4) stt_coupled = moveIt(stt,-1);
	else stt_coupled = moveIt(stt,+1);
	if (b4) {
		//convSiteType(stt, Cnew, stt->C2, ZZ);
		updateSites(stt, Cnew, stt->C2, 0);
		newSite = addSite(ZZ, Cnew->C1->C1, Cnew, stt);
	} else{
		//convSiteType(stt, stt->C1, Cnew, ZZ);
		updateSites(stt, stt->C1, Cnew, 0);
		newSite = addSite(ZZ, Cnew, Cnew->C2->C2, stt_coupled);
	}
	
	Spointer S1 = moveIt(sFE2, -1);
	Spointer S2 = moveIt(sFE2, +1);
	//sFE2 should already be the correct site type. It just needs to reshuffle the pointers.
	if (b4) {
		//updateSites(S1, S1->C1, sFE2->C1, 0);
		S1->C1 = S1->C1;
		S1->C2 = sFE2->C1;
		//convSiteType(sFE2, sFE2->C1, S2->C2, sFE2->type);
		sFE2->C1 = sFE2->C1;
		sFE2->C2 = S2->C2;
		removeSite(S2);
	} else{
		//convSiteType(S1, S1->C1, sFE2->C2, sFE2->type);
		sFE2->C1 = S1->C1;
		sFE2->C2 = sFE2->C2;
		removeSite(S1);
		//updateSites(S2, sFE2->C2, S2->C2, 0);
		S2->C1 = sFE2->C2;
		S2->C2 = S2->C2;
	}
	//removeSite(sFE2);
}

// ************************************************************
// ID66- Termination of ACR5 migration
// ************************************************************
void PAHProcess::proc_M5R_ACR5_termination_toR5(Spointer& stt, Cpointer C_1, Cpointer C_2, Spointer& sFE2, bool b4, int steps) {
	// The pentagon migrated N times and ended at the same position.
	if (sFE2 == stt) return;
	//Remove R5coords from m_pah->m_R5loc. This is done by starter function/
	//findR5internal(C_1->C2, C_2->C1);
	// First select carbons and sites affected.
	Cpointer CFE, CRem, CRem_next, CRem_before, CR5_otherside_1, CR5_otherside_2;
	Spointer checkR5_1, checkR5_2;
	if (b4) {
		//CFE = C_2->C1->C1;
		CFE = C_1->C2;
		//CR5_otherside_1 = C_2->C1;
		CR5_otherside_1 = C_1->C2->C2;
		//CR5_otherside_2 = C_2->C1->C1;
		CR5_otherside_2 = C_1->C2;
		checkR5_1 = moveIt(sFE2, -1);
		checkR5_2 = moveIt(sFE2, -2);
		CRem = sFE2->C1;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
	}
	else {
		//CFE = C_1->C2;
		CFE = C_2->C1->C1;
		//CR5_otherside_1 = C_1->C2;
		CR5_otherside_1 = C_2->C1->C1;
		//CR5_otherside_2 = C_1->C2->C2;
		CR5_otherside_2 = C_2->C1;
		checkR5_1 = moveIt(sFE2, +1);
		checkR5_2 = moveIt(sFE2, +2);
		CRem = sFE2->C2;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
	}
	//int end_site_type = (int)sFE2->type;

	//Add a new carbon between current R5 carbons of ACR5
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2;
	if (R5_dist < 2.5) dist2 = 1.4;
	else dist2 = R5_dist / 2.7 * 1.5;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir = get_vector(C_1->C2->coords,C_2->C1->coords);
	cpair normvec = (norm_vector(C_1->C2->coords, C_1->C2->C2->coords, C_1->C2->C2->C2->coords));
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, 1.4);
	updateA(Cnew, 'H', crossvec);
	
	//Remove carbon from end site 
	removeC(CRem, false);
	
	//First adjust starting site and add new site if needed.
	Spointer stt_coupled, newSite;
	if (steps == 0){
		if (b4) stt_coupled = moveIt(stt,+1);
		else stt_coupled = moveIt(stt,-1);
		if (b4) {
			int stype_diff = (int)stt->type - 101;
			Cpointer Cold = stt->C2;
			//convSiteType(stt, Cnew, stt->C2, ZZ);
			updateSites(stt, stt->C1,Cnew, -stype_diff);
			newSite = addSite((kmcSiteType)(1+stype_diff), Cnew, Cold, stt_coupled);
			Spointer prev_site = moveIt(newSite,+1);
			if ((int)prev_site->type>=500 && stype_diff>100) updateSites(newSite, newSite->C1, newSite->C2, +400);
		} else{
			int stype_diff = (int)stt->type - 101;
			Cpointer Cold = stt->C1;
			//convSiteType(stt, stt->C1, Cnew, ZZ);
			updateSites(stt, Cnew, stt->C2, -stype_diff);
			newSite = addSite((kmcSiteType)(1+stype_diff), Cold, Cnew, stt);
			Spointer prev_site = moveIt(newSite,-1);
			if ((int)prev_site->type>=500 && stype_diff>100) updateSites(newSite, newSite->C1, newSite->C2, +400);
		}
	} else{
		if (b4) stt_coupled = moveIt(stt,-1);
		else stt_coupled = moveIt(stt,+1);
		if (b4) {
			//convSiteType(stt, Cnew, stt->C2, ZZ);
			updateSites(stt, Cnew, stt->C2, 0);
			newSite = addSite(ZZ, Cnew->C1->C1, Cnew, stt);
		} else{
			//convSiteType(stt, stt->C1, Cnew, ZZ);
			updateSites(stt, stt->C1, Cnew, 0);
			newSite = addSite(ZZ, Cnew, Cnew->C2->C2, stt_coupled);
		}
	}
	
	
	Spointer S1 = moveIt(sFE2, -1);
	Spointer S2 = moveIt(sFE2, +1);
	//sFE2 should already be the correct site type. It just needs to reshuffle the pointers.
	if (b4) {
		S1->C1 = S1->C1;
		S1->C2 = sFE2->C2;
		removeSite(sFE2);
	} else{
		S2->C1 = sFE2->C1;
		S2->C2 = S2->C2;
		removeSite(sFE2);
	}
}

// ************************************************************
// ID66- Recursive embedded 5-member ring migration
// ************************************************************
void PAHProcess::proc_M5R_ACR5_multiple_sites(Spointer& stt, Cpointer C_1, Cpointer C_2, Spointer& sFE2, bool b4) {
	// The pentagon migrated N times and ended at the same position.
	if (sFE2 == stt) return;
	//Remove R5coords from m_pah->m_R5loc. This is done by starter function/
	//findR5internal(C_1->C2, C_2->C1);
	// First select carbons and sites affected.
	Cpointer CFE, CRem, CRem_next, CRem_before, CR5_otherside_1, CR5_otherside_2;
	Spointer checkR5_1, checkR5_2;
	if (b4) {
		CFE = C_2->C1->C1;
		CR5_otherside_1 = C_2->C1;
		CR5_otherside_2 = C_2->C1->C1;
		checkR5_1 = moveIt(sFE2, -1);
		checkR5_2 = moveIt(sFE2, -2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C1;
		else CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
	}
	else {
		CFE = C_1->C2;
		CR5_otherside_1 = C_1->C2;
		CR5_otherside_2 = C_1->C2->C2;
		checkR5_1 = moveIt(sFE2, +1);
		checkR5_2 = moveIt(sFE2, +2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C2;
		else CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
	}
	//int end_site_type = (int)sFE2->type;

	//Add a new carbon between current R5 carbons of ACR5
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2;
	if (R5_dist < 2.5) dist2 = 1.4;
	else dist2 = R5_dist / 2.7 * 1.5;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir = get_vector(C_1->C2->coords,C_2->C1->coords);
	cpair normvec = (norm_vector(C_1->C2->coords, C_1->C2->C2->coords, C_1->C2->C2->C2->coords));
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, 1.4);
	updateA(Cnew, 'H', crossvec);
	
	//Remove carbon from end site 
	removeC(CRem, false);
	OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol, 500);
	passbackPAH(mol);
	
	//First adjust starting site and add new site if needed.
	Spointer stt_coupled, newSite;
	if (b4) stt_coupled = moveIt(stt,-1);
	else stt_coupled = moveIt(stt,+1);
	if (b4) {
		convSiteType(stt, Cnew, stt->C2, ZZ);
		newSite = addSite(ZZ, Cnew->C1->C1, Cnew, stt);
	} else{
		convSiteType(stt, stt->C1, Cnew, ZZ);
		newSite = addSite(ZZ, Cnew, Cnew->C2->C2, stt_coupled);
	}
	
	Spointer S1 = moveIt(sFE2, -1);
	Spointer S2 = moveIt(sFE2, +1);
	
	//Adjust sites for ending position.
	//Migration outside.
	if ((int)checkR5_1->type == 0 && (int)sFE2->type == 0){
		if (b4) {
			Spointer sFE2_right = moveIt(sFE2, +1);
			updateSites(sFE2_right, CRem_before, sFE2_right->C2, +0);
			convSiteType(checkR5_1, CRem_next, CRem_before, R5);
			updateSites(checkR5_2, checkR5_2->C1, CRem_next, +0);
		}
		else {
			Spointer sFE2_left = moveIt(sFE2, -1);
			updateSites(sFE2_left, sFE2_left->C1, CRem_before, +0);
			convSiteType(checkR5_1, CRem_before, CRem_next, R5);
			updateSites(checkR5_2, CRem_next, checkR5_2->C2, +0);
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type == 0){
		//This means that the pentagon has migrated to the edge but will have a single carbon out of the structure.
		if (b4) {
			updateSites(S1, S1->C1, sFE2->C1, 0);
			updateSites(S2, sFE2->C1, S2->C2, 0);
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			updateSites(S2, sFE2->C2, S2->C2, 0);
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 1 && (int)sFE2->type <= 4){
		//This means that the pentagon has migrated to a basic site.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		Spointer S1 = moveIt(sFE2, -1);
		Spointer S2 = moveIt(sFE2, +1);
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){
		//This means that the pentagon has migrated to a bay containing an R5.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	
	//Update combined sites
	Spointer stt_1, stt_2, S3, S4;
	if (b4){
		stt_1 = moveIt(newSite, -1);
		stt_2 = moveIt(stt, +1);
	}else{
		stt_1 = moveIt(stt, -1);
		stt_2 = moveIt(newSite, +1);
	}
	S3 = moveIt(S1, -1);
	S4 = moveIt(S2, +1);
	updateCombinedSites(stt_1); updateCombinedSites(stt); updateCombinedSites(newSite); updateCombinedSites(stt_2);
	updateCombinedSites(S1); updateCombinedSites(S2); updateCombinedSites(S3); updateCombinedSites(S4);
}


// ************************************************************
// ID67- Recursive partially embedded 5-member ring migration
// ************************************************************
void PAHProcess::proc_M5R_R5R6_multiple_sites(Spointer& stt, Cpointer C_1, Cpointer C_2, Spointer& sFE2, bool b4) {
	// First select carbons and sites affected.
	Cpointer CFE, CRem, CRem_next, CRem_before, CR5_otherside_1, CR5_otherside_2;
	Spointer checkR5_1, checkR5_2;
	if (b4) {
		CFE = C_2->C1;
		CR5_otherside_1 = C_2->C2;
		CR5_otherside_2 = C_2->C1;
		checkR5_1 = moveIt(sFE2, -1);
		checkR5_2 = moveIt(sFE2, -2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C1;
		else CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
	}
	else {
		CFE = C_1;
		CR5_otherside_1 = C_1->C1;
		CR5_otherside_2 = C_1->C2;
		checkR5_1 = moveIt(sFE2, +1);
		checkR5_2 = moveIt(sFE2, +2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C2;
		else CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
	}

	//int end_site_type = (int)sFE2->type;
	//Add a new carbon between current R5 carbons of ACR5
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2 = 1.4;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir = get_vector(CFE->coords, CFE->C2->coords);
	cpair normvec;
	if (b4) normvec = invert_vector((norm_vector(CFE->coords, CFE->C2->coords, CFE->C2->C2->coords)));
	else normvec = (norm_vector(CFE->coords, CFE->C2->coords, CFE->C2->C2->coords));
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, 1.4);
	updateA(Cnew, 'H', crossvec);
	
	//Remove carbon from end site 
	removeC(CRem, false);
	//saveXYZ("KMC_DEBUG/before_optimisation");
	OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol, 500);
	passbackPAH(mol);
	//saveXYZ("KMC_DEBUG/after_optimisation");
	
	//First adjust starting site and add new site if needed.
	Spointer stt_coupled, newSite;
	if (b4) stt_coupled = moveIt(stt,+1);
	else stt_coupled = moveIt(stt,-1);
	if (b4) {
		convSiteType(stt, Cnew->C1->C1, Cnew, ZZ);
		newSite = addSite(FE, Cnew, stt_coupled->C1, stt_coupled);
	} else{
		convSiteType(stt, Cnew, Cnew->C2->C2, ZZ);
		newSite = addSite(FE, stt_coupled->C2, Cnew, stt);
	}
	if ((int)stt_coupled->type < 2000) updateSites(stt_coupled, stt_coupled->C1, stt_coupled->C2, 0);
	else updateSites(stt_coupled, stt_coupled->C1, stt_coupled->C2, 0);
	
	Spointer S1 = moveIt(sFE2, -1);
	Spointer S2 = moveIt(sFE2, +1);
	
	//Adjust sites for ending position.
	//Migration outside.
	if ((int)checkR5_1->type == 0 && (int)sFE2->type == 0){
		if (b4) {
			Spointer sFE2_right = moveIt(sFE2, +1);
			updateSites(sFE2_right, CRem_before, sFE2_right->C2, 0);
			convSiteType(checkR5_1, CRem_next, CRem_before, R5);
			updateSites(checkR5_2, checkR5_2->C1, CRem_next, 0);
		}
		else {
			Spointer sFE2_left = moveIt(sFE2, -1);
			updateSites(sFE2_left, sFE2_left->C1, CRem_before, 0);
			convSiteType(checkR5_1, CRem_before, CRem_next, R5);
			updateSites(checkR5_2, CRem_next, checkR5_2->C2, 0);
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type == 0){
		//This means that the pentagon has migrated to the edge but will have a single carbon out of the structure.
		if (b4) {
			updateSites(S1, S1->C1, sFE2->C1, 0);
			updateSites(S2, sFE2->C1, S2->C2, 0);
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			updateSites(S2, sFE2->C2, S2->C2, 0);
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 1 && (int)sFE2->type <= 4){
		//This means that the pentagon has migrated to a basic site.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		Spointer S1 = moveIt(sFE2, -1);
		Spointer S2 = moveIt(sFE2, +1);
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){
		//This means that the pentagon has migrated to a bay containing an R5.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	
	//Update combined sites
	Spointer stt_1, stt_2, S3, S4;
	if (b4){
		stt_1 = moveIt(newSite, -1);
		stt_2 = moveIt(stt, +1);
	}else{
		stt_1 = moveIt(stt, -1);
		stt_2 = moveIt(newSite, +1);
	}
	S3 = moveIt(S1, -1);
	S4 = moveIt(S2, +1);
	updateCombinedSites(stt_1); updateCombinedSites(stt); updateCombinedSites(newSite); updateCombinedSites(stt_2);
	updateCombinedSites(S1); updateCombinedSites(S2); updateCombinedSites(S3); updateCombinedSites(S4);
}

// ************************************************************
// ID68- Recursive embedded 5-member ring migration from FEACR5 site
// ************************************************************
void PAHProcess::proc_M5R_FEACR5_multiple_sites(Spointer& stt, Cpointer C_1, Cpointer C_2, Spointer& sFE2, bool b4) {
	// The pentagon migrated N times and ended at the same position.
	if (sFE2 == stt) return;
	// First select carbons and sites affected.
	Cpointer CFE, CRem, CRem_next, CRem_before, CR5_otherside_1, CR5_otherside_2;
	Spointer checkR5_1, checkR5_2;
	if (b4) {
		CFE = C_2->C1->C1->C1;
		CR5_otherside_1 = C_2->C1->C1;
		CR5_otherside_2 = C_2->C1->C1->C1;
		checkR5_1 = moveIt(sFE2, -1);
		checkR5_2 = moveIt(sFE2, -2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C1;
		else CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
	}
	else {
		CFE = C_1->C2->C2;
		CR5_otherside_1 = C_1->C2->C2;
		CR5_otherside_2 = C_1->C2->C2->C2;
		checkR5_1 = moveIt(sFE2, +1);
		checkR5_2 = moveIt(sFE2, +2);
		if (checkR5_1->type == FE && sFE2->type == FE) CRem = sFE2->C2;
		else CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
	}
	//int end_site_type = (int)sFE2->type;

	//Add a new carbon between current R5 carbons of ACR5
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2;
	if (R5_dist < 2.5) dist2 = 1.4;
	else dist2 = R5_dist / 2.7 * 1.5;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir, normvec;
	if (b4) {
		R5dir = get_vector(C_1->C2->coords,C_1->C2->C2->coords);
		normvec = (norm_vector(C_1->C2->coords, C_1->C2->C2->coords, C_1->C2->C2->C2->coords));
	}
	else {
		R5dir = get_vector(C_2->C1->C1->coords,C_2->C1->coords);
		normvec = (norm_vector(C_2->C1->C1->coords, C_2->C1->coords, C_1->coords));
	}
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, 1.4);
	updateA(Cnew, 'H', crossvec);
	
	//Remove carbon from end site 
	removeC(CRem, false);
	OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol, 500);
	passbackPAH(mol);
	
	//First adjust starting site and add new site if needed.
	Spointer stt_coupled, newSite;
	if (b4) stt_coupled = moveIt(stt,-1);
	else stt_coupled = moveIt(stt,+1);
	if (b4) {
		updateSites(stt, Cnew, stt->C2, 0);
		newSite = addSite(ZZ, Cnew->C1->C1, Cnew, stt);
	} else{
		updateSites(stt, stt->C1, Cnew, 0);
		newSite = addSite(ZZ, Cnew, Cnew->C2->C2, stt_coupled);
	}
	
	Spointer S1 = moveIt(sFE2, -1);
	Spointer S2 = moveIt(sFE2, +1);
	
	//Adjust sites for ending position.
	//Migration outside.
	if ((int)checkR5_1->type == 0 && (int)sFE2->type == 0){
		if (b4) {
			Spointer sFE2_right = moveIt(sFE2, +1);
			updateSites(sFE2_right, CRem_before, sFE2_right->C2, 0);
			convSiteType(checkR5_1, CRem_next, CRem_before, R5);
			updateSites(checkR5_2, checkR5_2->C1, CRem_next, 0);
		}
		else {
			Spointer sFE2_left = moveIt(sFE2, -1);
			updateSites(sFE2_left, sFE2_left->C1, CRem_before, 0);
			convSiteType(checkR5_1, CRem_before, CRem_next, R5);
			updateSites(checkR5_2, CRem_next, checkR5_2->C2, 0);
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type == 0){
		//This means that the pentagon has migrated to the edge but will have a single carbon out of the structure.
		if (b4) {
			updateSites(S1, S1->C1, sFE2->C1, 0);
			updateSites(S2, sFE2->C1, S2->C2, 0);
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			updateSites(S2, sFE2->C2, S2->C2, 0);
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 1 && (int)sFE2->type <= 4){
		//This means that the pentagon has migrated to a basic site.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){
		//This means that the pentagon has migrated to neighbour an edge R5 edge.
		Spointer S1 = moveIt(sFE2, -1);
		Spointer S2 = moveIt(sFE2, +1);
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){
		//This means that the pentagon has migrated to a bay containing an R5.
		if (b4) {
			updateSites(S2, sFE2->C1, S2->C2, 0);
			/*if ((int)S2->type<2000) {
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1, S4->C2,+500);
			}*/
		}
		else {
			updateSites(S1, S1->C1, sFE2->C2, 0);
			/*if ((int)S1->type<2000) {
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2,+500);
			}*/
		}
		removeSite(sFE2);
	}
	
	//Update combined sites
	Spointer stt_1, stt_2, S3, S4;
	if (b4){
		stt_1 = moveIt(newSite, -1);
		stt_2 = moveIt(stt, +1);
	}else{
		stt_1 = moveIt(stt, -1);
		stt_2 = moveIt(newSite, +1);
	}
	S3 = moveIt(S1, -1);
	S4 = moveIt(S2, +1);
	updateCombinedSites(stt_1); updateCombinedSites(stt); updateCombinedSites(newSite); updateCombinedSites(stt_2);
	updateCombinedSites(S1); updateCombinedSites(S2); updateCombinedSites(S3); updateCombinedSites(S4);
}

// ************************************************************
// ID68- Embedded 5-member ring migration from FEACR5 site
// ************************************************************
void PAHProcess::proc_M5R_FEACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	for (Cpointer Ccheck = C_1->C1; Ccheck != C_2; Ccheck = Ccheck->C2){
		if (getDistance_twoC(Ccheck, Ccheck->C2)>1.6 || getDistance_twoC(Ccheck, Ccheck->C2)<1.2 ){
			if (!m_pah->m_optimised){
				OpenBabel::OBMol mol = passPAH();
				mol = optimisePAH(mol);
				passbackPAH(mol);
			}
		}
	}

	//First check if R6 is to the left or the right of R5
	//Get R5 internal coordinates
	cpair R5coords;
	bool check_left = true;
	bool check_right = true;
	bool optimised = true;
	if (stt->type == ACR5) R5coords = findR5internal(stt->C1->C2, stt->C2->C1);
	else{
		if ( isR5internal(stt->C1->C2, stt->C1->C2->C2) ) {
			R5coords = findR5internal(stt->C1->C2, stt->C1->C2->C2);
			check_right = false;
		}
		else if ( isR5internal(stt->C2->C1->C1, stt->C2->C1) ) {
			R5coords = findR5internal(stt->C2->C1->C1, stt->C2->C1);
			check_left = false;
		}
		else {
			//R5 not found
			return;
		}
	}

	//Fundamental assumption: R5-R7 pairs cannot move away from each other!
	if (m_pah->m_R7loc.size()>=1){
		for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
			double distR5R7 = getDistance_twoC(*it, R5coords);
			if (distR5R7 < 3.1) {
				break;
			}
		}
	}

	//Check site to the left
	Spointer site_left = moveIt(stt,-1);
	if (check_left){
		//Check site to the left
		Spointer checkR5_1 = moveIt(site_left,-1);
		Spointer checkR5_2 = moveIt(site_left,-2);
		//Check for unsupported sites. This section heavily assumes that the Isolated Pentagon Rule is valid.
		if ((int)site_left->type == 0 && (int)checkR5_1->type == 0 && (int)checkR5_2->type == 0) check_left = false; // The result would be an indene, not supported. YET!
		if ((int)site_left->type == 100 || (int)site_left->type == 101 || (int)site_left->type == 501 || (int)site_left->type == 2002) check_left = false; // This would violate the IPR.
		if ((int)site_left->type == 9999 || (int)site_left->type == -1 || site_left->type == None) check_left = false;
		if ((int)site_left->type == 0){
			if ((int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 || (int)checkR5_1->type == 100) check_left = false;
			if ((int)checkR5_1->type >= 1002 && (int)checkR5_1->type <= 1004) check_left = false;
			if ((int)checkR5_1->type >= 2002 && (int)checkR5_1->type <= 2204) check_left = false;
			if ((int)checkR5_1->type >= 2204 && (int)checkR5_1->type <= 2205) check_left = false;
			if ((int)checkR5_1->type == 0){
				if ((int)checkR5_2->type == 101 || (int)checkR5_2->type == 501) check_left = false;
				if ((int)checkR5_2->type >= 1002 && (int)checkR5_2->type <= 1004) check_left = false;
				if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2205) check_left = false;
			}
		}
		Cpointer CR5_otherside_end = site_left->C2->C1;
		if (CR5_otherside_end->bridge) check_left = false;
		if (CR5_otherside_end->C2->bridge) check_left = false; 
		if (CR5_otherside_end->C1->bridge) check_left = false;
		if (CR5_otherside_end->C2->C2->bridge) check_left = false; 
		if (CR5_otherside_end->C1->C1->bridge) check_left = false;
		
		//Check for other side being valid
		Cpointer thirdC_after = findThirdC(CR5_otherside_end);
		if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
			Spointer opp_site_after = findSite(thirdC_after);
			//I have no clue how to flag an Spointer as error. This can cause seg faults.
			if (opp_site_after != m_pah->m_siteList.end()){
				int os_endtype = opp_site_after->type;
				if (os_endtype >= 200 && os_endtype <= 203) check_left = false;
				if (os_endtype == 101) check_left = false;
				if (os_endtype >= 600 && os_endtype <= 603) check_left = false;
				if (os_endtype >= 1000 && os_endtype <= 1003) check_left = false;
				if (os_endtype >= 500 && os_endtype <= 504) check_left = false;
				if (os_endtype >= 2000 && os_endtype <= 2205) check_left = false;
				if (os_endtype >= 2103 && os_endtype <= 2105) check_left = false;
				if (os_endtype >= 2204 && os_endtype <= 2205) check_left = false;
				if (os_endtype == 9999 || os_endtype == -1 || opp_site_after->type == None) check_left = false;
			}
		}
		
		//check that two pentagons (including internals) will not collide
		cpair R5coords_end = endposR5internal(CR5_otherside_end, CR5_otherside_end->C2);
		if (m_pah->m_R5loc.size()>=1){
			for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
				double distR5s = getDistance_twoC(*it, R5coords_end);
				if (distR5s < 2.8) {
					//This distance is a parameter of this jump process. Might need some more tuning. 
					//2.8 seems appropiate but may reject too many jumps.
					//Two pentagons will be next to each other violating the Isolated Pentagon Rule
					check_left = false;
					break;
				}
			}
		}
		
		//check that pentagon and heptagon (including internals) will not collide
		if (m_pah->m_R7loc.size()>=1){
			for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
				double distR5R7 = getDistance_twoC(*it, R5coords_end);
				if (distR5R7 < 2.6) {
					//This distance is a parameter of this jump process. Might need some more tuning. 
					//2.8 seems appropiate but may reject too many jumps.
					//Two pentagons will be next to each other violating the Isolated Pentagon Rule
					check_left = false;
					break;
				}
			}
		}
	}
	//Check site to the right
	Spointer site_right = moveIt(stt,+1);
	if (check_right){
		//Check site to the right
		Spointer checkR5_1 = moveIt(site_right,+1);
		Spointer checkR5_2 = moveIt(site_right,+2);
		//Check for unsupported sites. This section heavily assumes that the Isolated Pentagon Rule is valid.
		if ((int)site_right->type == 0 && (int)checkR5_1->type == 0 && (int)checkR5_2->type == 0) check_right = false; // The result would be an indene, not supported. YET!
		if ((int)site_right->type == 100 || (int)site_right->type == 101 || (int)site_right->type == 501 || (int)site_right->type == 2002) check_right = false; // This would violate the IPR.
		if ((int)site_right->type == 9999 || (int)site_right->type == -1 || site_right->type == None) check_right = false;
		if ((int)site_right->type == 0){
			if ((int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 || (int)checkR5_1->type == 100) check_right = false;
			if ((int)checkR5_1->type >= 1002 && (int)checkR5_1->type <= 1004) check_right = false;
			if ((int)checkR5_1->type >= 2002 && (int)checkR5_1->type <= 2204) check_right = false;
			if ((int)checkR5_1->type >= 2204 && (int)checkR5_1->type <= 2205) check_right = false;
			if ((int)checkR5_1->type == 0){
				if ((int)checkR5_2->type == 101 || (int)checkR5_2->type == 501) check_right = false;
				if ((int)checkR5_2->type >= 1002 && (int)checkR5_2->type <= 1004) check_right = false;
				if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2205) check_right = false;
			}
		}
		Cpointer CR5_otherside_end = site_right->C1->C2;
		if (CR5_otherside_end->bridge) check_right = false;
		if (CR5_otherside_end->C2->bridge) check_right = false; 
		if (CR5_otherside_end->C1->bridge) check_right = false;
		if (CR5_otherside_end->C2->C2->bridge) check_right = false; 
		if (CR5_otherside_end->C1->C1->bridge) check_right = false;
		
		//Check for other side being valid
		Cpointer thirdC_after = findThirdC(CR5_otherside_end);
		if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
			Spointer opp_site_after = findSite(thirdC_after);
			//I have no clue how to flag an Spointer as error. This can cause seg faults.
			if (opp_site_after != m_pah->m_siteList.end()){
				int os_endtype = opp_site_after->type;
				if (os_endtype >= 200 && os_endtype <= 203) check_right = false;
				if (os_endtype == 101) check_right = false;
				if (os_endtype >= 600 && os_endtype <= 603) check_right = false;
				if (os_endtype >= 1000 && os_endtype <= 1003) check_right = false;
				if (os_endtype >= 500 && os_endtype <= 504) check_right = false;
				if (os_endtype >= 2000 && os_endtype <= 2205) check_right = false;
				if (os_endtype >= 2103 && os_endtype <= 2105) check_right = false;
				if (os_endtype >= 2204 && os_endtype <= 2205) check_right = false;
				if (os_endtype == 9999 || os_endtype == -1 || opp_site_after->type == None) check_right = false;
			}
		}
		
		//check that two pentagons (including internals) will not collide
		cpair R5coords_end = endposR5internal(CR5_otherside_end->C1, CR5_otherside_end,true);
		if (m_pah->m_R5loc.size()>=1){
			for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
				double distR5s = getDistance_twoC(*it, R5coords_end);
				if (distR5s < 2.8) {
					//This distance is a parameter of this jump process. Might need some more tuning. 
					//2.8 seems appropiate but may reject too many jumps.
					//Two pentagons will be next to each other violating the Isolated Pentagon Rule
					check_right = false;
					break;
				}
			}
		}
		
		//check that pentagon and heptagon (including internals) will not collide
		if (m_pah->m_R7loc.size()>=1){
			for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
				double distR5R7 = getDistance_twoC(*it, R5coords_end);
				if (distR5R7 < 2.6) {
					//This distance is a parameter of this jump process. Might need some more tuning. 
					//2.8 seems appropiate but may reject too many jumps.
					//Two pentagons will be next to each other violating the Isolated Pentagon Rule
					check_right = false;
					break;
				}
			}
		}
	}
	if ( (check_left && check_right) || (!check_left && !check_right) ) {
		//This process should only return one site
		return;
	}
	bool b4 = check_left;
	Spointer sFE2, checkR5_1, checkR5_2;
	Cpointer CRem, CRem_before, CRem_next, CFE, CR5_otherside_1, CR5_otherside_2;
	if (b4) {
		sFE2 = moveIt(stt, -1);
		checkR5_1 = moveIt(stt, -2);
		checkR5_2 = moveIt(stt, -3);
		CFE = C_1->C2;
		CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
		CR5_otherside_1 = C_1->C2->C2;
		CR5_otherside_2 = C_1->C2;
	}
	else {
		sFE2 = moveIt(stt, 1);
		checkR5_1 = moveIt(stt, 2);
		checkR5_2 = moveIt(stt, 3);
		CFE = C_2->C1->C1;
		CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
		CR5_otherside_1 = C_2->C1->C1;
		CR5_otherside_2 = C_2->C1;
	}

	// check if ACR5 has an opposite site.
	Spointer opp_site, opp_site_second, opp_site_after;
	bool opp_site_bool = false; bool opp_site_bool_second = false; bool opp_site_bool_after = false;
	Cpointer thirdC = findThirdC(CR5_otherside_1);
	Cpointer thirdC2 = findThirdC(CR5_otherside_2);
	Cpointer thirdC_after = findThirdC(CRem_next);
	// Seven cases:
	// 1. One pentagon has one exposed edge  and migrates to a location where it will have one exposed edge. Normal migration.
	// 2. One pentagon has one exposed edge  and migrates to a location where it will have two exposed edges. 
	// 3. One pentagon has two exposed edges and migrates to a location where it will have two exposed edges. 
	// 4. One pentagon has two exposed edges and migrates to a location where it will have one exposed edge.
	// 5. One pentagon has two exposed edges and migrates to a location where it will have three exposed edges.
	// 6. One pentagon has three exposed edges and migrates to a location where it will have three exposed edges.
	// 7. One pentagon has three exposed edges and migrates to a location where it will have two exposed edges.
	
	if (thirdC != NULLC) {
		opp_site = findSite(thirdC);
		if (opp_site != m_pah->m_siteList.end()) opp_site_bool = true;
	}
	if (thirdC2 != NULLC) {
		opp_site_second = findSite(thirdC2);
		if (opp_site_second != m_pah->m_siteList.end() && opp_site_second!=opp_site) opp_site_bool_second = true;
	}
	if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
		opp_site_after = findSite(thirdC_after);
		if (opp_site_after != m_pah->m_siteList.end()) {
			opp_site_bool_after = true;
			int os_endtype = opp_site_after->type;
			if (os_endtype >= 200 && os_endtype <= 203) return;
			if (os_endtype == 101) return;
			if (os_endtype >= 600 && os_endtype <= 603) return;
			if (os_endtype >= 1000 && os_endtype <= 1003) return;
			if (os_endtype >= 500 && os_endtype <= 504) return;
			if (os_endtype >= 2000 && os_endtype <= 2205) return;
			if (os_endtype >= 2103 && os_endtype <= 2105) return;
			if (os_endtype >= 2204 && os_endtype <= 2205) return;
		}
	}

	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2;
	if (R5_dist < 2.5) dist2 = 1.4;
	else dist2 = R5_dist / 2.7 * 1.5;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir = get_vector(CFE->coords,CFE->C2->coords);
	cpair normvec = (norm_vector(CFE->coords, CFE->C2->coords, CFE->C2->C2->coords));
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, 1.4);
	updateA(Cnew, 'H', crossvec);
	removeC(CRem, false);
	//The ACR5 site is on the "short" side of an R5. 
	OpenBabel::OBMol newmol = passPAH();
	newmol = optimisePAH(newmol, 250); // Only 250 steps for this test. Although more are probably needed.
	passbackPAH(newmol);
	
	if ((int)checkR5_1->type == 0 && (int)sFE2->type == 0){
		Cpointer C1_new, C2_new;
		if (b4) {
			C1_new = checkR5_1->C1->C1;
			C2_new = C1_new->C2->C2->C2;
		}
		else {
			C2_new = checkR5_1->C2->C2;
			C1_new = C1_new->C1->C1->C1;
		}
		//R5 has moved to the edge and will now be free.
		//removeC(C1_new->C2, false); removeC(C1_new->C2, false);
		//addC(C1_new, normAngle(C1_new->C1->bondAngle1 - 60), normAngle(C1_new->C1->bondAngle1), 1, true);
		if (b4) {
			sFE2->C1 = Cnew->C1;
			sFE2->C2 = Cnew;
			updateSites(stt, sFE2->C2, C_2, -2001);
		}
		else {
			sFE2->C1 = CFE;
			sFE2->C2 = CFE->C2;
			updateSites(stt, C_1, sFE2->C1, -2001);
		}
		if (!optimised){
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FE);
			//convSiteType(checkR5_1, C1_new, C1_new->C2->C2, ZZ);
			updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, -1);
			m_pah->m_rings5_Lone--;
			redrawR5(checkR5_1, C1_new, C2_new);
			//proc_G5R_ZZ(checkR5_1, checkR5_1->C1, checkR5_1->C2);
		}
		else {
			if (b4) {
				convSiteType(sFE2, sFE2->C1->C1, sFE2->C2, RFE);
				convSiteType(checkR5_1, C1_new->C2, C2_new->C1, R5);
				updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
			}
			else {
				convSiteType(sFE2, sFE2->C1, sFE2->C2->C2, RFE);
				convSiteType(checkR5_1, C1_new->C2, C2_new->C1, R5);
				updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
			}
		}
	}
	else {
		//Reassign connectivity at next site
		if (b4) sFE2->C2 = CFE;
		else sFE2->C1 = Cnew;
	}
	if (!m_pah->m_optimised){
		if (b4) addR5internal(sFE2->C1, sFE2->C1->C2,true);
		else addR5internal(sFE2->C2->C1, sFE2->C2,true);
	}
	
	// edit sites. first identify the neighbouring sites of resulting RFE & R5
	Spointer S1, S2, S3, S4;
	if (b4) {
		S1 = moveIt(sFE2, -1);
		S2 = moveIt(stt, +1);
		if ((int)sFE2->type == 0 && (int)S1->type == 0){ //sFE2 is a FE
			//This means that the pentagon has migrated to the edge. This is handled above.
			//convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			//updateSites(S1, S1->C1, sFE2->C1, 500);
		}
		else if ((int)sFE2->type == 0 ){ 
			//This means that the pentagon has migrated to the edge but will have a carbon out of the structure.
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			updateSites(S1, S1->C1, sFE2->C1, 500);
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C2->C1);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
			if (!optimised) addR5internal(sFE2->C2->C1, sFE2->C2);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
			if (!optimised) addR5internal(sFE2->C2->C1, sFE2->C2);
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
			if (!optimised) addR5internal(sFE2->C2->C1, sFE2->C2);
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
			if (!optimised) addR5internal(sFE2->C2->C1, sFE2->C2);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
			if (!optimised) addR5internal(sFE2->C2->C1->C1, sFE2->C2->C1);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
			if (!optimised) addR5internal(sFE2->C2->C1->C1, sFE2->C2->C1);
		}
		updateSites(stt, sFE2->C2, C_2, -2001);
	}
	else {
		S1 = moveIt(stt, -1); 
		S2 = moveIt(sFE2, 1); 
		if ((int)sFE2->type == 0 && (int)S2->type == 0){ //sFE2 is a FE
			//This means that the pentagon has migrated to the edge. This is handled above.
			//convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			//updateSites(S2, sFE2->C2, S2->C2, 500);
		}
		else if ((int)sFE2->type == 0){
			//This means that the pentagon has migrated to the edge but will have a carbon out of the structure.
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			updateSites(S2, sFE2->C2, S2->C2, 500);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C2);
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C2->C1);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C1->C2->C2);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
			if (!optimised) addR5internal(sFE2->C1->C2, sFE2->C1->C2->C2);
		}
		updateSites(stt, stt->C1, sFE2->C1, -2001);
	}
	// update H atoms
	/*if (b4){
		updateA(S1->C2->C1, C_2, 'H');
	}
	else{
		updateA(C_1, S2->C1->C2, 'H');
	}*/
	if (opp_site_bool && !opp_site_bool_second && !opp_site_bool_after) {
		if ( (int)opp_site->type >= 2100) {
			Spointer S1_opp_site = moveIt(opp_site, -1);
			Spointer S2_opp_site = moveIt(opp_site, +1);
			if (S1_opp_site->type==R5 || S2_opp_site->type==R5){
				updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
			}
			else updateSites(opp_site, opp_site->C1, opp_site->C2, -100);
		}
		else updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		Spointer S1_opp_site = moveIt(opp_site, -1);
		Spointer S2_opp_site = moveIt(opp_site, +1);
		updateCombinedSites(opp_site); updateCombinedSites(S1_opp_site);  updateCombinedSites(S2_opp_site); 
	}
	else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		if (opp_site_after != S1 && opp_site_after != S2) {
			if ((int)opp_site_after->type >= 500 && (int)opp_site_after->type <= 700) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -400);
			else if ((int)opp_site_after->type >= 1000 && (int)opp_site_after->type <= 2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -800);
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			Spointer S1_opp_site = moveIt(opp_site_after, -1);
			Spointer S2_opp_site = moveIt(opp_site_after, +1);
			updateCombinedSites(opp_site_after); updateCombinedSites(S1_opp_site);  updateCombinedSites(S2_opp_site); 
		}
		updateCombinedSites(opp_site);
	}
	else if (opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
		//convSiteType(opp_site_second, opp_site_second->C1, opp_site_second->C2, (kmcSiteType)new_stype);
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, +1500);
		updateCombinedSites(opp_site);
		updateCombinedSites(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateCombinedSites(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, -1500);
		updateCombinedSites(opp_site_second);
		if (opp_site_after != S1 && opp_site_after != S2) {
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
			updateCombinedSites(opp_site_after);
		}
	}
	else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		if (opp_site_after != S1 && opp_site_after != S2) {
			if ( (int)opp_site_after->type >=2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +100);
			else updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			updateCombinedSites(opp_site_after);
		}
	}
	else if (opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
		updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
		updateCombinedSites(opp_site);
		updateCombinedSites(opp_site_second);
		updateCombinedSites(opp_site_after);
	}
	//printStruct();
	// update combined sites for all sites involved and their neighbours
	// (excluding new FE sites, since their combined site type will still be None)
	if (b4){
		S3 = moveIt(S1, -1);
	}
	else{
		S4 = moveIt(S2, 1);
	}
	updateCombinedSites(stt); updateCombinedSites(sFE2); 
	if (b4){
		updateCombinedSites(S1);
		updateCombinedSites(S3);
	}
	else{
		updateCombinedSites(S2);
		updateCombinedSites(S4); // neighbours
	}
}

// ************************************************************
// ID24- ACR5 migration without modifying typespace
// ************************************************************
void PAHProcess::proc_M5R_ACR5_ZZ_light(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//First check if R6 is to the left or the right of R5
	bool b4;
	int steps=-99999;
	unsigned int ii = 0;
	for (ii = 0;ii != m_pah->m_R5walker_sites.size();ii++){
		Spointer migr_site_start = std::get<0>(m_pah->m_R5walker_sites[ii]);
		Spointer site_check = moveIt(migr_site_start, std::get<2>(m_pah->m_R5walker_sites[ii]));
		Spointer migr_site_start_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
		Spointer site_check_2 = moveIt(migr_site_start_2, std::get<2>(m_pah->m_R5walker_sites[ii]));
		if (site_check == stt && site_check_2 == stt){
			steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
			break;
		}
	}
	//cpair R5coords;
	if (steps==-99999){
		std::cout << "Error on R5 migration in proc_M5R_ACR5_ZZ_light. Walker not found." << std::endl;
		printSitesMigration();
		return;
	}
	if (steps > 0) b4 = true;
	else if (steps<0) b4 = false;
	else{
		//We need to identify if we can move left or right.
		if ((int)stt->type>=2003){
			//An R5 was added to tell us the direction
			if (isR5internal(stt->C1->C2,stt->C1->C2->C2) ) b4 = true;
			else if (isR5internal(stt->C2->C1->C1,stt->C2->C1) ) b4 = false;
			else{
				//saveXYZ("KMC_DEBUG/FEACR5MIGR");
				std::cout << "Error on R5 migration in proc_M5R_ACR5_ZZ_light. Walker with needed R5loc not found." << std::endl;
				printSitesMigration();
				return;
			}
			cpair R5_coords_rem;
			if (b4) R5_coords_rem = findR5internal(stt->C1->C2,stt->C1->C2->C2);
			else R5_coords_rem = findR5internal(stt->C2->C1->C1,stt->C2->C1);
		}
		else{
			Spointer left_site = moveIt(stt,-1);
			bool check_left = checkSiteMigration (left_site,true);
			Spointer right_site = moveIt(stt,+1);
			bool check_right = checkSiteMigration (right_site,false);
			if (check_left && !check_right) b4 = true;
			else if (!check_left && check_right) b4 = false;
			else{
				//saveXYZ("KMC_DEBUG/FEACR5MIGR");
				std::cout << "Error on R5 migration in proc_M5R_ACR5_ZZ_light. Walker was allowed to migrate in both directions." << std::endl;
				printSitesMigration();
				return;
			}
		}
	}


	Spointer sFE2, checkR5_1, checkR5_2;
	Cpointer CRem, CRem_next, CFE, CR5_otherside_1, CR5_otherside_2;
	if (b4) {
		sFE2 = moveIt(stt, -1);
		checkR5_1 = moveIt(stt, -2);
		checkR5_2 = moveIt(stt, -3);
		CFE = C_1->C2; //Not modified but pointer created
		CRem = sFE2->C2; //Not modified but pointer created
		CRem_next = CRem->C1;
		if (steps == 0){
			CR5_otherside_1 = C_1->C2->C2;
			CR5_otherside_2 = C_1->C2;
			CRem_next = CRem->C1;
		} else {
			CR5_otherside_1 = C_1->C2;
			CR5_otherside_2 = C_1->C1;
			if (CR5_otherside_2->C1->A=='H') CRem_next = CR5_otherside_2->C1->C1;
			else CRem_next = CR5_otherside_2->C1;
		} 
	}
	else {
		sFE2 = moveIt(stt, 1);
		checkR5_1 = moveIt(stt, 2);
		checkR5_2 = moveIt(stt, 3);
		CFE = C_2->C1->C1; //Not modified but pointer created
		CRem = sFE2->C1; //Not modified but pointer created
		CRem_next = CRem->C2; 
		if (steps == 0){
			CR5_otherside_1 = C_2->C1->C1;
			CR5_otherside_2 = C_2->C1;
			CRem_next = CRem->C2;
		} else {
			CR5_otherside_1 = C_2->C1;
			CR5_otherside_2 = C_2->C2;
			if (CR5_otherside_2->C2->A=='H') CRem_next = CR5_otherside_2->C2->C2;
			else CRem_next = CR5_otherside_2->C2;
		} 
	}

	// check if ACR5 has an opposite site.
	Spointer opp_site, opp_site_second, opp_site_after;
	bool opp_site_bool = false; bool opp_site_bool_second = false; bool opp_site_bool_after = false;
	Cpointer thirdC = findThirdC(CR5_otherside_1);
	Cpointer thirdC2 = findThirdC(CR5_otherside_2);
	Cpointer thirdC_after = findThirdC(CRem_next);
	// Seven cases:
	// 1. One pentagon has one exposed edge  and migrates to a location where it will have one exposed edge. Normal migration.
	// 2. One pentagon has one exposed edge  and migrates to a location where it will have two exposed edges. 
	// 3. One pentagon has two exposed edges and migrates to a location where it will have two exposed edges. 
	// 4. One pentagon has two exposed edges and migrates to a location where it will have one exposed edge.
	// 5. One pentagon has two exposed edges and migrates to a location where it will have three exposed edges.
	// 6. One pentagon has three exposed edges and migrates to a location where it will have three exposed edges.
	// 7. One pentagon has three exposed edges and migrates to a location where it will have two exposed edges.
	
	if (thirdC != NULLC) {
		opp_site = findSite(thirdC);
		if (opp_site != m_pah->m_siteList.end()) opp_site_bool = true;
	}
	if (thirdC2 != NULLC) {
		opp_site_second = findSite(thirdC2);
		if (opp_site_bool && opp_site_second == opp_site) opp_site_bool_second = false;
		else if (opp_site_second != m_pah->m_siteList.end()) {
			opp_site_bool_second = true;
		}
	}
	bool end_site_allowed = true;
	if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
		opp_site_after = findSite(thirdC_after);
		if (opp_site_after != m_pah->m_siteList.end()) {
			opp_site_bool_after = true;
			int os_endtype = opp_site_after->type;
			if (os_endtype >= 200 && os_endtype <= 203) end_site_allowed = false;
			if (os_endtype == 101) end_site_allowed = false;
			if (os_endtype >= 600 && os_endtype <= 603) end_site_allowed = false;
			if (os_endtype >= 1000 && os_endtype <= 1003) end_site_allowed = false;
			if (os_endtype >= 500 && os_endtype <= 504) end_site_allowed = false;
			if (os_endtype >= 2000 && os_endtype <= 2205) end_site_allowed = false;
			if (os_endtype >= 2103 && os_endtype <= 2105) end_site_allowed = false;
			if (os_endtype >= 2204 && os_endtype <= 2205) end_site_allowed = false;
		}
	}
	if (!end_site_allowed){
		//m_pah->m_R5loc.push_back(R5coords);
		return;
	}

	//After this point we know that the process is accepted!
	//Check migration on the other side
	if (opp_site_bool && !opp_site_bool_second && !opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
	}else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
		addOppsiteR5Walker(opp_site_after,opp_site_after);
	}else if (opp_site_bool && opp_site_bool_second && !opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
	}else if (!opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		int jj = findWalker(opp_site_second);
		ii = remOppsiteR5Walker(ii, jj);
	}else if (!opp_site_bool && opp_site_bool_second && opp_site_bool_after){
		int jj = findWalker(opp_site_second);
		ii = remOppsiteR5Walker(ii, jj);
		addOppsiteR5Walker(opp_site_after,opp_site_after);
	}else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
			addOppsiteR5Walker(opp_site_after,opp_site_after);
	}else if (opp_site_bool && opp_site_bool_second && opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
		addOppsiteR5Walker(opp_site_after,opp_site_after);
	}
	
	if (b4) {
		std::get<2>(m_pah->m_R5walker_sites[ii])--;
	}
	else {
		std::get<2>(m_pah->m_R5walker_sites[ii])++;
	}
	checkR5Walkers(ii);
	bool opp_site_logic = false;
	
	if ((int)checkR5_1->type == 0 && (int)sFE2->type == 0){
		//R5 has moved to the edge and will now be free.
		//Since the FE2 will not move anymore it is probably a good idea to remove the walker.
		updateSites(stt, stt->C1, stt->C2, -1901);
		convSiteType(sFE2, sFE2->C1, sFE2->C2, R5);
		convSiteType(checkR5_1, checkR5_1->C1, checkR5_1->C2, R5);
		updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
		Spointer site_perf = std::get<0>(m_pah->m_R5walker_sites[ii]);
		Spointer site_perf_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
		if(site_perf==site_perf_2) proc_M5R_ACR5_termination_toR5(site_perf,site_perf->C1,site_perf->C2,sFE2,b4,steps);
		else proc_M5R_R5R6_multiple_sites(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		//sFE2 is deleted in the previous call so we must not call it again from here.
		//updateCombinedSitesMigration(stt); updateCombinedSitesMigration(checkR5_1); updateCombinedSitesMigration(checkR5_2);
		m_pah->m_R5walker_sites.erase(m_pah->m_R5walker_sites.begin()+ii);
		addR5internal(checkR5_1->C1,checkR5_1->C2);
		Spointer start_site, end_site, S1, S2, S3, S4, S5, S6;
		if (b4) {
			start_site = moveIt(stt,1);
			end_site = moveIt(start_site, -1);
			S1 = moveIt(end_site, -1);
			S2 = moveIt(start_site, +1);
			S3 = moveIt(S1, -1);
			S4 = moveIt(S2, 1);
			S5 = moveIt(S1, -2);
			S6 = moveIt(S2, 2);
			updateCombinedSitesMigration(start_site);
			updateCombinedSitesMigration(end_site); 
			updateCombinedSitesMigration(S1);
			updateCombinedSitesMigration(S3);
			updateCombinedSitesMigration(S2);
			updateCombinedSitesMigration(S4);
			updateCombinedSitesMigration(S5);
			updateCombinedSitesMigration(S6);
		}
		else {
			start_site = moveIt(stt,-1);
			end_site = moveIt(start_site, +1);
			S1 = moveIt(start_site, -1); 
			S2 = moveIt(end_site, 1);
			S3 = moveIt(S1, -1);
			S4 = moveIt(S2, 1);
			S5 = moveIt(S1, -2);
			S6 = moveIt(S2, 2);
			updateCombinedSitesMigration(start_site);
			updateCombinedSitesMigration(end_site); 
			updateCombinedSitesMigration(S2);
			updateCombinedSitesMigration(S4); // neighbours
			updateCombinedSitesMigration(S1);
			updateCombinedSitesMigration(S3);
			updateCombinedSitesMigration(S5);
			updateCombinedSitesMigration(S6);
		}
		//Need to check for opposite site logic before returning.
		opp_site_logic = true;
	}

	Spointer S1_check, S2_check;
	if (b4){
		S1_check = moveIt(sFE2, -1);
		S2_check = moveIt(stt, +1);
	} else{
		S1_check = moveIt(stt, -1); 
		S2_check = moveIt(sFE2, 1); 
	}
	if (opp_site_bool && !opp_site_bool_second && !opp_site_bool_after) {
		if ( (int)opp_site->type >= 2100) {
			Spointer S1_opp_site = moveIt(opp_site, -1);
			Spointer S2_opp_site = moveIt(opp_site, +1);
			if (S1_opp_site->type==R5 || S2_opp_site->type==R5){
				updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
			}
			else updateSites(opp_site, opp_site->C1, opp_site->C2, -100);
		}
		else updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		Spointer S1_opp_site = moveIt(opp_site, -1);
		Spointer S2_opp_site = moveIt(opp_site, +1);
		updateCombinedSitesMigration(opp_site); updateCombinedSitesMigration(S1_opp_site);  updateCombinedSitesMigration(S2_opp_site); 
	}
	else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		if (opp_site_after != S1_check && opp_site_after != S2_check) {
			if ((int)opp_site_after->type >= 500 && (int)opp_site_after->type <= 700) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -400);
			else if ((int)opp_site_after->type >= 1000 && (int)opp_site_after->type <= 2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -800);
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			Spointer S1_opp_site = moveIt(opp_site_after, -1);
			Spointer S2_opp_site = moveIt(opp_site_after, +1);
			updateCombinedSitesMigration(opp_site_after); updateCombinedSitesMigration(S1_opp_site);  updateCombinedSitesMigration(S2_opp_site); 
		}
		updateCombinedSitesMigration(opp_site);
	}
	else if (opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
		//convSiteType(opp_site_second, opp_site_second->C1, opp_site_second->C2, (kmcSiteType)new_stype);
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, +1500);
		updateCombinedSitesMigration(opp_site);
		updateCombinedSitesMigration(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateCombinedSitesMigration(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, -1500);
		updateCombinedSitesMigration(opp_site_second);
		if (opp_site_after != S1_check && opp_site_after != S2_check) {
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
			updateCombinedSitesMigration(opp_site_after);
		}
	}
	else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		if (opp_site_after != S1_check && opp_site_after != S2_check) {
			if ( (int)opp_site_after->type >=2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +100);
			else updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			updateCombinedSitesMigration(opp_site_after);
		}
	}
	else if (opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
		updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
		updateCombinedSitesMigration(opp_site);
		updateCombinedSitesMigration(opp_site_second);
		updateCombinedSitesMigration(opp_site_after);
	}

	if (opp_site_logic) return;
	
	// edit sites. first identify the neighbouring sites of resulting RFE & R5
	Spointer S1, S2, S3, S4, S5, S6;
	if (b4) {
		S1 = moveIt(sFE2, -1);
		S2 = moveIt(stt, +1);
		if ((int)sFE2->type == 0 && (int)S1->type == 0){ //sFE2 is a FE
			//This means that the pentagon has migrated to the edge. This is handled above.
			//convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			//updateSites(S1, S1->C1, sFE2->C1, 500);
		}
		else if ((int)sFE2->type == 0 ){ 
			//This means that the pentagon has migrated to a edge but will have a carbon out of the structure.
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			updateSites(S1, S1->C1, S1->C2, 500);
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
		}
		if ((int)stt->type>2100){
			if ((int)S2->type>=500) updateSites(stt, sFE2->C2, C_2, -1601);
			else updateSites(stt, sFE2->C2, C_2, -2001);
		} else updateSites(stt, sFE2->C2, C_2, -2001);
	}
	else {
		S1 = moveIt(stt, -1); 
		S2 = moveIt(sFE2, 1); 
		if ((int)sFE2->type == 0 && (int)S2->type == 0){ //sFE2 is a FE
			//This means that the pentagon has migrated to the edge. This is handled above.
			//convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			//updateSites(S2, sFE2->C2, S2->C2, 500);
		}
		else if ((int)sFE2->type == 0){
			//This means that the pentagon has migrated to the edge but will have a carbon out of the structure.
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			updateSites(S2, S2->C1, S2->C2, 500);
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
		}
		if ((int)stt->type>2100){
			if ((int)S1->type>=500) updateSites(stt, stt->C1, sFE2->C1, -1601);
			else updateSites(stt, stt->C1, sFE2->C1, -2001);
		} else updateSites(stt, stt->C1, sFE2->C1, -2001);
	}
	steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
	if((int)sFE2->type>=2003 && steps == 0) {
		if(b4) addR5internal(sFE2->C2->C1->C1,sFE2->C2->C1,true);
		else addR5internal(sFE2->C1->C2,sFE2->C1->C2->C2,true);
	}
	// update H atoms
	/*if (b4){
		updateA(S1->C2->C1, C_2, 'H');
	}
	else{
		updateA(C_1, S2->C1->C2, 'H');
	}*/
	
	//printStruct();
	// update combined sites for all sites involved and their neighbours
	// (excluding new FE sites, since their combined site type will still be None)
	if (b4){
		S3 = moveIt(S1, -1);
		S4 = moveIt(S2, 1);
		S5 = moveIt(S1, -2);
		S6 = moveIt(S2, 2);
	}
	else{
		S3 = moveIt(S1, -1);
		S4 = moveIt(S2, 1);
		S5 = moveIt(S1, -2);
		S6 = moveIt(S2, 2);
	}
	updateCombinedSitesMigration(stt);
	if (b4){
		updateCombinedSitesMigration(sFE2); 
		updateCombinedSitesMigration(S1);
		updateCombinedSitesMigration(S3);
		updateCombinedSitesMigration(S2);
		updateCombinedSitesMigration(S4);
		updateCombinedSitesMigration(S5);
		updateCombinedSitesMigration(S6);
	}
	else{
		updateCombinedSitesMigration(sFE2); 
		updateCombinedSitesMigration(S2);
		updateCombinedSitesMigration(S4); // neighbours
		updateCombinedSitesMigration(S1);
		updateCombinedSitesMigration(S3);
		updateCombinedSitesMigration(S5);
		updateCombinedSitesMigration(S6);
	}
}

// ************************************************************
// ID34- R5R6 migration without modifying typespace
// ************************************************************
void PAHProcess::proc_MR5_R6_light(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//First check if R6 is to the left or the right of R5
	bool b4;
	bool leave_edge = false;
	bool leave_corner = false;
	Spointer S1 = moveIt(stt,-1);
	Spointer S2 = moveIt(stt,+1);

	int steps=-99999;
	unsigned int ii = 0;
	for (ii = 0; ii != m_pah->m_R5walker_sites.size();ii++){
		Spointer migr_site_start = std::get<0>(m_pah->m_R5walker_sites[ii]);
		Spointer migr_site_start_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
		Spointer site_check = moveIt(migr_site_start, std::get<2>(m_pah->m_R5walker_sites[ii]));
		Spointer site_check2 = moveIt(migr_site_start_2, std::get<2>(m_pah->m_R5walker_sites[ii]));
		//Checking a walker that is on a corner
		if (migr_site_start != migr_site_start_2){
			//Walker is already at a corner
			if (site_check == stt){
				steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
				if (steps>0) {
					b4 = true;
					leave_corner = true;
				}
				else if (steps<0) {
					b4 = false;
					leave_corner = true;
				}
				else{
					if(migr_site_start == stt) {
						b4 = true;
						leave_corner = true;
					}
					else{
						//This kind of site should never have number ot steps = 0
						std::cout << "Error on R5 migration in proc_MR5_R6_light. " << std::endl;
						printSitesMigration();
						return;
					}
				}
				break;
			}
			else if(site_check2 == stt){
				steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
				if (steps>0) {
					b4 = true;
					leave_corner = true;
				}
				else if (steps<0) {
					b4 = false;
					leave_corner = true;
				}
				else{
					if (migr_site_start_2 == stt) {
						b4 = false;
						leave_corner = true;
					}
					else{
						//This kind of site should never have number ot steps = 0
						std::cout << "Error on R5 migration in proc_MR5_R6_light. " << std::endl;
						printSitesMigration();
						return;
					}
				}
				break;
			}
		} else{
			//Checking a walker that is NOT on a corner
			if (site_check == stt){
				steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
				if (steps>0){
					b4 = true;
					leave_edge = false;
				} 
				else if (steps<0) {
					b4 = false;
					leave_edge = false;
				}
				else{
					//The walker has an ACR5 but points to a corner. ERROR.
					std::cout << "Error on R5 migration in proc_MR5_R6_light. " << std::endl;
					printSitesMigration();
					return;
				}
				break;
			} else{
				if ((int)site_check->type>=501 && (int)site_check->type<=1004){
					//We need to check for the coupled site but only for corners.
					Spointer coupled_site = m_pah->m_siteList.end();
					if (std::get<2>(m_pah->m_R5walker_sites[ii])<0) coupled_site = moveIt(site_check,-1);
					else if (std::get<2>(m_pah->m_R5walker_sites[ii])>0) coupled_site = moveIt(site_check,+1);

					if(coupled_site == stt){
						steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
						if (steps>0){
							b4 = false;
							leave_edge = true;
						} 
						else if (steps<0) {
							b4 = true;
							leave_edge = true;
						}
						else{
							//The walker has an ACR5 but points to a corner. ERROR.
							std::cout << "Error on R5 migration in proc_MR5_R6_light. " << std::endl;
							printSitesMigration();
							return;
						}
						break;
					}
				}
			}
		}
	}
	if (steps==-99999){
		std::cout << "Error on R5 migration in proc_MR5_R6_light. Walker not found." << std::endl;
		printSitesMigration();
	}
	
	if (b4) {
		std::get<2>(m_pah->m_R5walker_sites[ii])--;
	}
	else {
		std::get<2>(m_pah->m_R5walker_sites[ii])++;
	}
	checkR5Walkers(ii);

	Spointer sFE2, checkR5_1, checkR5_2;
	Cpointer CRem, CRem_next, CFE, CR5_otherside_1, CR5_otherside_2;
	if (b4) {
		sFE2 = moveIt(stt, -1);
		checkR5_1 = moveIt (stt, -2);
		checkR5_2 = moveIt (stt, -3);
		CFE = C_1->C2;
		CRem = sFE2->C2;
		if (steps == 0){
			CRem_next = CRem->C1;
			CR5_otherside_1 = C_2->C2;
			CR5_otherside_2 = C_2->C1;
		} else if (steps > 0){
			CR5_otherside_1 = C_2->C2;
			CR5_otherside_2 = C_2->C1->C1;
			if (C_2->C1->C1->C1->A=='H') CRem_next = C_2->C1->C1->C1->C1;
			else CRem_next = C_2->C1->C1->C1;
		} else{
			CR5_otherside_1 = C_2->C2->C2;
			CR5_otherside_2 = C_2->C1;
			CRem_next = C_1->C1;
		}
	}
	else {
		sFE2 = moveIt(stt, 1);
		checkR5_1 = moveIt (stt, 2);
		checkR5_2 = moveIt (stt, 3);
		CFE = C_1;
		CRem = sFE2->C1;
		if (steps == 0){
			CRem_next = CRem->C2;
			CR5_otherside_1 = C_1->C1;
			CR5_otherside_2 = C_1->C2;
		} else if (steps > 0){
			CR5_otherside_1 = C_1->C1->C1;
			CR5_otherside_2 = C_1->C2;
			CRem_next = C_2->C2;
		} else{
			CR5_otherside_1 = C_1->C1;
			CR5_otherside_2 = C_1->C2->C2;
			if (C_1->C2->C2->C2->A=='H') CRem_next = C_1->C2->C2->C2->C2;
			else CRem_next = C_1->C2->C2->C2;
		}
	}

	// check if R5R6 has an opposite site.
	Spointer opp_site, opp_site_second, opp_site_after;
	bool opp_site_bool = false; bool opp_site_bool_second = false; bool opp_site_bool_after = false;
	Cpointer thirdC = findThirdC(CR5_otherside_1);
	Cpointer thirdC2 = findThirdC(CR5_otherside_2);
	Cpointer thirdC_after = findThirdC(CRem_next);
	// Seven cases:
	// 1. One pentagon has one exposed edge  and migrates to a location where it will have one exposed edge. Normal migration.
	// 2. One pentagon has one exposed edge  and migrates to a location where it will have two exposed edges. 
	// 3. One pentagon has two exposed edges and migrates to a location where it will have two exposed edges. 
	// 4. One pentagon has two exposed edges and migrates to a location where it will have one exposed edge.
	// 5. One pentagon has two exposed edges and migrates to a location where it will have three exposed edges.
	// 6. One pentagon has three exposed edges and migrates to a location where it will have three exposed edges.
	// 7. One pentagon has three exposed edges and migrates to a location where it will have two exposed edges.
	
	if (thirdC != NULLC) {
		opp_site = findSite(thirdC);
		if (opp_site != m_pah->m_siteList.end()) opp_site_bool = true;
	}
	if (thirdC2 != NULLC) {
		opp_site_second = findSite(thirdC2);
		if (opp_site_second != m_pah->m_siteList.end() && opp_site_second!=opp_site) opp_site_bool_second = true;
	}
	bool end_site_allowed = true;
	if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
		opp_site_after = findSite(thirdC_after);
		if (opp_site_after != m_pah->m_siteList.end()){
			opp_site_bool_after = true;
			int os_endtype = opp_site_after->type;
			if (os_endtype >= 200 && os_endtype <= 203) end_site_allowed = false;
			if (os_endtype == 101) end_site_allowed = false;
			if (os_endtype >= 600 && os_endtype <= 603) end_site_allowed = false;
			if (os_endtype >= 1000 && os_endtype <= 1003) end_site_allowed = false;
			if (os_endtype >= 500 && os_endtype <= 504) end_site_allowed = false;
			if (os_endtype >= 2000 && os_endtype <= 2205) end_site_allowed = false;
			if (os_endtype >= 2103 && os_endtype <= 2105) end_site_allowed = false;
			if (os_endtype >= 2204 && os_endtype <= 2205) end_site_allowed = false;
		}
	}
	if(!end_site_allowed){
		return;
	}

	//After this point we know that the process is accepted!
		//Check migration on the other side
	if (opp_site_bool && !opp_site_bool_second && !opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
	}else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
		addOppsiteR5Walker(opp_site_after,opp_site_after);
	}else if (opp_site_bool && opp_site_bool_second && !opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
	}else if (!opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		int jj = findWalker(opp_site_second);
		ii = remOppsiteR5Walker(ii, jj);
	}else if (!opp_site_bool && opp_site_bool_second && opp_site_bool_after){
		int jj = findWalker(opp_site_second);
		ii = remOppsiteR5Walker(ii, jj);
		addOppsiteR5Walker(opp_site_after,opp_site_after);
	}else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
			addOppsiteR5Walker(opp_site_after,opp_site_after);
	}else if (opp_site_bool && opp_site_bool_second && opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
		addOppsiteR5Walker(opp_site_after,opp_site_after);
	}

	Spointer S1_check, S2_check;
	if (b4){
		S1_check = moveIt(sFE2, -1);
		S2_check = moveIt(stt, +1);
	} else{
		S1_check = moveIt(stt, -1); 
		S2_check = moveIt(sFE2, 1); 
	}
	if (opp_site_bool && !opp_site_bool_second && !opp_site_bool_after) {
		if ( (int)opp_site->type >= 2100) {
			Spointer S1_opp_site = moveIt(opp_site, -1);
			Spointer S2_opp_site = moveIt(opp_site, +1);
			if (S1_opp_site->type==R5 || S2_opp_site->type==R5){
				updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
			}
			else updateSites(opp_site, opp_site->C1, opp_site->C2, -100);
		}
		else updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		Spointer S1_opp_site = moveIt(opp_site, -1);
		Spointer S2_opp_site = moveIt(opp_site, +1);
		updateCombinedSitesMigration(opp_site); updateCombinedSitesMigration(S1_opp_site); updateCombinedSitesMigration(S2_opp_site);
	}
	else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		if (opp_site_after != S1_check && opp_site_after != S2_check) {
			if ((int)opp_site_after->type >= 500 && (int)opp_site_after->type <= 700) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -400);
			else if ((int)opp_site_after->type >= 1000 && (int)opp_site_after->type <= 2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -800);
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			Spointer S1_opp_site = moveIt(opp_site_after, -1);
			Spointer S2_opp_site = moveIt(opp_site_after, +1);
			Spointer S3_opp_site = moveIt(opp_site_after, -2);
			Spointer S4_opp_site = moveIt(opp_site_after, +2);
			if ((int)opp_site_after->type < 2000){
				if (b4) updateSites(S1_opp_site, S1_opp_site->C1, S1_opp_site->C2, +500);
				else updateSites(S2_opp_site, S2_opp_site->C1, S2_opp_site->C2, +500);
				int kk = findWalker(opp_site_after);
				fixOppsiteR5Walker(kk);
			}
			updateCombinedSitesMigration(opp_site_after); updateCombinedSitesMigration(S1_opp_site); updateCombinedSitesMigration(S2_opp_site);
			updateCombinedSitesMigration(S3_opp_site); updateCombinedSitesMigration(S4_opp_site);
		}
		updateCombinedSitesMigration(opp_site);
	}
	else if (opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		if (opp_site != opp_site_second){
			updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
			//convSiteType(opp_site_second, opp_site_second->C1, opp_site_second->C2, (kmcSiteType)new_stype);
			updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, +1500);
			updateCombinedSitesMigration(opp_site);
			updateCombinedSitesMigration(opp_site_second);
		}
	}
	else if (!opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateCombinedSitesMigration(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, -1500);
		updateCombinedSitesMigration(opp_site_second);
		if (opp_site_after != S1_check && opp_site_after != S2_check) {
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
			updateCombinedSitesMigration(opp_site_after);
		}
	}
	else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		if (opp_site_after != S1_check && opp_site_after != S2_check) {
			if ( (int)opp_site_after->type >=2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +100);
			else updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			updateCombinedSitesMigration(opp_site_after);
		}
	}
	else if (opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		if ( (opp_site != opp_site_after) && (opp_site != opp_site_second) ){
			updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
			updateCombinedSitesMigration(opp_site);
			updateCombinedSitesMigration(opp_site_second);
			updateCombinedSitesMigration(opp_site_after);
		}
	}

	//Modify
	if (leave_edge){
		//The R5 will leave the current edge and move to other edge. 
		//This movement around the corner needs to correct the number of sites. 
		//Otherwise results would be incorrect.
		Spointer stt_coupled;
		if (b4) {
			stt_coupled = moveIt(stt,+1);
		}
		else {
			stt_coupled = moveIt(stt,-1);
		}
		updateSites(stt,stt->C1,stt->C2,-500);
		if((int)stt_coupled->type>2000) updateSites(stt_coupled,stt_coupled->C1,stt_coupled->C2,-101);
		else updateSites(stt_coupled,stt_coupled->C1,stt_coupled->C2,-501);
		Spointer start_site = std::get<0>(m_pah->m_R5walker_sites[ii]);
		Cpointer start_site_C1 = start_site->C1;
		Cpointer start_site_C2 = start_site->C2;

		proc_M5R_ACR5_around_corner(start_site,start_site_C1,start_site_C2,sFE2,b4,ii);
	}

	//After this point we know that the process is accepted!
	if (leave_corner){
		//The R5 will leave the current edge and move to other edge. 
		//This movement around the corner needs to correct the number of sites. 
		//Otherwise results would be incorrect.
		Spointer stt_coupled, start_site;
		if (b4) {
			stt_coupled = moveIt(stt,+1);
			start_site = std::get<0>(m_pah->m_R5walker_sites[ii]);
		}
		else {
			stt_coupled = moveIt(stt,-1);
			start_site = std::get<1>(m_pah->m_R5walker_sites[ii]);
		}
		updateSites(stt,stt->C1,stt->C2,-501);
		if((int)stt_coupled->type>2000) updateSites(stt_coupled,stt_coupled->C1,stt_coupled->C2,-100);
		else updateSites(stt_coupled,stt_coupled->C1,stt_coupled->C2,-500);
		Cpointer start_site_C1 = start_site->C1;
		Cpointer start_site_C2 = start_site->C2;

		proc_M5R_R5R6_out_of_corner(start_site,start_site_C1,start_site_C2,sFE2,b4,ii);
	}

	Spointer S3, S4, S5, S6;
	if (!leave_corner && !leave_edge){
		// edit sites. first identify the neighbouring sites of resulting RFE & R5
		if (b4) {
			S1 = moveIt(sFE2, -1);
			if ((int)sFE2->type == 0){ //sFE2 is a FE
				if ((int)S1->type == 0){ //S1 is a FE
					//R5 has moved to the edge and will now be free.
					//Since the FE2 will not move anymore it is probably a good idea to remove the walker.
					convSiteType(stt, stt->C1, stt->C2, RFE);
					convSiteType(sFE2, sFE2->C1, sFE2->C2, R5);
					convSiteType(checkR5_1, checkR5_1->C1, checkR5_1->C2, R5);
					updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
					Spointer site_perf = std::get<0>(m_pah->m_R5walker_sites[ii]);
					Spointer site_perf_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
					if(site_perf==site_perf_2) proc_M5R_ACR5_termination_toR5(site_perf,site_perf->C1,site_perf->C2,sFE2,b4,steps);
					else proc_M5R_R5R6_multiple_sites(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
					//sFE2 is deleted in the previous call so we must not call it again from here.
					//updateCombinedSitesMigration(stt); updateCombinedSitesMigration(checkR5_1); updateCombinedSitesMigration(checkR5_2);
					m_pah->m_R5walker_sites.erase(m_pah->m_R5walker_sites.begin()+ii);
					addR5internal(checkR5_1->C1,checkR5_1->C2);
					Spointer start_site, end_site, S1, S2, S3, S4, S5, S6;
					if (b4) {
						start_site = moveIt(stt,1);
						end_site = moveIt(start_site, -1);
						S1 = moveIt(end_site, -1);
						S2 = moveIt(start_site, +1);
						S3 = moveIt(S1, -1);
						S4 = moveIt(S2, 1);
						S5 = moveIt(S1, -2);
						S6 = moveIt(S2, 2);
						updateCombinedSitesMigration(start_site);
						updateCombinedSitesMigration(end_site); 
						updateCombinedSitesMigration(S1);
						updateCombinedSitesMigration(S2);
						updateCombinedSitesMigration(S3);
						updateCombinedSitesMigration(S4);
						updateCombinedSitesMigration(S5);
						updateCombinedSitesMigration(S6);
					}
					else {
						start_site = moveIt(stt,-1);
						end_site = moveIt(start_site, +1);
						S1 = moveIt(start_site, -1); 
						S1 = moveIt(start_site, -1); 
						S2 = moveIt(end_site, 1);
						S3 = moveIt(S1, -1);
						S4 = moveIt(S2, 1);
						S5 = moveIt(S1, -2);
						S6 = moveIt(S2, 2);
						updateCombinedSitesMigration(start_site);
						updateCombinedSitesMigration(end_site); 
						updateCombinedSitesMigration(S1);
						updateCombinedSitesMigration(S2);
						updateCombinedSitesMigration(S3);
						updateCombinedSitesMigration(S4);
						updateCombinedSitesMigration(S5);
						updateCombinedSitesMigration(S6);
					}
					return;
					//Need to check for opposite site logic before returning.			
				//Need to check for opposite site logic before returning.			
					//Need to check for opposite site logic before returning.			
				//Need to check for opposite site logic before returning.			
					//Need to check for opposite site logic before returning.			
				}
				else{
					convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
					int stype = S1->type;
					if (stype >= 2002 && stype <= 2204) {
						stype = stype + 100;
						convSiteType(S1, S1->C1, S1->C2, (kmcSiteType)stype);
					}
					else {
						stype = stype + 500;
						convSiteType(S1, S1->C1, S1->C2, (kmcSiteType)stype);
						//updateSites(S1, S1->C1, sFE2->C1, 5);
					}
				}
			}
			else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
				convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
			}
			else if ((int)sFE2->type == 2){ //sFE2 is a AC
				convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
			}
			else if ((int)sFE2->type == 3){ //sFE2 is a BY5
				convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
			}
			else if ((int)sFE2->type == 4){ //sFE2 is a BY6
				convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
			}
			else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
				updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
			}
			else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a R5 neighbour {
				updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
			}
			else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a R5 neighbour {
				updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
			}
			convSiteType(stt, stt->C2->C1, stt->C2, FE); //stt is originally the R5R6 site that will become the new FE site
			S2 = moveIt(stt, 1); // neighbour of stt
			int stype;
			if ( (int)S2->type >= 2103 && (int)S2->type <= 2205 ) stype = S2->type - 100;
			else stype = S2->type - 500;
			convSiteType(S2, S2->C1, S2->C2, (kmcSiteType)stype);
			/*if ((int)S2->type >= 10 && (int)S2->type <= 12){ //The site in the middle had 2R5s to each site, so lucky!
				updateSites(S2, stt->C2, S2->C2, -3); //convert the neighbour to its same version but next to an R5
			}
			else{
				int stype = S2->type - 55;
				convSiteType(S2, stt->C2, S2->C2, (kmcSiteType)stype);
				//updateSites(S2, stt->C2, S2->C2, -5); //convert the neighbour to its same version but next to an R5
			}
			*/
		}
		else {
			S2 = moveIt(sFE2, 1);
			if ((int)sFE2->type == 0){ //sFE2 is a FE
				if ((int)S2->type == 0){ //S2 is a FE
					convSiteType(stt, stt->C1, stt->C2, ZZ);
					convSiteType(sFE2, sFE2->C1, sFE2->C2, RFE);
					convSiteType(checkR5_1, checkR5_1->C1, checkR5_1->C2, R5);
					updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
				}
				else {
					convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
					int stype = S2->type;
					if (stype >= 2002 && stype <= 2204) {
						stype = stype + 100;
						convSiteType(S2, S2->C1, S2->C2, (kmcSiteType)stype);
					}
					else {
						stype = stype + 500;
						convSiteType(S2, S2->C1, S2->C2, (kmcSiteType)stype);
						//updateSites(S2, sFE2->C2, S2->C2, 5);
					}
				}
			}
			else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
				convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
			}
			else if ((int)sFE2->type == 2){ //sFE2 is a AC
				convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
			}
			else if ((int)sFE2->type == 3){ //sFE2 is a BY5
				convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
			}
			else if ((int)sFE2->type == 4){ //sFE2 is a BY6
				convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
			}
			else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
				updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
			}
			else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a BY5 {
				updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
			}
			else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a BY5 {
				updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
			}
			convSiteType(stt, stt->C1, stt->C1->C2, FE); //stt is originally the R5R6 site that will become the new FE site	
			S1 = moveIt(stt, -1); // neighbour of stt
			int stype;
			if ( (int)S1->type >= 2103 && (int)S1->type <= 2205 ) stype = S1->type - 100;
			else stype = S1->type - 500;
			convSiteType(S1, S1->C1, S1->C2, (kmcSiteType)stype);
			/*if ((int)S1->type >= 10 && (int)S1->type <= 12){ //The site in the middle had 2R5s to each site, so lucky!
				updateSites(S1, S1->C1, stt->C1, -3); //convert the neighbour to its same version but next to ONLY ONE R5
			}
			else{
				int stype = S1->type - 55;
				convSiteType(S1, S1->C1, stt->C1, (kmcSiteType)stype);
				//updateSites(S1, S1->C1, stt->C1, -5); //convert the neighbour to its same version but NOT next to an R5
			}*/
		}
		steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
		if((int)sFE2->type>=2003 && steps == 0) {
			if(b4) addR5internal(sFE2->C2->C1->C1,sFE2->C2->C1,true);
			else addR5internal(sFE2->C1->C2,sFE2->C1->C2->C2,true);
		}
	} //Migration without leaving corner or edge
	// update H atoms
	/*if (b4){
		updateA(S1->C2->C1, C_2, 'H');
	}
	else{
		updateA(C_1, S2->C1->C2, 'H');
	}*/
	//printStruct();
	// update combined sites for all sites involved and their neighbours
	// (excluding new FE sites, since their combined site type will still be None)
	if (!leave_corner && !leave_edge){
		if (b4){
			S3 = moveIt(S1, -1);
			S4 = moveIt(S2, 1);
			S5 = moveIt(S1, -2);
			S6 = moveIt(S2, 2);
		}
		else{
			S3 = moveIt(S1, -1);
			S4 = moveIt(S2, 1);
			S5 = moveIt(S1, -2);
			S6 = moveIt(S2, 2);
		}
		updateCombinedSitesMigration(stt); updateCombinedSitesMigration(sFE2); updateCombinedSitesMigration(S1); updateCombinedSitesMigration(S2);
		if (b4){
			updateCombinedSitesMigration(S3);
			updateCombinedSitesMigration(S4);
			updateCombinedSitesMigration(S5);
			updateCombinedSitesMigration(S6);
		}
		else{
			updateCombinedSitesMigration(S3);
			updateCombinedSitesMigration(S4); // neighbours
			updateCombinedSitesMigration(S5);
			updateCombinedSitesMigration(S6);
		}
	}
}

// ************************************************************
// ID23- ACR5 migration to double ZZ without modifying typespace
// ************************************************************
void PAHProcess::proc_M5R_ACR5_ZZ_ZZ_light(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	bool b4 = false;
	// Define a distribution that has two equally probably outcomes
	boost::bernoulli_distribution<> choiceDistrib;
	// Now build an object that will generate a sample using rng
	boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
	b4 = choiceGenerator(); // if FE3 on both sides, choose a random one
	/*if (moveIt(stt, -1)->comb == FE2 || moveIt(stt, 1)->comb == FE2) {
		if (moveIt(stt, -1)->comb == FE2 && moveIt(stt, 1)->comb == FE2) {
			// Define a distribution that has two equally probably outcomes
			boost::bernoulli_distribution<> choiceDistrib;
			// Now build an object that will generate a sample using rng
			boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
			b4 = choiceGenerator(); // if FE3 on both sides, choose a random one
		}
		if (moveIt(stt, -1)->comb == FE2){
			// Define a distribution that has two thirds to one third probably outcome
			boost::bernoulli_distribution<> choiceDistrib(1.0/3.0);
			// Now build an object that will generate a sample using rng
			boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
			b4 = choiceGenerator(); 
		}
	}
	else{
		// Define a distribution that has two equally probably outcomes
		boost::bernoulli_distribution<> choiceDistrib;
		// Now build an object that will generate a sample using rng
		boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
		b4 = choiceGenerator(); // if FE3 on both sides, choose a random one
	}*/

	int steps=-99999;
	unsigned int ii = 0;
	for (ii = 0;ii != m_pah->m_R5walker_sites.size();ii++){
		Spointer migr_site_start = std::get<0>(m_pah->m_R5walker_sites[ii]);
		Spointer migr_site_start_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
		Spointer site_check = moveIt(migr_site_start, std::get<2>(m_pah->m_R5walker_sites[ii]));
		Spointer site_check_2 = moveIt(migr_site_start_2, std::get<2>(m_pah->m_R5walker_sites[ii]));
		if (site_check == stt && site_check_2 == stt){
			steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
			break;
		}
	}
	if (steps==-99999){
		std::cout << "Error on R5 migration in proc_M5R_ACR5_ZZ_ZZ_light. Walker not found." << std::endl;
		printSitesMigration();
	}

	//Check for unsupported sites. This section heavily assumes that the Isolated Pentagon Rule is valid.
	/*bool allowed = true;
	if ((int)sFE2->type == 0 && (int)checkR5_1->type == 0 && (int)checkR5_2->type == 0) allowed = false; // The result would be an indene, not supported. YET!
	if ((int)sFE2->type == 101 || (int)sFE2->type == 501 || (int)sFE2->type == 2002) allowed = false; // This would violate the IPR.
	if ((int)sFE2->type == 0){
		if ((int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 || (int)checkR5_1->type == 100) allowed = false;
		//if ((int)checkR5_1->type >= 501 && (int)checkR5_1->type <= 65) return;
		if ((int)checkR5_1->type >= 1002 && (int)checkR5_1->type <= 1004) allowed = false;
		if ((int)checkR5_1->type >= 2002 && (int)checkR5_1->type <= 2204) allowed = false;
		if ((int)checkR5_1->type >= 2204 && (int)checkR5_1->type <= 2205) allowed = false;
		if ((int)checkR5_1->type == 0){
			//if ((int)checkR5_2->type >= 501 && (int)checkR5_2->type <= 65) return;
			if ((int)checkR5_2->type == 101 || (int)checkR5_2->type == 501) allowed = false;
			if ((int)checkR5_2->type >= 1002 && (int)checkR5_2->type <= 1004) allowed = false;
			if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2205) allowed = false;
		}
	}
	if (CRem_next->bridge) allowed = false;
	if (CRem_next->C2->bridge) allowed = false; if (CRem_next->C1->bridge) allowed = false;
	if (CRem_next->C2->C2->bridge) allowed = false; if (CRem_next->C1->C1->bridge) allowed = false;*/

	Spointer sFE2, checkR5_1, checkR5_2;
	Cpointer CRem, CFE, CRem_next, CR5_otherside_1, CR5_otherside_2;
	if (b4) {
		sFE2 = moveIt(stt, -1);
		checkR5_1 = moveIt(stt, -2);
		checkR5_2 = moveIt(stt, -3);
		CFE = C_1->C2; 
		CRem = sFE2->C2; 
		if (steps == 0){
			CR5_otherside_1 = C_1->C2->C2;
			CR5_otherside_2 = C_1->C2;
			CRem_next = CRem->C1;
		} else {
			CR5_otherside_1 = C_1->C2;
			CR5_otherside_2 = C_1->C1;
			if (CR5_otherside_2->C1->A=='H') CRem_next = CR5_otherside_2->C1->C1;
			else CRem_next = CR5_otherside_2->C1;
		} 
	}
	else {
		sFE2 = moveIt(stt, 1);
		checkR5_1 = moveIt(stt, 2);
		checkR5_2 = moveIt(stt, 3);
		CFE = C_2->C1->C1; //Not modified but pointer created
		CRem = sFE2->C1; //Not modified but pointer created
		if (steps == 0){
			CR5_otherside_1 = C_2->C1->C1;
			CR5_otherside_2 = C_2->C1;
			CRem_next = CRem->C2;
		} else {
			CR5_otherside_1 = C_2->C1;
			CR5_otherside_2 = C_2->C2;
			if (CR5_otherside_2->C2->A=='H') CRem_next = CR5_otherside_2->C2->C2;
			else CRem_next = CR5_otherside_2->C2;
		} 
	}
	
	// check if ACR5 has an opposite site.
	Spointer opp_site, opp_site_second, opp_site_after;
	bool opp_site_bool = false; bool opp_site_bool_second = false; bool opp_site_bool_after = false;
	Cpointer thirdC = findThirdC(CR5_otherside_1);
	Cpointer thirdC2 = findThirdC(CR5_otherside_2);
	Cpointer thirdC_after = findThirdC(CRem_next);
	// Seven cases:
	// 1. One pentagon has one exposed edge  and migrates to a location where it will have one exposed edge. Normal migration.
	// 2. One pentagon has one exposed edge  and migrates to a location where it will have two exposed edges. 
	// 3. One pentagon has two exposed edges and migrates to a location where it will have two exposed edges. 
	// 4. One pentagon has two exposed edges and migrates to a location where it will have one exposed edge.
	// 5. One pentagon has two exposed edges and migrates to a location where it will have three exposed edges.
	// 6. One pentagon has three exposed edges and migrates to a location where it will have three exposed edges.
	// 7. One pentagon has three exposed edges and migrates to a location where it will have two exposed edges.
	
	if (thirdC != NULLC) {
		opp_site = findSite(thirdC);
		if (opp_site != m_pah->m_siteList.end()) opp_site_bool = true;
	}
	if (thirdC2 != NULLC) {
		opp_site_second = findSite(thirdC2);
		if (opp_site_second != m_pah->m_siteList.end() && opp_site_second!=opp_site) opp_site_bool_second = true;
	}
	bool allowed = true;
	if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
		opp_site_after = findSite(thirdC_after);
		if (opp_site_after != m_pah->m_siteList.end()) {
			opp_site_bool_after = true;
			int os_endtype = opp_site_after->type;
			if (os_endtype >= 200 && os_endtype <= 203) allowed = false;
			if (os_endtype == 101) allowed = false;
			if (os_endtype >= 600 && os_endtype <= 603) allowed = false;
			if (os_endtype >= 1000 && os_endtype <= 1003) allowed = false;
			if (os_endtype >= 500 && os_endtype <= 504) allowed = false;
			if (os_endtype >= 2000 && os_endtype <= 2205) allowed = false;
			if (os_endtype >= 2103 && os_endtype <= 2105) allowed = false;
			if (os_endtype >= 2204 && os_endtype <= 2205) allowed = false;
		}
	}
	if (!allowed){
		return;
	}

	//After this point we know that the process is accepted!
		//Check migration on the other side
	if (opp_site_bool && !opp_site_bool_second && !opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
	}else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
		addOppsiteR5Walker(opp_site_after,opp_site_after);
	}else if (opp_site_bool && opp_site_bool_second && !opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
	}else if (!opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		int jj = findWalker(opp_site_second);
		ii = remOppsiteR5Walker(ii, jj);
	}else if (!opp_site_bool && opp_site_bool_second && opp_site_bool_after){
		int jj = findWalker(opp_site_second);
		ii = remOppsiteR5Walker(ii, jj);
		addOppsiteR5Walker(opp_site_after,opp_site_after);
	}else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
			addOppsiteR5Walker(opp_site_after,opp_site_after);
	}else if (opp_site_bool && opp_site_bool_second && opp_site_bool_after){
		int jj = findWalker(opp_site);
		ii = remOppsiteR5Walker(ii, jj);
		addOppsiteR5Walker(opp_site_after,opp_site_after);
	}
	if (b4) {
		std::get<2>(m_pah->m_R5walker_sites[ii])--;
	}
	else {
		std::get<2>(m_pah->m_R5walker_sites[ii])++;
	}
	checkR5Walkers(ii);
	bool opp_site_logic = false;
	
	if ((int)checkR5_1->type == 0 && (int)sFE2->type == 0){
		//R5 has moved to the edge and will now be free.
		//Since the FE2 will not move anymore it is probably a good idea to remove the walker.
		convSiteType(stt, stt->C1, stt->C2, RFE);
		convSiteType(sFE2, sFE2->C1, sFE2->C2, R5);
		convSiteType(checkR5_1, checkR5_1->C1, checkR5_1->C2, R5);
		updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
		Spointer site_perf = std::get<0>(m_pah->m_R5walker_sites[ii]);
		Spointer site_perf_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
		if(site_perf==site_perf_2) proc_M5R_ACR5_termination_toR5(site_perf,site_perf->C1,site_perf->C2,sFE2,b4,steps);
		else proc_M5R_R5R6_multiple_sites(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		//sFE2 is deleted in the previous call so we must not call it again from here.
		//updateCombinedSitesMigration(stt); updateCombinedSitesMigration(checkR5_1); updateCombinedSitesMigration(checkR5_2);
		m_pah->m_R5walker_sites.erase(m_pah->m_R5walker_sites.begin()+ii);
		addR5internal(checkR5_1->C1,checkR5_1->C2);
		Spointer start_site, end_site, S1, S2, S3, S4, S5, S6;
		if (b4) {
			start_site = moveIt(stt,1);
			end_site = moveIt(start_site, -1);
			S1 = moveIt(end_site, -1);
			S2 = moveIt(start_site, +1);
			S3 = moveIt(S1, -1);
			S4 = moveIt(S2, 1);
			S5 = moveIt(S1, -2);
			S6 = moveIt(S2, 2);
			updateCombinedSitesMigration(start_site);
			updateCombinedSitesMigration(end_site); 
			updateCombinedSitesMigration(S1);
			updateCombinedSitesMigration(S3);
			updateCombinedSitesMigration(S2);
			updateCombinedSitesMigration(S4);
			updateCombinedSitesMigration(S5);
			updateCombinedSitesMigration(S6);
		}
		else {
			start_site = moveIt(stt,-1);
			end_site = moveIt(start_site, +1);
			S1 = moveIt(start_site, -1); 
			S2 = moveIt(end_site, 1);
			S3 = moveIt(S1, -1);
			S4 = moveIt(S2, 1);
			S5 = moveIt(S1, -2);
			S6 = moveIt(S2, 2);
			updateCombinedSitesMigration(start_site);
			updateCombinedSitesMigration(end_site); 
			updateCombinedSitesMigration(S2);
			updateCombinedSitesMigration(S4); // neighbours
			updateCombinedSitesMigration(S1);
			updateCombinedSitesMigration(S3);
			updateCombinedSitesMigration(S5);
			updateCombinedSitesMigration(S6);
		}
		//Need to check for opposite site logic before returning.
		opp_site_logic = true;
	}

	Spointer S1_check, S2_check;
	if (b4){
		S1_check = moveIt(sFE2, -1);
		S2_check = moveIt(stt, +1);
	} else{
		S1_check = moveIt(stt, -1); 
		S2_check = moveIt(sFE2, 1); 
	}
	if (opp_site_bool && !opp_site_bool_second && !opp_site_bool_after) {
		if ( (int)opp_site->type >= 2100) {
			Spointer S1_opp_site = moveIt(opp_site, -1);
			Spointer S2_opp_site = moveIt(opp_site, +1);
			if (S1_opp_site->type==R5 || S2_opp_site->type==R5){
				updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
			}
			else updateSites(opp_site, opp_site->C1, opp_site->C2, -100);
		}
		else updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		Spointer S1_opp_site = moveIt(opp_site, -1);
		Spointer S2_opp_site = moveIt(opp_site, +1);
		updateCombinedSitesMigration(opp_site); updateCombinedSitesMigration(S1_opp_site);  updateCombinedSitesMigration(S2_opp_site); 
	}
	else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		if (opp_site_after != S1_check && opp_site_after != S2_check) {
			if ((int)opp_site_after->type >= 500 && (int)opp_site_after->type <= 700) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -400);
			else if ((int)opp_site_after->type >= 1000 && (int)opp_site_after->type <= 2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -800);
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			Spointer S1_opp_site = moveIt(opp_site_after, -1);
			Spointer S2_opp_site = moveIt(opp_site_after, +1);
			updateCombinedSitesMigration(opp_site_after); updateCombinedSitesMigration(S1_opp_site);  updateCombinedSitesMigration(S2_opp_site); 
		}
		updateCombinedSitesMigration(opp_site);
	}
	else if (opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
		//convSiteType(opp_site_second, opp_site_second->C1, opp_site_second->C2, (kmcSiteType)new_stype);
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, +1500);
		updateCombinedSitesMigration(opp_site);
		updateCombinedSitesMigration(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && !opp_site_bool_after) {
		updateCombinedSitesMigration(opp_site_second);
	}
	else if (!opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, -1500);
		updateCombinedSitesMigration(opp_site_second);
		if (opp_site_after != S1_check && opp_site_after != S2_check) {
			updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
			updateCombinedSitesMigration(opp_site_after);
		}
	}
	else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		if (opp_site_after != S1_check && opp_site_after != S2_check) {
			if ( (int)opp_site_after->type >=2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +100);
			else updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
			updateCombinedSitesMigration(opp_site_after);
		}
	}
	else if (opp_site_bool && opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
		updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
		updateCombinedSitesMigration(opp_site);
		updateCombinedSitesMigration(opp_site_second);
		updateCombinedSitesMigration(opp_site_after);
	}
	if(opp_site_logic) return;

	// edit sites. first identify the neighbouring sites of resulting RFE & R5
	Spointer S1, S2, S3, S4, S5, S6;
	if (b4) {
		S1 = moveIt(sFE2, -1);
		S2 = moveIt(stt, +1);
		if ((int)sFE2->type == 0 && (int)S1->type == 0){ //sFE2 is a FE
			//This means that the pentagon has migrated to the edge. This is handled above.
			//convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			//updateSites(S1, S1->C1, sFE2->C1, 500);
		}
		else if ((int)sFE2->type == 0 ){ 
			//This means that the pentagon has migrated to the edge but will have a carbon out of the structure.
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			updateSites(S1, S1->C1, S1->C2, 500);
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a R5 neighbour {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
		}
		convSiteType(stt, stt->C1, stt->C2, ZZ);
	}
	else {
		S1 = moveIt(stt, -1); 
		S2 = moveIt(sFE2, 1); 
		if ((int)sFE2->type == 0 && (int)S2->type == 0){ //sFE2 is a FE
			//This means that the pentagon has migrated to the edge. This is handled above.
			//convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			//updateSites(S2, sFE2->C2, S2->C2, 500);
		}
		else if ((int)sFE2->type == 0){
			//This means that the pentagon has migrated to the edge but will have a carbon out of the structure.
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			updateSites(S2, S2->C1, S2->C2, 500);
		}
		else if ((int)sFE2->type == 1){ //sFE2 is a ZZ
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACR5);
		}
		else if ((int)sFE2->type == 2){ //sFE2 is a AC
			convSiteType(sFE2, sFE2->C1, sFE2->C2, FEACR5);
		}
		else if ((int)sFE2->type == 3){ //sFE2 is a BY5
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ZZACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type == 4){ //sFE2 is a BY6
			convSiteType(sFE2, sFE2->C1, sFE2->C2, ACACR5); //BY6 with an R5 inside is treated same as a BY6
		}
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2115){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
		}
		else if ((int)sFE2->type >= 102 && (int)sFE2->type <= 104){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +2001);
		}
		else if ((int)sFE2->type >= 502 && (int)sFE2->type <= 504){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +1601);
		}
		convSiteType(stt, stt->C1, stt->C2, ZZ);
	}
	steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
	if((int)sFE2->type>=2003 && steps == 0) {
		if(b4) addR5internal(sFE2->C2->C1->C1,sFE2->C2->C1,true);
		else addR5internal(sFE2->C1->C2,sFE2->C1->C2->C2,true);
	}
	// update H atoms
	/*if (b4){
		updateA(S1->C2->C1, C_2, 'H');
	}
	else{
		updateA(C_1, S2->C1->C2, 'H');
	}*/
	
	//printStruct();
	// update combined sites for all sites involved and their neighbours
	// (excluding new FE sites, since their combined site type will still be None)
	if (b4){
		S3 = moveIt(S1, -1);
		S4 = moveIt(S2, 1);
		S5 = moveIt(S1, -2);
		S6 = moveIt(S2, 2);
	}
	else{
		S3 = moveIt(S1, -1);
		S4 = moveIt(S2, 1);
		S5 = moveIt(S1, -2);
		S6 = moveIt(S2, 2);
	}
	updateCombinedSitesMigration(stt); updateCombinedSitesMigration(sFE2); 
	if (b4){
		updateCombinedSitesMigration(S1);
		updateCombinedSitesMigration(S3);
		updateCombinedSitesMigration(S2);
		updateCombinedSitesMigration(S4);
		updateCombinedSitesMigration(S5);
		updateCombinedSitesMigration(S6); // neighbours
		if(S2->comb==FE2){
			Spointer S2_2 = moveIt(S2,+2);
			Spointer S2_3 = moveIt(S2,+3);
			updateCombinedSitesMigration(S2_2);
			updateCombinedSitesMigration(S2_3);
		}
	}
	else{
		updateCombinedSitesMigration(S2);
		updateCombinedSitesMigration(S4); // neighbours
		updateCombinedSitesMigration(S1);
		updateCombinedSitesMigration(S3);
		updateCombinedSitesMigration(S5);
		updateCombinedSitesMigration(S6); // neighbours
		if(S1->comb==FE2){
			Spointer S1_2 = moveIt(S1,-2);
			Spointer S1_3 = moveIt(S1,-3);
			updateCombinedSitesMigration(S1_2);
			updateCombinedSitesMigration(S1_3);
		}
	}

}

// ************************************************************
// ID71 - R5R7 edge healing
// ************************************************************
void PAHProcess::proc_MR5R7_edge(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	//printStruct();
	bool b4 = false;
	Cpointer CRem, CFE;
	if ( (isR7internal(stt->C1,stt->C1->C2) ) && (isR7internal(stt->C2->C1,stt->C2) )){
		//R5R7 to both sides.
		// Define a distribution that has two equally probably outcomes
		boost::bernoulli_distribution<> choiceDistrib;
		// Now build an object that will generate a sample using rng
		boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
		if(choiceGenerator()) {
			b4 = true;
			CRem = C_1;
			CFE = C_1->C2;
		} else {
			b4 = false;
			CRem = C_2;
			CFE = C_2->C1->C1;
		}
	}
	else if (isR7internal(stt->C1,stt->C1->C2) ){
		b4 = true;
		CRem = C_1;
		CFE = C_1->C2;
	}
	else {
		b4 = false;
		CRem = C_2;
		CFE = C_2->C1->C1;
	}

	// check if ACR5 has an opposite site.
	Spointer opp_site, opp_site_second;
	bool opp_site_bool = false; bool opp_site_bool_second = false;
	Cpointer thirdC = findThirdC(C_1->C2);
	Cpointer thirdC2 = findThirdC(C_2->C1);
	if (thirdC != NULLC) {
		opp_site = findSite(thirdC);
		if (opp_site != m_pah->m_siteList.end()) opp_site_bool = true;
	}
	if (thirdC2 != NULLC) {
		opp_site_second = findSite(thirdC2);
		if (opp_site_second != m_pah->m_siteList.end() && opp_site_second!=opp_site) opp_site_bool_second = true;
	}

	//Add a new carbon between current R5 carbons of ACR5
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	double dist2;
	if (R5_dist < 2.5) dist2 = 1.4;
	else dist2 = R5_dist / 2.7 * 1.5;
	double theta = asin(R5_dist/2.0/dist2);
	double magn = dist2 * cos(theta);
	cpair R5dir = get_vector(CFE->coords,CFE->C2->coords);
	cpair normvec;
	if (CFE->C2->C2->A=='C') normvec = invert_vector(norm_vector(CFE->coords, CFE->C2->coords, CFE->C2->C2->coords));
	else normvec = norm_vector(CFE->coords, CFE->C2->coords, CFE->C2->C2->coords);
	cpair crossvec = cross_vector(R5dir, normvec);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + magn * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + magn * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ magn * std::get<2>(crossvec));
	cpair Cnewdir = scale_vector(resultantvec);
	//cpair Cnewdir = get_vector(C_2->C1->coords,C_2->coords);
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, 1.4);
	updateA(Cnew, 'H', crossvec);
	
	//Remove carbon from end site 
	removeC(CRem, false);
	OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol);
	passbackPAH(mol);

	if (b4) {
		if( (int)stt->type>2000){
			//An ACR5 site or similar.
			updateSites(stt, Cnew, C_2, -2001);
			Spointer S_other = moveIt(stt,-1);
			updateSites(S_other, S_other->C1, Cnew, +1);
		}
		else{
			//An R5R6 site or similar
			convSiteType(stt, Cnew, C_2, FE);
			Spointer S_otherR7 = moveIt(stt,-1);
			updateSites(S_otherR7, S_otherR7->C1, Cnew, +1);
			Spointer S_otherR5 = moveIt(stt,+1);
			if((int)stt->type>2000) updateSites(S_otherR5, S_otherR5->C1, S_otherR5->C2, -100);
			else updateSites(S_otherR5, S_otherR5->C1, S_otherR5->C2, -500);
		}
	}
	else {
		if( (int)stt->type>2000){
			//An ACR5 site or similar.
			updateSites(stt, C_1, Cnew, -2001);
			Spointer S_other = moveIt(stt,+1);
			updateSites(S_other, Cnew, S_other->C2, +1);
		}
		else{
			//An R5R6 site or similar
			convSiteType(stt, C_1, Cnew, FE);
			Spointer S_otherR7 = moveIt(stt,+1);
			updateSites(S_otherR7, Cnew, S_otherR7->C2, +1);
			Spointer S_otherR5 = moveIt(stt,-1);
			if((int)stt->type>2000) updateSites(S_otherR5, S_otherR5->C1, S_otherR5->C2, -100);
			else updateSites(S_otherR5, S_otherR5->C1, S_otherR5->C2, -500);
		}
	}

	Spointer S1 = moveIt(stt,-1); 
	Spointer S2 = moveIt(stt,+1); 
	Spointer S3 = moveIt(stt,-2); 
	Spointer S4 = moveIt(stt,+2); 
	Spointer S5 = moveIt(stt,-3); 
	Spointer S6 = moveIt(stt,+3); 
	updateCombinedSites(stt); updateCombinedSites(S1); updateCombinedSites(S2); updateCombinedSites(S3); updateCombinedSites(S4); updateCombinedSites(S5); updateCombinedSites(S6);

	//Adjust opposite sites for starting and ending position.
	if (opp_site_bool && opp_site_bool_second) {
		if (opp_site == opp_site_second) opp_site_bool_second = false;
	}
	if (opp_site_bool) {
		if ( (int)opp_site->type >= 2100) {
			Spointer S1_opp_site = moveIt(opp_site, -1);
			Spointer S2_opp_site = moveIt(opp_site, +1);
			if (S1_opp_site->type==R5 || S2_opp_site->type==R5){
				updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
			}
			else updateSites(opp_site, opp_site->C1, opp_site->C2, -100);
		}
		else if  ((int)opp_site->type >= 2000) updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		else updateSites(opp_site, opp_site->C1, opp_site->C2, -500);
		updateCombinedSites(opp_site);
	}
	if (opp_site_bool_second) {
		if ( (int)opp_site_second->type >= 2100) {
			Spointer S1_opp_site = moveIt(opp_site_second, -1);
			Spointer S2_opp_site = moveIt(opp_site_second, +1);
			if (S1_opp_site->type==R5 || S2_opp_site->type==R5){
				updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, -2000);
			}
			else updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, -100);
		}
		else updateSites(opp_site_second, opp_site_second->C1, opp_site_second->C2, -2000);
		updateCombinedSites(opp_site_second);
	}
	//Ring transformation
	m_pah->m_rings5_Lone--;
	m_pah->m_rings7_Embedded--;
	m_pah->m_rings++;
	m_pah->m_rings++;
}


size_t PAHProcess::SiteListSize() const {
    return m_pah->m_siteList.size();
}

size_t PAHProcess::CarbonListSize() const {
    return m_pah->m_carbonList.size();
}

std::list<Site>& PAHProcess::SiteList() const {
    return m_pah->m_siteList;
}

//! obtains a vector of the PAH site list
std::vector<kmcSiteType> PAHProcess::SiteVector() const {
    std::vector<kmcSiteType> temp;
    for(Spointer i=SiteList().begin(); i!= SiteList().end(); ++i) {
        temp.push_back((*i).type);
    }
    return temp;
};

//! obtains a vector of the PAH site list
std::vector<std::string> PAHProcess::SiteVectorString() const {
    std::vector<std::string> temp;
    for(Spointer i=SiteList().begin(); i!= SiteList().end(); ++i) {
        temp.push_back(kmcSiteName((*i).type));
    }
    return temp;
};

//! obtains a vector of tuples from the PAH site list
std::vector<std::tuple<int, int, cpair, cpair>> PAHProcess::SiteVector_clone() const {
    std::vector<std::tuple<int, int, cpair, cpair>> temp;
    for(Spointer i=SiteList().begin(); i!= SiteList().end(); ++i) {
		std::tuple<int, int, cpair, cpair> i_site = std::make_tuple((int)(*i).type, (int)(*i).comb, (*i).C1->coords, (*i).C2->coords);
        temp.push_back(i_site);
    }
    return temp;
};

//! obtains a map of tuples from the PAH site map
std::map<int, std::vector<std::tuple<int, int, cpair, cpair>>> PAHProcess::SiteMap_clone() const {
    std::map<int, std::vector<std::tuple<int, int, cpair, cpair>>> temp;
	std::map<kmcSiteType, svector>::const_iterator map_it;
	for(map_it=m_pah->m_siteMap.begin(); map_it!=m_pah->m_siteMap.end(); map_it++){
		int sitetype = map_it->first;

		std::vector<Spointer> site_vector = map_it->second;
		std::vector<std::tuple<int, int, cpair, cpair>> site_vector_temp;

		for(unsigned int site_it=0; site_it!= site_vector.size(); site_it++) {
			Spointer S1 = site_vector[site_it];
			std::tuple<int, int, cpair, cpair> i_site = std::make_tuple((int)S1->type, (int)S1->comb, S1->C1->coords, S1->C2->coords);
        	site_vector_temp.push_back(i_site);
    	}
		temp[sitetype] = site_vector_temp;
	}
    return temp;
};

//! obtains a vector of tuples from the PAH Ccontainer
std::vector<std::tuple<cpair, int, int, cpair>> PAHProcess::EdgeCarbonVector_clone() const {
	std::vector<std::tuple<cpair, int, int, cpair>> temp;
	Cpointer Cnow = m_pah->m_cfirst;
	Cpointer Cprev = m_pah->m_cfirst;
	do{
		int brdige_val, H_val;
		if (Cnow-> bridge) brdige_val = 1;
		else brdige_val = 0;
		if (Cnow->A != 'C'){
			if (Cnow->A == 'H') H_val = 1;
			if (Cnow->A == 'M') H_val = 2;
		}
		else H_val = 0;
		std::tuple<cpair, int, int, cpair> i_carb = std::make_tuple(Cnow->coords, brdige_val, H_val, Cnow->growth_vector);
		temp.push_back(i_carb);
		if (Cnow-> bridge && Cprev != Cnow->C3){
			Cprev = Cnow;
			Cnow = Cnow->C3;
		}
		else{
			Cprev = Cnow;
			Cnow = Cnow->C2;
		}
	}while (temp.size() < m_pah->m_cpositions.size() + numberOfBridges()*2);
	if (Cnow!=m_pah->m_cfirst){
		std::cout << "Warning. EdgeCarbonVector_clone() did not finish at m_pah->m_cfirst." <<std::endl;
	}
	return temp;
}

//! obtains a vector of the carbons per site in PAH site list
std::vector<int> PAHProcess::SiteIntVector() const {
    std::vector<int> temp;
    for(Spointer i=SiteList().begin(); i!= SiteList().end(); ++i) {
		int site_number = (int)((*i).type);
		int carbon_number;
		if (site_number == 9999 || (*i).type == None){
			Cpointer Cspiral_1 = i->C1;
			Cpointer Cspiral_2 = i->C2;
			Cpointer Cspiral_bridge = NULLC;
			int carbon_spiral_counter = 0;
			do{
				carbon_spiral_counter++;
				if (Cspiral_1->bridge && Cspiral_1->C3 != Cspiral_bridge){
					Cspiral_bridge = Cspiral_1;
					Cspiral_1 = Cspiral_1->C3;
				}
				else {
					Cspiral_1 = Cspiral_1->C2;
				}
			}while (Cspiral_1 != Cspiral_2);
			carbon_number = carbon_spiral_counter;
		}
		else carbon_number = site_number % 10 + 1;
		temp.push_back(carbon_number);
    }
    return temp;
};

//! obtains a string containing the PAH site list
std::string PAHProcess::SiteString(char delimiter) const {
    std::ostringstream temp;
    std::vector<kmcSiteType> vec(SiteVector());
    temp << kmcSiteName(vec[0]);
    for(size_t i=1; i<vec.size(); i++) {
        temp << delimiter;
        temp << kmcSiteName(vec[i]);
    }
    return temp.str();
};

//! Returns how many carbons are in a site stt
int PAHProcess::SiteSize(Spointer& stt) const{
	Cpointer Cnow = stt->C1;
	Cpointer Cprev = Cnow;
	Cpointer Cend = stt->C2;
	int counter = 1;
	do{
		if (Cnow->bridge && !Cprev->bridge) {
			Cprev = Cnow;
			Cnow = Cnow->C3;
		}
		else {
			Cprev = Cnow;
			Cnow = Cnow->C2;
		}
		counter +=1;
	} while(Cnow != Cend);
	return counter;
}

//! Returns true if a site has the right number of carbons
bool PAHProcess::SiteRightSize(Spointer& stt) const{
	int stype = (int)stt->type;
	if (stype == 9999) return true; // An spiral can have any numberof carbons.
	if (stt->type==None) return true; // An error type can have any number of carbons.

	int counter = SiteSize(stt);
	stype = stype%10 + 2;
	if (counter != stype) return false;
	return true;
}

//! Called before migration process starts. Loops through the random walker sites and appends them.
void PAHProcess::startMigrationProcess(){
	std::vector<std::tuple<Spointer,Spointer,int>> migr_sites;
	std::vector<Spointer> migr_sites_appended; // Sites already included as walkers
	if(m_pah->m_siteMap[MIGR].size()>0){
		for(unsigned int ii=0;ii!=m_pah->m_siteMap[MIGR].size();ii++){
			Spointer st = (m_pah->m_siteMap[MIGR])[ii];
			std::tuple<Spointer,Spointer,int> migr_site_ii= std::make_tuple(st, st, 0);
			migr_sites.push_back(migr_site_ii);
			migr_sites_appended.push_back(st);
			cpair R5coords;
			if (st->type == ACR5) {
				R5coords = findR5internal(st->C1->C2, st->C2->C1);
				Cpointer C_check_other_side = st->C1->C2;
				Cpointer C_check_other_side2 = st->C2->C1;
				Cpointer C_other_side = findThirdC(C_check_other_side);
				Cpointer C_other_side2 = findThirdC(C_check_other_side2);
				if (C_other_side!=NULLC || C_other_side2!=NULLC) {
					Spointer opp_site;
					if (C_other_side!=NULLC) opp_site = findSite(C_other_side);
					else opp_site = findSite(C_other_side2);
					std::vector<Spointer>::iterator it;
					it = std::find(migr_sites_appended.begin(),migr_sites_appended.end(),opp_site);
					if(it==migr_sites_appended.end()) m_pah->m_R5loc.push_back(R5coords);
				}
			}
			/*else{
				if ( isR5internal(st->C1->C2, st->C1->C2->C2,false) || isR5internal(st->C1->C2, st->C1->C2->C2,true) ) {
					R5coords = findR5internal(st->C1->C2, st->C1->C2->C2);
				}
				else if ( isR5internal(st->C2->C1->C1, st->C2->C1,false) || isR5internal(st->C2->C1->C1, st->C2->C1,true) ) {
					R5coords = findR5internal(st->C2->C1->C1, st->C2->C1);
				}
				else {
					//R5 not found
					std::cout << "ERROR. R5 not found for " << st->type << " site in PAHProcess::startMigrationProcess." << std::endl;
				}
			}*/
		}
	}
	if(m_pah->m_siteMap[MIGR2].size()>0){
		for(unsigned int ii=0;ii!=m_pah->m_siteMap[MIGR2].size();ii++){
			Spointer st = (m_pah->m_siteMap[MIGR2])[ii];
			std::tuple<Spointer,Spointer,int> migr_site_ii= std::make_tuple(st, st, 0);
			migr_sites.push_back(migr_site_ii);
			migr_sites_appended.push_back(st);
			cpair R5coords;
			if (st->type == ACR5) {
				R5coords = findR5internal(st->C1->C2, st->C2->C1);
				Cpointer C_check_other_side = st->C1->C2;
				Cpointer C_check_other_side2 = st->C2->C1;
				Cpointer C_other_side = findThirdC(C_check_other_side);
				Cpointer C_other_side2 = findThirdC(C_check_other_side2);
				if (C_other_side!=NULLC || C_other_side2!=NULLC) {
					Spointer opp_site;
					if (C_other_side!=NULLC) opp_site = findSite(C_other_side);
					else opp_site = findSite(C_other_side2);
					std::vector<Spointer>::iterator it;
					it = std::find(migr_sites_appended.begin(),migr_sites_appended.end(),opp_site);
					if(it==migr_sites_appended.end()) m_pah->m_R5loc.push_back(R5coords);
				}
			}
			else{
				if ( isR5internal(st->C1->C2, st->C1->C2->C2) ) {
					R5coords = findR5internal(st->C1->C2, st->C1->C2->C2);
				}
				else if ( isR5internal(st->C2->C1->C1, st->C2->C1) ) {
					R5coords = findR5internal(st->C2->C1->C1, st->C2->C1);
				}
				else {
					//R5 not found
					std::cout << "ERROR. R5 not found for " << st->type << " site in PAHProcess::startMigrationProcess." << std::endl;
				}
			}
		}
	}
	if(m_pah->m_siteMap[R5R6_MIGR].size()>0){
		for(unsigned int ii=0;ii!=m_pah->m_siteMap[R5R6_MIGR].size();ii++){
			Spointer current_site = (m_pah->m_siteMap[R5R6_MIGR])[ii];
			//Spointer next_site;
			int coupled_site_dir = coupledSiteDirection(current_site);
			Spointer coupled_site = moveIt(current_site,coupled_site_dir);
			std::vector<Spointer>::iterator it, it2;
			it = std::find(migr_sites_appended.begin(),migr_sites_appended.end(),current_site);
			it2 = std::find(migr_sites_appended.begin(),migr_sites_appended.end(),coupled_site);
			if (it == migr_sites_appended.end() && it2 == migr_sites_appended.end()){
				int steps;
				bool b4;
				if (coupled_site_dir == -1){
					//if ( isR5internal(current_site->C1->C1,current_site->C1,true) || isR5internal(current_site->C1->C1,current_site->C1,false)){
					cpair R5coords = findR5internal(current_site->C1->C1,current_site->C1);
					//coupled_site = moveIt(current_site,-1);
					steps=0;
					b4 = true;
					//next_site = moveIt(current_site,+1);
					//Move the walker to its edge location. (Artificial move)
					//proc_M5R_R5R6_out_of_corner(current_site,current_site->C1,current_site->C2,next_site,false);
					//steps = -1;
				}else if (coupled_site_dir == 1){
					//else if(isR5internal(current_site->C2,current_site->C2->C2,true) || isR5internal(current_site->C2,current_site->C2->C2,false)){
					cpair R5coords = findR5internal(current_site->C2,current_site->C2->C2);
					//coupled_site = moveIt(current_site,+1);
					steps=0;
					b4 = false;
					//next_site = moveIt(current_site,-1);
					//Move the walker to its edge location. (Artificial move)
					//proc_M5R_R5R6_out_of_corner(current_site,current_site->C1,current_site->C2,next_site,true);
					//steps = 1;
				} else{
					std::cout << "Error. R5 not found in StartMigrationProcess for R5R6_MIGR site." << std::endl;
				}
				//std::tuple<Spointer,int> migr_site_ii = std::make_tuple(next_site, steps);
				std::tuple<Spointer,Spointer,int> migr_site_ii;
				
				if (b4) migr_site_ii = std::make_tuple(coupled_site, current_site, steps);
				else migr_site_ii = std::make_tuple(current_site, coupled_site, steps);
				migr_sites.push_back(migr_site_ii);
				migr_sites_appended.push_back(current_site);
				if (coupled_site->type != R5ACR5) migr_sites_appended.push_back(coupled_site);
			}
		}
	}

	//Append sites that are currently locked walkers
	std::vector<kmcSiteType> migr_site_types = {ACR5, FEACR5, ZZACR5, ACACR5, R5ACR5, R5R6, R5R6FER5R6};
	std::vector<Spointer> possible_locked_sites;
	for(unsigned int iii=0; iii!=migr_site_types.size();iii++){
		for(unsigned int ii=0;ii!=m_pah->m_siteMap[migr_site_types[iii]].size();ii++){
			Spointer st = (m_pah->m_siteMap[migr_site_types[iii]])[ii];
			possible_locked_sites.push_back(st);
		}
	}
	for(unsigned int ii=0;ii!=possible_locked_sites.size();ii++){
		Spointer st = possible_locked_sites[ii];
		std::vector<Spointer>::iterator it;
		it = std::find(migr_sites_appended.begin(),migr_sites_appended.end(),st);
		if(it==migr_sites_appended.end()){
			//The site st has not been appended. Needs to be checked.
			cpair R5coords;
			if (st->type == ACR5 || st->type == FEACR5 || st->type == ZZACR5 || st->type == R5ACR5){
				std::tuple<Spointer,Spointer,int> migr_site_ii= std::make_tuple(st, st, 0);
				migr_sites.push_back(migr_site_ii);
				migr_sites_appended.push_back(st);
				if (st->type == ACR5) {
					R5coords = findR5internal(st->C1->C2, st->C2->C1);
					Cpointer C_check_other_side = st->C1->C2;
					Cpointer C_check_other_side2 = st->C2->C1;
					Cpointer C_other_side = findThirdC(C_check_other_side);
					Cpointer C_other_side2 = findThirdC(C_check_other_side2);
					if (C_other_side!=NULLC || C_other_side2!=NULLC) {
						Spointer opp_site;
						if (C_other_side!=NULLC) opp_site = findSite(C_other_side);
						else opp_site = findSite(C_other_side2);
						std::vector<Spointer>::iterator it;
						it = std::find(migr_sites_appended.begin(),migr_sites_appended.end(),opp_site);
						if(it==migr_sites_appended.end()) m_pah->m_R5loc.push_back(R5coords);
					}
				}
				/*else {
					if ( isR5internal(st->C1->C2, st->C1->C2->C2,false) || isR5internal(st->C1->C2, st->C1->C2->C2,true) ) {
					R5coords = findR5internal(st->C1->C2, st->C1->C2->C2);
					}
					else if ( isR5internal(st->C2->C1->C1, st->C2->C1,false) || isR5internal(st->C2->C1->C1, st->C2->C1,true) ) {
						R5coords = findR5internal(st->C2->C1->C1, st->C2->C1);
					}
					else {
						//R5 not found
						std::cout << "ERROR. R5 not found for " << st->type << " site in PAHProcess::startMigrationProcess." << std::endl;
					}
				}*/
			}
			else{
				//site is R5R6 or R5R6FER5R6
				int direction_int = coupledSiteDirection(st);
				bool b4;
				Spointer coupled_site = moveIt(st,direction_int);
				std::vector<Spointer>::iterator it2;
				it2 = std::find(migr_sites_appended.begin(),migr_sites_appended.end(),coupled_site);
				if(it!=migr_sites_appended.end()){
					//coupled site has been added as migrator but not st
					std::cout<< "Error on R5 migration in startMigrationProcess(). Copuled site appears as walker but not current site." << std::endl;
					printSitesMigration();
				}
				else{
					if ( direction_int == -1){
						cpair R5coords = findR5internal(st->C1->C1,st->C1);
						b4 = true;
					}else if(direction_int == 1){
						cpair R5coords = findR5internal(st->C2,st->C2->C2);
						b4 = false;
					} else{
						if (st->type==R5R6) std::cout << "Error in startMigrationProcess(). R5 not found." << std::endl;
					}
					//std::tuple<Spointer,int> migr_site_ii = std::make_tuple(next_site, steps);
					std::tuple<Spointer,Spointer,int> migr_site_ii;
					
					if (b4) migr_site_ii = std::make_tuple(coupled_site, st, 0);
					else migr_site_ii = std::make_tuple(st, coupled_site, 0);
					migr_sites.push_back(migr_site_ii);
					migr_sites_appended.push_back(st);
					migr_sites_appended.push_back(coupled_site);
				}
			}
		}
	}

	m_pah->m_R5walker_sites = migr_sites;
}

//! Called after migration processes. Loops through the random walker sites and moves them to end locations.
void PAHProcess::performMigrationProcess(){
	for(unsigned int ii=0;ii!=m_pah->m_R5walker_sites.size();ii++){
		Spointer site_perf = std::get<0>(m_pah->m_R5walker_sites[ii]);
		Spointer site_perf_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
		int steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
		Spointer sFE2 = moveIt(site_perf,steps);
		bool b4;
		if(steps < 0) b4 = true;
		else b4 = false;
		if (site_perf->type == FE) proc_M5R_R5R6_multiple_sites(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		else if((int)site_perf->type <= 4) proc_M5R_ACR5_termination(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		else if((int)site_perf->type <= 104) proc_M5R_ACR5_termination(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		else if((int)site_perf->type <= 504) proc_M5R_ACR5_termination(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		else if((int)site_perf->type <= 1004) proc_M5R_ACR5_termination(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		else if ((int)site_perf->type>=2000 && (int)site_perf->type<=2100) proc_M5R_FEACR5_multiple_sites(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		else{
			//stt is already a termination site
			//std::cout << "Start site is termination site. Do nothing." << std::endl;
			if (steps!=0) {
				std::cout << "Error on R5 migration. Site was not moved in performMigrationProcess()." <<std::endl;
				printSitesMigration();
			}
		}
	}
	//Optimise once after all sites have been moved
	OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol, 500);
	passbackPAH(mol);
	//Clear the walkers vector.
	m_pah->m_R5walker_sites.clear();
}

//! Returns true if site is allowed for migration
bool PAHProcess::checkSiteMigration(Spointer stt, bool b4){
	Spointer checkR5_1, checkR5_2, checkR5_3;
	if (b4){
		checkR5_1 = moveIt(stt,-1);
		checkR5_2 = moveIt(stt,-2);
		checkR5_3 = moveIt(stt,-3);
	}
	else{
		checkR5_1 = moveIt(stt,+1);
		checkR5_2 = moveIt(stt,+2);
		checkR5_3 = moveIt(stt,+3);
	}

	//Check for unsupported sites. This section heavily assumes that the Isolated Pentagon Rule is valid.
	if ((int)stt->type == 0 && (int)checkR5_1->type == 0 && (int)checkR5_2->type == 0) return false; // The result would be an indene, not supported. YET!
	if ((int)stt->type == 100 || (int)stt->type == 101 || (int)stt->type == 501 || (int)stt->type == 2002) return false; // This would violate the IPR.
	//if ((int)stt->type >= 502 && (int)stt->type <= 502) return false; // This would violate the IPR.
	if ((int)stt->type >= 602 && (int)stt->type <= 604) return false; // This would violate the IPR.
	if ((int)stt->type >= 1002 && (int)stt->type <= 1004) return false; // This would violate the IPR.
	if ((int)stt->type == 9999 || (int)stt->type == -1 || stt->type == None) return false;
	if ((int)stt->type == 0){
		if ((int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 || (int)checkR5_1->type == 100) return false;
		if ((int)checkR5_1->type >= 1002 && (int)checkR5_1->type <= 1004) return false;
		if ((int)checkR5_1->type >= 2002 && (int)checkR5_1->type <= 2204) return false;
		if ((int)checkR5_1->type >= 2204 && (int)checkR5_1->type <= 2205) return false;
		if ((int)checkR5_1->type == 0){
			if ((int)checkR5_2->type == 101 || (int)checkR5_2->type == 501) return false;
			if ((int)checkR5_2->type >= 1002 && (int)checkR5_2->type <= 1004) return false;
			if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2205) return false;
		}
	}
	if ((int)stt->type == 1 && (int)checkR5_1->type == 0 && (int)checkR5_2->type == 0 && (int)checkR5_3->type == 2002) return false;
	if((int)stt->type == 2002 || (int)stt->type == 2003 || (int)stt->type == 2103 || (int)stt->type == 2104) return false;
	//Check for other sites
	Cpointer CR5_otherside_end;
	if (b4) CR5_otherside_end = stt->C2->C1;
	else CR5_otherside_end = stt->C1->C2;
	if (CR5_otherside_end->bridge) return false;
	if (CR5_otherside_end->C2->bridge) return false; 
	if (CR5_otherside_end->C1->bridge) return false;
	if (CR5_otherside_end->C2->C2->bridge) return false; 
	if (CR5_otherside_end->C1->C1->bridge) return false;
	
	//Check for other side being valid
	Cpointer thirdC_after = findThirdC(CR5_otherside_end);
	if (thirdC_after != NULLC && (int)checkR5_1->type%10 < 4){
		Spointer opp_site_after = findSite(thirdC_after);
		//I have no clue how to flag an Spointer as error. This can cause seg faults.
		if (opp_site_after != m_pah->m_siteList.end()){
			int os_endtype = opp_site_after->type;
			if (os_endtype >= 200 && os_endtype <= 203) return false;
			if (os_endtype == 101) return false;
			if (os_endtype >= 600 && os_endtype <= 603) return false;
			if (os_endtype >= 1000 && os_endtype <= 1003) return false;
			if (os_endtype >= 500 && os_endtype <= 504) return false;
			if (os_endtype >= 2000 && os_endtype <= 2205) return false;
			if (os_endtype >= 2103 && os_endtype <= 2105) return false;
			if (os_endtype >= 2204 && os_endtype <= 2205) return false;
			if (os_endtype == 9999 || os_endtype == -1 || opp_site_after->type == None) return false;
		}
	}
		
	//check that two pentagons (including internals) will not collide
	cpair R5coords_end;
	/*if (b4) R5coords_end = endposR5internal(CR5_otherside_end, CR5_otherside_end->C2);
	else R5coords_end = endposR5internal(CR5_otherside_end->C1, CR5_otherside_end,true);*/
	if (b4) {
		if (CR5_otherside_end->C2->A=='H') R5coords_end = endposR5internal(CR5_otherside_end, CR5_otherside_end->C2);
		else R5coords_end = endposR5internal(CR5_otherside_end, CR5_otherside_end->C2,true);
	}else{
		if (CR5_otherside_end->A=='H') R5coords_end = endposR5internal(CR5_otherside_end->C1, CR5_otherside_end);
		else R5coords_end = endposR5internal(CR5_otherside_end->C1, CR5_otherside_end,true);
	}
	if (m_pah->m_R5loc.size()>=1){
		for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
			double distR5s = getDistance_twoC(*it, R5coords_end);
			if (distR5s < 2.8) {
				//This distance is a parameter of this jump process. Might need some more tuning. 
				//2.8 seems appropiate but may reject too many jumps.
				//Two pentagons will be next to each other violating the Isolated Pentagon Rule
				return false;
			}
		}
	}
		
	//check that pentagon and heptagon (including internals) will not collide
	if (m_pah->m_R7loc.size()>=1){
		for (std::list<cpair>::iterator it = m_pah->m_R7loc.begin(); it!= m_pah->m_R7loc.end(); ++it){
			double distR5R7 = getDistance_twoC(*it, R5coords_end);
			if (distR5R7 < 2.6) {
				//This distance is a parameter of this jump process. Might need some more tuning. 
				//2.8 seems appropiate but may reject too many jumps.
				//Two pentagons will be next to each other violating the Isolated Pentagon Rule
				return false;
			}
		}
	}
	return true;
}

//! Returns -1 or +1 if the coupled site is to the left or right and 0 if there is no coupled site.
int PAHProcess::coupledSiteDirection(Spointer stt){
	if ((int)stt->type<500) return 0;
	if ((int)stt->type>2000 && (int)stt->type<2016) return 0;
	if ((int)stt->type>=2114 && (int)stt->type<=2115) return 0;
	Spointer check_left = stt;
	Spointer check_right = stt;
	bool left_bool = true;
	bool right_bool = true;
	for (int ii=0;ii!=7;ii++){
		check_left = moveIt(check_left,-1);
		check_right = moveIt(check_right,+1);
		while ((int)check_left->type>=1000 && (int)check_left->type<=1004){
			check_left = moveIt(check_left,-1);
		}
		while ((int)check_right->type>=1000 && (int)check_right->type<=1004){
			check_right = moveIt(check_right,+1);
		}
		if(check_left->type==R5FEACR5 || check_left->type==R5ACR5){
			Spointer left_site = moveIt(check_left,-1);
			if (left_site->type==R5) left_bool = false;
		}
		if(check_left->type==R5FEACR5 || check_left->type==R5ACR5){
			Spointer right_site = moveIt(check_right,+1);
			if (right_site->type==R5) right_bool = false;
		}
		
		if ((int)check_left->type<500) left_bool = false;
		if ((int)check_left->type>2000 && (int)check_left->type<2016) left_bool = false;
		if ((int)check_left->type>=2114 && (int)check_left->type<=2115) left_bool = false;
		if ((int)check_left->type==9999) left_bool = false;

		if ((int)check_right->type<500) right_bool = false;
		if ((int)check_right->type>2000 && (int)check_right->type<2016) right_bool = false;
		if ((int)check_right->type>=2114 && (int)check_right->type<=2115) right_bool = false;
		if ((int)check_right->type==9999) right_bool = false;

		bool invert;
		if (ii%2 == 0) invert = false;
		else invert = true;
		if(left_bool && !right_bool && !invert) return -1;
		if(!left_bool && right_bool && !invert) return 1;
		if(left_bool && !right_bool && invert) return 1;
		if(!left_bool && right_bool && invert) return -1;
	}
	return 0;
}

//! Modifies pointer for R5 walkers to avoid overlaps.
void PAHProcess::checkR5Walkers(){
	//Check that the current walker has a site available for addition of a C
	for (int ii=0; ii!=m_pah->m_R5walker_sites.size();ii++){
		Spointer start_site_ii = std::get<0>(m_pah->m_R5walker_sites[ii]);
		Spointer start_site_ii2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
		int ii_steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
		Spointer end_site_ii, end_site_ii2;
		end_site_ii = moveIt(start_site_ii,ii_steps);
		if (ii_steps!=0){
			for (int jj=0; jj!=m_pah->m_R5walker_sites.size();jj++){
				if (jj != ii){
					Spointer start_site_jj = std::get<0>(m_pah->m_R5walker_sites[jj]);
					Spointer start_site_jj2 = std::get<1>(m_pah->m_R5walker_sites[jj]);
					int jj_steps = std::get<2>(m_pah->m_R5walker_sites[jj]);
					Spointer end_site_jj = moveIt(start_site_jj,jj_steps);
					Spointer end_site_jj2 = moveIt(start_site_jj2,jj_steps);
					//Check if walker is now on a corner
					Spointer check_corner_site, check_site, check_site2;
					if (jj_steps<0) check_corner_site = moveIt(start_site_jj,jj_steps+1);
					else check_corner_site = moveIt(start_site_jj,jj_steps-1);
					if((int)check_corner_site->type>=501 && (int)check_corner_site->type<=1004 ){
						if (jj_steps<0) {
							check_site = moveIt(start_site_jj,jj_steps-1);
							check_site2 = moveIt(start_site_jj,jj_steps);
						} else{
							check_site = moveIt(start_site_jj,jj_steps);
							check_site2 = moveIt(start_site_jj,jj_steps+1);
						}
					}else{
						check_site = moveIt(start_site_jj,jj_steps);
						check_site2 = check_site;
					}
					if(check_site->type==FE || check_site2->type==FE){
						if (jj_steps>0) check_site2 = moveIt(check_site,+1);
						else check_site2 = moveIt(check_site,-1);
					}
					if (check_site == start_site_ii || check_site == start_site_ii2 || check_site2 == start_site_ii || check_site == start_site_ii2){
						//The next position of walker jj will become the start location of walker ii.
						//This will mess up the sites.
						//First - Move walker ii to current position and steps = 0
						//Second - Return the PAH to the program. It should handle walker jj now
						bool local_b4;
						if (jj_steps<0) {
							local_b4 = true;
						}
						else {
							local_b4 = false;
						}
						//saveXYZ("KMC_DEBUG/BEFORE_II_TERMINATION");
						proc_M5R_ACR5_termination(start_site_ii, start_site_ii->C1,start_site_ii->C2,end_site_ii,local_b4);
						if ((int)end_site_ii->type>500 && (int)end_site_ii->type<1100){
							if (ii_steps<0) {
								end_site_ii = moveIt(start_site_ii2,ii_steps-1);
								end_site_ii2 = moveIt(start_site_ii2,ii_steps);
							}
							else{
								end_site_ii = moveIt(start_site_ii2,ii_steps);
								end_site_ii2 = moveIt(start_site_ii2,ii_steps+1);
							}
						}
						else{
							end_site_ii = moveIt(start_site_ii,ii_steps);
							end_site_ii2 = moveIt(start_site_ii,ii_steps);
						}
						std::get<0>(m_pah->m_R5walker_sites[ii]) = end_site_ii;
						std::get<1>(m_pah->m_R5walker_sites[ii]) = end_site_ii2;
						std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
						//saveXYZ("KMC_DEBUG/AFTER_II_TERMINATION");
						//Reread the pointers
						/*end_site_jj = moveIt(start_site_jj,jj_steps);
						end_site_jj2 = moveIt(start_site_jj2,jj_steps);
						proc_M5R_ACR5_termination(start_site_jj, start_site_jj->C1,start_site_jj->C2,end_site_jj,local_b4);
						std::get<0>(m_pah->m_R5walker_sites[jj]) = end_site_jj;
						std::get<1>(m_pah->m_R5walker_sites[jj]) = end_site_jj2;
						if (local_b4) std::get<2>(m_pah->m_R5walker_sites[jj]) = -1;
						else std::get<2>(m_pah->m_R5walker_sites[jj]) = 1;*/
						//saveXYZ("KMC_DEBUG/AFTER_JJ_TERMINATION");
					}
				}
			}
		}
	}
}

//! Modifies pointer for R5 walkers to avoid overlaps. Overload that know walker jj is moving.
void PAHProcess::checkR5Walkers(int jj){
	//Check that the current walker with index jj has a site available for addition of a C
	for (int ii=0; ii!=m_pah->m_R5walker_sites.size();ii++){
		Spointer start_site_ii = std::get<0>(m_pah->m_R5walker_sites[ii]);
		Spointer start_site_ii2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
		int ii_steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
		Spointer end_site_ii, end_site_ii2;
		end_site_ii = moveIt(start_site_ii,ii_steps);
		if (ii_steps!=0){
			if (jj != ii){
				Spointer start_site_jj = std::get<0>(m_pah->m_R5walker_sites[jj]);
				Spointer start_site_jj2 = std::get<1>(m_pah->m_R5walker_sites[jj]);
				int jj_steps = std::get<2>(m_pah->m_R5walker_sites[jj]);
				Spointer end_site_jj = moveIt(start_site_jj,jj_steps);
				Spointer end_site_jj2 = moveIt(start_site_jj2,jj_steps);
				//Check if walker is now on a corner
				Spointer check_corner_site, check_site, check_site2;
				if (jj_steps<0) check_corner_site = moveIt(start_site_jj,jj_steps+1);
				else check_corner_site = moveIt(start_site_jj,jj_steps-1);
				if((int)check_corner_site->type>=501 && (int)check_corner_site->type<=1004 ){
					if (jj_steps<0) {
						check_site = moveIt(start_site_jj,jj_steps-1);
						check_site2 = moveIt(start_site_jj,jj_steps);
					} else{
						check_site = moveIt(start_site_jj,jj_steps);
						check_site2 = moveIt(start_site_jj,jj_steps+1);
					}
				}else{
					check_site = moveIt(start_site_jj,jj_steps);
					check_site2 = check_site;
				}
				if(check_site->type==FE || check_site2->type==FE){
					if(check_site->type==FE){
						if (jj_steps>0) check_site2 = moveIt(check_site,+1);
						else check_site2 = moveIt(check_site,-1);
					}else{
						if (jj_steps>0) check_site = moveIt(check_site2,+1);
						else check_site = moveIt(check_site2,-1);
					}
				}
				if (check_site == start_site_ii || check_site == start_site_ii2 || check_site2 == start_site_ii || check_site == start_site_ii2){
					//The next position of walker jj will become the start location of walker ii.
					//This will mess up the sites.
					//First - Move walker ii to current position and steps = 0
					//Second - Return the PAH to the program. It should handle walker jj now
					bool local_b4;
					if (jj_steps<0) {
						local_b4 = true;
					}
					else {
						local_b4 = false;
					}
					//saveXYZ("KMC_DEBUG/BEFORE_II_TERMINATION");
					proc_M5R_ACR5_termination(start_site_ii, start_site_ii->C1,start_site_ii->C2,end_site_ii,local_b4);
					if ((int)end_site_ii->type>500 && (int)end_site_ii->type<1100){
						if (ii_steps<0) {
							end_site_ii = moveIt(start_site_ii2,ii_steps-1);
							end_site_ii2 = moveIt(start_site_ii2,ii_steps);
						}
						else{
							end_site_ii = moveIt(start_site_ii2,ii_steps);
							end_site_ii2 = moveIt(start_site_ii2,ii_steps+1);
						}
					}
					else{
						end_site_ii = moveIt(start_site_ii,ii_steps);
						end_site_ii2 = moveIt(start_site_ii,ii_steps);
					}
					std::get<0>(m_pah->m_R5walker_sites[ii]) = end_site_ii;
					std::get<1>(m_pah->m_R5walker_sites[ii]) = end_site_ii2;
					std::get<2>(m_pah->m_R5walker_sites[ii]) = 0;
					if((int)end_site_ii->type>=2003 && local_b4==false) {
						if(ii_steps < 0) addR5internal(end_site_ii->C2->C1->C1,end_site_ii->C2->C1,true);
						else addR5internal(end_site_ii->C1->C2,end_site_ii->C1->C2->C2,true);
					} else if ((int)end_site_ii2->type>=2003 && local_b4==true) {
						if(ii_steps < 0) addR5internal(end_site_ii2->C2->C1->C1,end_site_ii2->C2->C1,true);
						else addR5internal(end_site_ii2->C1->C2,end_site_ii2->C1->C2->C2,true);
					}
					//saveXYZ("KMC_DEBUG/AFTER_II_TERMINATION");
					//Reread the pointers
					/*end_site_jj = moveIt(start_site_jj,jj_steps);
					end_site_jj2 = moveIt(start_site_jj2,jj_steps);
					proc_M5R_ACR5_termination(start_site_jj, start_site_jj->C1,start_site_jj->C2,end_site_jj,local_b4);
					std::get<0>(m_pah->m_R5walker_sites[jj]) = end_site_jj;
					std::get<1>(m_pah->m_R5walker_sites[jj]) = end_site_jj2;
					if (local_b4) std::get<2>(m_pah->m_R5walker_sites[jj]) = -1;
					else std::get<2>(m_pah->m_R5walker_sites[jj]) = 1;*/
					//saveXYZ("KMC_DEBUG/AFTER_JJ_TERMINATION");
				}
			}
		}
	}
}

//! Modifies pointer for R5 walker that moves into one of the coupled sites of other walker.
void PAHProcess::checkRemR5Walkers(int jj, bool b4, Spointer sFE2){
	for (int ii=0; ii!=m_pah->m_R5walker_sites.size();ii++){
		if (ii!=jj){
			Spointer check_site = std::get<0>(m_pah->m_R5walker_sites[ii]);
			Spointer check_site2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
			if (sFE2 == check_site){
				if (b4) std::get<0>(m_pah->m_R5walker_sites[ii]) = moveIt(sFE2,+1);
				else std::get<0>(m_pah->m_R5walker_sites[ii]) = moveIt(sFE2,-1);
			} 
			if (sFE2 == check_site2) {
				if (b4) std::get<1>(m_pah->m_R5walker_sites[ii]) = moveIt(sFE2,+1);
				else std::get<1>(m_pah->m_R5walker_sites[ii]) = moveIt(sFE2,-1);
			}
		}
	}
}

//! Find the walker associated with a site .
int PAHProcess::findWalker(Spointer current_site){
	unsigned int ii = -9999;
	for (unsigned int jj=0;jj!=m_pah->m_R5walker_sites.size();jj++){
		Spointer start_site = std::get<0>(m_pah->m_R5walker_sites[jj]);
		Spointer start_site_2 = std::get<1>(m_pah->m_R5walker_sites[jj]);
		int steps = std::get<2>(m_pah->m_R5walker_sites[jj]);
		Spointer check_site = moveIt(start_site, steps);
		Spointer check_site_2 = moveIt(start_site_2, steps);
		if (check_site == current_site || check_site_2 == current_site){
			ii = jj;
			return ii;
		}
	}
	std::cout << "Error. findWalker could not find associated walker." << std::endl;
	std::ostringstream msg;
            msg << "Error. findWalker could not find associated walker." << std::endl;
	printSitesMigration();
	throw std::runtime_error(msg.str());
    assert(false);
}

//! Adds walker when an opposite side site can now migrate.
void PAHProcess::addOppsiteR5Walker(Spointer opp_site, Spointer opp_site_coupled){
	std::tuple<Spointer,Spointer,int> new_walker = std::make_tuple(opp_site,opp_site_coupled,0);
	m_pah->m_R5walker_sites.push_back(new_walker);
}

//! Fixes opposite side walker if needed.
void PAHProcess::fixOppsiteR5Walker(int ii){
	Spointer start_site = std::get<0>(m_pah->m_R5walker_sites[ii]);
	Spointer start_site_2 = std::get<1>(m_pah->m_R5walker_sites[ii]);
	if (start_site==start_site_2 && start_site->type==R5R6){
		int dir = coupledSiteDirection(start_site);
		start_site_2 = moveIt(start_site, dir);
		if (dir==-1){
			std::get<0>(m_pah->m_R5walker_sites[ii]) = start_site_2;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = start_site;
		} else{
			std::get<0>(m_pah->m_R5walker_sites[ii]) = start_site;
			std::get<1>(m_pah->m_R5walker_sites[ii]) = start_site_2;
		}
	}
}

//! Remove opposite side site walker and modifies structure if needed. Removes walker jj and modifies walker ii. Returns the new position of ii.
int PAHProcess::remOppsiteR5Walker(int ii, int jj){
	//ii is the walker that is moving.
	//jj is the walker on the opposite side.
	int ii_steps = std::get<2>(m_pah->m_R5walker_sites[ii]);
	if (ii_steps == 0){
		//The current walker has 0 steps and an opposite side site. This means that the walker comes from other edge.
		//Move walker jj to where it jumped edges.
		Spointer site_perf = std::get<0>(m_pah->m_R5walker_sites[jj]);
		Spointer site_perf_2 = std::get<1>(m_pah->m_R5walker_sites[jj]);
		int steps = std::get<2>(m_pah->m_R5walker_sites[jj]);
		Spointer sFE2 = moveIt(site_perf,steps);
		bool b4;
		if(steps < 0) b4 = true;
		else b4 = false;
		if (site_perf->type == FE) proc_M5R_R5R6_multiple_sites(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		else if((int)site_perf->type < 2000) proc_M5R_ACR5_termination(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		else if ((int)site_perf->type>=2000 && (int)site_perf->type<=2100) proc_M5R_FEACR5_multiple_sites(site_perf,site_perf->C1,site_perf->C2,sFE2,b4);
		//Delete walker jj
		m_pah->m_R5walker_sites.erase(m_pah->m_R5walker_sites.begin() + jj);
	} else{
		//The current walker has non-zero steps and an opposite side site. This means that the walker is on this edge.
		//Remove walker jj.
		m_pah->m_R5walker_sites.erase(m_pah->m_R5walker_sites.begin() + jj);
	}
	if (ii>jj) ii--;
	return ii;
}