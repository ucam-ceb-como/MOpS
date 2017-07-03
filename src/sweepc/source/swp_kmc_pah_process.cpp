/*!
  * \author     Zakwan Zainuddin (zz260)
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
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>

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
    std::vector<kmcSiteType> sites = SiteVector();
    p.createPAH(sites, m_pah->m_rings, m_pah->m_rings5_Lone, m_pah->m_rings5_Embedded, m_pah->numofC(), m_pah->numofH());
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
        /**
         * If this condition is true then this PAH is a benzene ring.
         * The current rate of benzene removal appears to be too high so this in fact disables the process by setting the rate to 0.
         */
        if(getCHCount().first == 6) return 0;

        return (unsigned int) (m_pah->m_siteMap[st].size());
    }
    return (unsigned int) m_pah->m_siteMap[st].size();
}
//! Get Ring Counts
intpair PAHProcess::getRingsCount() const {
    return intpair(m_pah->m_rings, m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded);
}

int PAHProcess::getRings5EmbeddedCount() const {
    return m_pah->m_rings5_Embedded;
}

//! Get number of bridges
int PAHProcess::numberOfBridges() const {
    int num = 0;
	//NICK TO DO - Calculate the number of bridges based on new data type
    //for(Ccontainer::const_iterator i=m_pah->m_carbonList.begin(); i!=m_pah->m_carbonList.end(); i++) {
    //    if((*i)->bridge) num++;
    //}
    //num /= 2;
    return num;
}
//! Print Structure in console
//void PAHProcess::printStruct() const{
//    // iterator set to first C in structure
//    Cpointer prev = m_pah->m_cfirst;
//    Cpointer now = m_pah->m_cfirst->C2;
//    // counts the number of C at edge
//    unsigned int count = 0;
//    // a string to display if current C is a bridge
//    std::string btxt;
//    angletype angle;
//    do {
//        count++;
//        // checks if bridge
//        if(prev->bridge) {
//            btxt = " bridge";
//            angle = prev->bondAngle2;
//        } else {
//            btxt = "       ";
//            angle = prev->bondAngle1;
//        }
//        // outputs text in the form of "X: (C-H or C-C)    bond angle with next C"
//        std::cout << "coords ("<<prev->coords.first<<","<< prev->coords.second<<")"<< btxt << ": C-" << prev->A << "\t" << angle << '\n';
//        // moves iterator to next C
//        Cpointer oldnow = now;
//        now = moveCPointer(prev, now);
//        prev = oldnow;
//    } while 
//    (!(count!=1&&prev->coords.first == m_pah->m_cfirst->coords.first
//   &&prev->coords.second == m_pah->m_cfirst->coords.second));
//    // displays C and H counts
//    std::cout << "C Count: " << m_pah->m_counts.first << '\n';
//    std::cout << "H Count: " << m_pah->m_counts.second << '\n';
//}
//
//bool PAHProcess::havebridgeC(){
//    for(Ccontainer::const_iterator i=m_pah->m_carbonList.begin(); i!=m_pah->m_carbonList.end(); i++) {
//        if((*i)->bridge) return true;
//    }
//    return false;
//}
//! Print Structure in console, with arrow pointing at current C
//void PAHProcess::printStruct(Cpointer c) const{
//    // iterator set to first C in structure
//    Cpointer prev = m_pah->m_cfirst;
//    Cpointer now = m_pah->m_cfirst->C2;
//    // counts the number of C at edge
//    unsigned int count = 0;
//    // a string to display if current C is a bridge
//    std::string btxt;
//    angletype angle;
//    do {
//        count++;
//        // checks if bridge
//        if(prev->bridge) {
//            btxt = " bridge";
//            angle = prev->bondAngle2;
//        } else {
//            btxt = "       ";
//            angle = prev->bondAngle1;
//        }
//        // outputs text in the form of "X: (C-H or C-C)    bond angle with next C"
//        std::cout << prev << btxt << ": C-" << prev->A << "\t" << angle;
//        if(prev == c) std::cout << "\t<---";
//        std::cout << '\n';
//        // moves iterator to next C
//        Cpointer oldnow = now;
//        now = moveCPointer(prev, now);
//        prev = oldnow;
//    } while (prev != m_pah->m_cfirst);
//    // displays C and H counts
//    std::cout << "C Count: " << m_pah->m_counts.first << '\n';
//    std::cout << "H Count: " << m_pah->m_counts.second << '\n';
//}
//! Print Sites in console
//void PAHProcess::printSites() const{
//    Spointer i;
//    std::string st;
//    cout << "*******************\n";
//    cout << "Sites List:\n_____\n";
//    // displays total site count
//    cout << "Total Site Count: " << m_pah->m_siteList.size() << '\n';
//    for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
//        // convert site type into string
//        st = kmcSiteName(i->type);
//        // displays site type
//        cout << st << '\n';
//    }
//    cout << "********************\n";
//}
//! Print Sites in console, with an arrow pointing at site stt
void PAHProcess::printSites(Spointer& stt) const{
    Spointer i;
    std::string st;
    string pointer = "  <-  ";
    cout << "*******************\n";
    cout << "Sites List:\n_____\n";
    // displays total site count
    cout << "Total Site Count: " << m_pah->m_siteList.size() << '\n';
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
//! Print sites & site members in console, with arrow pointing at site stt
//NICK TO DO - Fix this
//void PAHProcess::printSitesMemb(Spointer& stt) const{
//	Spointer i;
//	std::string st;
//	string pointer = "  <-  ";
//	cout << "*******************\n";
//	cout << "Sites List:\n_____\n";
//	// displays total site count
//	cout << "Total Site Count: " << m_pah->m_siteList.size() << '\n';
//	for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
//		// convert site type into string
//		cout<< kmcSiteName(i->type);
//		cout<<"\t";
//		cout<<kmcSiteName(i->comb);
//		cout<<'\t'<<i->C1<<'\t'<<i->C2;
//		// checks if site i equivalent to stt, if yes put an arrow
//		if(i == stt) cout<<pointer;
//		// displays site type (and arrow)
//		cout <<'\n';
//	}
//	cout << "********************\n";
//}
//void PAHProcess::printSitesMemb() const{
//    Spointer i;
//    std::string st;
//    string pointer = "  <-  ";
//    cout << "*******************\n";
//    cout << "Sites List:\n_____\n";
//    // displays total site count
//    cout << "Total Site Count: " << m_pah->m_siteList.size() << '\n';
//    for(i=m_pah->m_siteList.begin(); i!=m_pah->m_siteList.end(); i++) {
//        // convert site type into string
//        cout<< kmcSiteName(i->type);
//        cout<<"\t";
//        cout<<kmcSiteName(i->comb);
//        cout<<'\t'<<i->C1<<'\t'<<i->C2;
//        // checks if site i equivalent to stt, if yes put an arrow
//        //if(i == stt) cout<<pointer;
//        // displays site type (and arrow)
//        cout <<'\n';
//    }
//    cout << "********************\n";
//}

//! Store results in external file, returns success/failure
// Stores structure into a DOT file for Graphviz
//bool PAHProcess::saveDOT(const std::string &filename) const {
//    bool success = saveDOT(filename, std::string("none"));
//    return success;
//}
// Stores structure into a DOT file for Graphviz with a title on graph
//bool PAHProcess::saveDOT(const std::string &filename, const std::string &title) const {
//    ofstream dotfile;
//    // creates or opens filename
//    //cout << "Attempting to open " << filename << "...\n";//++++
//    dotfile.open(filename.c_str());
//    if(dotfile.is_open()) { // if successful
//        // vector to store indices of C bonded with H atom
//        std::vector<int> CH_index;
//        // set iterator to starting position
//        Cpointer now = m_pah->m_cfirst->C2;
//        Cpointer prev = m_pah->m_cfirst;
//        //cout << filename << " had been opened, prepare writing..\n";//++++
//        // draw undirected graph
//        dotfile<< "Graph G {\n";
//        // sets max size for drawing structure (in inches)
//        dotfile<< "\tgraph [size=15];\n";
//        // sets style for C nodes
//        dotfile<< "\tnode [shape=circle, style=filled, color=black, height=0.2, width=0.2];\n";
//        // sets style for bonds
//        dotfile<< "\tedge [weight=2, style=\"setlinewidth(10)\"];\n";
//        // sets caption for graph
//        if(title!="none") {
//            dotfile<< "\tgraph [label=\""<<title<<"\"];\n";
//        }
//        // get coordinates of current C
//        coordtype x=prev->coords.first;
//        coordtype y=prev->coords.second;
//        int c=0; // index of current C
//        int count=0;
//        float CC_CHratio = 10; // length ratio between C-C and C-H bonds
//        do {
//            //if(now->bridge) {
//                //cout<<"Bridge at "<<now<<'\n';
//            //}
//            c++;
//            // write position coordinates of C atom with index c, with ! suffixed to pin node
//            dotfile<<"\tC"<<c<<" [pos = \""<<x<<','<<y<<"!\"";
//            if(c==1) dotfile<<", color=red";
//            dotfile<<"];\n";
//            if(prev->A=='H') { // if bonded with H
//                // write position coordinates of H atom with bonded C index, with ! suffixed to pin node
//                dotfile<<"\tH"<<c<<" [pos = \""<<(x+x_inc(prev->bondAngle2)/CC_CHratio);
//                // set style for H node
//                dotfile<<','<<(y+y_inc(prev->bondAngle2)/CC_CHratio)<<"!\", color=skyblue, fontcolor=skyblue, height=0.1, width=0.1];\n";
//                // stores index of current C
//                CH_index.push_back(c);
//            }
//            // move iterator to next C
//            Cpointer oldnow = now;
//            now = moveCPointer(prev, now);
//            prev = oldnow;
//            // obtain coordinates
//            x=prev->coords.first;
//            y=prev->coords.second;
//        }
//        while
//        (!(count!=1&&prev->coords.first == m_pah->m_cfirst->coords.first
//   &&prev->coords.second == m_pah->m_cfirst->coords.second));
//        // connects each C atom
//        for(int i=1; i<c; i++) {
//            dotfile<<"\tC"<<i<<" -- "<<"C"<<(i+1)<<";\n";
//        }
//        dotfile<<"\tC"<<c<<" -- "<<"C1;\n"; // connects last C with first C
//        // connects each H atom with corresponding C atoms
//        for(int i=0; i!=(int)CH_index.size(); i++) {
//            dotfile<<"\tH"<<CH_index[i]<<" -- "<<"C"<<CH_index[i]<<";\n";
//        }
//        // close dot file
//        dotfile<<"}";
//        dotfile.close();
//        //cout<<".DOT file writing complete.\n";//++++
//        return true;
//    }
//    // if unable to open/create dot file
//    else cout << "unable to open file\n";
//    return false;
//}
// Protected Read Processes
//! Get other member of the site a particular C atom is a member of
//Cpointer PAHProcess::getPair(const Cpointer Carb, bool after) const {/*
//    if(after) {
//        return Carb->S2->C2;
//    }else {
//        return Carb->S1->C1;
//    }*/
//    return NULL;
//}
//! Returns the next carbon atom after current, coming from the previous
//Cpointer PAHProcess::moveCPointer(Cpointer &previous, Cpointer &current) const {
//    //cout << "Moving from " << current << " with " << previous << " preceding..\n";
//    if(current == 0x0) {
//        cout << "ERROR: Moving from NULL C Pointer\n\n";
//        return NULL;
//    }
//    // check for bridge
//    if(!(current->bridge)) {
//    // if not bridge continue to next C
//        //cout << "current is not bridge, going to " << current->C2 << '\n';
//        return current->C2;
//        //cout << "moved on to " << current << '\n';
//    } else {
//        if(previous == current->C1) {
//            // if coming from main PAH, move iterator to bridged PAH
//            return current->C3;
//        }else if(previous == current->C3) {
//            // if coming from bridged PAH, move to next on main PAH
//            return current->C2;
//        }
//    }
//    return current->C2;
//    //cout << "moved on to " << current << '\n';
//}
//! Check if process is allowed [returns true for now]
bool PAHProcess::allowed(const Spointer& st, StructureProc proc) const {
    return true;
}
//! Jump to a position coordinate given starting position and angle towards new position
//cpair PAHProcess::jumpToPos(const cpair& starting, const angletype& direction) const{
//    cpair temp;
//    // calculate new coordinates
//    temp.first = starting.first + x_inc(direction);
//    temp.second = starting.second + y_inc(direction);
//    return temp;
//}
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
//bool PAHProcess::checkHindrance(const Spointer& st) const {
//    // position to check, start at C1 of site
//    cpair mpos;
//    switch(st->type) {
//    case FE:
//        // move position to O1, O2, O3, and O4
//        mpos = jumpToPos(st->C1->coords, st->C1->bondAngle2); // position O1
//        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//        mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2-60)); // position O2
//        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//        mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2-120)); // position O3
//        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//        mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2-180)); // position O4
//        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//        break;
//    case RFE:
//    case ZZ:
//        // move position to O1, O2 and O3
//        mpos = jumpToPos(st->C1->coords, normAngle(st->C1->bondAngle1+120)); // position O1
//        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//        mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle1+60)); // position O2
//        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//        mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle1)); // position O3
//        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//        break;
//    case RZZ:
//    case AC:
//        // move position to O1 and O2
//        mpos = jumpToPos(st->C1->coords, st->C1->bondAngle2); // position O1
//        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//        mpos = jumpToPos(mpos, normAngle(st->C1->bondAngle2-60)); // position O2
//        if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//        break;
//    default:
//        return true;
//    }
//    return false;
//}
/*! Check steric hindrance for phenyl addition reactions
  ! O -> positions to check
  !                                                                           
  ! -C        O2 -- O3                                                         
  !   \       /      \                                                         
  !    C -- O1       O4                                                          
  !   /       \      /                                                        
  ! -C        O6 -- O5                                                           
*/
//bool PAHProcess::checkHindrancePhenyl(const Cpointer C_1) const {
//    // position to check, start at O1
//    cpair mpos = jumpToPos(C_1->coords, C_1->bondAngle2); // position O1
//    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//    mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2+60)); // position O2
//    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//    mpos = jumpToPos(mpos, C_1->bondAngle2); // position O3
//    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//    mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2-60)); //position O4
//    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//    mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2-120)); // position O5
//    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//    mpos = jumpToPos(mpos, normAngle(C_1->bondAngle2-180)); // position O6
//    if(m_pah->m_cpositions.count(mpos)) return true; // check if position occupied
//    return false;
//}
//! Choose random site of a site type st, returns iterator to site
Spointer PAHProcess::chooseRandomSite(kmcSiteType st, rng_type &rng) {
    // to choose any principal site
    if(st == any) {
        // choose any site index from m_pah->m_siteList
        // There does not seem to be any handling of the case that there are no sites.

        // Set up an object to generate an integer uniformly distributed on [0, size - 1]
        typedef boost::uniform_smallint<unsigned int> site_index_distrib;
        site_index_distrib siteIndexDistrib(0,  m_pah->m_siteList.size()-1);
        boost::variate_generator<rng_type &, site_index_distrib> siteIndexGenerator(rng, siteIndexDistrib);

        // move iterator to site index and return iterator
        return moveIt(m_pah->m_siteList.begin(), siteIndexGenerator());
    } else if(st == benz) { //to choose sites for phenyl addition
        return chooseRandomSite(PHsites, rng);
    } else {
        //cout << "~~Choosing from " << m_pah->m_siteMap[st].size() << " sites...\n";
        // choose site index from site vector associated with site type st
        int sz = ((int) m_pah->m_siteMap[st].size())-1; // size - 1
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
//Cpointer PAHProcess::addC() {
//    Cpointer cb;
//    // Create new carbon atom at memory pointed by h
//    cb = new Carbon;
//    if(!m_pah->m_carbonList.insert(cb).second)
//        std::cout<<"ERROR: ADDING SAME CARBON POINTER TO SET\n";
//    //m_pah->m_cpositions.insert(cb->coords); // store coordinates
//    addCount(1,0); // add a C count
//    return cb;
//}
/*! Create a new carbon atom attached next to C1
  ! C1
  !    \ <----------- angle1
  !      C __ (C1->)C2
  !        angle 2
*/
//Cpointer PAHProcess::addC(Cpointer C_1, angletype angle1, angletype angle2, bool bulk) {
//    Cpointer cb;
//    // Create new carbon atom
//    cb = new Carbon;
//    // Set details of new Carbon
//    cb->C1 = C_1;
//    cb->C2 = C_1->C2;
//    cb->bondAngle1 = normAngle(angle2); // convert angle to +ve/-ve form
//    // set new coordinates and store
//    //cb->coords = jumpToPos(C_1->coords, angle1);
//    if(!m_pah->m_carbonList.insert(cb).second)
//        std::cout<<"ERROR: ADDING SAME CARBON POINTER TO SET\n";
//    //m_pah->m_cpositions.insert(cb->coords);
//    // Edit details of connected carbon(s)
//    if(C_1->C2 != NULL) {
//        // change member pointer of original neighbour of C_1
//        C_1->C2->C1 = cb;
//    }
//    C_1->bondAngle1 = normAngle(angle1);
//    C_1->C2 = cb;
//    if(!bulk) addCount(1,0);
//    if(C_1 == m_pah->m_clast) m_pah->m_clast = cb;
//    return cb;
//}
/*! Create a carbon atom bridging next to C_1
  !              newC
  !                /  <---- angle (from C1)
  !          R - C1
  !                \
  !                 R
*/
//Cpointer PAHProcess::bridgeC(Cpointer C_1) {
//    Cpointer cb;
//    // Create new carbon atom
//    cb = new Carbon;
//    // Set details of new Carbon
//    cb->C3 = C_1;
//    cb->bondAngle2 = normAngle(C_1->bondAngle2 - 180); // opposite direction of C_1->bondAngle2
//    cb->bridge = true;
//    // set new coordinates and store
//    cb->coords = jumpToPos(C_1->coords, C_1->bondAngle2);
//    if(!m_pah->m_carbonList.insert(cb).second)
//        std::cout<<"ERROR: ADDING SAME CARBON POINTER TO SET\n";
//    m_pah->m_cpositions.insert(cb->coords);
//    // Set details of C_1
//    C_1->bridge = true;
//    C_1->C3 = cb;
//    C_1->A = 'C';
//    addCount(1,0);
//    return cb;
//}/*
//! Creates a bulk carbon atom connected to C_1
//          C2        C2                                                  
//         /           \                                                  
//  C1 - C              C - C1                                            
//         \           /                                                  
//          C3        C3                                                  
//Cpointer PAHProcess::addBC(Cpointer C_1) {
//    Cpointer cb;
//    // Create new carbon atom
//    cb = new Carbon;
//    // determine if C_1 is C1, C2 or C3
//    cb->C1 = NULL; cb->C2 = NULL; cb->C3 = NULL;
//    if(C_1->bondAngle2 == 0 || C_1->bondAngle2 == 180 || C_1->bondAngle2 == -180) 
//        cb->C1 = C_1;
//    else if(C_1->bondAngle2 == 60 || C_1->bondAngle2 == 120)
//        cb->C2 = C_1;
//    else if(C_1->bondAngle2 == -60 || C_1->bondAngle2 == -120)
//        cb->C3 = C_1;
//    else std::cout<<"ERROR: addBC(): Edge bondAngle2 invalid!!\n";
//    // Set details of new Carbon
//    cb->bondAngle1 = C_1->bondAngle2; // angle C1 coming from
//    cb->bondAngle2 = 0;
//    cb->bridge = false;
//    cb->A = 'C';
//    // set new coordinates and store
//    cb->coords = jumpToPos(C_1->coords, C_1->bondAngle2);
//    m_pah->m_cpositions.insert(cb->coords);
//    // Set details of C_1
//    C_1->C3 = cb;
//    C_1->A = 'C';
//    addCount(1,0);
//    return cb;
//}*/
//! Connects a carbon atom to another carbon (to close loop)
//void PAHProcess::connectToC(Cpointer C_1, Cpointer C_2) {
//    C_1->C2 = C_2;
//    C_2->C1 = C_1;
//}
//! Removes a carbon atom from perimeter taking into account if the atom becomes bulk carbon within the PAH
//void PAHProcess::removeC(Cpointer C_1, bool bulk) {
//    if(C_1 == m_pah->m_cfirst) {
//        m_pah->m_cfirst = C_1->C2;
//    }else if(C_1 == m_pah->m_clast) {
//        m_pah->m_clast = C_1->C1;
//    }
//    /**
//     * In order to model curved PAHs the code simply tracks the list of site types which makes up the edge of the PAH.
//     * Therefore, the coordinates of the edge carbon atoms (coords) have been made redundant.
//     * This check which depends on the coordinates of the edge carbon atom can no longer be applied.
//     */
//    //if(m_pah->m_cpositions.find(C_1->coords)== m_pah->m_cpositions.end()) {
//    //    cout<<"ERROR: removeC: coordinates ("<<C_1->coords.first<<','<<C_1->coords.second<<") not in m_pah->m_cpositions!\n";
//    //    cout<<"Coordinates of nearby 5 C atoms:\n";
//    //    Cpointer now = C_1->C1->C1->C1->C1->C1;
//    //    for(int i=0; i!=11; i++) {
//    //        cout<<'('<<now->coords.first<<','<<now->coords.second<<")";
//    //        if(now == C_1) cout<<" - ";
//    //        cout<<'\n';
//    //        now = now->C2;
//    //    }
//    //    saveDOT("COORD_PROB.dot");
//    //}
//    
//    // Change details of neighbouring Carbons
//    if(C_1->C1 != NULL) { // check prev C atom
//        C_1->C1->C2 = C_1->C2; // connect prev C atom to next C atom
//    }
//    if(C_1->C2 != NULL) { // check next C atom
//        C_1->C2->C1 = C_1->C1; // connect next C atom to prev C atom
//    }
//    if(C_1->bridge) { // change details of C atom bridging to it
//        //C_1->C3->bondAngle2;// = 0;
//        C_1->C3->C3 = NULL;
//        C_1->C3->bridge = false;
//    }
//    // Remove coordinates of C from m_pah->m_cpositions
//    //m_pah->m_cpositions.erase(m_pah->m_cpositions.find(C_1->coords));
//    // delete Carbon object
//    delete C_1;
//    if(m_pah->m_carbonList.erase(C_1) == 0)
//        std::cout<<"ERROR: removeC: NOT ERASING ANY POINTERS!\n";
//    if(!bulk) { // if desorption, decrease C count accordingly
//        addCount(-1, 0); 
//    }
//}
//! Adds a site before site b4site
Spointer PAHProcess::addSite(kmcSiteType stype, Spointer& b4site) {
    // Create new Site
    Site st;
    // Set details of site
    st.type = stype;
	//st.C1 = C_1;
	//st.C2 = C_2;
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
//! Adds a site with members at end of SiteList
Spointer PAHProcess::addSite(kmcSiteType stype) {
    // Create new Site
    Site st;
    // Set details of site
    st.type = stype;
    //st.C1 = C_1;
    //st.C2 = C_2;
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
void PAHProcess::addR5toSite(Spointer& st) {
    int stype = (int) st->type;
    if(stype<5) // i.e. uncombined principal sites (eg FE)
        stype += 6;
    else if(stype < 10) // i.e. principal sites with 5R at one side (eg RFE)
        stype += 4;
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
    //st->C1 = Carb1;
    //st->C2 = Carb2;
}
//! Changes site type into combined site without R5 (e.g. RFE -> FE)
void PAHProcess::remR5fromSite(Spointer& st) {
    int stype = (int) st->type;
    if(stype<13 && stype>9) // i.e. principal sites with 5R at both sides (eg RFER)
        stype -= 4;
    else if(stype<10 && stype>5) //i.e. principal sites with 5R at one side (eg RFE)
        stype -= 6;
    else {
        std::cout<<"ERROR: rem5RfromSite: illegal site type to remove 5R from.\n";
        return;}
    // removes site from m_pah->m_siteMap (principal site)
    delSiteFromMap(st->type, st);
    // change site type
    st->type = (kmcSiteType) stype;
    // add site to m_pah->m_siteMap
    m_pah->m_siteMap[st->type].push_back(st);
    // update member C
    //st->C1 = Carb1;
    //st->C2 = Carb2;
}
//! Changes site type into another site type
void PAHProcess::convSiteType(Spointer& st, kmcSiteType t) {
    // removes site from m_pah->m_siteMap (principal site)
    delSiteFromMap(st->type, st);
    // change site type
    st->type = t;
    // add site to m_pah->m_siteMap
    m_pah->m_siteMap[t].push_back(st);
    // update member C
    //st->C1 = Carb1;
    //st->C2 = Carb2;
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
//void PAHProcess::updateA(Cpointer C, char sp) {
//    if(C->bridge) {
//        C->A = 'C'; // a bridge is obviously connected to 3 other C atoms
//        C->bondAngle2=normAngle(C->bondAngle1+120);
//    }else {
//        // Check if C is a reactive surface carbon by finding angle it makes with neighbours
//        angletype x = normAngle(C->bondAngle1 - C->C1->bondAngle1);
//        // Sets A if C is a reactive surface carbon
//        if (x<0) {
//            C->A = sp;
//            C->bondAngle2=normAngle(C->bondAngle1+120);
//        } else {
//            C->A = 'C';
//            C->bondAngle2=normAngle(C->bondAngle1-120);
//        }
//    }
//}
//! Overload function, updateA for all C from C_1 to C_2 inclusive
//void PAHProcess::updateA(Cpointer C_1, Cpointer C_2, char spc) {
//    // Iterate through each Carbon atom
//    Cpointer now = C_1;
//    Cpointer prev;
//    // iterator will not cross bridges if C_1 is a bridge
//    if(C_1->bridge) prev = C_1->C3;
//    else prev = C_1->C1;
//    // a count to check if updateA is looping uncontrollably
//    unsigned int count = 0;
//    do {
//        updateA(now, spc); 
//        Cpointer oldnow = now;
//        now = moveCPointer(prev, now);
//        prev = oldnow;
//        count++;
//        if(count == 1000)
//            cout << "WARNING: Loop count has reached 1000 (Sweep::KMC_ARS::PAHProcess::updateA)\n";
//        else if(count > 5000){
//            cout << "ERROR: Loop count has reached more than 5000. Stopping simulation now.. (Sweep::KMC_ARS::PAHProcess::updateA)\n";
//            std::ostringstream msg;
//            msg << "ERROR: Possibility of function looping indefinitely."
//                << " (Sweep::KMC_ARS::PAHProcess::updateA)";
//            throw std::runtime_error(msg.str());
//            assert(false);
//        }
//
//    }
//    while (prev != C_2);
//}
//! Update all principal sites
//void PAHProcess::updateSites() {
//    //cout<<"Clearing Sites..\n";
//    m_pah->m_siteList.clear(); //cout << "m_pah->m_siteList cleared!\n";
//    m_pah->m_siteMap.clear();//cout << "m_pah->m_siteMap cleared!\n";
//    //check if m_pah->m_cfirst is not bonded to C, if not iterate until reach one with H
//    for(int i=0; i!=200;i++) {
//        if (m_pah->m_cfirst->A != 'C') {
//            break;
//        } else {
//            m_pah->m_clast = m_pah->m_cfirst;
//            m_pah->m_cfirst = m_pah->m_cfirst->C2;
//        }
//        if (i==199) {
//            std::cout << "ERROR: updateSites(): 200th loop reached without finding C-H\n";
//            std::ostringstream msg;
//            msg << "ERROR: Possibility of function looping indefinitely."
//                << " (Sweep::KMC_ARS::PAHProcess::updateSites)";
//            throw std::runtime_error(msg.str());
//            assert(false);
//            return;
//        }
//    }
//    Cpointer prev = m_pah->m_cfirst;
//    Cpointer now = m_pah->m_cfirst->C2;
//    Cpointer siteC1=prev;
//    kmcSiteType sitetype;
//    unsigned int bulk = 0;
//    do {
//        if(now->A != 'C') {
//            switch (bulk) {
//            case 0:
//                sitetype = FE; break;
//            case 1:
//                sitetype = ZZ; break;
//            case 2:
//                sitetype = AC; break;
//            case 3:
//                sitetype = BY5; break;
//            case 4:
//                sitetype = BY6; break;
//            default:
//                std::cout << "ERROR: updateSites(): more than 4 bulk carbon atoms detected\n"; break;
//            }
//            Spointer end_of_siteList = m_pah->m_siteList.end();
//            addSite(sitetype, siteC1, now, end_of_siteList);
//            bulk = 0;
//            siteC1 = now;
//        }else {
//            bulk++;
//        }
//        Cpointer oldnow = now;
//        now = moveCPointer(prev, now);
//        prev = oldnow;
//    } while (prev != m_pah->m_cfirst);
//    // Stores iterators into vectors
//    //m_pah->m_siteMap.clear();
//    //Spointer it;
//    //for(it=m_pah->m_siteList.begin(); it!=m_pah->m_siteList.end(); it++) {
//        //m_pah->m_siteMap[it->type].push_back(it);
//    //}
//    setCount(m_pah->m_counts.first, (int) m_pah->m_siteList.size());
//    //stericHindrance();
//    //cout << "Principal Sites Updated, H count: "<< m_pah->m_counts[1] << "..\n";
//}
//! Updates particular site
void PAHProcess::updateSites(Spointer& st, // site to be updated
//                               ,Cpointer Carb1, Cpointer Carb2, // new C members
                               int bulkCchange) { // addition to number of bulk C in site
//    // check if site type change is valid (as long as site still principal site)
      int stype = (int) st->type;
//    if(stype < 6) {
//        if (!(stype == 5 && bulkCchange == 0)) {
//    if((stype + bulkCchange) < 0 || (stype + bulkCchange) > 4) {
//        //printSites(st);//++++
//        cout << "ERROR: updateSites: Bulk C change invalid (Principal)\n";
//        std::ostringstream msg;
//            msg << "ERROR: Bulk C change invalid (Principal). Trying to add "
//                << bulkCchange << " bulk C to a " << kmcSiteName(st->type)
//                << " (Sweep::KMC_ARS::PAHProcess::updateSites)";
//            throw std::runtime_error(msg.str());
//            assert(false);
//        return;
//            }
//        }
//    }
//    if(stype < 13 && stype > 5) {
//        if((stype + bulkCchange) < 6 || (stype + bulkCchange) > 12) {
//            //printSites(st);//++++
//            cout << "ERROR: updateSites: Bulk C change invalid (with R5)\n";
//            std::ostringstream msg;
//            msg << "ERROR: Bulk C change invalid (with R5)."
//                << " (Sweep::KMC_ARS::PAHProcess::updateSites)";
//            throw std::runtime_error(msg.str());
//            assert(false);
//            return;
//        }
//    }
//    // removes site from m_pah->m_siteMap (principal site)
    delSiteFromMap(st->type, st);
	if (stype > 0){
		st->type = (kmcSiteType)((int)st->type + bulkCchange);
	}
	else{
		st->type = (kmcSiteType)((int)st->type - bulkCchange);
	}

    // add site to m_pah->m_siteMap
	  //NICK TO DO - Figure out logic for updating site type. This may already work correctly
    m_pah->m_siteMap[st->type].push_back(st);
//    // update member C
//    st->C1 = Carb1;
//    st->C2 = Carb2;
}
void PAHProcess::MergeSites(PAHProcess& rhs, rng_type &rng) {
	Spointer Sp1, Sp2;
	bool fail = false;
	int guard = 0;
	bool sel1 = false, sel2 = false;

	if (m_pah->m_siteMap[FE2].size() > 0){
		while (guard < 10 && !sel1){
			guard++;
			Sp1 = chooseRandomSite(FE2, rng);
			if (moveIt(Sp1, 1)->comb == FE2){
				Spointer Sp3 = moveIt(Sp1, 1);
				Sp1 = Sp3;
			}
			//if (Sp1->type != NFE && moveIt(Sp1, -1)->type != NFE){
				sel1 = true;
			//}
		}

	}
	else{
		fail = true;
	}
	if (!sel1){
		fail = true;
	}
	if (!fail){
		guard = 0;
		if (rhs.m_pah->m_siteMap[FE2].size() > 0){
			while (guard < 10 && !sel2){
				Sp2 = rhs.chooseRandomSite(FE2, rng);
				if (rhs.moveIt(Sp2, 1)->comb == FE2){
					Sp2 = rhs.moveIt(Sp2, 1);
				}
				//if (Sp2->type != NFE && rhs.moveIt(Sp2, -1)->type != NFE){
					sel2 = true;
				//}
			}
		}
		else{
			fail = true;
		}
		if (!sel2){
			fail = true;
		}
	}

	if (!fail){
		//convSiteType(moveIt(Sp1,-1), ACBR);
		convSiteType(moveIt(Sp1, -1), AC);
		Spointer st3;
		for (st3 = rhs.moveIt(Sp2,1); st3 != rhs.m_pah->m_siteList.end(); st3++){
			addSite(st3->type, Sp1);
		}
		for (st3 = rhs.m_pah->m_siteList.begin(); st3 != rhs.moveIt(Sp2,-1); st3++){
			addSite(st3->type, Sp1);
		}
		//convSiteType(Sp1, ACBL);
		convSiteType(Sp1, AC);
		Sp1 = m_pah->m_siteList.begin();
		updateCombinedSites();
		updateHinderedSites();
	}

}
void PAHProcess::updateHinderedSites() {
	
	bool hindered;
	int Sides;
	int count;
	double Angle, target;
	Spointer st1;
	SpointerRev st2;
	//Now, re-assign hindered site types
	for (Spointer st = m_pah->m_siteList.begin(); st != m_pah->m_siteList.end(); st++) {
		Sides = 0;
		Angle = 0.0;
		count = 0;
		hindered = false;
		switch (st->type) {
		case FE:
		case NFE:
		case AC:
		case NAC:
		case NACBL:
		case NACBR:
			for (st1 = moveIt(st, 1); st1 != m_pah->m_siteList.end(); st1++) {
				count++;
				Sides += SiteSides(st1);
				Angle += SiteAngle(st1);
				if (Sides < 5 || (Sides >5 && Sides < 8)){
					target = -100.0;
				}
				else if ((Sides+1)%2 == 0){
					target = ((double)Sides - 1.0)*180.0 -240.0 - ((double) count - 1.0)*240;
				}
				else{
					target = ((double)Sides)*180.0 - 480.0 - ((double)count - 1.0) * 240;
				}
				if (Angle == target) {
					hindered = true;
				}
				if (hindered) break;
			}
			if (hindered) break;
			for (st1 = m_pah->m_siteList.begin(); st1 != st; st1++) {
				count++;
				Sides += SiteSides(st1);
				Angle += SiteAngle(st1);
				if (Sides < 5 || (Sides >5 && Sides < 8)){
					target = -100.0;
				}
				else if ((Sides + 1) % 2 == 0){
					target = ((double)Sides - 1.0)*180.0 - 240.0 - ((double)count - 1.0) * 240;
				}
				else{
					target = ((double)Sides)*180.0 - 480.0 - ((double)count - 1.0) * 240;
				}
				if (Angle == target) {
					hindered = true;
				}
				if (hindered) break;
			}
			if (hindered) break;
			for (st1 = moveIt(st, 1); st1 != m_pah->m_siteList.end(); st1++) {
				count--;
				Sides -= SiteSides(st1);
				Angle -= SiteAngle(st1);
				if (Sides < 5 || (Sides >5 && Sides < 8)){
					target = -100.0;
				}
				else if ((Sides + 1) % 2 == 0){
					target = ((double)Sides - 1.0)*180.0 - 240.0 - ((double)count - 1.0) * 240;
				}
				else{
					target = ((double)Sides)*180.0 - 480.0 - ((double)count - 1.0) * 240;
				}
				if (Angle == target) {
					hindered = true;
					break;
				}
				if (hindered) break;
			}
			if (hindered) break;
			for (st1 = m_pah->m_siteList.begin(); st1 != st; st1++) {
				count--;
				Sides -= SiteSides(st1);
				Angle -= SiteAngle(st1);
				if (Sides < 5 || (Sides >5 && Sides < 8)){
					target = -100.0;
				}
				else if ((Sides + 1) % 2 == 0){
					target = ((double)Sides - 1.0)*180.0 - 240.0 - ((double)count - 1.0) * 240;
				}
				else{
					target = ((double)Sides)*180.0 - 480.0 - ((double)count - 1.0) * 240;
				}
				if (Angle == target) {
					hindered = true;
					break;
				}
				if (hindered) break;
			}
			if (hindered) break;

		}

		if (hindered){
			if (st->type > 0){
				convSiteType(st, (kmcSiteType) (0 - (int)st->type));
			}
		}
		else{
			if (st->type < 0){
				convSiteType(st, (kmcSiteType)(0 - (int)st->type));
			}
		}
	}
}

int PAHProcess::SiteSides(Spointer& stt) {
	switch (stt->type){
	case FE:
	case NFE:
	case ERFE:
	case ERFEER:
		return 1;
	case ZZ:
	case NZZ:
	case ERZZ:
	case ERZZER:
		return 2;
	case AC:
	case NAC:
	case ERAC:
	case ERACER:
	case ACBL:
	case ACBR:
	case NACBL:
	case NACBR:
		return 3;
	case BY5:
	case ERBY5:
	case ER5:
		return 4;
	case BY6:
		return 5;
	}
}

double PAHProcess::SiteAngle(Spointer& stt) {
	switch (stt->type){
	case FE:
	case NFE:
	case ERFE:
	case ERFEER:
		return 0.0;
	case ZZ:
	case NZZ:
	case ERZZ:
	case ERZZER:
		return 120.0;
	case AC:
	case NAC:
	case ERAC:
	case ERACER:
	case ACBL:
	case ACBR:
	case NACBL:
	case NACBR:
		return 240.0;
	case BY5:
	case ERBY5:
	case ER5:
		return 360.0;
	case BY6:
		return 480.0;
	}
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
		if (i->type == FE || i->type == NFE) {
			indicesOfFE.push_back(loopIndex);
		}
		else {
			indicesOfNonFE.push_back(loopIndex);
		}
		loopIndex++;
    }

	Spointer S1;
    S1 = m_pah->m_siteList.begin();

	for (int i = 0 ; i < indicesOfFE.size(); ++i) {
        Spointer S2 = moveIt(S1,indicesOfFE[i]);
		updateCombinedSites(S2);
	}

	for (int i = 0 ; i < indicesOfNonFE.size(); ++i) {
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
    switch(st->type) {
	case FE: 
	case NFE:
        // Check for FE3 (if there's FE on each side of the FE)
		if ((moveIt(st, 1)->type == FE || moveIt(st, 1)->type == NFE) && (moveIt(st, -1)->type == FE || moveIt(st, -1)->type == NFE)) {
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
		else if ((moveIt(st, 1)->type == FE || moveIt(st, 1)->type == NFE) || (moveIt(st, -1)->type == FE || moveIt(st, -1)->type == NFE)){
            Spointer S1,S2;
            S1 = moveIt(st,-1); S2 = moveIt(st, 1);
            // Check if that FE is not a FE3
			if ((S2->type == FE || S2->type == NFE) && moveIt(S2, 1)->type != FE && moveIt(S2, 1)->type != NFE) {
                st->comb = FE2;
                m_pah->m_siteMap[FE2].push_back(st);
		//
                // An FE2 site is a combined site where an FE site has an FE site only 
                // For example, for ZZ - FE - FE - ZZ, both of the FE sites has a combi
                //
                // Zig-zag oxidation reactions are based on the number of side-by-side 
                // If these reactions were based on the number of FE2 sites, we would o
                // So we can either calculate the rate based on the number of FE2 sites
                // or - as has been done here - remove half of the FE2 sites.
                //
                if(S2->comb == FE2) delSiteFromMap(S2->comb, st);
                //
                if(S2->comb != FE2) updateCombinedSites(S2);
			}
			else if ((S1->type == FE || S1->type == NFE) && moveIt(S1, -1)->type != FE && moveIt(S1, -1)->type != NFE) {
                st->comb = FE2;
                m_pah->m_siteMap[FE2].push_back(st);
                if(S1->comb == FE2) delSiteFromMap(S1->comb, st);
                if(S1->comb != FE2) updateCombinedSites(S1);
            } else
                st->comb = None;
            break;
        }
        // Check for FE_HACA
        else if(moveIt(st,1)->type != FE && moveIt(st,-1)->type != FE && moveIt(st,1)->type != RFE && moveIt(st,-1)->type != RFE
			&& moveIt(st, 1)->type != NFE && moveIt(st, -1)->type != NFE) {
            //if(st->C1->C1->bridge || st->C2->C2->bridge)
            st->comb = FE_HACA;
            m_pah->m_siteMap[FE_HACA].push_back(st);
            break;
        }
        else st->comb = None;
        break;
    case AC:
        // Check for AC_FE3
		//NICK TO DO - Properly check for bridges. Might not be needed because bridges will be new principal site type
        if((moveIt(st,2)->comb==FE3 || moveIt(st,-2)->comb==FE3) && (moveIt(st,3)->comb!=FE3 || moveIt(st,-3)->comb!=FE3)) {
            st->comb = AC_FE3;
            m_pah->m_siteMap[AC_FE3].push_back(st);
            break;
        }else st->comb = None;
        break;
    case BY5:
    // Check for BY5_FE3
		//NICK TO DO - Properly check for bridges. Might not be needed because bridges will be new principal site type
        if((moveIt(st,2)->comb==FE3 || moveIt(st,-2)->comb==FE3) && (moveIt(st,3)->comb!=FE3 || moveIt(st,-3)->comb!=FE3)) {
            st->comb = BY5_FE3;
            m_pah->m_siteMap[BY5_FE3].push_back(st);
            break;
        }else st->comb = None;
        break;
  //  case eR5:
  //  // Check for eR5_FE3
		////NICK TO DO - Properly check for bridges. Might not be needed because bridges will be new principal site type
		//if ((moveIt(st, -1)->type == eRFE && moveIt(st, -2)->comb == FE2) || (moveIt(st, 1)->type == eRFE && moveIt(st, 2)->comb == FE2)) {
  //          st->comb = eR5_FE3;
  //          m_pah->m_siteMap[BY5_FE3].push_back(st);
  //          break;
  //      }else st->comb = None;
  //      break;
    case RAC:
    // Check for RAC_FE3
		//NICK TO DO - Properly check for bridges. Might not be needed because bridges will be new principal site type
        if((moveIt(st,2)->comb==FE3 || moveIt(st,-2)->comb==FE3) && (moveIt(st,3)->comb!=FE3 || moveIt(st,-3)->comb!=FE3)) {
            st->comb = RAC_FE3;
            m_pah->m_siteMap[RAC_FE3].push_back(st);
            break;
        }else st->comb = None;
        break;
    default:
        st->comb = None;
        break;
    }
    if((ori_comb == FE3) && (moveIt(st,1)->comb!=FE3 || moveIt(st,-1)->comb!=FE3)) {
        Spointer n1, n2;
        n1 = moveIt(st, -2);
        n2 = moveIt(st, 2);
        if((n1->type == AC || n1->type == BY5)) updateCombinedSites(n1); //NICK TO DO - Add in check for RAC?
        if((n2->type == AC || n2->type == BY5)) updateCombinedSites(n2); //NICK TO DO - Add in check for RAC?
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
    intpair rings;
    intpair CH;
    int rings_Embedded;

    // Structure for Benzene
    std::string BENZENE_Sites = "FE,FE,FE,FE,FE,FE";
    intpair BENZENE_Rings(1,0);
    intpair BENZENE_CH(6,6);
    int BENZENE_RINGS_EMBEDDED = 0;
    // Structure for Naphthalene
    std::string NAPHTHALENE_Sites = "FE,FE,FE,ZZ,FE,FE,FE,ZZ";
    intpair NAPHTHALENE_Rings(2,0);
    intpair NAPHTHALENE_CH(10,0);
    int NAPHTHALENE_RINGS_EMBEDDED = 0;
    // Structure for Pyrene
    std::string PYRENE_Sites = "ZZ,FE,FE,ZZ,FE,ZZ,FE,FE,ZZ,FE";
    intpair PYRENE_Rings(4,0);
    intpair PYRENE_CH(16,10);
    int PYRENE_RINGS_EMBEDDED = 0;
    // Structure for BENZOPYRENE
    std::string BENZOPYRENE_Sites = "ZZ,FE,FE,ZZ,FE,ZZ,ZZ,FE,FE,FE,AC,FE";
    intpair BENZOPYRENE_Rings(5,0);
    intpair BENZOPYRENE_CH(20,12);
    int BENZOPYRENE_RINGS_EMBEDDED = 0;
    // Structure for Coronene
    std::string CORONENE_Sites = "FE,ZZ,FE,ZZ,FE,ZZ,FE,ZZ,FE,ZZ,FE,ZZ";
    intpair CORONENE_Rings(7,0);
    intpair CORONENE_CH(24,12);
    int CORONENE_RINGS_EMBEDDED = 0;
    // Test Structure
    std::string TEST_Sites = "FE,AC,FE,ZZ,RFE,R5,RAC,RFE,R5,RFE,FE,AC,FE,FE,FE,FE,BY5,AC,FE,BY5,FE,ZZ,FE,FE";
    intpair TEST_Rings(8,2);
    intpair TEST_CH(18,10);
    int TEST_RINGS_EMBEDDED = 0;
    // Choose structure
    switch(ss) {
    case BENZENE_C:
        chosen = BENZENE_Sites;
        rings = BENZENE_Rings;
        CH = BENZENE_CH;
        rings_Embedded = BENZENE_RINGS_EMBEDDED;
        break;
    case NAPHTHALENE_C:
        chosen = NAPHTHALENE_Sites;
        rings = NAPHTHALENE_Rings;
        CH = NAPHTHALENE_CH;
        rings_Embedded = NAPHTHALENE_RINGS_EMBEDDED;
        break;
    case PYRENE_C:
        chosen = PYRENE_Sites;
        rings = PYRENE_Rings;
        CH = PYRENE_CH;
        rings_Embedded = PYRENE_RINGS_EMBEDDED;
        break;
    case BENZOPYRENE_C:
        chosen = BENZOPYRENE_Sites;
        rings = BENZOPYRENE_Rings;
        CH = BENZOPYRENE_CH;
        rings_Embedded = BENZOPYRENE_RINGS_EMBEDDED;
        break;
    case CORONENE_C:
        chosen = CORONENE_Sites;
        rings = CORONENE_Rings;
        CH = CORONENE_CH;
        rings_Embedded = CORONENE_RINGS_EMBEDDED;
        break;
    case TEST_STRUCT:
        chosen = TEST_Sites;
        rings = TEST_Rings;
        CH = TEST_CH;
        rings_Embedded = TEST_RINGS_EMBEDDED;
        break;
    default: 
            std::cout<<"ERROR: Starting Structure undefined.. (PAHProcess::initialise)\n\n";
            assert(false);
            abort();
    }
    // Create Structure
    return initialise(chosen, rings.first, rings.second, rings_Embedded, CH.first, CH.second);
    //printSites(m_pah->m_siteList.begin());
}

/*!
 * @param[in]    ss    Initialise the PAH structure corresponding to the starting structure, ss
 *
 * @return       Initialized PAH structure
 */
PAHStructure& PAHProcess::initialise(StartingStructure ss){
	if (m_pah == NULL) {
		PAHStructure* pah = new PAHStructure();
		m_pah = pah;
	}
    //}else if(m_pah->m_cfirst != NULL)
    //     m_pah->clear();
    switch(ss) {
        //Cpointer newC;
        
    case BENZENE_C:
        //cout << "newC pointer created\n";
        // add first C atom
        //m_pah->m_cfirst = addC();
        ////cout << "m_cfirst is " << m_cfirst << '\n';
        //// adds next C atoms according to structure
        //newC = addC(m_pah->m_cfirst, 0, 0);
        ////cout << "newC is " << newC << '\n';
        //newC = addC(newC, -60, 0);
        ////cout << "newC is " << newC << '\n';
        //newC = addC(newC, -120, 0);
        ////cout << "newC is " << newC << '\n';
        //newC = addC(newC, -180, 0);
        ////cout << "newC is " << newC << '\n';
        //// adds the last C atom, with bond angle towards m_cfirst
        //m_pah->m_clast = addC(newC, 120, 60);
        ////cout << "m_clast is " << m_clast << '\n';
        //// closes structure
        //connectToC(m_pah->m_clast, m_pah->m_cfirst);
        ////cout << "m_clast and m_cfirst connected, updating A\n";
        //// update H atoms
        //updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(BENZENE_C, BENZENE_H);
        // set ring counts
        m_pah->m_rings = 1;
        m_pah->m_rings5_Lone = 0;
        m_pah->m_rings5_Embedded = 0;
        // update all sites and combined sites
        //updateSites();
		addSite(FE);
		addSite(FE);
		addSite(FE);
		addSite(FE);
		addSite(FE);
		addSite(FE);
        updateCombinedSites();
        //cout << "Benzene Initialised!\n";
        break;
    case NAPHTHALENE_C:
        // add first C atom
        //m_pah->m_cfirst = addC();
        //// adds next C atoms according to structure
        //newC = addC(m_pah->m_cfirst, 0, 0);
        //newC = addC(newC, 60, 0);
        //newC = addC(newC, 0, 0);
        //newC = addC(newC, -60, 0);
        //newC = addC(newC, -120, 0);
        //newC = addC(newC, -180, 0);
        //newC = addC(newC, -120, 0);
        //newC = addC(newC, -180, 0);
        //// adds the last C atom, with bond angle towards m_cfirst
        //m_pah->m_clast = addC(newC, 120, 60);
        //// closes structure
        //connectToC(m_pah->m_clast, m_pah->m_cfirst);
        //// update H atoms
        //updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(10, 8);
        // set ring counts
        m_pah->m_rings = 2;
        m_pah->m_rings5_Lone = 0;
        m_pah->m_rings5_Embedded = 0;
        // update all sites and combined sites
		//updateSites();
		addSite(FE);
		addSite(FE);
		addSite(FE);
		addSite(ZZ);
		addSite(FE);
		addSite(FE);
		addSite(FE);
		addSite(ZZ);
        updateCombinedSites();
        //cout << "Naphthalene Initialised!\n";
        break;
    case PYRENE_C:
        // add first C atom
        //m_pah->m_cfirst = addC();
        //// adds next C atoms according to structure
        //newC = addC(m_pah->m_cfirst, 0, 0);
        //newC = addC(newC, -60, 0);
        //newC = addC(newC, 0, 0);
        //newC = addC(newC, -60, 0);
        //newC = addC(newC, -120, 0);
        //newC = addC(newC, -180, 0);
        //newC = addC(newC, -120, 0);
        //newC = addC(newC, -180, 0);
        //newC = addC(newC, 120, 0);
        //newC = addC(newC, 180, 0);
        //newC = addC(newC, 120, 0);
        //newC = addC(newC, 60, 0);
        //// adds the last C atom, with bond angle towards m_cfirst
        //m_pah->m_clast = addC(newC, 0, 60);
        //// closes structure
        //connectToC(m_pah->m_clast, m_pah->m_cfirst);
        //// update H atoms
        //updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(PYRENE_C, PYRENE_H);
        // set ring counts
        m_pah->m_rings = 4;
        m_pah->m_rings5_Lone = 0;
        m_pah->m_rings5_Embedded = 0;
        // update all sites and combined sites
		//updateSites();
		addSite(FE);
		addSite(ZZ);
		addSite(FE);
		addSite(FE);
		addSite(ZZ);
		addSite(FE);
		addSite(ZZ);
		addSite(FE);
		addSite(FE);
		addSite(ZZ);
        updateCombinedSites();
        //cout << "Pyrene Initialised!\n";
        break;
    case BENZOPYRENE_C:
        // add first C atom
  //      m_pah->m_cfirst = addC();
  //      // adds next C atoms according to structure
  //      newC = addC(m_pah->m_cfirst, 0, 0);
  //      newC = addC(newC, -60, 0);
  //      newC = addC(newC, 0, 0);
  //      newC = addC(newC, -60, 0);
  //      newC = addC(newC, -120, 0);
  //      newC = addC(newC, -180, 0);
  //      newC = addC(newC, -120, 0);
  //      newC = addC(newC, -180, 0);
  //      newC = addC(newC, 120, 0);
  //      newC = addC(newC, 180, 0);
  //      newC = addC(newC, 120, 0);
		//newC = addC(newC, 180, 0);
		//newC = addC(newC, 120, 0);
		//newC = addC(newC, 60, 0);
		//newC = addC(newC, 0, 0);
		//newC = addC(newC, -60, 0);
  //      // adds the last C atom, with bond angle towards m_cfirst
  //      m_pah->m_clast = addC(newC, 0, 60);
  //      // closes structure
  //      connectToC(m_pah->m_clast, m_pah->m_cfirst);
  //      // update H atoms
  //      updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');
        // set C & H counts
        setCount(BENZOPYRENE_C, BENZOPYRENE_H);
        // set ring counts
        m_pah->m_rings = 5;
        m_pah->m_rings5_Lone = 0;
        m_pah->m_rings5_Embedded = 0;
        // update all sites and combined sites
		//updateSites();
		addSite(FE);
		addSite(FE);
		addSite(AC);
		addSite(FE);
		addSite(ZZ);
		addSite(FE);
		addSite(FE);
		addSite(ZZ);
		addSite(FE);
		addSite(ZZ);
		addSite(ZZ);
		addSite(FE);

        updateCombinedSites();
        //cout << "Benzopyrene Initialised!\n";
        break;
    case BY5_C:
        /**
         * Add the first carbon atom.
         */ 
        //m_pah->m_cfirst = addC();

        ///**
        // * Keep adding carbon atoms to describe the PAH.
        // */ 
        //newC = addC(m_pah->m_cfirst, 0, 0);
        //newC = addC(newC, -60, 0);
        //newC = addC(newC, 0, 0);
        //newC = addC(newC, -60, 0);
        //newC = addC(newC, -120, 0);
        //newC = addC(newC, -180, 0);
        //newC = addC(newC, 120, 0);
        //newC = addC(newC, -180, 0);
        //newC = addC(newC, -120, 0);
        //newC = addC(newC, -60, 0);
        //newC = addC(newC, -120, 0);
        //newC = addC(newC, -180, 0);
        //newC = addC(newC, 120, 0);
        //newC = addC(newC, 60, 0);
        //newC = addC(newC, 120, 0);
        //newC = addC(newC, 60, 0);

        ///**
        // * Add the last carbon atom.
        // */ 
        //m_pah->m_clast = addC(newC, 0, 60);

        ///**
        // * Connect the last and first carbon atoms together to close the PAH structure.
        // */ 
        //connectToC(m_pah->m_clast, m_pah->m_cfirst);

        ///**
        // * Loop through all the edge carbon atoms and update their status, i.e., whether they are connected to a hydrogen atom.
        // */ 
        //updateA(m_pah->m_cfirst, m_pah->m_clast, 'H');

        /**
         * Set the total number of carbon (BY5_C) and hydrogen atoms (BY5_H).
         */ 
        setCount(BY5_C, BY5_H);

        /**
         * Set the number of 6 and 5-member aromatic rings. m_rings and m_rings5, respectively.
         */ 
        m_pah->m_rings = 4;
        m_pah->m_rings5_Lone = 0;
        m_pah->m_rings5_Embedded = 0;

        /**
         * Update the basic site types (updateSites) and the combined site types (updateCombinedSites).
         */ 
		//updateSites();
		//NICK TO DO - Must add in sites into site vectors
        updateCombinedSites();

        //cout << "5-member bay site initialised!\n";
        break;
     default: 
            std::cout<<"ERROR: Starting Structure undefined.. (PAHProcess::initialise)\n\n";
            assert(false);
            abort();
    }
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
PAHStructure& PAHProcess::initialise(std::string siteList_str, int R6_num, int R5_num_Lone, int R5_num_Embedded, int numC, int numH){
	if (m_pah == NULL) {
		PAHStructure* pah = new PAHStructure();
		m_pah = pah;
	}
    //}else if(m_pah->m_cfirst != NULL)
    //    m_pah->clear();
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
                throw std::runtime_error(msg.str());
                assert(false);
        }
        siteList_vec.push_back(temp);
    }
    createPAH(siteList_vec, R6_num, R5_num_Lone, R5_num_Embedded, numC, numH);
    return *m_pah;
}

// Create Structure from vector of site types
void PAHProcess::createPAH(std::vector<kmcSiteType>& vec, int R6, int R5_Lone, int R5_Embedded, int numC, int numH) {
    // current C, bondangle and coordinates
    //Cpointer newC=addC();
    //m_pah->m_cfirst = newC;
    //m_pah->m_clast = NULLC;
    // number of bulk C to be added
    //int bulkC;
    // type of site; if type 0, basic site types (FE - BY6); if type 1, R5 and basic sites with
    // a R5 at one side (RFE - RBY5); if type 2, basic sites wit R5 at each side (RFER - RACR)
    unsigned short int site_t;
    // start drawing..
    
    /**
     * A means of tracking whether the code has reached the last element of the site vector, vec. 
     */
    size_t final_iter = vec.size();
    --final_iter;
    
    for(size_t i=0; i<vec.size(); i++) {
        //Cpointer S_C1 = newC;
        //// get number of bulk C to be added and site type
        //if((int)vec[i] <= 4) {
        //    bulkC = (int) vec[i]; site_t = 0;
        //}else if((int)vec[i] >= 5 && (int)vec[i] <= 9) {
        //    bulkC = (int) vec[i] - 5; site_t = 1;
        //}else if((int)vec[i] >= 10 && (int)vec[i] <= 12) {
        //    bulkC = (int) vec[i] - 8; site_t = 2;
        //}
        ///**
        // * If this condition is true, vec[i] is the site tye ACR5.
        // * It is essentially a special type of armchair site (a basic site type and associated with two bulk carbon atoms).
        // */
        //else if((int)vec[i] == 18) {
        //    bulkC = (int) vec[i] - 16; site_t = 0;

        //}else {
        //    cout << "createPAH: Combined site types in list of sites. Please use only\n"
        //        << "principal site types\n";
        //    std::ostringstream msg;
        //    msg << "ERROR: Combined site types found in list of sites."
        //        << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
        //        throw std::runtime_error(msg.str());
        //        assert(false);
        //    return;
        //}
        //kmcSiteType prevType;
        //if(i==0) prevType = vec.back();
        //else prevType = vec[i-1];
        //switch(site_t) {
        //case 0:
        //    newC = drawType0Site(newC, bulkC, i == final_iter); break;
        //case 1:
        //    newC = drawType1Site(newC, bulkC, i == final_iter); break;
        //case 2:
        //    newC = drawType2Site(newC, bulkC, i == final_iter); break;
        //default:
        //    cout << "createPAH: Invalid site_t number...\n";
        //    std::ostringstream msg;
        //    msg << "ERROR: invalid site classification."
        //        << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
        //        throw std::runtime_error(msg.str());
        //        assert(false);
        //    return;
        //}
        addSite(vec[i]);
    }
    // check if PAH closes correctly
	// NICK TO DO - Figure out new logic to see if PAH closes properly
    //if(m_pah->m_clast == NULLC || newC != m_pah->m_cfirst) {
    //    // PAH did not close properly. invalid structure
    //    cout << "createPAH: PAH did not close properly. Could be problem "
    //        <<"with site list input...\n";
    //    std::ostringstream msg;
    //    msg << "ERROR: PAH did not close properly.."
    //        << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
    //    //saveDOT("KMC_DEBUG/KMC_PAH_X_CLOSE.dot");
    //    throw std::runtime_error(msg.str());
    //    assert(false);
    //    return;
    //}
    m_pah->m_rings = R6;
    m_pah->m_rings5_Lone = R5_Lone;
    m_pah->m_rings5_Embedded = R5_Embedded;
    //for(Ccontainer::iterator i = m_pah->m_carbonList.begin();
    //    i != m_pah->m_carbonList.end(); i++)
    //    updateA(*i, 'H');
    //int totalC_num = 2*m_pah->m_rings + (CarbonListSize()+m_pah->m_rings5_Lone+m_pah->m_rings5_Embedded)/2 + numberOfBridges() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded + 1;
	int totalC_num = numC; 
	int totalH_num = numH; 
    m_pah->setnumofC(totalC_num);
	m_pah->setnumofH(totalH_num);
    updateCombinedSites();
}

/** Draw a basic site type which does not have a 5-member aromatic ring on either side.
 *
 * @param[in,out]   Cnow       Pointer to the current edge carbon atom.
 * @param[in]       bulkC      The number of bulk carbon atoms associated with the site.
 * @param[in]       lastSite   Boolean to indicate whether the code is drawing the last site in the site vector.
 *
 * @return Pointer to the current edge carbon atom.
 */
//Cpointer PAHProcess::drawType0Site(Cpointer Cnow, int bulkC, bool lastSite) {
//    // draw site
//    angletype angle = normAngle(Cnow->bondAngle1-60);
//    for(int c=0; c<=bulkC; ++c) {
//        /**
//         * In order to model curved PAHs the code simply tracks the list of site types which makes up the edge of the PAH.
//         * Therefore, the coordinates of the edge carbon atoms (coords) have been made redundant.
//         * This check which depends on the coordinates of the edge carbon atom can no longer be applied.
//         */
//        // check if adding on existing C atom (bridge)
//        //cpair pos = jumpToPos(Cnow->coords, angle);
//        //if(m_pah->m_cpositions.count(pos)) {
//        //    // this coordinate is filled
//        //    Cpointer Cpos = findC(pos);
//        //    if(Cpos != m_pah->m_cfirst) { // it is a bridged C atom
//        //        Cpointer Cbridge = Cpos->C1;
//        //        Cnow->bondAngle1 = angle;
//        //        Cpos->bridge = true; Cbridge->bridge = true;
//        //        Cpos->C3 = Cbridge; Cbridge->C3 = Cpos;
//        //        connectToC(Cnow, Cpos);
//        //        Cbridge->C2 = NULLC;
//        //        angle = normAngle(angle + 120);
//        //        Cbridge->bondAngle1 = angle;
//        //        Cnow = Cbridge;
//        //        --bulkC;
//        //    }else { // reached end of PAH
//        //        Cnow->bondAngle1 = angle;
//        //        connectToC(Cnow, Cpos);
//        //        m_pah->m_clast = Cnow;
//        //        return Cpos;
//        //    }
//        //}else {
//        //}
//
//        /**
//         * If this is the last edge carbon atom to be added and the code is drawing the last site type
//         * close the PAH structure by connecting the first carbon atom and this last carbon atom to be added.
//         * Otherwise, simply add another edge carbon atom.
//         */
//        if(c == bulkC && lastSite){
//                Cnow->bondAngle1 = angle;
//            connectToC(Cnow, m_pah->m_cfirst);
//                m_pah->m_clast = Cnow;
//            return m_pah->m_cfirst;
//            }
//        else {
//            Cnow = addC(Cnow, angle, angle, false);
//            angle = normAngle(angle+60);
//        }
//    }
//    return Cnow;
//}

/** Draw a basic site type which has a 5-member aromatic ring on one side.
 *
 * @param[in,out]   Cnow       Pointer to the current edge carbon atom.
 * @param[in]       bulkC      The number of bulk carbon atoms associated with the site.
 * @param[in]       lastSite   Boolean to indicate whether the code is drawing the last site in the site vector.
 *
 * @return Pointer to the current edge carbon atom.
 */
//Cpointer PAHProcess::drawType1Site(Cpointer Cnow, int bulkC, bool lastSite) {
//    //draw R5 site
//    angletype angle = Cnow->bondAngle1-60;
//    if(bulkC == 0) {
//        angle = normAngle(angle-30);
//
//        /**
//         * In order to model curved PAHs the code simply tracks the list of site types which makes up the edge of the PAH.
//         * Therefore, the coordinates of the edge carbon atoms (coords) have been made redundant.
//         * This check which depends on the coordinates of the edge carbon atom can no longer be applied.
//         * 
//         * Instead the code checks whether this is the last site type to be added.
//         * If so, close the PAH structure by connecting the first carbon atom and this last carbon atom to be added.
//         * Otherwise, simply add another edge carbon atom.
//         */
//        //cpair pos = jumpToPos(Cnow->coords, angle);
//        if(lastSite) {
//            //Cpointer Cpos = findC(pos);
//            Cnow->bondAngle1 = angle;
//            connectToC(Cnow, m_pah->m_cfirst);
//            m_pah->m_clast = Cnow;
//            return m_pah->m_cfirst;
//        }
//        else Cnow = addC(Cnow, angle, normAngle(angle-30), false);
//
//        return Cnow;
//    }else { //draw RXX site
//            return drawType0Site(Cnow, bulkC, lastSite);
//    }
//}

/** Draw a basic site type which has a 5-member aromatic ring either side.
 *
 * @param[in,out]   Cnow       Pointer to the current edge carbon atom.
 * @param[in]       bulkC      The number of bulk carbon atoms associated with the site.
 * @param[in]       lastSite   Boolean to indicate whether the code is drawing the last site in the site vector.
 *
 * @return Pointer to the current edge carbon atom.
 */
//Cpointer PAHProcess::drawType2Site(Cpointer Cnow, int bulkC, bool lastSite) {
//    //angletype angle = Cnow->bondAngle1;
//    //Cnow->bondAngle1 = normAngle(angle-60);
//    return drawType0Site(Cnow, bulkC, lastSite);
//}

//! Finds C atom with specific coordinates
//Cpointer PAHProcess::findC(cpair coordinates) {
//    for(Ccontainer::iterator i=m_pah->m_carbonList.begin();
//        i != m_pah->m_carbonList.end(); ++i) {
//            if((*i)->coords == coordinates)
//                return (*i);
//    }
//    return NULLC;
//}

// Check to validate if coordinates of C matches bond angles
//bool PAHProcess::checkCoordinates() const{
//    // start at first C in site list
//    Cpointer start = m_pah->m_siteList.begin()->C1;
//    Cpointer now = start;
//    Cpointer next = now->C2;
//    unsigned int count=0;
//    do{
//        count++;
//        Cpointer oldnext = next;
//        cpair corr_coords;
//        // first check if next is a valid Carbon pointer
//        if(next == NULL) {
//            cout<<"checkCoordinates() failed at "<<count<<"th C atom..\n"
//                <<"Message: This atom is not a valid pointer to Carbon -- "
//                <<next<<"\n";
//            //printStruct(next);
//            return false;
//        }
//        if(next->bridge) {// if next C is a bridge
//            if(now == next->C1) { // check if current C is next->C1
//                corr_coords = jumpToPos(now->coords, now->bondAngle1);
//                next = next->C3; // cross bridge
//            }else if(now == next->C3) { // check if current C is from bridge
//                corr_coords = jumpToPos(now->coords, now->bondAngle2);
//                next = next->C2;
//            }
//        }else { // if next C is not a bridge
//            corr_coords = jumpToPos(now->coords, now->bondAngle1);
//            next = next->C2;
//        }
//        now = oldnext;
//        // check coordinates
//        if(now->coords != corr_coords) {
//            cout<<"checkCoordinates() failed at "<<count<<"th C atom..\n"
//                <<"Coordinates: ("<<now->coords.first<<','<<now->coords.second<<")\n"
//                <<"Correct      : ("<<corr_coords.first<<','<<corr_coords.second<<")\n";
//            saveDOT(std::string("KMC_DEBUG/error_struct.dot"));
//            return false;
//        }
//    }while(now != start);
//    return true;
//}

// Check to see if all sites are connected to each other
bool PAHProcess::checkSiteContinuity() const {
    std::list<Site>::const_iterator i=m_pah->m_siteList.begin();
	//NICK TO DO - Figure out logic to see if the sites actually connect and PAH closes
    //for(unsigned int k = 0; k!=(unsigned int)m_pah->m_siteList.size();k++) {
    //    Cpointer lhs = i->C2; i++;
    //    if(i == m_pah->m_siteList.end()) i = m_pah->m_siteList.begin();
    //    Cpointer rhs = i->C1;
    //    if(lhs != rhs) 
    //        return false;
    //}
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
				if (i->type != FE && i->type != NFE) {
                    error = true;
                    error_stype_comb = FE3;
                    error_stype = i->type;
                }
                break;
            case FE2:
				if (i->type != FE && i->type != NFE) {
                    error = true;
                    error_stype_comb = FE2;
                    error_stype = i->type;
                }
                break;
            case AC_FE3:
				if (i->type != AC && i->type != NAC) {
                    error = true;
                    error_stype_comb = AC_FE3;
                    error_stype = i->type;
                }
                break;
            case FE_HACA:
				if (i->type != FE && i->type != NFE) {
                    error = true;
                    error_stype_comb = FE_HACA;
                    error_stype = i->type;
                }
                break;
            case BY5_FE3:
                //if((i->type != BY5) && (i->type != eBY5)) {
                if(i->type != BY5) {
                    error = true;
                    error_stype_comb = BY5_FE3;
                    error_stype = i->type;
                }
                break;
            case RAC_FE3:
                if(i->type != RAC) {
                    error = true;
                    error_stype_comb = RAC_FE3;
                    error_stype = i->type;
                }
                break;
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

//! Structure processes: returns success or failure
bool PAHProcess::performProcess(const JumpProcess& jp, rng_type &rng, int PAH_ID)
{
    //printStruct();
    //cout << "Start Performing Process..\n";
    kmcSiteType stp = jp.getSiteType();
    int id = jp.getID();

	//if (PAH_ID == 100320){
	//	cout << "ID is: " << id << endl;
	//}

    // choose random site of type stp to perform process
    Spointer site_perf = chooseRandomSite(stp, rng); //cout<<"[random site chosen..]\n";
    // stores pointers to site Carbon members
    //Cpointer site_C1 = site_perf->C1;
    //Cpointer site_C2 = site_perf->C2;
    //cout << jp.getName() << '\n';
    //printSites(site_perf);
    //cout<<'\t'<<kmcSiteName(site_perf->type)<<' '<<site_C1<<' '<<site_C2<<'\n';
    // find structure change function
    std::ostringstream dotname, dotname2;
    switch(id) {
        case 1:
            proc_G6R_AC(site_perf); break;
        case 2:
            proc_G6R_FE(site_perf); break;
        case 3:
            //dotname << "KMC_DEBUG/" << site_perf->C1 << "_1.dot";
            //dotname2 << "KMC_DEBUG/" << site_perf->C1 << "_2.dot";
            //saveDOT(dotname.str());
            proc_L6_BY6(site_perf);
            //saveDOT(dotname2.str());
            break;
        case 4:
            //proc_PH_benz(site_perf, site_C1, site_C2, rng); break;
        case 5:
            proc_D6R_FE3(site_perf); break;
        case 6:
            proc_O6R_FE3_O2(site_perf); break;
        case 7:
            proc_O6R_FE3_OH(site_perf); break;
        case 8:
            //proc_O6R_FE_HACA_O2(site_perf); break;
        case 9:
            //proc_O6R_FE_HACA_OH(site_perf); break;
        case 10:
            proc_G5R_ZZ(site_perf); break;
        case 11:
            proc_D5R_R5(site_perf); break;
        case 12:
            proc_C6R_AC_FE3(site_perf, rng); break;
        case 13:
            proc_C5R_RFE(site_perf); break;
        case 14:
            proc_C5R_RAC(site_perf); break;
        case 15:
            proc_M5R_RZZ(site_perf); break;
        case 16:
            proc_C6R_BY5_FE3(site_perf, rng); break;
        case 17:
            proc_C6R_BY5_FE3violi(site_perf, rng); break;
        case 18:
            proc_L5R_BY5(site_perf); break;
        case 19:
            proc_M6R_BY5_FE3(site_perf, rng); break;
        case 20:
            proc_O6R_FE2(site_perf); break;
        case 21:
            proc_O6R_FE2(site_perf); break;
        case 22:
            proc_B6R_ACR5(site_perf); break;
        case 23:
            proc_M5R_eR5_FE3_ZZ(site_perf, rng); break;
        case 24:
            proc_G6R_RZZ(site_perf); break;
        case 25:
            proc_G6R_RFER(site_perf); break;
        case 26:
            proc_G6R_R5(site_perf); break;
        case 27:
            proc_L6_RBY5(site_perf); break;
        case 28:
            proc_L6_RACR(site_perf); break;
        case 29:
            proc_G5R_RFE(site_perf); break;
        case 30:
            proc_C6R_BY5_FE3(site_perf, rng); break;
        case 31:
            proc_C6R_BY5_FE3(site_perf, rng); break;
        case 32:
            proc_C6R_RAC_FE3(site_perf, rng); break;
		case 33:
			proc_G6R_ACBR(site_perf); break;
		case 34:
			proc_G6R_ACBL(site_perf); break;
        default:
            cout<<"ERROR: PAHProcess::performProcess: Process not found\n";
            return false;
    }/*
    if(m_pah->m_siteList.begin()->C1 != m_pah->m_cfirst)
        cout<<"WARNING: C1 of Site 1 does not correspond to m_cfirst..\n";
    if(m_pah->m_cfirst->A != 'H')
        cout<<"WARNING: A of m_cfirst is not H..\n";
    //cout<<"----PROCESS PERFORMED!-----\n";*/
	//Redetermine any sites that are non-reactive (due to hinderances)

	if (PAH_ID == 262 && m_pah->m_siteList.size() > 0){
		updateHinderedSites();
	}
	else{
		updateHinderedSites();
	}

    Spointer S1,S2,S3,S4;
    S1 = moveIt(site_perf, -1); S2 = moveIt(site_perf, 1);
    S3 = moveIt(site_perf, -2); S4 = moveIt(site_perf, 2);
    if(!checkCombinedSiteType(site_perf) || !checkCombinedSiteType(S1)
        || !checkCombinedSiteType(S2) || !checkCombinedSiteType(S3)
        || !checkCombinedSiteType(S4)) {
        std::cout<<"Structure produced invalid combined site type after performing process "
            << "ID"<<id<<" on PAH ID: "<<m_pah->m_parent->ID()<<"...\n"
            <<"*************\nAfter performing process --\n";
        printSites(site_perf);
        std::ostringstream msg;
        msg << "ERROR: Structure produced invalid combined site type after performing process "
            << "ID"<<id<<" on PAH ID: "<<m_pah->m_parent->ID()<<"..."
            << " (Sweep::KMC_ARS::PAHProcess::performProcess)";
        throw std::runtime_error(msg.str());
        assert(false);
        abort();
    }
	//NICK TO DO - Other way to make sure PAH is still valid
    //int calc_total = 2*m_pah->m_rings + (CarbonListSize()+m_pah->m_rings5_Lone+m_pah->m_rings5_Embedded)/2 + numberOfBridges() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded + 1;
    //if(calc_total != getCHCount().first) {
    //    //saveDOT("KMC_DEBUG/KMC_C_Counts_ERROR.dot");
    //    cout<<"Calculated total did not tally with double C counts!\n";
    //    cout<<"Last performed process: "<<jp.getName()<<", ID"<<jp.getID()<<'\n';
    //    std::ostringstream msg;
    //    msg << "\nCalculated total: "<<calc_total<<'\n'
    //        << "Real total: " << getCHCount().first << '\n';
    //    throw std::runtime_error(msg.str());
    //}
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
//NICK TO DO - FOR ALL PROCESS, MUST ADD CARBONS THAT ARE GAINED OR LOST BY EACH PROCESS!!!!!!!!!
void PAHProcess::proc_G6R_AC(Spointer& stt) {
    //printSitesMemb(stt);
    //printStruct();//++++
    //Cpointer newC1;
    //Cpointer newC2;

    /**
     * In order to model curved PAHs the code simply tracks the list of site types which makes up the edge of the PAH.
     * Therefore, the coordinates of the edge carbon atoms (coords) have been made redundant.
     * The code checks for hindrance by checking whether the proposed carbon atoms to be added would occupy a position which is already present.
     * This check which depends on the coordinates of the edge carbon atom can no longer be applied.
     *
     */
    //if(checkHindrance(stt)) {
    //    /*cout<<"Site hindered, process not performed.\n"*/ return;}

/*    if(!(C_1->C2->bridge) || !(C_2->C1->bridge)) {*/ // check if bulk C in AC site is a bridge
        // Add and remove C
        //if(C_1->C2->C3 != NULL) C_1->C2->C3->C3 = NULL;
        //if(C_2->C1->C3 != NULL) C_2->C1->C3->C3 = NULL;
        //removeC(C_1->C2, true);
        //removeC(C_2->C1, true);
        //newC1 = addC(C_1, normAngle(C_1->bondAngle1+120), 0);
        //newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1-60), normAngle(newC1->C1->bondAngle1-120));
    //}else {
        // update bridges info, both are no longer bridges
        //C_1->C2->C1 = C_1->C2->C3;// update neighbour
        //C_2->C1->C2 = C_2->C1->C3;
        //C_1->C2->bridge = false;
        //C_2->C1->bridge = false;
        //angletype a = C_2->C1->bondAngle1;
        //C_2->C1->bondAngle1 = C_2->C1->bondAngle2;
        //C_2->C1->bondAngle2 = a;
        //C_1->C2->C3 = NULL;
        //C_2->C1->C3 = NULL;
        // connect C_1 and C_2
    //    connectToC(C_1, C_2);
    //    // Add C
    //    newC1 = addC(C_1, normAngle(C_1->bondAngle1+120), 0);
    //    newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1-60), normAngle(newC1->C1->bondAngle1-120));
    //}
    //printStruct();
    // Add and remove H
    //updateA(C_1->C1, C_2->C2, 'H');

	//Check if this is growth on a bridge. If so, change other side of bridge to regular AC
	//if (stt->type == ACBR){
	//	int count = 1;
	//	Spointer Sptt = stt;
	//	while (count > 0){
	//		Sptt = moveIt(Sptt, 1);
	//		if (Sptt->type == ACBL) count--;
	//		if (Sptt->type == ACBR) count++;
	//	}
	//	convSiteType(Sptt, AC);
	//}

	//if (stt->type == ACBL){
	//	int count = 1;
	//	Spointer Sptt = stt;
	//	while (count > 0){
	//		Sptt = moveIt(Sptt, 1);
	//		if (Sptt->type == ACBR) count--;
	//		if (Sptt->type == ACBL) count++;
	//	}
	//	convSiteType(Sptt, AC);
	//}

    // neighbouring sites:
    Spointer S1 = moveIt(stt, -1); 
    Spointer S2 = moveIt(stt, 1);
    // Update Site and neighbours
	convSiteType(stt, FE);

	//if ((S1t == eRZZ && S3t == eR5 && S5t == eRFE) || (S2t == eRZZ && S4t == eR5 && S6t == eRFE)){
	//	//// Convert to a BY6 site
	//	bool b4 = true;
	//	Spointer Srem1, Srem2, Srem3;
	//	if (S1t == eRZZ){
	//		Srem1 = moveIt(stt, -1);
	//		Srem2 = moveIt(stt, -2);
	//		Srem3 = moveIt(stt, -3);
	//	}
	//	else{
	//		Srem1 = moveIt(stt, 1);
	//		Srem2 = moveIt(stt, 2);
	//		Srem3 = moveIt(stt, 3);
	//		b4 = false;
	//	}
	//	//// Remove sites and combine the neighbouring sites into BY6. 
	//	//// First remove all three from site map. Elementary site types first..
	//	delSiteFromMap(Srem1->type, Srem1);
	//	delSiteFromMap(Srem2->type, Srem2);
	//	delSiteFromMap(Srem3->type, Srem3);
	//	//// then for combined site types..
	//	delSiteFromMap(Srem1->comb, Srem1);
	//	delSiteFromMap(Srem2->comb, Srem2);
	//	delSiteFromMap(Srem3->comb, Srem3);
	//	//// remove the sites
	//	removeSite(Srem1);
	//	removeSite(Srem2);
	//	removeSite(Srem3);

	//	////add in BY6 site
	//	if (b4){
	//		addSite(BY6, stt);
	//		updateSites(S2, 1);
	//	}
	//	else{
	//		addSite(BY6, moveIt(stt, 1));
	//		updateSites(S1, 1);
	//	}

	//}
	//else if ((S1t == eRFE && S3t == eR5 && S5t == eRZZ) || (S2t == eRFE && S4t == eR5 && S6t == eRZZ)){
	//	//// Convert to a BY6 site
	//	bool b4 = true;
	//	Spointer Srem1, Srem2, Srem3;
	//	if (S1t == eRFE){
	//		Srem1 = moveIt(stt, -1);
	//		Srem2 = moveIt(stt, -2);
	//		Srem3 = moveIt(stt, -3);
	//	}
	//	else{
	//		Srem1 = moveIt(stt, 1);
	//		Srem2 = moveIt(stt, 2);
	//		Srem3 = moveIt(stt, 3);
	//		b4 = false;
	//	}

	//	//// Remove sites and combine the neighbouring sites into BY6. 
	//	//// First remove all three from site map. Elementary site types first..
	//	delSiteFromMap(Srem1->type, Srem1);
	//	delSiteFromMap(Srem2->type, Srem2);
	//	delSiteFromMap(Srem3->type, Srem3);
	//	//// then for combined site types..
	//	delSiteFromMap(Srem1->comb, Srem1);
	//	delSiteFromMap(Srem2->comb, Srem2);
	//	delSiteFromMap(Srem3->comb, Srem3);
	//	//// remove the sites
	//	removeSite(Srem1);
	//	removeSite(Srem2);
	//	removeSite(Srem3);

	//	////add in BY6 site
	//	if (b4){
	//		addSite(BY6, stt);
	//		updateSites(S2, 1);
	//	}
	//	else{
	//		addSite(BY6, moveIt(stt, 1));
	//		updateSites(S1, 1);
	//	}

	//}
	//else {
	//	if (S1->type == ACR5){
	//		convSiteType(S1, eRZZ);
	//		addSite(eRFE, S1);
	//		addSite(eR5, S1);
	//	}
	//	else {
	//		updateSites(S1, 1);
	//	}

	//	if (S2->type == ACR5){
	//		Spointer S3 = moveIt(S2, 1);
	//		convSiteType(S2, eRZZ);
	//		addSite(eR5, S3);
	//		addSite(eRFE, S3);

	//	}
	//	else {
	//		updateSites(S2, 1);
	//	}
	//}
	updateSites(S1, 1);
	updateSites(S2, 1);
    // Update combined site for Site and neighbours
    Spointer S3, S4;
	S1 = moveIt(stt, -1); S2 = moveIt(stt, 1);
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt);
    updateCombinedSites(S1); updateCombinedSites(S2);
    updateCombinedSites(S3); updateCombinedSites(S4);
    // add ring counts
    m_pah->m_rings++;
	addCount(2, 0);
    //printSites(stt);

	//S1t = (int)moveIt(stt, -1)->type;
	//S2t = (int)moveIt(stt, 1)->type;
	//S3t = (int)moveIt(stt, -2)->type;
	//S4t = (int)moveIt(stt, 2)->type;
	//S5t = (int)moveIt(stt, -3)->type;
	//S6t = (int)moveIt(stt, 3)->type;
	//St = (int)stt->type;
	//cout << "AC End" << endl << St << endl << S1t << endl << S2t << endl << S3t << endl << S4t << endl << S5t << endl << S6t << endl << "Done" << endl;
}
// 
// ************************************************************
// ID2- R6 growth on FE (AR2 on Matlab)
// ************************************************************
void PAHProcess::proc_G6R_FE(Spointer& stt) {
//    printSites(stt);
    //Cpointer newC1, newC2, newC3, newC4;

    /**
     * In order to model curved PAHs the code simply tracks the list of site types which makes up the edge of the PAH.
     * Therefore, the coordinates of the edge carbon atoms (coords) have been made redundant.
     * The code checks for hindrance by checking whether the proposed carbon atoms to be added would occupy a position which is already present.
     * This check which depends on the coordinates of the edge carbon atom can no longer be applied.
     *
     */
    //if(checkHindrance(stt)) {
    //    /*cout<<"Site hindered, process not performed.\n"*/ return;}

    // Add C
    //newC1 = addC(C_1, normAngle(C_1->bondAngle1+120),0);
    //newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1 - 60),0);
    //newC3 = addC(newC2, normAngle(newC2->C1->bondAngle1 - 60),0);
    //newC4 = addC(newC3, normAngle(newC3->C1->bondAngle1 - 60), normAngle(newC3->C1->bondAngle1 - 120));
    //C_1->C3 = C_2; C_2->C3 = C_1;
    // Add and remove H
    //updateA(C_1->C1, C_2->C2, 'H');
    // neighbouring sites:
    Spointer S1 = moveIt(stt, -1); 
    Spointer S2 = moveIt(stt, 1);
    updateSites(stt, 0);
	//if ((S1t == eRZZ && S3t == eR5 && S5t == eRFE) || (S2t == eRZZ && S4t == eR5 && S6t == eRFE)){
	//	//// Convert to a BY6 site
	//	bool b4 = true;
	//	Spointer Srem1, Srem2, Srem3;
	//	if (S1t == eRZZ){
	//		Srem1 = moveIt(stt, -1);
	//		Srem2 = moveIt(stt, -2);
	//		Srem3 = moveIt(stt, -3);
	//	}
	//	else{
	//		Srem1 = moveIt(stt, 1);
	//		Srem2 = moveIt(stt, 2);
	//		Srem3 = moveIt(stt, 3);
	//		b4 = false;
	//	}
	//	//// Remove sites and combine the neighbouring sites into BY6. 
	//	//// First remove all three from site map. Elementary site types first..
	//	delSiteFromMap(Srem1->type, Srem1);
	//	delSiteFromMap(Srem2->type, Srem2);
	//	delSiteFromMap(Srem3->type, Srem3);
	//	//// then for combined site types..
	//	delSiteFromMap(Srem1->comb, Srem1);
	//	delSiteFromMap(Srem2->comb, Srem2);
	//	delSiteFromMap(Srem3->comb, Srem3);
	//	//// remove the sites
	//	removeSite(Srem1);
	//	removeSite(Srem2);
	//	removeSite(Srem3);

	//	////add in BY6 site
	//	if (b4){
	//		addSite(BY6, stt);
	//		updateSites(S2, 1);
	//	}
	//	else{
	//		addSite(BY6, moveIt(stt,1));
	//		updateSites(S1, 1);
	//	}

	//}
	//else if((S1t == eRFE && S3t == eR5 && S5t == eRZZ) || (S2t == eRFE && S4t == eR5 && S6t == eRZZ)){
	//	//// Convert to a BY6 site
	//	bool b4 = true;
	//	Spointer Srem1, Srem2, Srem3;
	//	if (S1t == eRFE){
	//		Srem1 = moveIt(stt, -1);
	//		Srem2 = moveIt(stt, -2);
	//		Srem3 = moveIt(stt, -3);
	//	}
	//	else{
	//		Srem1 = moveIt(stt, 1);
	//		Srem2 = moveIt(stt, 2);
	//		Srem3 = moveIt(stt, 3);
	//		b4 = false;
	//	}

	//	//// Remove sites and combine the neighbouring sites into BY6. 
	//	//// First remove all three from site map. Elementary site types first..
	//	delSiteFromMap(Srem1->type, Srem1);
	//	delSiteFromMap(Srem2->type, Srem2);
	//	delSiteFromMap(Srem3->type, Srem3);
	//	//// then for combined site types..
	//	delSiteFromMap(Srem1->comb, Srem1);
	//	delSiteFromMap(Srem2->comb, Srem2);
	//	delSiteFromMap(Srem3->comb, Srem3);
	//	//// remove the sites
	//	removeSite(Srem1);
	//	removeSite(Srem2);
	//	removeSite(Srem3);

	//	////add in BY6 site
	//	if (b4){
	//		addSite(BY6, stt);
	//		updateSites(S2, 1);
	//	}
	//	else{
	//		addSite(BY6, moveIt(stt, 1));
	//		updateSites(S1, 1);
	//	}
	//
	//}
	//else {
	//	if (S1->type == ACR5){
	//		convSiteType(S1, eRZZ);
	//		addSite(eRFE, S1);
	//		addSite(eR5, S1);
	//	}
	//	else {
	//		updateSites(S1, 1);
	//	}

	//	if (S2->type == ACR5){
	//		Spointer S3 = moveIt(S2, 1);
	//		convSiteType(S2, eRZZ);
	//		addSite(eR5, S3);
	//		addSite(eRFE, S3);

	//	}
	//	else {
	//		updateSites(S2, 1);
	//	}
	//}
	updateSites(S1, 1);
	updateSites(S2, 1);
    // Add new Sites
    Spointer newS1 = addSite(FE, stt);
	Spointer newS2 = addSite(FE, moveIt(stt, 1));
    // Update combined sites for all new sites and original neighbours
    Spointer S3, S4, S5, S6;
	S1 = moveIt(stt, -1); S2 = moveIt(stt, 1);
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
	S5 = moveIt(S3, -1); S6 = moveIt(S4, 1);
    updateCombinedSites(stt);
    updateCombinedSites(newS1); updateCombinedSites(newS2); // new sites
    updateCombinedSites(S1); updateCombinedSites(S2); // original neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours of neighbours
	updateCombinedSites(S5); updateCombinedSites(S6); // neighbours of neighbours of neighbours
    // Add H count
    addCount(4, 2);
    // add ring counts
    m_pah->m_rings++;

	//S1t = (int)moveIt(stt, -1)->type;
	//S2t = (int)moveIt(stt, 1)->type;
	//S3t = (int)moveIt(stt, -2)->type;
	//S4t = (int)moveIt(stt, 2)->type;
	//S5t = (int)moveIt(stt, -3)->type;
	//S6t = (int)moveIt(stt, 3)->typ;e
	//St = (int)stt->type;
	//cout << "FE End" << endl << St << endl << S1t << endl << S2t << endl << S3t << endl << S4t << endl << S5t << endl << S6t << endl << "Done" << endl;
}
// 
// ************************************************************
// ID3- BY6 closure reaction (AR14 on Matlab)
// ************************************************************
void PAHProcess::proc_L6_BY6(Spointer& stt) {
    //printSites(stt);
    // Remove C

    //Cpointer now = C_1->C2;
    //do{
    //    Cpointer next;
    //    if(!now->bridge) {
    //        next = now->C2;
    //        //if(now->C3 != NULL) now->C3->C3 = NULL;
    //        removeC(now, true);
    //        now = next;
    //    }
    //    else {
    //        next = now->C3->C2; //
    //        // bridged bulk C will not be removed from edge
    //        Cpointer b = now->C3; // C bridged to now
    //        b->C3 = NULL; now->C3 = NULL;
    //        b->bridge = false; now->bridge = false;
    //        b->bondAngle1 = b->bondAngle2;
    //        now->bondAngle2 = 0; b->bondAngle2 = 0;
    //        // connect C_1 to next and the two bulk atoms still remaining
    //        connectToC(C_1, next);
    //        connectToC(b, now);
    //        now = next;
    //    }
    //}while(now!=C_2);
    //// Change bond angle between C_1 and C_2
    //C_1->bondAngle1 = normAngle(C_1->bondAngle1+120);
    //// Add and remove H
    //updateA(C_1->C1, C_2->C2, 'H');

    //// Remove BY6 site and combine the neighbouring sites. 
    //// First remove all three from site map. Elementary site types first..
    delSiteFromMap(moveIt(stt, -1)->type, moveIt(stt, -1));
    delSiteFromMap(moveIt(stt, 1)->type, moveIt(stt, 1));
    delSiteFromMap(stt->type, stt);
    //// then for combined site types..
    delSiteFromMap(moveIt(stt, -1)->comb, moveIt(stt, -1));
    delSiteFromMap(moveIt(stt, 1)->comb, moveIt(stt, 1));
    delSiteFromMap(stt->comb, stt);
    // Convert the BY6 site into the resulting site after reaction,
    // finding resulting site type:
    int ntype1 = abs( (int) moveIt(stt, -1)->type);
    int ntype2 = abs((int) moveIt(stt, 1)->type);
//    if(ntype1 < 6 && ntype2 < 6) {
        int newType = (ntype1+ntype2+1);
        // convert site
   //     if(newType>5) {
   //         //saveDOT(std::string("BY6ClosureProblem.dot"));
   //         std::cerr<<"ERROR: newType is > 5 (PAHProcess::proc_L6_BY6)\n";
			//cout << ntype1 << endl;
			//cout << ntype2 << endl;
   //     }
		convSiteType(stt, (kmcSiteType)newType);
//    }
/*    else {
        int newType = 0; 
        if((ntype1-6)>=0) {
            ntype1 -= 6;
            newType = 6;
        }
        if((ntype2-6)>=0) {
            ntype2 -= 6;
            if(newType>0) newType = 10;
            else newType = 6;
        }
        newType += ntype1+ntype2+2;
		if (newType == 20) newType = 4;
        convSiteType(stt, (kmcSiteType) newType);
    }  */  
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
}
// 
// ************************************************************
// ID4- phenyl addition (AR15 in Matlab)
// ************************************************************
//NICK TO DO - Figure this one out
//void PAHProcess::proc_PH_benz(Spointer& stt, rng_type &rng) {
//    Cpointer chosen;
//    bool before; // true if C1 of site is chosen, false if C2
//    // choose one of the C atoms if site type is FE/AC/ZZ
//    if(stt->type == FE || stt->type == AC || stt->type == ZZ) {
//        // Define a distribution that has two equally probably outcomes
//        boost::bernoulli_distribution<> choiceDistrib;
//        // Now build an object that will generate a sample using rng
//        boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);
//
//        if(choiceGenerator()) {
//            chosen = C_1;
//            before=true;
//        }
//        else {
//            chosen = C_2;
//            before=false;
//        }
//    }
//    else { //if not, choose the C which is not part of a 5-membered ring
//        if(moveIt(stt,-1)->type == R5) {
//            chosen = C_2;
//            before=false;
//        }
//        else {
//            chosen = C_1;
//            before=true;
//        }
//    }
//    // neighbouring site to be updated:
//    Spointer neighbour;
//    if(before) {
//        neighbour = moveIt(stt,-1);
//    }
//    else {
//        neighbour = moveIt(stt,1);
//    }
//    // check hindrance
//    if(checkHindrancePhenyl(chosen)) return;
//    // add C atoms
//    Cpointer newC;
//    newC = bridgeC(chosen);
//    newC = addC(newC, normAngle(chosen->bondAngle2+60), 0);
//    newC = addC(newC, normAngle(newC->C1->bondAngle1-60),0);
//    newC = addC(newC, normAngle(newC->C1->bondAngle1-60),0);
//    newC = addC(newC, normAngle(newC->C1->bondAngle1-60),0);
//    newC = addC(newC, normAngle(newC->C1->bondAngle1-60),normAngle(newC->C1->bondAngle1-120));
//    connectToC(newC, chosen->C3);
//    // which site do we add the new sites before?
//    Spointer addBefore;
//    if(before) addBefore = stt;
//    else addBefore = neighbour;
//    // Add new sites (4 new FE sites)
//    Cpointer Cnow = chosen->C3->C2;
//    Spointer nS1 = addSite(FE, Cnow, Cnow->C2, addBefore); Cnow = Cnow->C2;
//    Spointer nS2 = addSite(FE, Cnow, Cnow->C2, addBefore); Cnow = Cnow->C2;
//    Spointer nS3 = addSite(FE, Cnow, Cnow->C2, addBefore); Cnow = Cnow->C2;
//    Spointer nS4 = addSite(FE, Cnow, Cnow->C2, addBefore); Cnow = Cnow->C2;
//    // add and remove H
//    updateA(chosen->C1, chosen->C2, 'H');
//    // update sites and neighbour
//    // new member C for stt and neighbour
//    Cpointer s_C1, s_C2, n_C1, n_C2;
//    Spointer n1, n2; // new neighbours for updated sites
//    if(before) {
//        s_C1 = moveIt(stt, -1)->C2;
//        s_C2 = stt->C2;
//        n_C1 = neighbour->C1;
//        n_C2 = moveIt(neighbour, 1)->C1;
//        n1 = moveIt(stt, 1);
//        n2 = moveIt(neighbour, -1);
//    }else {
//        s_C1 = stt->C1;
//        s_C2 = moveIt(stt, 1)->C1;
//        n_C1 = moveIt(neighbour, -1)->C2;
//        n_C2 = neighbour->C2;
//        n1 = moveIt(stt, -1);
//        n2 = moveIt(neighbour, 1);
//    }
//    updateSites(stt, s_C1, s_C2, 2);
//    updateSites(neighbour, n_C1, n_C2, 2);
//    // update combined sites for all new sites and neighbours (and their neighbours)
//    updateCombinedSites(stt); updateCombinedSites(neighbour);
//    updateCombinedSites(nS1); updateCombinedSites(nS2); updateCombinedSites(nS3); updateCombinedSites(nS4);
//    updateCombinedSites(n1); updateCombinedSites(n2);
//    // update H count
//    addCount(0,4);
//    // add ring counts
//    m_pah->m_rings++;
//}
// 
// ************************************************************
// ID5- R6 desorption at FE (AR8 in Matlab)
// ************************************************************
void PAHProcess::proc_D6R_FE3(Spointer& stt) {
    //printSites(stt);
    // cannot happen if R6 is next to bridge
    //if(stt->C1->C1->C1->bridge || stt->C2->C2->C2->bridge) return;
	//NICK TO DO - Check for bridges
    //// member C atoms of resulting FE site
    //Cpointer C1_new, C2_new;
    //C1_new = C_1->C1->C1;
    //C2_new = C_2->C2->C2;
    // delete neighbouring FE sites from site map. Identify them first:
    Spointer S1 = moveIt(stt,-1);
    Spointer S2 = moveIt(stt,1);
    // then remove them
    removeSite(S1);
    removeSite(S2);
    // update site stt and new neighbours
    // new neighbours:
    S1 = moveIt(stt,-1); S2 = moveIt(stt,1);
    updateSites(stt, 0); // only change C1 and C2, site still FE
    updateSites(S1, -1);
    updateSites(S2, -1);
    // update combined sites
    Spointer S3, S4;
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt); // update resulting site
    updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
    // remove C
    //for(int i=0; i<4; i++) {
    //    removeC(C1_new->C2, false);
    //}
    //C1_new->bondAngle1 = normAngle(C1_new->bondAngle1-120);
    //// add H atoms
    //updateA(C1_new->C1, C2_new->C2, 'H');
    // update H count
    addCount(-4,-2);
    // add ring counts
    m_pah->m_rings--;
}
// ************************************************************
// ID6- R6 oxidation at FE by O2 (AR10 in Matlab)
// ************************************************************
void PAHProcess::proc_O6R_FE3_O2(Spointer& stt) {
    proc_D6R_FE3(stt);
}
// ************************************************************
// ID7- R6 oxidation at FE by OH (AR11 in Matlab)
// ************************************************************
void PAHProcess::proc_O6R_FE3_OH(Spointer& stt) {
    proc_D6R_FE3(stt);
}
// ************************************************************
// ID8- R6 oxidation at AC by O2 (AR12 in Matlab)
// ************************************************************
//NICK TO DO - Figure this one out
//void PAHProcess::proc_O6R_FE_HACA_O2(Spointer& stt, Cpointer C_1, Cpointer C_2) {
//    //printSites(stt);
//    // member C atoms of resulting AC site
//    Cpointer C1_res, C2_res;
//    C1_res = C_1->C1;
//    C2_res = C_2->C2;
//    // check if process will result in a bridge
//    bool bridge = false;
//    cpair pos = jumpToPos(C1_res->coords, normAngle(C1_res->bondAngle1-120));
//    if(m_pah->m_cpositions.count(pos)) return;//bridge = true;
//    // check if site is next to a bridge
//    bridge = (C1_res->bridge || C2_res->bridge);
//    if(bridge) return;
//    // remove C
//    removeC(C_1, false);
//    removeC(C_2, false);
//    //if(!bridge) {
//        addC(C1_res, normAngle(C1_res->bondAngle1-120), 0, true);
//        addC(C1_res->C2, normAngle(C1_res->bondAngle1+60), normAngle(C1_res->bondAngle1+120), true);
//    //}
//    // update H
//    updateA(C1_res, C2_res, 'H');
//    // update sites and neighbours
//    Spointer S1, S2, S3, S4;
//    S1 = moveIt(stt,-1); S2 = moveIt(stt,1);
//    updateSites(stt, C1_res, C2_res, 2);
//    updateSites(S1, S1->C1, C1_res, -1);
//    updateSites(S2, C2_res, S2->C2, -1);
//    // update combined sites
//    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
//    updateCombinedSites(stt); // update resulting site
//    updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
//    updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
//    // add ring counts
//    m_pah->m_rings--;
//}
// ************************************************************
// ID9- R6 oxidation at AC by OH (AR13 in Matlab)
// ************************************************************
//void PAHProcess::proc_O6R_FE_HACA_OH(Spointer& stt, Cpointer C_1, Cpointer C_2) {
//    proc_O6R_FE_HACA_O2(stt, C_1, C_2);
//}
// ************************************************************
// ID10- R5 growth on ZZ (AR3 in Matlab)
// ************************************************************
void PAHProcess::proc_G5R_ZZ(Spointer& stt) {
    //printSites(stt);
    
    /**
     * In order to model curved PAHs the code simply tracks the list of site types which makes up the edge of the PAH.
     * Therefore, the coordinates of the edge carbon atoms (coords) have been made redundant.
     * The code checks for hindrance by checking whether the proposed carbon atoms to be added would occupy a position which is already present.
     * This check which depends on the coordinates of the edge carbon atom can no longer be applied.
     *
     */
    //if(checkHindrance(stt)) {
    //   /* cout<<"Site hindered, process not performed.\n";*/ return;}
    
    // member C atoms of resulting R5 site
    //removeC(C_2->C1, true);
    //Cpointer C1_res, C2_res;
    //C1_res = addC(C_1, normAngle(C_1->bondAngle1+120), 0);
    //C2_res = addC(C1_res, normAngle(C_1->bondAngle1-90), normAngle(C_1->bondAngle1-180));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    //updateA(C_1->C1, C_2->C2, 'H');
    // update sites and neighbours
    convSiteType(stt, R5);
    Spointer S1, S2, S3, S4;
    // neighbours
    S1 = moveIt(stt,-1); 
    S2 = moveIt(stt,1);
    addR5toSite(S1);
    addR5toSite(S2);
    // update combined sites
    S3 = moveIt(S1, -1); 
    S4 = moveIt(S2, 1);
    updateCombinedSites(stt); // update resulting site
    updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
    // add ring counts
    m_pah->m_rings5_Lone++;
	addCount(2, 0);
}
// ************************************************************
// ID11- R5 desorption (AR7 in Matlab)
// ************************************************************
void PAHProcess::proc_D5R_R5(Spointer& stt) {
    //printSites(stt);
    // member C atoms of resulting R5 site
    //Cpointer C1_res, C2_res;
    //C1_res = C_1->C1;
    //C2_res = C_2->C2;
    //removeC(C_1, false);
    //removeC(C_2, false);
    //addC(C1_res, normAngle(C1_res->C1->bondAngle1-60), normAngle(C1_res->C1->bondAngle1), true);
    //updateA(C1_res->C1, C2_res->C2, 'H');
    // update sites and neighbours
    Spointer S1, S2, S3, S4;
    // neighbours
    S1 = moveIt(stt,-1); S2 = moveIt(stt,1);
    convSiteType(stt, ZZ);
    remR5fromSite(S1);
    remR5fromSite(S2);
    // update combined sites
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt); // update resulting site
    updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
    // add ring counts
    m_pah->m_rings5_Lone--;
	addCount(-2, 0);
}
// ************************************************************
// ID12- R6 conversion to R5 (AR9 in Matlab)
// ************************************************************
void PAHProcess::proc_C6R_AC_FE3(Spointer& stt, rng_type &rng) {
    //printSites(stt);
    
    //if(checkHindrance(stt)) {
    //    /*cout<<"Site hinderzed, process not performed.\n"*/ return;}
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

    //Cpointer C1_res, C2_res, C1_R5, C2_R5;
    Spointer FE_res;
    // removing R6 first. save the two member C for resulting FE
    if(b4) {
        FE_res = moveIt(stt, -2); // resulting FE (the FE3)
        //C1_res = FE_res->C1->C1->C1;
        //C2_res = C_1->C2;
    }else {
        FE_res = moveIt(stt, 2);
        //C1_res = C_2->C1;
        //C2_res = FE_res->C2->C2->C2;
    }
    // remove C atoms
    //for(int i=0; i!=4; i++) removeC(C1_res->C2, false);
    //C1_res->bondAngle1 = normAngle(C1_res->bondAngle1-120);
    //// now add R5 on resulting ZZ site (used to be AC)
    //Cpointer Cstart;
    //if(b4) Cstart = C2_res;
    //else Cstart = C_1;
    //removeC(Cstart->C2, true);
    //C1_R5 = addC(Cstart, normAngle(Cstart->bondAngle1+120), 0, false);
    //C2_R5 = addC(C1_R5, normAngle(Cstart->bondAngle1-90), normAngle(Cstart->bondAngle1-180), false);
    //// update H atoms
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
        convSiteType(stt, R5); // convert the AC site to R5
        addR5toSite(FE_res); // convert the FE site to RFE
        addR5toSite(S2); // update neighbour of resulting R5
        updateSites(S1, -1); // update neighbour of resulting RFE
    }else {
        S1 = moveIt(stt, -1);
        S2 = moveIt(FE_res, 1);
        convSiteType(stt, R5); // convert the AC site to R5
        addR5toSite(FE_res); // convert the FE site to RFE
        addR5toSite(S1); // update neighbour of resulting R5
        updateSites(S2, -1); // update neighbour of resulting RFE
    }
    // update combined sites for all the sites involved and their neighbours
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt); updateCombinedSites(FE_res);// update resulting sites
    updateCombinedSites(S1); updateCombinedSites(S2); // update neighbours
    updateCombinedSites(S3); updateCombinedSites(S4); // update neighbours of neighbours
    // update H count
    addCount(-2, -2);
    // add ring counts
    m_pah->m_rings--;
    m_pah->m_rings5_Lone++;
}
// ************************************************************
// ID13- R5 conversion to R6 on FE (AR5 in Matlab)
// ************************************************************
void PAHProcess::proc_C5R_RFE(Spointer& stt) {
    //printSites(stt);

    /**
     * In order to model curved PAHs the code simply tracks the list of site types which makes up the edge of the PAH.
     * Therefore, the coordinates of the edge carbon atoms (coords) have been made redundant.
     * The code checks for hindrance by checking whether the proposed carbon atoms to be added would occupy a position which is already present.
     * This check which depends on the coordinates of the edge carbon atom can no longer be applied.
     *
     */
    //if(checkHindrance(stt)) {
    //    /*cout<<"Site hindered, process not performed.\n"*/ return;}
    
    // check if R5 is before or after the site
    bool b4=false;
    if(moveIt(stt,-1)->type == R5) b4 = true;
    // R5 will be removed and replaced by ZZ. first identify the C member of the
    // resulting AC site, C_AC which is from bulk. Identify where the R5 site is too.
    Spointer sR5;
    //Cpointer C_AC;
    if(b4) {
        sR5 = moveIt(stt, -1);
        //C_AC = sR5->C1->C1;
    }else {
        sR5 = moveIt(stt, 1);
        //C_AC = sR5->C2->C2;
    }
    // remove R5 first, leaving a ZZ site (adding another C atom after removing C)
    //Cpointer Cstart;
    //if(b4) Cstart = C_AC; 
    //else Cstart = C_2->C1;
    //for(int i=0; i!=2; i++) removeC(Cstart->C2, false);
    //addC(Cstart, normAngle(Cstart->bondAngle1-120), normAngle(Cstart->bondAngle1-60), true);
    //// this new C atom is irrelevant. Next add a R6 on the resulting FE (from RFE)
    //Cpointer C1_new, C2_new, C3_new, C4_new; // save all new C atoms
    //if(b4) 
    //    Cstart = C_2->C1;
    //else 
    //    Cstart = C_1;
    //C1_new = addC(Cstart, normAngle(Cstart->bondAngle1+120), 0);
    //C2_new = addC(C1_new, normAngle(C1_new->C1->bondAngle1-60), 0);
    //C3_new = addC(C2_new, normAngle(C2_new->C1->bondAngle1-60), 0);
    //C4_new = addC(C3_new, normAngle(C3_new->C1->bondAngle1-60), normAngle(C3_new->C1->bondAngle1-120));
    // edit sites. first identify the neighbouring sites of resulting AC & FE3
    Spointer S1, S2, S3, S4;
    if(b4) {
        S1 = moveIt(sR5, -1); // neighbour of R5
        S2 = moveIt(stt, 1); // neighbour of RFE (stt)
        convSiteType(sR5, AC); // convert R5 to AC
        remR5fromSite(S1); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        remR5fromSite(stt); // convert RFE to FE
        addSite(FE, stt); // add new FE site before the original RFE
        addSite(FE, S2); // add new FE site after the original RFE
        updateSites(S2, +1); // update resulting FE3 neighbour
    }else {
        S1 = moveIt(stt, -1); // neighbour of RFE (stt)
        S2 = moveIt(sR5, 1); // neighbour of R5
        convSiteType(sR5, AC); // convert R5 to AC
        remR5fromSite(S2); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        remR5fromSite(stt); // convert RFE to FE
        addSite(FE, stt); // add new FE site before the original RFE
        addSite(FE, sR5); // add new FE site after the original RFE
        updateSites(S1, +1); // update resulting FE3 neighbour
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
    addCount(2, 2);
    // add ring counts
    m_pah->m_rings++;
    m_pah->m_rings5_Lone--;
}
// ************************************************************
// ID14- R5 conversion to R6 on AC (AR4 in Matlab)
// ************************************************************
void PAHProcess::proc_C5R_RAC(Spointer& stt) {
    //printSites(stt);

    // check if R5 is before or after RAC
    bool b4 = false;
    if(moveIt(stt,-1)->type == R5) b4 = true;
    // R5 will be removed and replaced by AC. first identify the C member of the
    // resulting AC site, C_AC. Identify where the R5 site is too.
    Spointer sR5;
    //Cpointer C_AC;
    if(b4) {
        sR5 = moveIt(stt, -1);
        //C_AC = sR5->C1->C1;
    }else {
        sR5 = moveIt(stt, 1);
        //C_AC = sR5->C2->C2;
    }
    // check if there's a bridge in the BY5
	//NICK TO DO - CHECK if there is a bridge
    bool bridge = false;
    //if(b4) bridge = C_2->C1->bridge;
    //else bridge = C_1->C2->bridge;
    // remove R5 first, leaving a ZZ site (adding another C atom after removing C)
    //Cpointer Cstart;
    //if(b4) Cstart = C_AC; 
    //else Cstart = C_2->C1;
    //for(int i=0; i!=2; i++) removeC(Cstart->C2, false);
    //addC(Cstart, normAngle(Cstart->bondAngle1-120), normAngle(Cstart->bondAngle1-60), true);
    //// this new C atom is irrelevant. Next add a R6 on the resulting AC (from RAC)
    //Cpointer C1_new, C2_new; // save all new C atoms
    //if(b4) Cstart = C_AC->C2->C2;
    //else Cstart = C_1;
    //if(!bridge){
    //    for(int i=0; i!=2; i++) removeC(Cstart->C2, true);
    //} else {//else if there's a bridge, convert the bridge atoms to normal edge atoms
    //    Cpointer b1 = Cstart->C2;
    //    Cpointer b2 = b1->C3;
    //    connectToC(Cstart, b2->C2);
    //    b1->C3 = NULL;
    //    b2->C3 = NULL;
    //    b1->bridge = false;
    //    b2->bridge = false;
    //    b1->C1 = b2;
    //    b2->C2 = b1;
    //    angletype a = b2->bondAngle1;
    //    b2->bondAngle1 = b2->bondAngle2;
    //    b2->bondAngle2 = a;
    //}
    //C1_new = addC(Cstart, normAngle(Cstart->bondAngle1+120), 0);
    //C2_new = addC(C1_new, normAngle(C1_new->C1->bondAngle1-60), normAngle(C1_new->C1->bondAngle1-120));
    // edit sites. first identify the neighbouring sites of resulting AC & FE3
    Spointer S1, S2, S3, S4;
    if(b4) {
        S1 = moveIt(sR5, -1); // neighbour of R5
        S2 = moveIt(stt, 1); // neighbour of RAC (stt)
        convSiteType(sR5, AC); // convert R5 to AC
        remR5fromSite(S1); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, FE); // convert RAC to FE
        updateSites(S2, +1); // update resulting FE neighbour
    } else {
        S1 = moveIt(stt, -1); // neighbour of RAC (stt)
        S2 = moveIt(sR5, 1); // neighbour of R5
        convSiteType(sR5, AC); // convert R5 to AC
        remR5fromSite(S2); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, FE); // convert RAC to FE
        updateSites(S1, +1); // update resulting FE neighbour
    }
    // update H atoms
    //if(b4) updateA(C_AC->C1, C_2->C2, 'H');
    //else updateA(C_1->C1, C_AC->C2, 'H');
    // update combined sites for all sites involved and their neighbours
    // (excluding new FE sites, since their combined site type will still be None)
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
void PAHProcess::proc_M5R_RZZ(Spointer& stt) {
    //printSites(stt);
    
    /**
     * In order to model curved PAHs the code simply tracks the list of site types which makes up the edge of the PAH.
     * Therefore, the coordinates of the edge carbon atoms (coords) have been made redundant.
     * The code checks for hindrance by checking whether the proposed carbon atoms to be added would occupy a position which is already present.
     * This check which depends on the coordinates of the edge carbon atom can no longer be applied.
     *
     */
    //if(checkHindrance(stt)) {
    //    /*cout<<"Site hindered, process not performed.\n"*/ return;}
    
    // check if R5 is before or after RZZ
    bool b4 = false;
    if(moveIt(stt,-1)->type == R5) b4 = true;
    // R5 will be removed and replaced by ZZ. first identify the C member of the
    // resulting RZZ site, C_RZZ. Identify where the R5 site is too.
    Spointer sR5;
    //Cpointer C_RZZ;
    if(b4) {
        sR5 = moveIt(stt, -1);
        //C_RZZ = sR5->C1->C1;
    }else {
        sR5 = moveIt(stt, 1);
        //C_RZZ = sR5->C2->C2;
    }
    // remove R5 first, leaving a ZZ site (adding another C atom after removing C)
    //Cpointer Cstart;
    //if(b4) Cstart = C_RZZ; 
    //else Cstart = C_2->C1;
    //for(int i=0; i!=2; i++) removeC(Cstart->C2, false);
    //addC(Cstart, normAngle(Cstart->bondAngle1-120), normAngle(Cstart->bondAngle1-60),true);
    //// Next add a R5 on the neighbouring ZZ
    //Cpointer C1_R5, C2_R5; // save all new C atoms
    //if(b4) Cstart = C_RZZ->C2->C2;
    //else Cstart = C_1;
    //removeC(Cstart->C2, true);
    //C1_R5 = addC(Cstart, normAngle(Cstart->bondAngle1+120), 0);
    //C2_R5 = addC(C1_R5, normAngle(C1_R5->C1->bondAngle1-90), normAngle(C1_R5->C1->bondAngle1-180));
    // edit sites. first identify the neighbouring sites of resulting RZZ & R5
    Spointer S1, S2, S3, S4;
    if(b4) {
        S1 = moveIt(sR5, -1); // neighbour of R5
        S2 = moveIt(stt, 1); // neighbour of RAC (stt)
        convSiteType(sR5, RZZ); // convert R5 to RZZ
        remR5fromSite(S1); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, R5); // convert RZZ to R5
        addR5toSite(S2); // update S2
    } else {
        S1 = moveIt(stt, -1); // neighbour of RAC (stt)
        S2 = moveIt(sR5, 1); // neighbour of R5
        convSiteType(sR5, RZZ); // convert R5 to RZZ
        addR5toSite(S1); // convert S1 from RS1 to S1 (e.g RFE -> FE)
        convSiteType(stt, R5); // convert RZZ to R5
        remR5fromSite(S2); // update S2
    }
    // update H atoms
    //if(b4) updateA(C_RZZ->C1, C_2->C2, 'H');
    //else updateA(C_1->C1, C_RZZ->C2, 'H');
    // update combined sites for all sites involved and their neighbours
    // (excluding new FE sites, since their combined site type will still be None)
    S3 = moveIt(S1, -1);
    S4 = moveIt(S2, 1);
    updateCombinedSites(stt); updateCombinedSites(sR5); // new FE and AC
    updateCombinedSites(S1); updateCombinedSites(S2); 
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours
    //printSites(stt);
    //cout<<sp.None;
}

// ************************************************************
// ID16- R6 migration & conversion to R5 at BY5 (pyrene+R5; pathway 1; AR22 in Matlab)
// ************************************************************
void PAHProcess::proc_C6R_BY5_FE3(Spointer& stt, rng_type &rng) {
    //printSites(stt);

    // check if there are any bridges in the BY5, cancel process is yes
	//NICK TO DO - Check for bridges
    //Cpointer now=C_1->C2;
    //for(int i=0; i!=3; i++) {
    //    if(now->bridge) return;
    //    else now = now->C2;
    //}
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
    // first the R6 will be changed to R5 by removing a C. locate where the FE3 is and which
    // C is to be removed.
    Spointer sFE3;
    //Cpointer CRem, CFE;// C to be removed and remaining C from original FE3 site
    if(b4) {
        sFE3 = moveIt(stt, -2); 
        //CFE = sFE3->C2;
        //CRem = sFE3->C1;
    }else {
        sFE3 = moveIt(stt, 2);
        //CFE = sFE3->C1;
        //CRem = sFE3->C2;
    }
    //CRem->C1->bondAngle1 = normAngle(CRem->C1->bondAngle1-30);
    //removeC(CRem, false);
    // now close BY5 to form R6. First remove all bulk C
    //for(int i=0; i!=3; i++) removeC(C_1->C2, true);
    //// add a C atom
    //addC(C_1, normAngle(C_1->bondAngle1+120), normAngle(C_1->bondAngle1+60));
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
        convSiteType(sFE3, R5); // convert FE3 to R5
        convSiteType(stt, RFE); // convert BY5 to RFE
        updateSites(S2, +1); // add a bulk C to S2
        updateSites(S1, -1);
        addR5toSite(S1); // remove a bulk C from S1 and add R5
    } else {
        S2 = moveIt(sFE3, 1); // neighbour of FE3
        S1 = moveIt(stt, -1); // neighbour of BY5 (stt)
        convSiteType(sFE3, R5); // convert FE3 to R5
        convSiteType(stt, RFE); // convert BY5 to RFE
        updateSites(S1, +1); // add a bulk C to S1
        updateSites(S2, -1);
        addR5toSite(S2); // remove a bulk C from S2 and add R5
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
}
// ************************************************************
// ID17- R6 migration & conversion to R5 at BY5 (pyrene+R5; pathway 2-violi; AR24 in Matlab)
// ************************************************************
void PAHProcess::proc_C6R_BY5_FE3violi(Spointer& stt, rng_type &rng) {
    proc_C6R_BY5_FE3(stt, rng);
}
// ************************************************************
// ID18- BY5 closure (AR16 in Matlab)
// ************************************************************
void PAHProcess::proc_L5R_BY5(Spointer& stt) {
    //printSites(stt);

    // Remove C
    //Cpointer now = C_1->C2;
    //do{
    //    Cpointer next;
    //    if(!now->bridge) {
    //        next = now->C2;
    //        //if(now->C3 != NULL) now->C3->C3 = NULL;
    //        removeC(now, true);
    //        now = next;
    //    }
    //    else {
    //        next = now->C3->C2; //
    //        // bridged bulk C will not be removed from edge
    //        Cpointer b = now->C3; // C bridged to now
    //        b->C3 = NULL; now->C3 = NULL;
    //        b->bridge = false; now->bridge = false;
    //        b->bondAngle1 = b->bondAngle2;
    //        now->bondAngle2 = 0; b->bondAngle2 = 0;
    //        // connect C_1 to next and the two bulk atoms still remaining
    //        connectToC(C_1, next);
    //        connectToC(b, now);
    //        now = next;
    //    }
    //}while(now!=C_2);
    //// Change bond angle between C_1 and C_2
    //C_1->bondAngle1 = normAngle(C_1->bondAngle1+90);
    //// Add and remove H
    //updateA(C_1->C1, C_2->C2, 'H');

    // Convert the BY5 site into the resulting site after reaction,
	Spointer S1 = moveIt(stt, -1);
	Spointer S2 = moveIt(stt, 1);
	int S1type1 = abs((int)S1->type);
	int S1type2 = abs((int)S2->type);
	int newType1, newType2;
	
	convSiteType(stt, ER5);

	//Convert S1 site to site with ER5 beside it
	//Check if site is already beside an ER5
	if (S1type1 > 10 && S1type1 % 30 < 10){
		newType1 = S1type1 + 20;
	}
	else{
		newType1 = S1type1 + 30;
	}
	convSiteType(S1, (kmcSiteType) newType1);

	//Convert S1 site to site with ER5 beside it
	//Check if site is already beside an ER5
	if (S1type2 > 10 && S1type2 % 30 < 10){
		newType2 = S1type2 + 20;
	}
	else{
		newType2 = S1type2 + 30;
	}
	convSiteType(S2, (kmcSiteType)newType2);

	//Update combined sites
	Spointer S3 = moveIt(S1, -2); Spointer S4 = moveIt(S2, 2);
	updateCombinedSites(S3); updateCombinedSites(S4);
	updateCombinedSites(stt);
	updateCombinedSites(S1); updateCombinedSites(S2);

 
    //printSites(stt);
    // update H count
    addCount(0,-2);
    // add ring counts
    m_pah->m_rings5_Lone++;
    //cout<<"WARNING: BY5 closure called. Process not specified yet.\n";
}

// ************************************************************
// ID19- R6 desorption at bay -> pyrene (AR21 in Matlab)
// ************************************************************
void PAHProcess::proc_M6R_BY5_FE3(Spointer& stt, rng_type &rng) {
    //printSites(stt);

    // check if there are any bridges in the BY5, cancel process if yes
    //Cpointer now=C_1->C2;
    //for(int i=0; i!=3; i++) {
    //    if(now->bridge) return;
    //    else now = now->C2;
    //}
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
    if(b4) {
        sFE3 = moveIt(stt, -2);
    }else {
        sFE3 = moveIt(stt, 2);
    }
    // then remove all bulk C in the BY5
    //for(int i=0; i!=3; i++) removeC(C_1->C2, true);
    // then add a C to convert it to R6
    //addC(C_1, normAngle(C_1->bondAngle1+120), normAngle(C_1->bondAngle1+60));
    //// next remove R6 neighbouring with the original R5. removing bulk C:
    //Cpointer Cstart;
    //if(b4) Cstart = C_1->C1->C1->C1->C1;
    //else Cstart = C_2;
    //for(int i=0; i!=3; i++) removeC(Cstart->C2, false);
    //// then add a C
    //addC(Cstart, normAngle(Cstart->bondAngle1-120), normAngle(Cstart->bondAngle1-60), true);
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
        updateSites(sFE3, +1); // convert FE3 (FE) to ZZ
        convSiteType(stt, FE); // convert BY5 to FE
        updateSites(S1, -1); // remove a bulk C from S1
        updateSites(S2, +1); // add a bulk C to S2
    } else {
        S1 = moveIt(stt, -1); // neighbour of BY5 (stt)
        S2 = moveIt(sFE3, 1); // neighbour of FE3
        updateSites(sFE3, +1); // convert FE3 (FE) to ZZ
        convSiteType(stt, FE); // convert BY5 to FE
        updateSites(S1, +1); // add a bulk C to S1
        updateSites(S2, -1); // remove a bulk C from S2
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

    addCount(-2, -2);
    //printSites(stt);
   // cout<<sp.None;
}
// ************************************************************
// ID20 & ID21- R6 oxidation at ZZ site
// ************************************************************
void PAHProcess::proc_O6R_FE2(Spointer& stt) {
    //printSites(stt);
    // identify if the other FE site is before or after this site
    //std::ostringstream dotname, dotname2;
    //dotname << "KMC_DEBUG/" << stt->C1 << "_1.dot";
    //dotname2 << "KMC_DEBUG/" << stt->C1 << "_2.dot";
    //saveDOT(dotname.str());
    bool b4 = false; // <-- position of other FE site
    Spointer other = moveIt(stt, -1);
    Spointer S1, S2;
    if(other->type == FE) {
        b4 = true;
        S1 = moveIt(other, -1);
        S2 = moveIt(stt, 1);
    }
    else {
        other = moveIt(stt, 1);
        S1 = moveIt(other, 1);
        S2 = moveIt(stt, -1);
    }
    // Check if there are any R5 nearby, if yes, don't perform process
    // then check for nearby bridges
    if(b4) {
        if((int) moveIt(other, -1)->type > 4 || (int) moveIt(stt, 1)->type > 4)
            return;
    } else {
        if((int) moveIt(other, 1)->type > 4 || (int) moveIt(stt, -1)->type > 4)
            return;
    }
    //Cpointer C_bulk;
    //for(C_bulk = S1->C1; C_bulk != S1->C2; C_bulk=C_bulk->C2) {
    //    if(C_bulk->bridge) return;
    //}
    //for(C_bulk = S2->C1; C_bulk != S2->C2; C_bulk=C_bulk->C2) {
    //    if(C_bulk->bridge) return;
    //}
    // Identify C1 & C2 for new ZZ site
    //Cpointer C_start, C_end;
    //if(b4) {
    //    C_start = stt->C1->C1->C1;
    //    C_end = stt->C2->C2;
    //} else {
    //    C_start = stt->C1->C1;
    //    C_end = stt->C2->C2->C2;
    //}
    //// Remove the 3 C atoms after C_start
    //removeC(C_start->C2, false);
    //removeC(C_start->C2, false);
    //removeC(C_start->C2, false);
    // add a C atom after C_start (bulk in ZZ)
    //addC(C_start, normAngle(C_start->bondAngle1-120), normAngle(C_start->bondAngle1-60), true);
    //// update H
    //updateA(C_start, C_end, 'H');
    // Remove one of the FE in FE2
    removeSite(other);
    // Update Sites and neighbouring sites
    Spointer S3, S4;
    S1 = moveIt(stt, -1); S2 = moveIt(stt, 1);
    updateSites(stt, +1); // FE --> ZZ
    updateSites(S1, -1); // S1 --> reduce 1
    updateSites(S2, -1); // S2 --> reduce 1
    // update combined sites for all sites and their neighbours
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt);
    updateCombinedSites(S1); updateCombinedSites(S2);
    updateCombinedSites(S3); updateCombinedSites(S4);
    addCount(-2,-1);
    m_pah->m_rings--;
    //saveDOT(dotname2.str());
}
//
// ************************************************************
// ID22- Bay-capping
// ************************************************************
void PAHProcess::proc_B6R_ACR5(Spointer& stt) {
    //printSitesMemb(stt);
    //printStruct();//++++
    //Cpointer newC1;
    //Cpointer newC2;

	int S1t = (int)moveIt(stt, -1)->type;
	int S2t = (int)moveIt(stt, 1)->type;
	int S3t = (int)moveIt(stt, -2)->type;
	int S4t = (int)moveIt(stt, 2)->type;
	int S5t = (int)moveIt(stt, -3)->type;
	int S6t = (int)moveIt(stt, 3)->type;
	int S7t = (int)moveIt(stt, -4)->type;
	int S8t = (int)moveIt(stt, 4)->type;
	int St = (int)stt->type;
	//cout << "AC Begin" << endl << St << endl << S1t << endl << S2t << endl << S3t << endl << S4t << endl << S5t << endl << S6t << endl << "Done" << endl;

	if (S1t == BY6 || S2t == BY6){
		return; //This should not happen
	}

	if ((S1t == AC && S3t == BY6) || (S2t == AC && S4t == BY6)){
		return; //This should not happen
	}

	if ((S1t == FE && S3t == BY6 && S5t == AC) || (S2t == FE && S4t == BY6 && S6t == AC)){
		return; //This should not happen
	}

	if ((S1t == FE && S3t == BY6 && S5t == FE && S7t == BY6) || (S2t == FE && S4t == BY6 && S6t == FE && S8t == BY6)){
		return; //This should not happen
	}

	if ((S1t == FE && S3t == BY5) || (S2t == FE && S4t == BY5)){
		return; //This should not happen
	}

	if ((S1t == AC && S3t == BY5) || (S2t == AC && S4t == BY5)){
		return; //This should not happen
	}

	if ((S1t == FE && S3t == ERZZ) || (S2t == FE && S4t == ERZZ)){
		return; //This should not happen
	}

	if ((S1t == ACR5 && S3t != FE) || (S2t == ACR5 && S4t != FE)){
		return; //This should not happen
	}

	if ((S1t == ZZ && S3t == BY6 && S5t == ZZ) || (S2t == ZZ && S4t == BY6 && S6t == ZZ)){
		return; //This should not happen
	}

    //if(checkHindrance(stt)) {
    //    /*cout<<"Site hindered, process not performed.\n"*/ return;}
    //if(!(C_1->C2->bridge) || !(C_2->C1->bridge)) { // check if bulk C in AC site is a bridge
    //    // Add and remove C
    //    //if(C_1->C2->C3 != NULL) C_1->C2->C3->C3 = NULL;
    //    //if(C_2->C1->C3 != NULL) C_2->C1->C3->C3 = NULL;
    //    removeC(C_1->C2, true);
    //    removeC(C_2->C1, true);
    //    newC1 = addC(C_1, normAngle(C_1->bondAngle1+120), 0);
    //    newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1-60), normAngle(newC1->C1->bondAngle1-120));
    //}else {
    //    // update bridges info, both are no longer bridges
    //    C_1->C2->C1 = C_1->C2->C3;// update neighbour
    //    C_2->C1->C2 = C_2->C1->C3;
    //    C_1->C2->bridge = false;
    //    C_2->C1->bridge = false;
    //    angletype a = C_2->C1->bondAngle1;
    //    C_2->C1->bondAngle1 = C_2->C1->bondAngle2;
    //    C_2->C1->bondAngle2 = a;
    //    //C_1->C2->C3 = NULL;
    //    //C_2->C1->C3 = NULL;
    //    // connect C_1 and C_2
    //    connectToC(C_1, C_2);
    //    // Add C
    //    newC1 = addC(C_1, normAngle(C_1->bondAngle1+120), 0);
    //    newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1-60), normAngle(newC1->C1->bondAngle1-120));
    //}
    ////printStruct();
    //// Add and remove H
    //updateA(C_1->C1, C_2->C2, 'H');
    // neighbouring sites:
    Spointer S1 = moveIt(stt, -1); 
    Spointer S2 = moveIt(stt, 1);
    // Update Site and neighbours
    convSiteType(stt, (kmcSiteType) 0);

	//if ((S1t == eRZZ && S3t == eR5 && S5t == eRFE) || (S2t == eRZZ && S4t == eR5 && S6t == eRFE)){
	//	//// Convert to a BY6 site
	//	bool b4 = true;
	//	Spointer Srem1, Srem2, Srem3;
	//	if (S1t == eRZZ){
	//		Srem1 = moveIt(stt, -1);
	//		Srem2 = moveIt(stt, -2);
	//		Srem3 = moveIt(stt, -3);
	//	}
	//	else{
	//		Srem1 = moveIt(stt, 1);
	//		Srem2 = moveIt(stt, 2);
	//		Srem3 = moveIt(stt, 3);
	//		b4 = false;
	//	}
	//	//// Remove sites and combine the neighbouring sites into BY6. 
	//	//// First remove all three from site map. Elementary site types first..
	//	delSiteFromMap(Srem1->type, Srem1);
	//	delSiteFromMap(Srem2->type, Srem2);
	//	delSiteFromMap(Srem3->type, Srem3);
	//	//// then for combined site types..
	//	delSiteFromMap(Srem1->comb, Srem1);
	//	delSiteFromMap(Srem2->comb, Srem2);
	//	delSiteFromMap(Srem3->comb, Srem3);
	//	//// remove the sites
	//	removeSite(Srem1);
	//	removeSite(Srem2);
	//	removeSite(Srem3);

	//	////add in BY6 site
	//	if (b4){
	//		addSite(BY6, stt);
	//		updateSites(S2, 1);
	//	}
	//	else{
	//		addSite(BY6, moveIt(stt, 1));
	//		updateSites(S1, 1);
	//	}

	//}
	//else if ((S1t == eRFE && S3t == eR5 && S5t == eRZZ) || (S2t == eRFE && S4t == eR5 && S6t == eRZZ)){
	//	//// Convert to a BY6 site
	//	bool b4 = true;
	//	Spointer Srem1, Srem2, Srem3;
	//	if (S1t == eRFE){
	//		Srem1 = moveIt(stt, -1);
	//		Srem2 = moveIt(stt, -2);
	//		Srem3 = moveIt(stt, -3);
	//	}
	//	else{
	//		Srem1 = moveIt(stt, 1);
	//		Srem2 = moveIt(stt, 2);
	//		Srem3 = moveIt(stt, 3);
	//		b4 = false;
	//	}

	//	//// Remove sites and combine the neighbouring sites into BY6. 
	//	//// First remove all three from site map. Elementary site types first..
	//	delSiteFromMap(Srem1->type, Srem1);
	//	delSiteFromMap(Srem2->type, Srem2);
	//	delSiteFromMap(Srem3->type, Srem3);
	//	//// then for combined site types..
	//	delSiteFromMap(Srem1->comb, Srem1);
	//	delSiteFromMap(Srem2->comb, Srem2);
	//	delSiteFromMap(Srem3->comb, Srem3);
	//	//// remove the sites
	//	removeSite(Srem1);
	//	removeSite(Srem2);
	//	removeSite(Srem3);

	//	////add in BY6 site
	//	if (b4){
	//		addSite(BY6, stt);
	//		updateSites(S2, 1);
	//	}
	//	else{
	//		addSite(BY6, moveIt(stt, 1));
	//		updateSites(S1, 1);
	//	}

	//}
	//else {
	//	if (S1->type == ACR5){
	//		convSiteType(S1, eRZZ);
	//		addSite(eRFE, S1);
	//		addSite(eR5, S1);
	//	}
	//	else {
	//		updateSites(S1, 1);
	//	}

	//	if (S2->type == ACR5){
	//		Spointer S3 = moveIt(S2, 1);
	//		convSiteType(S2, eRZZ);
	//		addSite(eR5, S3);
	//		addSite(eRFE, S3);

	//	}
	//	else {
	//		updateSites(S2, 1);
	//	}
	//}

	updateSites(S1, 1);
	updateSites(S2, 1);

    // Update combined site for Site and neighbours
    Spointer S3, S4;
	S1 = moveIt(stt, -1); S2 = moveIt(stt, 1);
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt);
    updateCombinedSites(S1); updateCombinedSites(S2);
    updateCombinedSites(S3); updateCombinedSites(S4);
    // add ring counts
    m_pah->m_rings++;
    m_pah->m_rings5_Lone--;
    m_pah->m_rings5_Embedded++;
	addCount(2, 0);
    //printSites(stt);
}

// ************************************************************
// ID23- Embedded 5-member ring migration to ZZ
// ************************************************************
void PAHProcess::proc_M5R_eR5_FE3_ZZ(Spointer& stt, rng_type &rng) {
    bool b4 = false;
    if(moveIt(stt, -1)->comb == FE2 || moveIt(stt, 1)->comb == FE2) {
        if(moveIt(stt,-1)->comb == FE2 && moveIt(stt,1)->comb == FE2) {
            // Define a distribution that has two equally probably outcomes
            boost::bernoulli_distribution<> choiceDistrib;
            // Now build an object that will generate a sample using rng
            boost::variate_generator<rng_type&, boost::bernoulli_distribution<> > choiceGenerator(rng, choiceDistrib);

            b4= choiceGenerator(); // if FE3 on both sides, choose a random one
        }
        else {
            if(moveIt(stt,-1)->comb == FE2)
                b4 = true;
        }    
    }else{
        return;    
    }
    
    Spointer sFE2;
    //Cpointer CRem, CFE;
    if(b4) {
        sFE2 = moveIt(stt, -1); 
        //CFE = sFE2->C2;
        //CRem = sFE2->C1;
    }else {
        sFE2 = moveIt(stt, 1);
        //CFE = sFE2->C1;
        //CRem = sFE2->C2;
    }
    //CRem->C1->bondAngle1 = normAngle(CRem->C1->bondAngle1-30);
    //removeC(CRem, false);
    //// add a C atom
    //addC(C_1, normAngle(C_1->bondAngle1+30), normAngle(C_1->bondAngle1-30));
    // delete neighbouring (to sFE3) FE sites from site map. Identify them first:
    // then remove them
    if(b4){
        Spointer Srem1 = moveIt(sFE2,-1);
        removeSite(Srem1);
    }else{
        Spointer Srem2 = moveIt(sFE2,1);
        removeSite(Srem2);
    }
    // edit sites. first identify the neighbouring sites of resulting RFE & R5
    Spointer S1, S2, S3, S4;
    //Cpointer C21 = C_2->C1;
    //Cpointer C22 = C_2->C2;
    //Cpointer C11 = C_1->C1;
    //Cpointer C12 = C_1->C2;
    if(b4) {
        S1 = moveIt(sFE2, -1); // neighbour of FE3
        Spointer stt2;
        stt2 = moveIt(stt, 1);
        convSiteType(sFE2, R5); // convert FE3 to R5
        convSiteType(stt, RFE); // convert eR5 to RFE
        addSite(ZZ, stt2);        
        updateSites(S1, -1);
        addR5toSite(S1); // remove a bulk C from S1 and add R5
    } else {
        S2 = moveIt(sFE2, 1); // neighbour of FE3
        convSiteType(sFE2, R5); // convert FE3 to R5
        convSiteType(stt, RFE); // convert eR5 to RFE
        addSite(ZZ, stt);
        updateSites(S2, -1);
        addR5toSite(S2); // remove a bulk C from S2 and add R5
    }
    // update H atoms
    //if(b4){
    //    updateA(S1->C2->C1, C_2, 'H');
    //}else{
    //    updateA(C_1, S2->C1->C2, 'H');
    //}
    
    // update combined sites for all sites involved and their neighbours
    // (excluding new FE sites, since their combined site type will still be None)
    if(b4){
        S3 = moveIt(S1, -1);
    }else{
        S4 = moveIt(S2, 1);
    }
    updateCombinedSites(stt); updateCombinedSites(sFE2); // new R5 and RFE
    if(b4){
        updateCombinedSites(S1); updateCombinedSites(S3);
    }else{
        updateCombinedSites(S2); updateCombinedSites(S4); // neighbours
    }
}

// ************************************************************
// ID24 - R6 growth on RZZ 
// ************************************************************
void PAHProcess::proc_G6R_RZZ(Spointer& stt) {
    proc_G6R_AC(stt);
}

// ************************************************************
// ID25 - R6 growth on RFER 
// ************************************************************
void PAHProcess::proc_G6R_RFER(Spointer& stt) {
    proc_G6R_AC(stt);
}

// ************************************************************
// ID26 - R6 growth on R5
// ************************************************************
void PAHProcess::proc_G6R_R5(Spointer& stt) {
    proc_G6R_FE(stt);
}

// ************************************************************
// ID27 - RBY5 closure reaction
// ************************************************************
void PAHProcess::proc_L6_RBY5(Spointer& stt) {
    proc_L6_BY6(stt);
}

// ************************************************************
// ID28 - RACR closure reaction
// ************************************************************
void PAHProcess::proc_L6_RACR(Spointer& stt) {
    proc_L6_BY6(stt);
}

// ************************************************************
// ID29 - R5 growth on RFE 
// ************************************************************
void PAHProcess::proc_G5R_RFE(Spointer& stt) {
    proc_G5R_ZZ(stt);
}

// ************************************************************
// ID30 - R6 migration & conversion to R5 at RAC
// ************************************************************
void PAHProcess::proc_C6R_RAC_FE3(Spointer& stt, rng_type &rng) {
    proc_C6R_BY5_FE3(stt, rng);
}

// ************************************************************
// ID31 - R6 migration & conversion to R5 at RAC
// ************************************************************
void PAHProcess::proc_C6R_RAC_FE3violi(Spointer& stt, rng_type &rng) {
    proc_C6R_BY5_FE3(stt, rng);
}

// ************************************************************
// ID32 - R6 desorption at RAC -> pyrene
// ************************************************************
void PAHProcess::proc_M6R_RAC_FE3(Spointer& stt, rng_type &rng) {
    proc_M6R_BY5_FE3(stt, rng);
}

// ************************************************************
// ID33 - R6 growth on AC
// ************************************************************
void PAHProcess::proc_G6R_ACBR(Spointer& stt) {
	proc_G6R_AC(stt);
}

// ************************************************************
// ID34 - R6 growth on AC
// ************************************************************
void PAHProcess::proc_G6R_ACBL(Spointer& stt) {
	proc_G6R_AC(stt);
}

size_t PAHProcess::SiteListSize() const {
    return m_pah->m_siteList.size();
}

//size_t PAHProcess::CarbonListSize() const {
//    return m_pah->m_carbonList.size();
//}

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
