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


/*!
 * Save the PAH's structure as a DOT file. These can be converted into images
 * using the program Graphviz. E.g. to create EPS figures:
 * neato filename.dot -O -Teps
 *
 * @param filename   The filename to write the DOT file for
 */
void PAHStructure::SaveDOT(const string& filename, const string& title="") const {
	// Create and open the output stream
    ofstream dotfile;
    dotfile.open(filename.c_str());

    if(dotfile.good()) { // if successful
        // vector to store indices of C bonded with H atom
        std::vector<int> CH_index;
        // set iterator to starting position
        Cpointer now = m_cfirst->C2;
        Cpointer prev = m_cfirst;

        // draw undirected graph
        dotfile<< "Graph G {\n";
        // sets max size for drawing structure (in inches)
        dotfile<< "\tgraph [size=15];\n";
        // sets style for C nodes
        dotfile<< "\tnode [shape=circle, style=filled, color=black, height=0.2, width=0.2];\n";
        // sets style for bonds
        dotfile<< "\tedge [weight=2, style=\"setlinewidth(10)\"];\n";
        // sets caption for graph
        if(title!="none") {
            dotfile<< "\tgraph [label=\""<<title<<"\"];\n";
        }
        // get coordinates of current C
        coordtype x=prev->coords.first;
        coordtype y=prev->coords.second;
        int c=0; // index of current C
        int count=0;
        float CC_CHratio = 10; // length ratio between C-C and C-H bonds
        do {
            //if(now->bridge) {
                //cout<<"Bridge at "<<now<<'\n';
            //}
            c++;
            // write position coordinates of C atom with index c, with ! suffixed to pin node
            dotfile<<"\tC"<<c<<" [pos = \""<<x<<','<<y<<"!\"";
            if(c==1) dotfile<<", color=red";
            dotfile<<"];\n";
            if(prev->A=='H') { // if bonded with H
                // write position coordinates of H atom with bonded C index, with ! suffixed to pin node
                dotfile<<"\tH"<<c<<" [pos = \""<<(x+x_inc(prev->bondAngle2)/CC_CHratio);
                // set style for H node
                dotfile<<','<<(y+y_inc(prev->bondAngle2)/CC_CHratio)<<"!\", color=skyblue, fontcolor=skyblue, height=0.1, width=0.1];\n";
                // stores index of current C
                CH_index.push_back(c);
            }
            // move iterator to next C
            Cpointer oldnow = now;
            now = Carbon::MoveCPointer(prev, now);
            prev = oldnow;
            // obtain coordinates
            x=prev->coords.first;
            y=prev->coords.second;
        }
        while
        (!(count !=1 && prev->coords.first == m_cfirst->coords.first && prev->coords.second == m_cfirst->coords.second));
        // connects each C atom
        for(int i=1; i<c; i++) {
            dotfile<<"\tC"<<i<<" -- "<<"C"<<(i+1)<<";\n";
        }
        dotfile<<"\tC"<<c<<" -- "<<"C1;\n"; // connects last C with first C
        // connects each H atom with corresponding C atoms
        for(int i=0; i!=(int)CH_index.size(); i++) {
            dotfile<<"\tH"<<CH_index[i]<<" -- "<<"C"<<CH_index[i]<<";\n";
        }
        // close dot file
        dotfile<<"}";
        dotfile.close();
    }
}


void PAHStructure::saveDOTperLoop(int PAH_ID, int i)
{

    string filename = "KMC_DEBUG/ID_";
    filename.append(Strings::cstr(PAH_ID));
    filename.append("_");
    filename.append(Strings::cstr(i));
    filename.append(".dot");
    SaveDOT(filename, Strings::cstr(PAH_ID));
}

/*!
 * Print the PAH structure to console. Used for debugging purposes.
 */
void PAHStructure::PrintStruct() const {
    // iterator set to first C in structure
    Cpointer prev = m_cfirst;
    Cpointer now = m_cfirst->C2;
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
        std::cout << "coords ("<<prev->coords.first<<","<< prev->coords.second<<")"<< btxt << ": C-" << prev->A << "\t" << angle << '\n';
        // moves iterator to next C
        Cpointer oldnow = now;
        now = Carbon::MoveCPointer(prev, now);
        prev = oldnow;
    } while
    (!(count!=1&&prev->coords.first == m_cfirst->coords.first
   &&prev->coords.second == m_cfirst->coords.second));
    // displays C and H counts
    std::cout << "C Count: " << m_counts.first << '\n';
    std::cout << "H Count: " << m_counts.second << '\n';
}

/*!
 * Print the site list to console. Used for debugging only.
 */
void PAHStructure::PrintSites() const {
    std::string st;
    cout << "*******************\n";
    cout << "Sites List:\n_____\n";
    // displays total site count
    cout << "Total Site Count: " << m_siteList.size() << '\n';
    for(std::list<Site>::const_iterator i=m_siteList.begin(); i!=m_siteList.end(); i++) {
        // convert site type into string
        st = kmcSiteName(i->type);
        // displays site type
        cout << st << '\n';
    }
    cout << "********************\n";
}

/*!
 * Print sites and site members to console. Used for debugging only.
 */
void PAHStructure::PrintSitesMemb() const {
    std::string st;
    string pointer = "  <-  ";
    cout << "*******************\n";
    cout << "Sites List:\n_____\n";
    // displays total site count
    cout << "Total Site Count: " << m_siteList.size() << '\n';
    for(std::list<Site>::const_iterator i=m_siteList.begin(); i!=m_siteList.end(); i++) {
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
    double val=0.0;

    std::set<cpair>::iterator itEnd=m_cpositions.end();
    for (std::set<cpair>::iterator it=m_cpositions.begin();it!=itEnd;++it)
    {
        val=(*it).first;
        out.write((char*)&val, sizeof(val));
        val=(*it).second;
        out.write((char*)&val, sizeof(val));
    }
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


