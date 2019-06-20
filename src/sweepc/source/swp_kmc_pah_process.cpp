/*!
  * \author     Zakwan Zainuddin (zz260) Gustavo Leon (gl413)
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
    std::vector<kmcSiteType> sites = SiteVector();
	p.createPAH(sites, m_pah->m_rings, m_pah->m_rings5_Lone, m_pah->m_rings5_Embedded, m_pah->m_rings7_Lone, m_pah->m_rings7_Embedded, m_pah->m_InternalCarbons);
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
        if(getCHCount().first == 6) return 0;
        return (unsigned int) (m_pah->m_siteMap[st].size());
    }
    return (unsigned int) m_pah->m_siteMap[st].size();
}
//! Get Ring Counts
std::tuple <int, int, int> PAHProcess::getRingsCount() const {
	std::tuple <int, int, int> rings_tuple = std::make_tuple(m_pah->m_rings, m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded, m_pah->m_rings7_Lone + m_pah->m_rings7_Embedded);
	return rings_tuple;
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
//! Print Structure in console
void PAHProcess::printStruct() const{
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
		std::cout << "coords (" << std::get<0>(prev->coords) << "," << std::get<1>(prev->coords) << "," << std::get<2>(prev->coords) << ")" << btxt << ": C-" << prev->A << "\t" << angle << '\n';
        // moves iterator to next C
        Cpointer oldnow = now;
        now = moveCPointer(prev, now);
        prev = oldnow;
    } while 
		(!(count != 1 && std::get<0>(prev->coords) == std::get<0>(m_pah->m_cfirst->coords)
		&& std::get<1>(prev->coords) == std::get<1>(m_pah->m_cfirst->coords)
		&& std::get<2>(prev->coords) == std::get<2>(m_pah->m_cfirst->coords)));
    // displays C and H counts
    std::cout << "C Count: " << m_pah->m_counts.first << '\n';
    std::cout << "H Count: " << m_pah->m_counts.second << '\n';
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
    cout << "Total Site Count: " << m_pah->m_siteList.size() << '\n';
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
    cout << "Total Site Count: " << m_pah->m_siteList.size() << '\n';
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
    cout << "Total Site Count: " << m_pah->m_siteList.size() << '\n';
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
		mol = optimisePAH(mol, 3000, "Ghemical");
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
void PAHProcess::save_trajectory_xyz(const double timer, const std::string &filename, bool optimise) {
	OpenBabel::OBMol mol = passPAH(optimise);
	if (optimise){
		mol = optimisePAH(mol, 3000, "mmff94");
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
	//check if vector is unitary
	cpair temp2 = scale_vector(temp);
    return temp2;
}
//! Rescales a vector.
cpair PAHProcess::scale_vector(cpair vec) const{
	//check if vector is unitary
	double tol = 1e-3;
	cpair temp = vec;
	if (abs(std::get<0>(vec)*std::get<0>(vec) + std::get<1>(vec)*std::get<1>(vec) + std::get<2>(vec)*std::get<2>(vec) - 1.0)  > tol)
	{
		double magnitude = sqrt(std::get<0>(vec)*std::get<0>(vec) + std::get<1>(vec)*std::get<1>(vec) + std::get<2>(vec)*std::get<2>(vec));
		temp = std::make_tuple(std::get<0>(vec)/magnitude, std::get<1>(vec)/magnitude, std::get<2>(vec)/magnitude);
		/*std::get<0>(temp) = std::get<0>(temp)/magnitude;
		std::get<1>(temp) = std::get<1>(temp)/magnitude;
		std::get<2>(temp) = std::get<2>(temp)/magnitude;*/
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
	if (theta <= 180.0){
		std::get<0>(temp) = std::get<0>(temp) * -1.0;
		std::get<1>(temp) = std::get<1>(temp) * -1.0;
		std::get<2>(temp) = std::get<2>(temp) * -1.0;
	}
	cpair temp2 = scale_vector(temp);
	return temp2;
}
//! Returns the cross product of two vectors. If v1 is the surface normal vector and v2 goes from C->newC v1xv2 redefines the next growth vector.
cpair PAHProcess::cross_vector (cpair vec1, cpair vec2) const{
	//check if vector is unitary
	cpair vec1_adj = scale_vector(vec1);
	cpair vec2_adj = scale_vector(vec2);
	cpair temp = std::make_tuple(std::get<1>(vec1_adj) * std::get<2>(vec2_adj) - std::get<2>(vec1_adj) * std::get<1>(vec2_adj), std::get<2>(vec1_adj) * std::get<0>(vec2_adj) - std::get<0>(vec1_adj) * std::get<2>(vec2_adj), std::get<0>(vec1_adj) * std::get<1>(vec2_adj) - std::get<1>(vec1_adj) * std::get<0>(vec2_adj));
	cpair temp2 = scale_vector(temp);
	return temp2;
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
		if (abs(std::get<0>(coords) - std::get<0>(*it)) < tol && abs(std::get<1>(coords) - std::get<1>(*it)) < tol && abs(std::get<2>(coords) - std::get<2>(*it)) < tol){
			return true;
		}
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
	double tol = 1e0;
	cpair mpos = jumpToPos(C_1->coords, C_1->growth_vector, 1.4);
	for (std::set<cpair>::iterator it = m_pah->m_cpositions.begin(); it != m_pah->m_cpositions.end(); ++it) {
		if (abs(std::get<0>(mpos) - std::get<0>(*it)) < tol && abs(std::get<1>(mpos) - std::get<1>(*it)) < tol && abs(std::get<2>(mpos) - std::get<2>(*it)) < tol){
			return true;
		}
	}
	return false;
}

bool PAHProcess::checkHindrance_twoC(const Cpointer C_1, const Cpointer C_2) const {
	double tol = 1e0;
	if (abs(std::get<0>(C_1->coords) - std::get<0>(C_2->coords)) < tol && abs(std::get<1>(C_1->coords) - std::get<1>(C_2->coords)) < tol && abs(std::get<2>(C_1->coords) - std::get<2>(C_2->coords)) < tol){
		return false;
	}
	else return true;
}

double PAHProcess::getDistance_twoC(const Cpointer C_1, const Cpointer C_2) const {
	double xdist, ydist, zdist;
	xdist = std::get<0>(C_1->coords) - std::get<0>(C_2->coords);
	ydist = std::get<1>(C_1->coords) - std::get<1>(C_2->coords);
	zdist = std::get<2>(C_1->coords) - std::get<2>(C_2->coords);
	return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
}

double PAHProcess::getDistance_twoC(const cpair C_1, const cpair C_2) const {
	double xdist, ydist, zdist;
	xdist = std::get<0>(C_1) - std::get<0>(C_2);
	ydist = std::get<1>(C_1) - std::get<1>(C_2);
	zdist = std::get<2>(C_1) - std::get<2>(C_2);
	return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
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

Spointer PAHProcess::RZZchooseRandomSite(kmcSiteType st, rng_type &rng) {
	// to choose any principal site
		//cout << "~~Choosing from " << m_pah->m_siteMap[st].size() << " sites...\n";
	int R5_unavail_sites = (int)m_pah->m_siteMap[R5R6].size() + (int)m_pah->m_siteMap[R5R6ZZ].size() + (int)m_pah->m_siteMap[R5R6AC].size() + (int)m_pah->m_siteMap[R5R6BY5].size();
	if (R5_unavail_sites > 0){
		Spointer S1, S2;
		std::vector<int> list;
		for (int i = 0; i < (int)m_pah->m_siteMap[st].size(); i++) {
			S1 = moveIt(m_pah->m_siteMap[st][i], -1);
			S2 = moveIt(m_pah->m_siteMap[st][i], +1);
			if (S1->type != R5R6 && S2->type != R5R6 && S1->type != R5R6ZZ && S2->type != R5R6ZZ && S1->type != R5R6AC && S2->type != R5R6AC && S1->type != R5R6BY5 && S2->type != R5R6BY5){
				//RZZ site available for migrations
				list.push_back(i);
			}
		}
		int sz = ((int)m_pah->m_siteMap[st].size()) - R5_unavail_sites - 1; // size - 1
		if (sz >= 0) {
			typedef boost::uniform_smallint<unsigned int> site_index_distrib;
			site_index_distrib siteIndexDistrib(0, static_cast<unsigned int>(sz));
			boost::variate_generator<rng_type &, site_index_distrib> siteIndexGenerator(rng, siteIndexDistrib);

			const unsigned r = siteIndexGenerator();
			return m_pah->m_siteMap[st][list[r]];
		}
		else{
			Spointer normal_site = chooseRandomSite(st, rng);
			return normal_site;
		}
	}
	else{
		Spointer normal_site = chooseRandomSite(st, rng);
		return normal_site;
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
    return cb;
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
//! Adds an R5 to the list of R5s and R7s
void PAHProcess::addR5internal(Cpointer C_1, Cpointer C_2) {
	double R5_dist = getDistance_twoC(C_1, C_2);
	cpair R5dir = get_vector(C_1->coords,C_2->coords);
	cpair normvec = norm_vector(C_1->coords, C_2->coords, C_2->C2->coords);
	cpair crossvec = cross_vector(R5dir, normvec);
	double theta = atan(R5_dist/2.0/0.7);
	double magn = 0.7 / cos(theta);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + 0.7 * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + 0.7 * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ 0.7 * std::get<2>(crossvec));
	cpair intdir = scale_vector(resultantvec);
	cpair mpos = jumpToPos(C_1->coords, intdir, magn);
	m_pah->m_R5loc.push_back(mpos);
}
//! Removes an R5 from the list of R5s and R7s
void PAHProcess::removeR5internal(Cpointer C_1, Cpointer C_2) {
	std::list<cpair>::iterator it1, it2;
	double minimal_dist = 1e3;
	for (it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		double dist_x = std::get<0>(C_1->coords) + std::get<0>(C_2->coords) - 2*std::get<0>(*it1);
		double dist_y = std::get<1>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<1>(*it1);
		double dist_z = std::get<2>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<2>(*it1);
		double dist = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
		if (dist < minimal_dist){
			minimal_dist = dist;
			it2 = it1;
		}
	}
	cpair temp = *it2;
	m_pah->m_R5loc.erase(it2);
}
//! Removes an R5 from the list of R5s and R7s
void PAHProcess::removeR7internal(Cpointer C_1, Cpointer C_2) {
	std::list<cpair>::iterator it1, it2;
	double minimal_dist = 1e3;
	for (it1 = m_pah->m_R7loc.begin(); it1 != m_pah->m_R7loc.end(); ++it1) {
		double dist_x = std::get<0>(C_1->coords) + std::get<0>(C_2->coords) - 2*std::get<0>(*it1);
		double dist_y = std::get<1>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<1>(*it1);
		double dist_z = std::get<2>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<2>(*it1);
		double dist = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
		if (dist < minimal_dist){
			minimal_dist = dist;
			it2 = it1;
		}
	}
	cpair temp = *it2;
	m_pah->m_R7loc.erase(it2);
}
//! Return internal R5 associated to two carbons and deletes it from R5 list.
cpair PAHProcess::findR5internal(Cpointer C_1, Cpointer C_2) {
	std::list<cpair>::iterator it1, it2;
	double minimal_dist = 1e3;
	for (it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		double dist_x = std::get<0>(C_1->coords) + std::get<0>(C_2->coords) - 2*std::get<0>(*it1);
		double dist_y = std::get<1>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<1>(*it1);
		double dist_z = std::get<2>(C_1->coords) + std::get<1>(C_2->coords) - 2*std::get<2>(*it1);
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
//! Are the two carbon atoms members of an R5 with coordinates in R5Internal??
bool PAHProcess::isR5internal(Cpointer C_1, Cpointer C_2) {
	cpair R5_pos_loc = endposR5internal(C_1, C_2);
	std::list<cpair>::iterator it1;
	double minimal_dist = 0.5;
	for (it1 = m_pah->m_R5loc.begin(); it1 != m_pah->m_R5loc.end(); ++it1) {
		double dist = getDistance_twoC(R5_pos_loc, *it1);
		if (dist <= minimal_dist) return true;
	}
	return false;
}
//! Return coords of final position of an internal R5 based on two carbons
cpair PAHProcess::endposR5internal(Cpointer C_1, Cpointer C_2) {
	//C_1 is the carbon that stays in the R5, C_2 the R6 that becomes part of the R5
	double R5_dist = getDistance_twoC(C_1, C_2);
	cpair R5dir = get_vector(C_1->coords,C_2->coords);
	cpair normvec = invert_vector(norm_vector(C_1->coords, C_2->coords, C_2->C2->coords));
	cpair crossvec = cross_vector(R5dir, normvec);
	double theta = atan(R5_dist/2.0/0.7);
	double magn = 0.7 / cos(theta);
	cpair resultantvec = std::make_tuple(R5_dist/2.0 * std::get<0>(R5dir) + 0.7 * std::get<0>(crossvec), R5_dist/2.0 * std::get<1>(R5dir) + 0.7 * std::get<1>(crossvec), R5_dist/2.0 * std::get<2>(R5dir)+ 0.7 * std::get<2>(crossvec));
	cpair intdir = scale_vector(resultantvec);
	cpair mpos = jumpToPos(C_1->coords, intdir, magn);
	return mpos;
}

//! Passes a PAH from MOpS to OpenBabel. Returns a mol object.
OpenBabel::OBMol PAHProcess::passPAH(bool detectBonds) {
	//R6 Bay detection
	std::list<int> R6pairs1, R6pairs2;
	std::list<Cpointer>CR6_pair1, CR6_pair2;
	std::list<Cpointer>::iterator resR6, resR62;
	//Second neighbours bond detection
	std::list<int> C_intlist, first_neighbour, second_neighbour;
	std::list<Cpointer>C_list, C_first_neighbour, C_second_neighbour;
	std::list<int>::iterator sn_iter, sn_iter1, sn_iter2;
	if (detectBonds){
		for (std::list<Site>::iterator site_it = m_pah->m_siteList.begin(); site_it != m_pah->m_siteList.end(); site_it++) {
			if ( (int)site_it->type % 10 == 4 ){
				Cpointer CR6_1 = site_it->C1;
				Cpointer CR6_2 = site_it->C2;
				if (getDistance_twoC (CR6_1,CR6_2) <=2.17){
					CR6_pair1.push_back(CR6_1);
					CR6_pair2.push_back(CR6_2);
				}
			}
		}
		Ccontainer::iterator itCsn;
		for (itCsn = m_pah->m_carbonList.begin(); itCsn != m_pah->m_carbonList.end(); ++itCsn) {
			Cpointer C_check = *itCsn;
			C_list.push_back(C_check);
			C_first_neighbour.push_back(C_check->C2);
			C_second_neighbour.push_back(C_check->C2->C2);
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
			cout << "NaN in coordinates while passing them to OB.\n";
		} 
		//atom->SetHyb(2);
		
		//BY6 bond detection.
		if (detectBonds){
			resR6 = std::find(std::begin(CR6_pair1), std::end(CR6_pair1), C_change);
			if (resR6 != std::end(CR6_pair1)){
				R6pairs1.push_back(atom->GetIdx());
				R6C1 = true;
			}
			resR62 = std::find(std::begin(CR6_pair2), std::end(CR6_pair2), C_change);
			if (resR62 != std::end(CR6_pair2)){
				R6pairs2.push_back(atom->GetIdx());
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
			cout << "NaN in internal coordinates";
		} 
		//atom->SetHyb(2);
		counter ++;
	}
	
	if (detectBonds){
		//Bond generation
		mol.ConnectTheDots();
		//connectPAH(mol);
		//Deletes bonds longer than 1.75A. 
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
		//Adds Bonds within C and its C2. Deletes bonds between C and C->C2->C2
		std::list<int>::iterator sn_iter;
		std::list<int>::iterator sn_iter1 = first_neighbour.begin();
		std::list<int>::iterator sn_iter2 = second_neighbour.begin();
		for(sn_iter = C_intlist.begin(); sn_iter != C_intlist.end(); ++sn_iter){
			OpenBabel::OBBond* my_bond = mol.GetBond(*sn_iter, *sn_iter1);
			OpenBabel::OBBond* my_bond2 = mol.GetBond(*sn_iter, *sn_iter2);
			mol.AddBond(*sn_iter, *sn_iter1,5);
			if (my_bond2 != NULL) mol.DeleteBond(my_bond2);
			++sn_iter1;
			++sn_iter2;
		}
		
		//Deletes bonds within unclosed BY6.
		std::list<int>::iterator it_R6pairs1, it_R6pairs2;
		it_R6pairs2 = R6pairs2.begin();
		for(it_R6pairs1 = R6pairs1.begin(); it_R6pairs1 != R6pairs1.end(); ++it_R6pairs1){
			OpenBabel::OBBond* my_bond = mol.GetBond(*it_R6pairs1, *it_R6pairs2);
			if (my_bond != NULL) mol.DeleteBond(my_bond);
			++it_R6pairs2;
		}
		//Assign aromatic bond orders to all bonds.
		/*for (OpenBabel::OBBondIterator bond_iter=mol.BeginBonds(); bond_iter != mol.EndBonds(); bond_iter++){
			OpenBabel::OBBond* my_bond = *bond_iter;
			int BOrder = my_bond->GetBO();
			int a1 = my_bond->GetBeginAtomIdx(); int a2 = my_bond->GetEndAtomIdx(); 
			cout << a1 << "-" << a2 << " Bond order = " << BOrder << "\n";
		}*/
		
		//PerceiveBondOrders calls several routines that try to identify aromatic and unaromatic parts of a molecule. This is very nice but expensive.
		//If this is not called, OB recognises bonds as single bonds.
		mol.PerceiveBondOrders();
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
	return mol;
}
//Needed for connect the dots.
bool SortAtomZ(const std::pair<OpenBabel::OBAtom*,double> &a, const std::pair<OpenBabel::OBAtom*,double> &b){
    return (a.second < b.second);
}

//! Connects the atoms in a PAH using OpenBabel routines. Equivalent to OpenBabel::OBMol::ConnectTheDots();
void PAHProcess::connectPAH(OpenBabel::OBMol my_mol) {
    int j,k,max;
    bool unset = false;
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

//! Passes a PAH from MOpS to OpenBabel. Returns a mol object.
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
	vector<OpenBabel::OBRing*> *rlist = (vector<OpenBabel::OBRing*>*)mol.GetData("RingList");
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
}
int forcefield_error_counter = 0;
//! Minimisation of a PAH
OpenBabel::OBMol PAHProcess::optimisePAH(OpenBabel::OBMol mol, int nsteps, std::string forcefield) {
	mol.BeginModify();
	//Defines a forcefield object
	//OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField("Ghemical");
	OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(forcefield);
	if (!pFF) {
    cerr << ": could not find forcefield MMFF94s." <<endl;
    exit (-1);
	}
	/*pFF->SetLogFile(&cerr);
	pFF->SetLogLevel(OBFF_LOGLVL_LOW);
	pFF->SetVDWCutOff(6.0);
	pFF->SetElectrostaticCutOff(10.0);
	pFF->SetUpdateFrequency(10);
	pFF->EnableCutOff(false);*/
	
	//Initialise minimisation
	if (!pFF->Setup(mol)) {
      cout << "Error: could not setup force field.\n" << endl;
	  cout << "Sites before calling optimiser:\n";
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
	bool done = true;
	pFF->SteepestDescentInitialize(nsteps, 1e-6);
	
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
	}
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
	stype += 101;
	if ((kmcSiteType)stype == None) {
		std::cout << "ERROR: add5RtoSite: illegal site type to add 5R to.\n";
		printSites(st);
		return;
	}
	// removes site from m_pah->m_siteMap (principal site)
	delSiteFromMap(st->type, st);
	// change site type
	st->type = (kmcSiteType)stype;
	// add site to m_pah->m_siteMap
	m_pah->m_siteMap[st->type].push_back(st);
	// update member C
	st->C1 = Carb1;
	st->C2 = Carb2;

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
    int stype = (int) st->type;
	stype -= 101;
	if ((kmcSiteType)stype == None) {
		std::cout << "ERROR: rem5RfromSite: illegal site type to remove 5R from.\n";
		printSites(st);
		return;
	}
	// removes site from m_pah->m_siteMap (principal site)
	delSiteFromMap(st->type, st);
	// change site type
	st->type = (kmcSiteType)stype;
	// add site to m_pah->m_siteMap
	m_pah->m_siteMap[st->type].push_back(st);
	// update member C
	st->C1 = Carb1;
	st->C2 = Carb2;
	/*
    if(stype<13 && stype>9) // i.e. principal sites with 5R at both sides (eg RFER)
        stype -= 4;
	else if (stype < 10 && stype>5) //i.e. principal sites with 5R at one side (eg RFE)
		stype -= 6;
	else if (stype > 20 && stype < 24) //i.e. principal site is R5R6ZZ, R5R6AC or R5R6BY5 and are losing one carbon atom from the structure only
		stype -= 1;
	else if (stype >= 25 && stype <= 27)
		stype -= 5;
	else if (stype == 35 || stype == 36 || stype == 40)
		stype -= 5;
	else if (stype >= 62 && stype <= 65)
		stype -= 6;
    else {
        std::cout<<"ERROR: rem5RfromSite: illegal site type to remove 5R from.\n";
		saveDOT("KMC_DEBUG/Error.dot");
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
		addR5internal(C1_res,C2_res);
	}
}
int convSiteType_error_counter =0;
//! Changes site type into another site type
void PAHProcess::convSiteType(Spointer& st, Cpointer Carb1, Cpointer Carb2, kmcSiteType t) {
    // removes site from m_pah->m_siteMap (principal site)
    delSiteFromMap(st->type, st);
    // change site type
    st->type = t;
	if (!checkSiteValid(st)) {
		cout << "Invalid site convSiteType. This may be alright if the edge is unreactive but it may also be an error. \n";
		printSites(st);
		ifstream  src("KMC_DEBUG/BEFORE.xyz");
		std::string filename = "KMC_DEBUG/BEFORE_convSiteType";
		filename.append(std::to_string(convSiteType_error_counter));
		filename.append(".xyz");
		ofstream dst(filename);
		dst << src.rdbuf();
		std::string filename2 = "KMC_DEBUG/KMC_PAH_convSiteType_error_";
		filename2.append(std::to_string(convSiteType_error_counter));
		saveXYZ(filename2);
		cout<<"Saving file: "<< filename<<"\n";
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
    setCount(m_pah->m_counts.first, (int) m_pah->m_siteList.size());
    //stericHindrance();
    //cout << "Principal Sites Updated, H count: "<< m_pah->m_counts[1] << "..\n";
}

int updatesites_error_counter = 0;
//! Updates particular site
void PAHProcess::updateSites(Spointer& st, // site to be updated
                               Cpointer Carb1, Cpointer Carb2, // new C members
                               int bulkCchange) { // addition to number of bulk C in site
    // check if site type change is valid (as long as site still principal site)
	int stype = (int)st->type;
	if (!checkSiteValid(stype)){
		cout << "ERROR: updateSites: Invalid site type before update\n";
		std::ostringstream msg;
		msg << "ERROR: updateSites: Invalid site type before update\n";
		//saveDOT("KMC_DEBUG/KMC_PAH_X_UPDATE_prev.dot");
		ifstream  src("KMC_DEBUG/BEFORE.xyz");
		std::string filename = "KMC_DEBUG/BEFORE_";
		filename.append(std::to_string(updatesites_error_counter));
		filename.append(".xyz");
		//filename.append(std::to_string(this->m_pah->m_parent->ID()));
		ofstream dst(filename);
		dst << src.rdbuf();
		std::string filename2 = "KMC_DEBUG/KMC_PAH_X_UPDATE_prev_";
		filename2.append(std::to_string(updatesites_error_counter));
		//filename2.append("_");
		//filename2.append(std::to_string(this->m_pah->m_parent->ID()));
		saveXYZ(filename2);
		cout<<"Saving file: "<< filename<<"\n";
		cout<<"Saving file: "<< filename2<<".xyz\n";
		//throw std::runtime_error(msg.str());
		//assert(false);
		printSites(st);
		++updatesites_error_counter;
		return;
	}
	if ((stype + bulkCchange) < 0){
		cout << "ERROR: updateSites: Bulk C change invalid (Principal)\n";
		std::ostringstream msg;
		msg << "ERROR: Bulk C change invalid (Principal). Trying to add "
			<< bulkCchange << " bulk C to a " << kmcSiteName(st->type)
			<< " (Sweep::KMC_ARS::PAHProcess::updateSites)";
		//saveDOT("KMC_DEBUG/KMC_PAH_X_UPDATE.dot");
		ifstream  src("KMC_DEBUG/BEFORE.xyz");
		std::string filename = "KMC_DEBUG/BEFORE_";
		filename.append(std::to_string(updatesites_error_counter));
		filename.append(".xyz");
		//filename.append(std::to_string(this->m_pah->m_parent->ID()));
		ofstream dst(filename);
		dst << src.rdbuf();
		std::string filename2 = "KMC_DEBUG/KMC_PAH_X_UPDATE_";
		filename2.append(std::to_string(updatesites_error_counter));
		saveXYZ(filename2);
		cout<<"Saving file: "<< filename<<"\n";
		cout<<"Saving file: "<< filename2<<".xyz\n";
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
		removeR5internal(Carb1,Carb2);
	}
	if (stype + bulkCchange == 2102) {
		stype = 1002; bulkCchange = 0;
	}
	if (stype + bulkCchange == 2004 || stype + bulkCchange == 2014) {
		//There are two possible sites ZZACR5 and FEACR5FE. Decide which one.
		Cpointer Ccheck = st->C1->C2; Cpointer Ccheck2 = Ccheck->C2; Cpointer Ccheck3 = st->C2->C1;	Cpointer Ccheck4 = Ccheck3->C2;
		bool decide_PAH = false;
		if (isR5internal(Ccheck, Ccheck2) || isR5internal(Ccheck3, Ccheck4)) decide_PAH = true;
		if(decide_PAH == true){
			stype = 2004; bulkCchange = 0;
		}
		else {
			stype = 2014; bulkCchange = 0;
		}
	}
	if (stype + bulkCchange == 2005 || stype + bulkCchange == 2015) {
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
		//There are two possible sites ZZACR5 and FEACR5FE. Decide which one.
		Cpointer Ccheck = st->C1->C2;
		Cpointer Ccheck2 = Ccheck->C2;
		Cpointer Ccheck3 = st->C2->C1;
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
		cout << "Added" << bulkCchange << "carbons to a None site.\n";
		bulkCchange = 0;
		ifstream  src("KMC_DEBUG/BEFORE.xyz");
		std::string filename = "KMC_DEBUG/BEFORE_";
		filename.append(std::to_string(updatesites_error_counter));
		filename.append(".xyz");
		ofstream dst(filename);
		dst << src.rdbuf();
		std::string filename2 = "KMC_DEBUG/KMC_PAH_X_UPDATE_None_site_";
		filename2.append(std::to_string(updatesites_error_counter));
		//filename2.append(std::to_string(this->m_pah->m_parent->ID()));
		saveXYZ(filename2);
		cout<<"Saving file: "<< filename<<"\n";
		cout<<"Saving file: "<< filename2<<".xyz\n";
		//throw std::runtime_error(msg.str());
		//assert(false);
		printSites(st);
		++updatesites_error_counter;
	}
	if (!checkSiteValid(stype + bulkCchange)){
		cout << "ERROR: updateSites: Created undefined site:\n";
		cout << "Type = " << stype << "\n";
		cout << "BulkC change = " << bulkCchange << "\n";
		/*cout << "ERROR: updateSites: Bulk C + Type gave invalid site\n";
		std::ostringstream msg;
		msg << "ERROR: Bulk C + Type gave invalid site. Trying to add "
			<< bulkCchange << " bulk C to a " << kmcSiteName(st->type)
			<< " (Sweep::KMC_ARS::PAHProcess::updateSites)";*/
		//saveDOT("KMC_DEBUG/KMC_PAH_X_UPDATE.dot");
		ifstream  src("KMC_DEBUG/BEFORE.xyz");
		std::string filename = "KMC_DEBUG/BEFORE_";
		filename.append(std::to_string(updatesites_error_counter));
		filename.append(".xyz");
		//filename.append(std::to_string(this->m_pah->m_parent->ID()));
		ofstream dst(filename);
		dst << src.rdbuf();
		std::string filename2 = "KMC_DEBUG/KMC_PAH_X_UPDATE_BulkplusC_";
		filename2.append(std::to_string(updatesites_error_counter));
		//filename2.append("_");
		//filename2.append(std::to_string(this->m_pah->m_parent->ID()));
		saveXYZ(filename2);
		cout<<"Saving file: "<< filename<<"\n";
		cout<<"Saving file: "<< filename2<<".xyz\n";
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
        // Check for FE3 (if there's FE on each side of the FE)
        if(moveIt(st,1)->type == FE && moveIt(st,-1)->type == FE) {
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
        else if(moveIt(st,1)->type == FE || moveIt(st,-1)->type == FE){
            Spointer S1,S2;
            S1 = moveIt(st,-1); S2 = moveIt(st, 1);
            // Check if that FE is not a FE3
            if(S2->type == FE && moveIt(S2,1)->type != FE) {
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
            } else if(S1->type == FE && moveIt(S1,-1)->type != FE) {
                st->comb = FE2;
                m_pah->m_siteMap[FE2].push_back(st);
                if(S1->comb == FE2) delSiteFromMap(S1->comb, st);
                if(S1->comb != FE2) updateCombinedSites(S1);
            } else
                st->comb = None;
            break;
        }
        // Check for FE_HACA
        else if(moveIt(st,1)->type != FE && moveIt(st,-1)->type != FE && moveIt(st,1)->type != RFE && moveIt(st,-1)->type != RFE) {
            //if(st->C1->C1->bridge || st->C2->C2->bridge)
            st->comb = FE_HACA;
            m_pah->m_siteMap[FE_HACA].push_back(st);
            break;
        }
        else st->comb = None;
        break;
    case AC:
        // Check for AC_FE3
        if((moveIt(st,2)->comb==FE3 || moveIt(st,-2)->comb==FE3) && (moveIt(st,3)->comb!=FE3 || moveIt(st,-3)->comb!=FE3) && !st->C1->C2->bridge) {
            st->comb = AC_FE3;
            m_pah->m_siteMap[AC_FE3].push_back(st);
            break;
        }else st->comb = None;
        break;
    case BY5:
    // Check for BY5_FE3
        if((moveIt(st,2)->comb==FE3 || moveIt(st,-2)->comb==FE3) && (moveIt(st,3)->comb!=FE3 || moveIt(st,-3)->comb!=FE3) && !st->C1->C2->bridge) {
            st->comb = BY5_FE3;
            m_pah->m_siteMap[BY5_FE3].push_back(st);
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
	int rings_Embedded;
	std::list <cpair> IntCarbons;
    // Structure for Benzene
	std::string BENZENE_Sites = "FE,FE,FE,FE,FE,FE";
	auto BENZENE_Rings = std::make_tuple (1, 0, 0) ;
	intpair BENZENE_CH(6, 6);
	int BENZENE_RINGS_EMBEDDED = 0;
	// Structure for Naphthalene
	std::string NAPHTHALENE_Sites = "FE,FE,FE,ZZ,FE,FE,FE,ZZ";
	auto NAPHTHALENE_Rings = std::make_tuple(2, 0, 0);
	intpair NAPHTHALENE_CH(10, 0);
	int NAPHTHALENE_RINGS_EMBEDDED = 0;
	// Structure for Pyrene
	std::string PYRENE_Sites = "ZZ,FE,FE,ZZ,FE,ZZ,FE,FE,ZZ,FE";
	auto PYRENE_Rings = std::make_tuple (4, 0, 0);
	intpair PYRENE_CH(16, 10);
	int PYRENE_RINGS_EMBEDDED = 0;
	std::list <cpair> PYRENE_intCarbons;
	PYRENE_intCarbons.push_back(std::make_tuple(1.4*cos(-60.0*M_PI/180.0), 1.4*sin(-60.0*M_PI/180.0), 0.0));
	PYRENE_intCarbons.push_back(std::make_tuple(1.4*cos(-60.0*M_PI/180.0) + 1.4, 1.4*sin(-60.0*M_PI/180.0), 0.0));
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
	int TEST_RINGS_EMBEDDED = 0;
	//Internal Carbons missing
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
		IntCarbons = PYRENE_intCarbons;
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
	return initialise(chosen, std::get<0>(rings), std::get<1>(rings), rings_Embedded, std::get<2>(rings), rings_Embedded, IntCarbons);
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
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();
		//Set Internal Carbons list
		P_intCarbons.push_back(std::make_tuple(0.0, -2*1.4*cos(30.0*M_PI/180.0), 0.0));
		P_intCarbons.push_back(std::make_tuple(1.4, -2*1.4*cos(30.0*M_PI/180.0), 0.0));
		m_pah->m_InternalCarbons = P_intCarbons;
		//cout << "Pyrene Initialised!\n";
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
        // update all sites and combined sites
        updateSites();
        updateCombinedSites();
        //cout << "Benzopyrene Initialised!\n";
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
PAHStructure& PAHProcess::initialise(std::string siteList_str, int R6_num, int R5_num_Lone, int R5_num_Embedded, int R7_num_Lone, int R7_num_Embedded, std::list<cpair> internalCarbons){
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
                throw std::runtime_error(msg.str());
                assert(false);
        }
        siteList_vec.push_back(temp);
    }
	createPAH(siteList_vec, R6_num, R5_num_Lone, R5_num_Embedded, R7_num_Lone, R7_num_Embedded, internalCarbons);
    return *m_pah;
}

// Create Structure from vector of site types
void PAHProcess::createPAH(std::vector<kmcSiteType>& vec, int R6, int R5_Lone, int R5_Embedded, int R7_Lone, int R7_Embedded, std::list<cpair> inCarbs) {
    // current C, bondangle and coordinates
    Cpointer newC=addC();
    m_pah->m_cfirst = newC;
    m_pah->m_clast = NULLC;
    // number of bulk C to be added
    int bulkC;
    // type of site; if type 0, basic site types (FE - BY6); if type 1, R5 and basic sites with
    // a R5 at one side (RFE - RBY5); if type 2, basic sites wit R5 at each side (RFER - RACR)
    unsigned short int site_t;
    // start drawing..
    for(size_t i=0; i<vec.size(); i++) {
        Cpointer S_C1 = newC;
        // get number of bulk C to be added and site type
        if((int)vec[i] <= 4) {
            bulkC = (int) vec[i]; site_t = 0;
        }else if((int)vec[i] >= 5 && (int)vec[i] <= 9) {
            bulkC = (int) vec[i] - 5; site_t = 1;
        }else if((int)vec[i] >= 10 && (int)vec[i] <= 12) {
            bulkC = (int) vec[i] - 8; site_t = 2;
		}
		/**
		* If this condition is true, vec[i] is the site tye ACR5.
		* It is essentially a special type of armchair site (a basic site type and associated with two bulk carbon atoms).
		*/
		else if ((int)vec[i] == 30) {
			bulkC = (int)vec[i] - 16; site_t = 0;

		}
		else {
            cout << "createPAH: Combined site types in list of sites. Please use only\n"
                << "principal site types\n";
            std::ostringstream msg;
            msg << "ERROR: Combined site types found in list of sites."
                << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
                throw std::runtime_error(msg.str());
                assert(false);
            return;
        }
        kmcSiteType prevType;
        if(i==0) prevType = vec.back();
        else prevType = vec[i-1];
        switch(site_t) {
        case 0:
            newC = drawType0Site(newC, bulkC); break;
        case 1:
            newC = drawType1Site(newC, bulkC, prevType); break;
        case 2:
            newC = drawType2Site(newC, bulkC); break;
        default:
            cout << "createPAH: Invalid site_t number...\n";
            std::ostringstream msg;
            msg << "ERROR: invalid site classification."
                << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
                throw std::runtime_error(msg.str());
                assert(false);
            return;
        }
        addSite(vec[i], S_C1, newC);
    }
    // check if PAH closes correctly
	if (m_pah->m_clast == NULLC || checkHindrance_twoC(newC, m_pah->m_cfirst)) {
        // PAH did not close properly. invalid structure
        cout << "createPAH: PAH did not close properly. Could be problem "
            <<"with site list input...\n";
        std::ostringstream msg;
        msg << "ERROR: PAH did not close properly.."
            << " (Sweep::KMC_ARS::PAHProcess::createPAH)";
        //saveDOT("KMC_DEBUG/KMC_PAH_X_CLOSE.dot");
        throw std::runtime_error(msg.str());
        assert(false);
        return;
    }
    m_pah->m_rings = R6;
	m_pah->m_rings5_Lone = R5_Lone;
	m_pah->m_rings5_Embedded = 0;
	m_pah->m_rings7_Lone = 0;
	m_pah->m_rings7_Embedded = 0;
	m_pah->m_InternalCarbons = inCarbs;
    for(Ccontainer::iterator i = m_pah->m_carbonList.begin();
        i != m_pah->m_carbonList.end(); i++)
        updateA(*i, 'H');
	//int totalC_num = 2 * m_pah->m_rings + (CarbonListSize() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded) / 2 + numberOfBridges() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded + 1;
	int totalC_num = 2 * m_pah->m_rings + (CarbonListSize() + 3 * m_pah->m_rings5_Lone + 3 * m_pah->m_rings5_Embedded + 5 * m_pah->m_rings7_Lone + 5 * m_pah->m_rings7_Embedded) / 2 + numberOfBridges() + 1;
	m_pah->setnumofC(totalC_num);
    m_pah->setnumofH((int)vec.size());
    updateCombinedSites();
}

//! For createPAH function: drawing type 0 sites
Cpointer PAHProcess::drawType0Site(Cpointer Cnow, int bulkC) {
    // draw site
    angletype angle = normAngle(Cnow->bondAngle1-60);
    for(int c=0; c<=bulkC; ++c) {
        // check if adding on existing C atom (bridge)
        cpair pos = jumpToPos(Cnow->coords, angle, 0, 1.4);
        if(m_pah->m_cpositions.count(pos)) {
            // this coordinate is filled
            Cpointer Cpos = findC(pos);
            if(Cpos != m_pah->m_cfirst) { // it is a bridged C atom
                Cpointer Cbridge = Cpos->C1;
                Cnow->bondAngle1 = angle;
                Cpos->bridge = true; Cbridge->bridge = true;
                Cpos->C3 = Cbridge; Cbridge->C3 = Cpos;
                connectToC(Cnow, Cpos);
                Cbridge->C2 = NULLC;
                angle = normAngle(angle + 120);
                Cbridge->bondAngle1 = angle;
                Cnow = Cbridge;
                --bulkC;
            }else { // reached end of PAH
                Cnow->bondAngle1 = angle;
                connectToC(Cnow, Cpos);
                m_pah->m_clast = Cnow;
                return Cpos;
            }
        }else {
            Cnow = addC(Cnow, angle, angle, 1.4, false);
            angle = normAngle(angle+60);
        }
    }
    return Cnow;
}

//! For createPAH function: drawing type 1 sites
Cpointer PAHProcess::drawType1Site(Cpointer Cnow, int bulkC, kmcSiteType prevType) {
    //draw R5 site
    angletype angle = Cnow->bondAngle1-60;
    if(bulkC == 0) {
        angle = normAngle(angle-30);
        cpair pos = jumpToPos(Cnow->coords, angle, 0, 1.4);
        if(m_pah->m_cpositions.count(pos)) { // reached end of PAH
            Cpointer Cpos = findC(pos);
            Cnow->bondAngle1 = angle;
            connectToC(Cnow, Cpos);
            m_pah->m_clast = Cnow;
            return Cpos; // m_cfirst
        }
        else Cnow = addC(Cnow, angle, normAngle(angle-30), 1.4, false);
        return Cnow;
    }else { //draw RXX site
            return drawType0Site(Cnow, bulkC);
    }
}

//! For createPAH function: drawing type 2 sites
Cpointer PAHProcess::drawType2Site(Cpointer Cnow, int bulkC) {
    //angletype angle = Cnow->bondAngle1;
    //Cnow->bondAngle1 = normAngle(angle-60);
    return drawType0Site(Cnow, bulkC);
}

//! Finds C atom with specific coordinates
Cpointer PAHProcess::findC(cpair coordinates) {
	double tol = 1e-1;
    for(Ccontainer::iterator i=m_pah->m_carbonList.begin();
        i != m_pah->m_carbonList.end(); ++i) {
		if (abs(std::get<0>((*i)->coords) - std::get<0>(coordinates)) < tol && abs(std::get<1>((*i)->coords) - std::get<1>(coordinates)) < tol && abs(std::get<2>((*i)->coords) - std::get<2>(coordinates)) < tol){
			return (*i);
		}
    }
    return NULLC;
}

//! Find Site to the other side of a bridge
Spointer PAHProcess::findSite(Cpointer C_1) {
	Spointer temp;
	Cpointer Carb1, Carb2;
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
}

//! Find Site to the other side of a bridge
Cpointer PAHProcess::findThirdC(Cpointer C_1) {
	for(Ccontainer::iterator i=m_pah->m_carbonList.begin(); i != m_pah->m_carbonList.end(); ++i) {
		Cpointer C_check = *i;
		if (C_check != C_1 && C_check != C_1->C1 && C_check != C_1->C1->C1 && C_check != C_1->C2 && C_check != C_1->C2->C2){
			double dist = getDistance_twoC(C_check, C_1);
			if (dist < 1.8) return C_check;
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
		case R5ZZACR5: return true;
		case R5ACR5R5: return true;
		case ACR5RFER: return true;
		case RAC_FE3: return true;
		case None: return true;
		case Inv: return true;
		case any: return true;
		case benz: return true;
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
//! Structure processes: returns success or failure
bool PAHProcess::performProcess(const JumpProcess& jp, rng_type &rng, int PAH_ID)
{
    //printStruct();
    //cout << "Start Performing Process..\n";
    kmcSiteType stp = jp.getSiteType();
	Spointer site_perf;
    int id = jp.getID();
    
    // choose random site of type stp to perform process
	if (id == 15 || id == 14){
		site_perf = RZZchooseRandomSite(stp, rng); //cout<<"[random site chosen..]\n";
		//site_perf = chooseRandomSite(stp, rng); //cout<<"[random site chosen..]\n";
	}
	else {
		site_perf = chooseRandomSite(stp, rng); //cout<<"[random site chosen..]\n";
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
	saveXYZ("KMC_DEBUG/BEFORE");
	//Copy site list before performing process
	std::list<std::string> Sitelist_before = copySites(site_perf);
	///////
    switch(id) {
        case 1:
            proc_G6R_AC(site_perf, site_C1, site_C2); break;
        case 2:
            proc_G6R_FE(site_perf, site_C1, site_C2); break;
        case 3:
            //dotname << "KMC_DEBUG/" << site_perf->C1 << "_1.dot";
            //dotname2 << "KMC_DEBUG/" << site_perf->C1 << "_2.dot";
            //saveDOT(dotname.str());
            proc_L6_BY6(site_perf, site_C1, site_C2);
            //saveDOT(dotname2.str());
            break;
        case 4:
            proc_PH_benz(site_perf, site_C1, site_C2, rng); break;
        case 5:
            proc_D6R_FE3(site_perf, site_C1, site_C2); break;
        case 6:
            proc_O6R_FE3_O2(site_perf, site_C1, site_C2); break;
        case 7:
            proc_O6R_FE3_OH(site_perf, site_C1, site_C2); break;
        case 8:
            proc_O6R_FE_HACA_O2(site_perf, site_C1, site_C2); break;
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
            proc_O6R_FE2(site_perf, site_C1, site_C2);
            break;
        case 21:
            proc_O6R_FE2(site_perf, site_C1, site_C2); break;
		case 22:
			proc_D6R_FE_AC(site_perf, site_C1, site_C2); break;
		case 23:
			proc_B6R_ACR5(site_perf, site_C1, site_C2); break;
		case 24:
			proc_M5R_ACR5_ZZ(site_perf, site_C1, site_C2, rng); break;
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
			proc_MR5_R6(site_perf, site_C1, site_C2, rng); break;
		case 35:
			proc_GR7_R5R6AC(site_perf, site_C1, site_C2); break;
		case 36:
			proc_GR7_FEACR5(site_perf, site_C1, site_C2); break;
		case 37:
			proc_G6R_R5R6ZZ(site_perf, site_C1, site_C2); break;
		case 38:
			proc_L7_ACACR5(site_perf, site_C1, site_C2); break;
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
			<< "ID" << id << " on PAH ID: " << PAH_ID << "...\n"
			<< "*************\nAfter performing process --\n";
		printBeforeSites(Sitelist_before);
		printSites(site_perf);
		ifstream  src("KMC_DEBUG/BEFORE.xyz");
		std::string filename = "KMC_DEBUG/BEFORE_performProcess_";
		filename.append(std::to_string(perform_process_error_counter));
		filename.append(".xyz");
		ofstream dst(filename);
		dst << src.rdbuf();
		std::string filename2 = "KMC_DEBUG/KMC_PAH_performProcess_checkSiteValid_";
		filename2.append(std::to_string(perform_process_error_counter));
		saveXYZ(filename2);
		std::ostringstream msg;
		msg << "ERROR: Structure produced invalid combined site type after performing process "
			<< "ID" << id << " on PAH ID: " << PAH_ID << "..."
			<< " (Sweep::KMC_ARS::PAHProcess::performProcess)";
		//throw std::runtime_error(msg.str());
		//assert(false);
		//abort();
		cout<<"Saving file: "<< filename<<"\n";
		cout<<"Saving file: "<< filename2<<".xyz\n";
		++perform_process_error_counter;
	}
    if(!checkCombinedSiteType(site_perf) || !checkCombinedSiteType(S1)
        || !checkCombinedSiteType(S2) || !checkCombinedSiteType(S3)
        || !checkCombinedSiteType(S4)) {
        std::cout<<"ERROR. Structure produced invalid combined site type after performing process "
            << "ID"<<id<<" on PAH ID: "<< PAH_ID <<"...\n"
            <<"*************\nAfter performing process --\n";
        printBeforeSites(Sitelist_before);
		printSites(site_perf);
		ifstream  src("KMC_DEBUG/BEFORE.xyz");
		std::string filename = "KMC_DEBUG/BEFORE_performProcess_";
		filename.append(std::to_string(perform_process_error_counter));
		filename.append(".xyz");
		ofstream dst(filename);
		dst << src.rdbuf();
		std::string filename2 = "KMC_DEBUG/KMC_PAH_performProcess_checkCombinedSite_";
		filename2.append(std::to_string(perform_process_error_counter));
		saveXYZ(filename2);
        std::ostringstream msg;
        msg << "ERROR: Structure produced invalid combined site type after performing process "
            << "ID"<<id<<" on PAH ID: "<< PAH_ID <<"..."
            << " (Sweep::KMC_ARS::PAHProcess::performProcess)";
		cout<<"Saving file: "<< filename<<"\n";
		cout<<"Saving file: "<< filename2<<".xyz\n";
		++perform_process_error_counter;	
        //throw std::runtime_error(msg.str());
        //assert(false);
        //abort();
    }
	//int calc_total = 2 * m_pah->m_rings + (CarbonListSize() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded) / 2 + numberOfBridges() + m_pah->m_rings5_Lone + m_pah->m_rings5_Embedded + 1;
	int calc_total = 2 * m_pah->m_rings + (CarbonListSize() + 3 * m_pah->m_rings5_Lone + 3 * m_pah->m_rings5_Embedded + 5 * m_pah->m_rings7_Lone + 5 * m_pah->m_rings7_Embedded) / 2 + numberOfBridges() + 1;
	if(calc_total != getCHCount().first) {
        //saveDOT("KMC_DEBUG/KMC_C_Counts_ERROR.dot");
		ifstream  src("KMC_DEBUG/BEFORE.xyz");
		std::string filename = "KMC_DEBUG/BEFORE_performProcess_";
		filename.append(std::to_string(perform_process_error_counter));
		filename.append(".xyz");
		ofstream dst(filename);
		dst << src.rdbuf();
		std::string filename2 = "KMC_DEBUG/KMC_PAH_performProcess_C_Count_";
		filename2.append(std::to_string(perform_process_error_counter));
		saveXYZ(filename2);
        cout<<"ERROR. Calculated total did not tally with double C counts!\n";
        cout<<"Last performed process: "<<jp.getName()<<", ID = "<<jp.getID()<<"PAH ID = "<<PAH_ID<<'\n';
		cout<<"R6s = "<<m_pah->m_rings<<", R5s = "<<m_pah->m_rings5_Lone<<" alone + "<<m_pah->m_rings5_Embedded<<" embedded, R7s = "<<m_pah->m_rings7_Embedded<<"\n";
        std::ostringstream msg;
        msg << "\nCalculated total: "<<calc_total<<'\n'
            << "Real total: " << getCHCount().first << '\n';
        printBeforeSites(Sitelist_before);
		printSites(site_perf);
		//throw std::runtime_error(msg.str());
		cout<<"Saving file: "<< filename<<"\n";
		cout<<"Saving file: "<< filename2<<".xyz\n";
		++perform_process_error_counter;
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
	if (double dist = getDistance_twoC(C_1,C_2) < 2.6 ) {
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
    Cpointer newC1;
    Cpointer newC2;
    if(checkHindrance_newC(C_1) || checkHindrance_newC(C_2)) {
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
    Spointer S1 = moveIt(stt, -1); 
    Spointer S2 = moveIt(stt, 1);
	Spointer S3 = moveIt(S1, -1);
	Spointer S4 = moveIt(S2, 1);
    // Update Site and neighbours
	if ( stt->type == RFER ){
		convSiteType(stt, newC1, newC2, FE);
		convSiteType(S1, S1->C1, newC1, R5R6); // neighbours
		convSiteType(S2, newC2, S2->C2, R5R6);
		updateSites(S3, S3->C1, S3->C2, +400); // neighbours of neighbours
		updateSites(S4, S4->C1, S4->C2, +400);
	}
	else if ( stt->type == R5R6FER ){
		convSiteType(stt, newC1, newC2, FE);
		if (S1->type==R5) {
			convSiteType(S1, S1->C1, newC1, R5R6); // neighbours
			updateSites(S2, newC2, S2->C2, +1501);
		}
		else if (S2->type==R5){
			convSiteType(S2, newC2, S2->C2, R5R6);
			updateSites(S1, S1->C1, newC1, +1501);
		}
		else {
			cout<<"R5R6FER not next to R5. Error.\n";
			saveXYZ("R5R6FER_error");
		}
	}
	else if ( stt->type == R5R6FER5R6 ){
		convSiteType(stt, newC1, newC2, FE);
		updateSites(S2, newC2, S2->C2, +1501);
		updateSites(S1, S1->C1, newC1, +1501);
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
	if (double dist = getDistance_twoC(newC2,C_2) > 1.7 ) {
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
    if(checkHindrance_newC(C_1) || checkHindrance_newC(C_2)) {
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
	if (double dist = getDistance_twoC(newC4, C_2) > 1.65){
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
    Cpointer now = C_1->C2;
    do{
        Cpointer next;
        if(!now->bridge) {
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
	int ntype_site = (int)stt->type;
	if (ntype_site == 2005) convSiteType(stt, stt->C1, stt->C2, (kmcSiteType)4);
    int ntype1 = (int) moveIt(stt, -1)->type;
    int ntype2 = (int) moveIt(stt, 1)->type;
	int newType;
    if(ntype1 < 5 && ntype2 < 5) {
        newType = (ntype1+ntype2+2);
        // convert site
        if(newType>4) {
            //saveDOT(std::string("BY6ClosureProblem.dot"));
            std::cerr<<"ERROR: newType is > 4 (PAHProcess::proc_L6_BY6)\n";
        }
        updateSites(stt, moveIt(stt,-1)->C1, moveIt(stt,1)->C2,(newType-4));
    }
    else {
		int new_point = 0;
		if (ntype_site == 4) newType = (ntype1 + ntype2 + 2);
		else if (ntype_site == 104) {
			new_point = 502;
			if (ntype1 == 100) {
				newType = (new_point + ntype2);
				Spointer S1 = moveIt(stt, -2);
				int stype1 = (int)S1->type + 400;
				convSiteType(S1, S1->C1, S1->C2, (kmcSiteType)stype1);
			}
			else {
				newType = (new_point + ntype1);
				Spointer S2 = moveIt(stt, 2);
				int stype2 = (int)S2->type + 400;
				convSiteType(S2, S2->C1, S2->C2, (kmcSiteType)stype2);
			}
		}
		else if (ntype_site == 204) {
			newType = 1002;
			Spointer S1 = moveIt(stt, -2); Spointer S2 = moveIt(stt, 2);
			int stype1 = (int)S1->type + 400; int stype2 = (int)S2->type + 400;
			convSiteType(S1, S1->C1, S1->C2, (kmcSiteType)stype1);
			convSiteType(S2, S2->C1, S2->C2, (kmcSiteType)stype2);
			}
		else if (ntype_site == 504) {
			new_point = 2003;
			if (ntype1 >= 501 && ntype1 <= 504) ntype1 -= 501;
			else if (ntype2 >= 501 && ntype1 <= 504) ntype2 -= 501;
			newType = (new_point + ntype1 + ntype2);
		}
		else if (ntype_site == 604) {
			new_point = 2103;
			if (ntype1 == 100) {
				ntype1 = 0;
				ntype2 -= 501;
			}
			else {
				ntype2 = 0;
				ntype1 -= 501;
			}
			newType = (new_point + ntype1 + ntype2);
		}
		else if (ntype_site == 504) {
			new_point = 2002;
			ntype1 = ntype1 % 10;
			ntype2 = ntype2 % 10;
			newType = (new_point + ntype1 + ntype2);
			/*if (ntype1 >= 501 && ntype1 <= 504) ntype1 -= 501;
			else if (ntype2 >= 501 && ntype1 <= 504) ntype2 -= 501;*/
		}
		else if (ntype_site == 2014 || ntype_site == 2004) {
			new_point = 2;
			ntype1 = ntype1 % 10;
			ntype2 = ntype2 % 10;
			newType = (new_point + ntype1 + ntype2);
			/*if (ntype1 >= 501 && ntype1 <= 504) ntype1 -= 501;
			else if (ntype2 >= 501 && ntype1 <= 504) ntype2 -= 501;*/
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
	if (double dist = getDistance_twoC(C_1,C_2) > 1.6 ) {
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
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
	cpair ivec1 = invert_vector(vec1);
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
	if ( (int)checkR5_1->type == 101 || (int)checkR5_1->type == 501 ) return;
	if ( (int)checkR5_2->type == 101 || (int)checkR5_2->type == 501 ) return;
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
	if (S1->type >= 2002 && S1->type <= 2104 && S2->type >= 2002 && S2->type <= 2104){
		//Needs debugging!!! DONE?
		removeR5internal(C1_new, C2_new);
		convSiteType(stt, C1_new->C1, C2_new->C2, ZZ);
		updateSites(S1, S1->C1, C1_new->C1, -2002);
		updateSites(S2, C2_new->C2, S2->C2, -2002);
		m_pah->m_rings5_Embedded--;
		redrawR5(stt, C1_new->C1, C2_new->C2);
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
// ID8- R6 oxidation at AC by O2 (AR12 in Matlab)
// ************************************************************
void PAHProcess::proc_O6R_FE_HACA_O2(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSites(stt);
    // member C atoms of resulting AC site
    Cpointer C1_res, C2_res;
	//printStruct(C_1);
    C1_res = C_1->C1;
    C2_res = C_2->C2;
	cpair Cdir = get_vector(C_2->coords,C_2->C2->coords);
	cpair FEdir = get_vector(C_1->coords,C_2->coords);
	cpair Hdir1 = C_2->growth_vector;
	cpair Hdir2 = C_1->growth_vector;
    // check if process will result in a bridge
	cpair pos = jumpToPos(C1_res->coords, Cdir, 1.4);
    //cpair pos = jumpToPos(C1_res->coords, normAngle(C1_res->bondAngle1-120), 0, 1.4);
	bool bridge = checkHindrance_C_PAH((pos));
	if (bridge) {
		//cout<<"Oxidation generated a bridge\n";
		//saveXYZ("KMC_DEBUG/Oxidation_to_bridge_before");
	}
	bool hept_bool = false;
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
            throw std::runtime_error(msg.str());
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
		Cpointer Cnew2 = addC(C1_res->C2, FEdir, 1.4, true);
		updateA(C2_res,'H', Hdir2);
		//saveXYZ("KMC_DEBUG/oxidation_issue");
		double dist = getDistance_twoC(C1_res,C2_res);
		if (dist > 3.2 && m_pah->m_rings7_Embedded >= 1){
			//Assume that this site was in an heptagon.
			removeR7internal(C1_res, C2_res);
			hept_bool = true;
			Cpointer Cnew3 = addC(Cnew2, FEdir, 1.4, true);
		}
		//cout<<"Distance = "<<dist<<", R7s = "<<m_pah->m_rings7_Embedded<<"\n";
		//Cpointer Cnew1 = addC(C1_res, normAngle(C1_res->bondAngle1 - 120), 0, 1.4, true);
		//Cpointer Cnew2 = addC(C1_res->C2, normAngle(C1_res->bondAngle1 + 60), normAngle(C1_res->bondAngle1 + 120), 1.4, true);
		// update H
		//updateA(C1_res, C2_res, 'H');
    }
    // update sites and neighbours
    Spointer S1, S2, S3, S4;
    S1 = moveIt(stt,-1); S2 = moveIt(stt,1);
	if (!hept_bool){
		if (S1->type == R5R6 && (S2->type == R5R6)) {
			S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
			updateSites(stt, C1_res, C2_res, 202);
			updateSites(S1, S1->C1, C1_res, -1);
			updateSites(S2, C2_res, S2->C2, -1);
			if ((int)S3->type > 500 && (int)S3->type < 2000) updateSites(S3, S3->C1, S3->C2, -400);
			if ((int)S4->type > 500 && (int)S4->type < 2000) updateSites(S4, S4->C1, S4->C2, -400);
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}	
		else if (S1->type == R5R6 || (S2->type == R5R6)) {
			if (S1->type == R5R6){
				S3 = moveIt(S1, -1);
				updateSites(stt, C1_res, C2_res, 102);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
				if ((int)S3->type <2000)  updateSites(S3, S3->C1, S3->C2, -400);
			}
			else {
				S4 = moveIt(S2, 1);
				updateSites(stt, C1_res, C2_res, 102);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
				if ((int)S4->type <2000) updateSites(S4, S4->C1, S4->C2, -400);
			}
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		else if (S1->type == R5ACR5 || (S2->type == R5ACR5)) {
			if (S1->type == R5ACR5 && S2->type == R5ACR5){
				updateSites(stt, C1_res, C2_res, 1002);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else if (S1->type == R5ACR5){
				updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else {
				S4 = moveIt(S2, 1);
				updateSites(stt, C1_res, C2_res, 102);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
				updateSites(S4, S4->C1, S4->C2, -400);
			}
			OpenBabel::OBMol mol = passPAH();
			mol = optimisePAH(mol);
			passbackPAH(mol);
		}
		else if (S1->type == ACR5 && S2->type == ACR5) {
			updateSites(stt, C1_res, C2_res, 1003);
			updateSites(S1, S1->C1, C1_res, -1501);
			updateSites(S2, C2_res, S2->C2, -1501);
		}
		else if (S1->type == ACR5 || S2->type == ACR5) {
			updateSites(stt, C1_res, C2_res, 502);
			updateSites(S1, S1->C1, C1_res, -1);
			updateSites(S2, C2_res, S2->C2, -1);
		}
		else if (((int)S1->type >= 2003 && (int)S1->type <= 2005) && ((int)S2->type >= 2003 && (int)S2->type <= 2005)){
			cpair R5check1 = endposR5internal(C1_res, C1_res->C1);
			cpair R5check2 = endposR5internal(C2_res, C2_res->C2);
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
			cpair R5check1 = endposR5internal(C1_res, C1_res->C1);
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
			cpair R5check1 = endposR5internal(C1_res, C1_res->C1);
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
				m_pah->m_rings5_Embedded--;*/
			}
		}
		else {
			updateSites(stt, C1_res, C2_res, 2);
			updateSites(S1, S1->C1, C1_res, -1);
			updateSites(S2, C2_res, S2->C2, -1);
		}
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
		if (S1->type == ACR5 && S2->type == ACR5) {
			updateSites(stt, C1_res, C2_res, 2204);
			updateSites(S1, S1->C1, C1_res, -1501);
			updateSites(S2, C2_res, S2->C2, -1501);
		}
		else if (S1->type == ACR5 || S2->type == ACR5) {
			updateSites(stt, C1_res, C2_res, 2103);
			updateSites(S1, S1->C1, C1_res, -1);
			updateSites(S2, C2_res, S2->C2, -1);
		}
		else if (((int)S1->type >= 2003 && (int)S1->type <= 2005) && ((int)S2->type >= 2003 && (int)S2->type <= 2005)){
			double dist1 = getDistance_twoC(C1_res, C1_res->C1);
			double dist2 = getDistance_twoC(C2_res, C2_res->C2);
			if (dist1 > 1.1 && dist2 > 1.1){
				updateSites(stt, C1_res, C2_res, 2002);
				updateSites(S1, S1->C1, C1_res, -1);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else if (dist1 > 1.1){
				updateSites(stt, C1_res, C2_res, 502);
				updateSites(S1, S1->C1, C1_res, -1501);
				updateSites(S2, C2_res, S2->C2, -1);
			}
			else if (dist2 > 1.1){
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
			double dist2 = getDistance_twoC(C1_res, C1_res->C1);
			if (dist2 > 1.1){
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
			double dist2 = getDistance_twoC(C2_res, C2_res->C2);
			if (dist2 > 1.1){
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
			updateSites(stt, C1_res, C2_res, 2003);
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
    proc_O6R_FE_HACA_O2(stt, C_1, C_2);
}
// ************************************************************
// ID10- R5 growth on ZZ (AR3 in Matlab)
// ************************************************************
void PAHProcess::proc_G5R_ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2) {
    //printSites(stt);
	if (double dist = getDistance_twoC(C_1,C_2) < 2.1 ) {
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
    if(checkHindrance_newC(C_1) || checkHindrance_newC(C_2)) {
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
	addR5internal(C1_res,C2_res);
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
    removeC(C_1, false);
    removeC(C_2, false);
	double Cdist = getDistance_twoC(C1_res, C2_res);
	Cpointer newC = addC(C1_res, Cdir, Cdist/2.35*1.4, true);
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
    if(checkHindrance_newC(C_1) || checkHindrance_newC(C_2)) {
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
	if (b4 && (int)moveIt(stt, -4)->type >= 2000 && (int)moveIt(stt, -4)->type <= 2204) return; // Two pentagons will collide if this is the case
	if (b4 && (int)moveIt(stt, +4)->type == 100) return;
	if (b4 && (int)moveIt(stt, +4)->type >= 500 && (int)moveIt(stt, +4)->type <= 504) return; // Two pentagons will collide if this is the case
	if (b4 && (int)moveIt(stt, +4)->type >= 2000 && (int)moveIt(stt, +4)->type <= 2204) return; // Two pentagons will collide if this is the case
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
	OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol);
	passbackPAH(mol);
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
    if(checkHindrance_newC(C_1) || checkHindrance_newC(C_2)) {
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
		starting_direction = C_1->growth_vector;
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
		C3_new = addC(C2_new, invert_vector(Cdir), dist);
		updateA(C3_new, 'H', opp_vec);
		C4_new = addC(C3_new, invert_vector(starting_direction), dist);
		updateA(C4_new, 'H', invert_vector(Cdir));
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
	cpair starting_direction, Hdir1, Hdir2, FEdir, Hdir;
    if(b4) {
        sR5 = moveIt(stt, -1);
        C_AC = sR5->C1->C1;
		starting_direction = get_vector(C_2->C1->coords, C_2->coords);
		Hdir1 = get_vector(C_1->C2->C2->coords, C_1->C2->coords);
		FEdir = get_vector(C_1->C2->coords, C_2->coords);
		Hdir2 = get_vector(C_2->C1->coords, C_2->coords);
		Hdir = starting_direction;
    }else {
        sR5 = moveIt(stt, 1);
        C_AC = sR5->C2->C2;
		starting_direction = get_vector(C_1->C2->coords, C_1->coords);
		Hdir1 = get_vector(C_1->C2->coords, C_1->coords);
		FEdir = get_vector(C_1->coords, C_2->C1->coords);
		Hdir2 = get_vector(C_2->C1->C1->coords, C_2->C1->coords);
		Hdir = Hdir1;
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
	cpair CZZdir = invert_vector(Hdir1);
    for(int i=0; i!=2; i++) removeC(Cstart->C2, false);
	addC(Cstart, CZZdir, 1.4, true);
	updateA(C_AC, 'H', Hdir);
    //addC(Cstart, normAngle(Cstart->bondAngle1-120), normAngle(Cstart->bondAngle1-60), 1.4, true);
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
    if(checkHindrance_newC(C_1) || checkHindrance_newC(C_2)) {
        /*cout<<"Site hindered, process not performed.\n"*/ return;}
    
    // check if R5 is before or after RZZ
    bool b4 = false;
    if(moveIt(stt,-1)->type == R5) b4 = true;
    // R5 will be removed and replaced by ZZ. first identify the C member of the
    // resulting RZZ site, C_RZZ. Identify where the R5 site is too.
    Spointer sR5;
	Cpointer C_RZZ, Cstart, Ccheck;
	cpair ZZCdir, Hdir1, Hdir2, Hdir;
	double ZZdist;
    if(b4) {
        sR5 = moveIt(stt, -1);
		Hdir1 = sR5->C1->growth_vector;
		Hdir2 = sR5->C2->growth_vector;
        C_RZZ = sR5->C1->C1;
		ZZCdir = get_vector(C_1->C2->coords, C_1->C2->C2->coords );
		Hdir = C_2->growth_vector;
		ZZdist = getDistance_twoC(C_1->C2, C_2);
    }else {
        sR5 = moveIt(stt, 1);
		Hdir1 = sR5->C1->growth_vector;
		Hdir2 = sR5->C2->growth_vector;
        C_RZZ = sR5->C2->C2;
		ZZCdir = get_vector(C_1->coords, C_1->C2->coords );
		Hdir = C_1->growth_vector;
		ZZdist = getDistance_twoC(C_1, C_2->C1);
    }
	cpair starting_direction = get_vector(sR5->C1->C1->coords, sR5->C1->coords);
	cpair R5vec = get_vector(sR5->C1->coords, sR5->C2->coords);
	//Check for ZZ site
	Ccheck = sR5->C1->C1;
	cpair mpos = jumpToPos(Ccheck->coords, ZZCdir, 1.4);
	//cpair mpos = jumpToPos(Ccheck->coords, normAngle(Ccheck->C1->bondAngle1 - 60), 0, 1.4);
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
		// update H atoms
		//if (b4) updateA(C_RZZ->C1, C_2->C2, 'H');
		//else updateA(C_1->C1, C_RZZ->C2, 'H');
	}
	else { // migration over bridge.
		if (b4){
			//Adjust angle of previous carbon
			//Cstart->bondAngle1 = normAngle(Cstart->C1->bondAngle1 - 60);
			//Connect carbon to existing ZZ carbon
			saveXYZ("KMC_DEBUG/Migration_over_bridge");
			Cpointer C_newbridge, C_sharedbridge, C_oldbridge;
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
	}
    
    // edit sites. first identify the neighbouring sites of resulting RZZ & R5
    Spointer S1, S2, S3, S4;
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
    updateCombinedSites(stt); updateCombinedSites(sR5); // new FE and AC
    updateCombinedSites(S1); updateCombinedSites(S2); 
    updateCombinedSites(S3); updateCombinedSites(S4); // neighbours
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
    Cpointer CRem, CFE;// C to be removed and remaining C from original FE3 site
	cpair Cdir = get_vector(C_1->C2->coords, C_1->coords);
    if(b4) {
        sFE3 = moveIt(stt, -2); 
        CFE = sFE3->C2;
        CRem = sFE3->C1;
    }else {
        sFE3 = moveIt(stt, 2);
        CFE = sFE3->C1;
        CRem = sFE3->C2;
    }
	// close BY5 to form R6. First remove all bulk C
	for (int i = 0; i != 3; i++) removeC(C_1->C2, true);
	// add a C atom
	Cpointer Cnew;
	Cnew = addC(C_1, C_1->growth_vector, 1.4);
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
		updateSites(S1, S1->C1, CR5_1, -1);
        addR5toSite(S1, S1->C1, S1->C2); // remove a bulk C from S1 and add R5
    } else {
        S2 = moveIt(sFE3, 1); // neighbour of FE3
        S1 = moveIt(stt, -1); // neighbour of BY5 (stt)
		convSiteType(sFE3, CR5_1, CR5_1->C2, R5); // convert FE3 to R5
		convSiteType(stt, C_1->C2, CR5_1, RFE); // convert BY5 to RFE
        updateSites(S1, S1->C1, C_1->C2, +1); // add a bulk C to S1
		updateSites(S2, CR5_1->C2, S2->C2, -1);
        addR5toSite(S2, S2->C1, S2->C2); // remove a bulk C from S2 and add R5
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
   addR5internal(CR5_1, CR5_1->C2);
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
	Cpointer now = C_1->C2, Cbridge1, Cbridge2;
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
			Cbridge1 = b;
			Cbridge2 = now;
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
	int ntype_site = (int)stt->type;
	int ntype1 = (int)moveIt(stt, -1)->type;
	int ntype2 = (int)moveIt(stt, 1)->type;
	int newType;
	if (ntype1 < 5 && ntype2 < 5) {
		newType = (2002 + ntype1 + ntype2);
		// convert site
		if (newType>2005) {
			//saveDOT(std::string("BY5ClosureProblem.dot"));
			std::cerr << "ERROR: newType is > 65 (PAHProcess::proc_L5R_BY5)\n";
		}
		Spointer Srem1 = moveIt(stt, -1);
		Spointer Srem2 = moveIt(stt, 1);
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
	else if ( (ntype1 >= 2002 && ntype1 <= 2003) || (ntype2 >= 2002 && ntype2 <= 2003) ){
		ntype1 = ntype1 % 2002;
		ntype2 = ntype2 % 2002;
		newType = 2104 + ntype1 + ntype2;
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}
	else {
		int new_point = 0;
		newType = (ntype1 + ntype2 + 2002);
		if ((kmcSiteType)newType == None) {
			//saveDOT(std::string("BY5ClosureProblem.dot"));
			std::cerr << "ERROR: newType is None (PAHProcess::proc_L5R_BY5)\n";
		}
		convSiteType(stt, moveIt(stt, -1)->C1, moveIt(stt, 1)->C2, (kmcSiteType)newType);
	}
	
	if (bridged_before){
		if (ostype >= 100) {
			ostype = ostype % 10;
			ostype += 2100;
		}
		else ostype += 2000;
		convSiteType(other_side, other_side->C1, other_side->C2, (kmcSiteType)ostype);
		updateCombinedSites(other_side);
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
	m_pah->m_rings5_Embedded++;
	//Optimise if needed
	if (double dist = getDistance_twoC(C_1, C_2) > 2.6){
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
	}
	else addR5internal(C_1,C_2);
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
	Cpointer Cres = addC(Cstart, CZZdir, 1.4, true);
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
		if( (int)S1->type >= 501 && (int) S1->type <= 504) updateSites(sFE3, Cstart, C_1, +500); // convert FE3 to R5R6
    } else {
        S1 = moveIt(stt, -1); // neighbour of BY5 (stt)
        S2 = moveIt(sFE3, 1); // neighbour of FE3
        updateSites(sFE3, Cstart, C_2->C2->C2, +1); // convert FE3 (FE) to ZZ
        convSiteType(stt, C_1->C2, C_2, FE); // convert BY5 to FE
        updateSites(S1, S1->C1, C_1->C2, +1); // add a bulk C to S1
        updateSites(S2, sFE3->C2, S2->C2, -1); // remove a bulk C from S2
		if((int)S2->type >= 501 && (int)S2->type <= 504) updateSites(sFE3, Cstart, C_2->C2->C2, +500); // convert FE3 to R5R6
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
}
// ************************************************************
// ID20 & ID21- R6 oxidation at ZZ site
// ************************************************************
void PAHProcess::proc_O6R_FE2(Spointer& stt, Cpointer C_1, Cpointer C_2) {
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
    Cpointer C_bulk;
    for(C_bulk = S1->C1; C_bulk != S1->C2; C_bulk=C_bulk->C2) {
        if(C_bulk->bridge) return;
    }
    for(C_bulk = S2->C1; C_bulk != S2->C2; C_bulk=C_bulk->C2) {
        if(C_bulk->bridge) return;
    }
    // Identify C1 & C2 for new ZZ site
    Cpointer C_start, C_end;
    if(b4) {
        C_start = stt->C1->C1->C1;
        C_end = stt->C2->C2;
    } else {
        C_start = stt->C1->C1;
        C_end = stt->C2->C2->C2;
    }
	cpair CZZdir = invert_vector(C_start->C2->growth_vector);
	cpair Cdir = get_vector(C_start->coords, C_start->C2->coords);
    // Remove the 3 C atoms after C_start
    removeC(C_start->C2, false);
    removeC(C_start->C2, false);
    removeC(C_start->C2, false);
    // add a C atom after C_start (bulk in ZZ)
	addC(C_start, CZZdir, 1.4, true);
	updateA(C_start, 'H', Cdir);
	updateA(C_end, 'H', Cdir);
    //addC(C_start, normAngle(C_start->bondAngle1-120), normAngle(C_start->bondAngle1-60), 1.4, true);
    // update H
    //updateA(C_start, C_end, 'H');
    // Remove one of the FE in FE2
    removeSite(other);
    // Update Sites and neighbouring sites
    Spointer S3, S4;
    S1 = moveIt(stt, -1); S2 = moveIt(stt, 1);
    updateSites(stt, C_start, C_end, +1); // FE --> ZZ
    updateSites(S1, S1->C1, C_start, -1); // S1 --> reduce 1
    updateSites(S2, C_end, S2->C2, -1); // S2 --> reduce 1
    // update combined sites for all sites and their neighbours
    S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
    updateCombinedSites(stt);
    updateCombinedSites(S1); updateCombinedSites(S2);
    updateCombinedSites(S3); updateCombinedSites(S4);
    addCount(0,-1);
    m_pah->m_rings--;
    //saveDOT(dotname2.str());
}

// ************************************************************
// ID22- R6 desorption from FE to form AC (Added by GLC)
// ************************************************************
void PAHProcess::proc_D6R_FE_AC(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_O6R_FE_HACA_O2(stt, C_1, C_2);
}

//
// ************************************************************
// ID22- Bay-capping
// ************************************************************
void PAHProcess::proc_B6R_ACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//printSitesMemb(stt);
	//printStruct();//++++
	
	OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol);
	passbackPAH(mol);
	//saveXYZ("After_1stmin");
	Cpointer newC1;
	Cpointer newC2;
	cpair CR5_1, CR5_2;
	cpair vec1, vec2, Hdir1, Hdir2;

	//if(checkHindrance(stt)) {
	//    /*cout<<"Site hindered, process not performed.\n"*/ return;}
	if (!(C_1->C2->bridge) || !(C_2->C1->bridge)) { // check if bulk C in AC site is a bridge
		//Getting vectors
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
	m_pah->m_rings++;
	m_pah->m_rings5_Lone--;
	m_pah->m_rings5_Embedded++;
	//printSites(stt);
	//saveXYZ("Before_second_min");
	OpenBabel::OBMol newmol = passPAH();
	newmol = optimisePAH(newmol);
	passbackPAH(newmol);
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
			if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2204) return;
		}
	}
	if (CRem_next->bridge) return;
	if (CRem_next->C2->bridge) return; if (CRem_next->C1->bridge) return;
	if (CRem_next->C2->C2->bridge) return; if (CRem_next->C1->C1->bridge) return;
	
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
		opp_site_bool = true;
	}
	if (thirdC2 != NULLC) {
		opp_site_second = findSite(thirdC2);
		opp_site_bool_second = true;
	}
	if (thirdC_after != NULLC){
		opp_site_after = findSite(thirdC_after);
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
	//There are two main cases. The R5 of the ACR5 is shown on a short or a long side.
	double bond_distance = 1.4;
	double R5_dist = getDistance_twoC(CFE, CFE->C2);
	bool optimised = false;
	if (R5_dist < 1.6){
		//The ACR5 site is on the "short" side of an R5. 
		optimised = true;
		OpenBabel::OBMol mol = passPAH();
		mol = optimisePAH(mol);
		passbackPAH(mol);
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
	cpair R5coords_end = endposR5internal(CRem_before, CRem_next);
	if (m_pah->m_R5loc.size()>=1){
		for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
			double distR5s = getDistance_twoC(*it, R5coords_end);
			if (distR5s < 2.4) {
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
	cpair Cdir = C_1->growth_vector;
	// add a C atom
	Cpointer Cnew = addC(CFE, Cnewdir, bond_distance);
	updateA(Cnew, 'H', crossvec);
	removeC(CRem, false);
	if (optimised){
		//The ACR5 site is on the "short" side of an R5. 
		OpenBabel::OBMol newmol = passPAH();
		newmol = optimisePAH(newmol, 1000);
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
			m_pah->m_rings5_Embedded--;
			redrawR5(checkR5_1, C1_new, C2_new);
			//proc_G5R_ZZ(checkR5_1, checkR5_1->C1, checkR5_1->C2);
		}
		else {
			convSiteType(sFE2, sFE2->C1, sFE2->C2->C2, RFE);
			convSiteType(checkR5_1, C1_new->C2, C2_new->C1, R5);
			updateSites(checkR5_2, checkR5_2->C1, checkR5_2->C2, +100);
		}
	}

	
	//Reassign connectivity at next site
	if (b4) sFE2->C2 = C_2->C1->C1;
	else sFE2->C1 = C_1->C2->C2;

	// edit sites. first identify the neighbouring sites of resulting RFE & R5
	Spointer S1, S2, S3, S4;
	Cpointer C21 = C_2->C1;
	Cpointer C22 = C_2->C2;
	Cpointer C11 = C_1->C1;
	Cpointer C12 = C_1->C2;
	if (b4) {
		S1 = moveIt(sFE2, -1);
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
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2005){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
			if (!optimised) addR5internal(sFE2->C2->C1, sFE2->C2);
		}
		convSiteType(stt, sFE2->C2, stt->C2, ZZ);

	}
	else {
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
		else if ((int)sFE2->type >= 2003 && (int)sFE2->type <= 2005){ //sFE2 is a BY5 {
			updateSites(sFE2, sFE2->C1, sFE2->C2, +101);
			if (!optimised) addR5internal(sFE2->C1, sFE2->C1->C2);
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
		if ( (int)opp_site->type >= 2100) updateSites(opp_site, opp_site->C1, opp_site->C2, -100);
		else updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		Spointer S1_opp_site = moveIt(opp_site, -1);
		Spointer S2_opp_site = moveIt(opp_site, +1);
		/*if ((int)S1_opp_site->type >= 501 && (int)S1_opp_site->type <= 504 ) {
			updateSites(opp_site, opp_site->C1, opp_site->C2, +400);
		}
		else if ((int)S2_opp_site->type >= 501 && (int)S2_opp_site->type <= 504 ) {
			updateSites(opp_site, opp_site->C1, opp_site->C2, +400);
		}*/
		updateCombinedSites(opp_site);
	}
	else if (opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		updateSites(opp_site, opp_site->C1, opp_site->C2, -2000);
		if ((int)opp_site_after->type >= 500 && (int)opp_site_after->type <= 700) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -400);
		else if ((int)opp_site_after->type >= 1000 && (int)opp_site_after->type <= 2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, -800);
		updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
		updateCombinedSites(opp_site);
		updateCombinedSites(opp_site_after);
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
		updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +500);
		updateCombinedSites(opp_site_second);
		updateCombinedSites(opp_site_after);
	}
	else if (!opp_site_bool && !opp_site_bool_second && opp_site_bool_after) {
		if ( (int)opp_site_after->type >=2000) updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +100);
		else updateSites(opp_site_after, opp_site_after->C1, opp_site_after->C2, +2000);
		updateCombinedSites(opp_site_after);
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
	if (checkHindrance_newC(C_1)|| checkHindrance_newC(C_2)) {
		/*cout<<"Site hindered, process not performed.\n"*/ return;
	}
	Spointer S1 = moveIt(stt, -1);
	Spointer S2 = moveIt(stt, 1);
	Spointer S1check = moveIt(stt, -1);
	Spointer S2check = moveIt(stt, 1);
	Spointer SR5;
	bool b4 = false;
	if (stt->type == R5R6ZZ) { //G6R_R5R6ZZ
		if ( (((int) S1->type >= 501 && (int) S1->type <= 504) || ((int) S1->type >= 1002 && (int) S1->type <= 1004) || ((int) S1->type >= 2103 && (int) S1->type <= 2204 )) && (((int) S2->type >= 501 && (int) S2->type <= 504) || ((int) S2->type >= 1002 && (int) S2->type <= 1004) || ((int) S2->type >= 2103 && (int) S2->type <= 2204 ) ) ) {
			for (int i = 1; i<=3; i++){
				//Both sides of the R5R6ZZ seem to be possible. Proceed to identify them.
				int S1_before = (int)moveIt(S1check, -1)->type;
				int S2_after = (int)moveIt(S2check, +1)->type;
				if ( ((S1_before >= 501 && S1_before <= 504) || (S1_before >= 1002 && S1_before <= 1004) || (S1_before >= 2103 && S1_before <= 2204 ) ) && ( (S2_after >= 501 && S2_after <= 504) || (S2_after >= 1002 && S2_after <= 1004) || (S2_after >= 2103 && S2_after <= 2204 ) ) ){
					//Both second neighbour sites could be the coupled site to the R5R6ZZ. Move to next neighbours.
				}
				else if ( ((S1_before >= 501 && S1_before <= 504) || (S1_before >= 1002 && S1_before <= 1004) || (S1_before >= 2103 && S1_before <= 2204 ) ) && !( (S2_after >= 501 && S2_after <= 504) || (S2_after >= 1002 && S2_after <= 1004) || (S2_after >= 2103 && S2_after <= 2204 ) ) ){
					SR5 = S2;
					S2 = moveIt(SR5, +1);
					b4 = false;
					break;
				}
				else if ( !((S1_before >= 501 && S1_before <= 504) || (S1_before >= 1002 && S1_before <= 1004) || (S1_before >= 2103 && S1_before <= 2204 ) ) && ( (S2_after >= 501 && S2_after <= 504) || (S2_after >= 1002 && S2_after <= 1004) || (S2_after >= 2103 && S2_after <= 2204 ) ) ){
					SR5 = S1;
					S1 = moveIt(SR5, -1);
					b4 = true;
					break;
				}
				else {
					//No associated site. Odd. Throw error.
					cout << "G6R_R5R6ZZ has a site with no correct neighbours. Error.";
					std::list<std::string> Sitelist_before = copySites(stt);
					printBeforeSites(Sitelist_before);
					printSites(stt);
					std::string filename = "KMC_DEBUG/RZZ_neighbour_error_";
					filename.append(std::to_string(G6RRZZ_error_counter));
					saveXYZ(filename);
					G6RRZZ_error_counter++;
					return;
				}
				S1check = moveIt(S1check, -1);
				S2check = moveIt(S2check, +1);
			}
		}
		else {
			if ( ((int) S1->type >= 501 && (int) S1->type <= 504 ) || ((int) S1->type >= 602 && (int) S1->type <= 604 ) || ((int) S1->type >= 1002 && (int) S1->type <= 1004) || ((int) S1->type >= 2103 && (int) S1->type <= 2204 )) {
				SR5 = S1;
				S1 = moveIt(SR5, -1);
				b4 = true;
			}
			else if ( ((int) S2->type >= 501 && (int) S2->type <= 504) || ((int) S2->type >= 602 && (int) S2->type <= 604) || ((int) S2->type >= 1002 && (int) S2->type <= 1004 ) || ((int) S2->type >= 2103 && (int) S2->type <= 2204 )) {
				SR5 = S2;
				S2 = moveIt(SR5, +1);
				b4 = false;
			}
			else {
				cout << "R5R6ZZ not neighbouring an R5R6 site. Error.\n";
				std::list<std::string> Sitelist_before = copySites(stt);
				printBeforeSites(Sitelist_before);
				printSites(stt);
				std::string filename = "KMC_DEBUG/RZZ_neighbour_error_";
				filename.append(std::to_string(G6RRZZ_error_counter));
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
			else updateSites(S1, S1->C1, stt->C1, +1);
		}
		else if ((int) SR5->type >= 501 && (int) SR5->type <= 504 ){
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
			else updateSites(S2, stt->C2, S2->C2, +1);
		}
		else if ((int) SR5->type >= 501 && (int) SR5->type <= 504 ){
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
	if (getDistance_twoC(stt->C1->C1, stt->C2->C2) >= 3.0){
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
	updateCombinedSites(newS1); updateCombinedSites(newS2); // new sites
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
void PAHProcess::proc_MR5_R6(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng) {
	//printStruct();
	//First check if R6 is to the left or the right of R5
	bool b4 = false;
	Spointer sFE2, checksite, checkR5_1, checkR5_2;
	Cpointer CRem, CRem_before, CRem_next, CFE;
	if ((int)moveIt(stt, -1)->type > 4 && (int)moveIt(stt, +1)->type > 4){
		//Pentagons to both sides, JP not allowed
		return;
	}
	checksite = moveIt(stt, -1);
	if ((int)checksite->type < 4) b4 = true ;
	cpair Cdir, Hdir;
	if (b4) {
		sFE2 = moveIt(stt, -1);
		checkR5_1 = moveIt (stt, -2);
		checkR5_2 = moveIt (stt, -3);
		CFE = C_1->C2;
		CRem = sFE2->C2;
		CRem_next = CRem->C1;
		CRem_before = CRem->C2;
		Cdir = C_2->growth_vector;
		Hdir = C_1->growth_vector;
	}
	else {
		sFE2 = moveIt(stt, 1);
		checkR5_1 = moveIt (stt, 2);
		checkR5_2 = moveIt (stt, 3);
		CFE = C_1;
		CRem = sFE2->C1;
		CRem_next = CRem->C2;
		CRem_before = CRem->C1;
		Cdir = get_vector(C_1->C2->coords, C_2->coords);
		Hdir = C_2->growth_vector;
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
			if ((int)checkR5_2->type >= 2002 && (int)checkR5_2->type <= 2204) return;
		}
	}
	if (CRem_next->bridge) return;
	if (CRem_next->C2->bridge) return; if (CRem_next->C1->bridge) return;
	if (CRem_next->C2->C2->bridge) return; if (CRem_next->C1->C1->bridge) return;
	//check that two pentagons (including internals) will not collide
	cpair R5coords = endposR5internal(CRem_before, CRem_next);
	if (m_pah->m_R5loc.size()>=2){
		for (std::list<cpair>::iterator it = m_pah->m_R5loc.begin(); it!= m_pah->m_R5loc.end(); ++it){
			double distR5s = getDistance_twoC(*it, R5coords);
			if (distR5s > 0.1 && distR5s <2.1) {
				return;
			}
		}
	}
	// add a C atom
	Cpointer newC = addC(CFE, Cdir, 1.4);
	updateA(newC, 'H', Hdir);
	//addC(CFE, normAngle(CFE->bondAngle1 + 30), normAngle(CFE->bondAngle1 - 30), 1.4);
	//CRem->C1->bondAngle1 = normAngle(CRem->C1->bondAngle1 - 30);
	//printStruct(CRem);
	removeC(CRem, false);
	if (b4) sFE2->C2 = C_2->C1;
	else sFE2->C1 = C_1->C2;

	// edit sites. first identify the neighbouring sites of resulting RFE & R5
	Spointer S1, S2, S3, S4;
	Cpointer C21 = C_2->C1;
	Cpointer C22 = C_2->C2;
	Cpointer C11 = C_1->C1;
	Cpointer C12 = C_1->C2;
	if (b4) {
		S1 = moveIt(sFE2, -1);
		if ((int)sFE2->type == 0){ //sFE2 is a FE
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			int stype = S1->type;
			if (stype >= 2002 && stype <= 2104) {
				stype = stype + 100;
				convSiteType(S1, S1->C1, sFE2->C1, (kmcSiteType)stype);
			}
			else if (stype == 0){
				convSiteType(S1, S1->C1->C1, S1->C2->C2, ZZ);
				convSiteType(sFE2, S1->C2, sFE2->C2, FE);
				Spointer S3 = moveIt(S1, -1);
				updateSites(S3, S3->C1, S3->C2->C1, -1);
				m_pah->m_rings5_Embedded--;
				removeR5internal(S1->C1->C2, S1->C2->C1); redrawR5(S1, S1->C1, S1->C2);
			}
			else {
				stype = stype + 500;
				convSiteType(S1, S1->C1, sFE2->C1, (kmcSiteType)stype);
				//updateSites(S1, S1->C1, sFE2->C1, 5);
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
			convSiteType(sFE2, sFE2->C1, sFE2->C2, R5R6);
			int stype = S2->type;
			if (stype >= 2002 && stype <= 2104) {
				stype = stype + 100;
				convSiteType(S2, sFE2->C2, S2->C2, (kmcSiteType)stype);
			}
			else if (stype == 0){
				convSiteType(S2, S2->C1->C1, S2->C2->C2, ZZ);
				convSiteType(sFE2, sFE2->C1, S2->C1, FE);
				Spointer S4 = moveIt(S2, +1);
				updateSites(S4, S4->C1->C2, S4->C2, -1);
				removeR5internal(S2->C1->C2, S2->C2->C1); redrawR5(S2, S2->C1, S2->C2);
				m_pah->m_rings5_Embedded--;
			}
			else {
				stype = stype + 500;
				convSiteType(S2, sFE2->C2, S2->C2, (kmcSiteType)stype);
				//updateSites(S2, sFE2->C2, S2->C2, 5);
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

// ************************************************************
// ID35- R7 growth on embedded-obstructed R5
// ************************************************************
void PAHProcess::proc_GR7_R5R6AC(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	//printSitesMemb(stt);
	//printStruct();//++++
	OpenBabel::OBMol mol = passPAH();
	mol = optimisePAH(mol);
	passbackPAH(mol);
	Cpointer newC1;
	Cpointer newC2;
	cpair Hdir1 = C_2->growth_vector;
	cpair Hdir2 = C_1->growth_vector;
	cpair FEdir = get_vector(C_1->coords,C_2->coords);
	if(checkHindrance_newC(C_1) || checkHindrance_newC(C_2)) {
		/*cout<<"Site hindered, process not performed.\n"*/ return;
	}
	if (!(C_1->C2->bridge) && !(C_2->C1->bridge)) { // check if bulk C in AC site is a bridge
		// Add and remove C
		//if(C_1->C2->C3 != NULL) C_1->C2->C3->C3 = NULL;
		//if(C_2->C1->C3 != NULL) C_2->C1->C3->C3 = NULL;
		
		removeC(C_1->C2, true);
		removeC(C_1->C2, true);
		removeC(C_2->C1, true);
		newC1 = addC(C_1, C_1->growth_vector, 1.4);
		updateA(C_1,'C', C_1->growth_vector);
		updateA(newC1, 'H', Hdir1);
		newC2 = addC(newC1, FEdir, 1.4);
		updateA(C_2,'C', C_2->growth_vector);
		updateA(newC2, 'H', Hdir2);
		
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
		//newC2 = addC(newC1, normAngle(newC1->C1->bondAngle1 - 60), normAngle(newC1->C1->bondAngle1 - 90), 1.4*1.5);
		newC1 = addC(C_1, C_1->growth_vector, 1.1);
		updateA(C_1,'C', C_1->growth_vector);
		updateA(newC1, 'H', Hdir1);
		newC2 = addC(newC1, FEdir, 1.4);
		updateA(C_2,'C', C_2->growth_vector);
		updateA(newC2, 'H', Hdir2);
	}
	//printStruct();
	// Add and remove H
	//updateA(C_1->C1, C_2->C2, 'H');
	// neighbouring sites:
	Spointer S1 = moveIt(stt, -1);
	Spointer S2 = moveIt(stt, 1);
	// Update Site and neighbours
	convSiteType (stt, newC1, newC2, FE); //FE of an heptagon
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
	else if ( (int)S1->type >= 501 && (int)S1->type <= 504){ 
		int stype = (int)S1->type;
		stype = (stype % 10) + 2001;
		convSiteType(S1, S1->C1, newC1, (kmcSiteType)stype); // neighbours
	}	
	else updateSites(S1, S1->C1, newC1, 1); // neighbours
	
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
	else if ( (int)S2->type >= 501 && (int)S2->type <= 504){ 
		int stype = (int)S2->type;
		stype = (stype % 10) + 2001;
		convSiteType(S2, newC2, S2->C2, (kmcSiteType)stype); // neighbours
	}	
	else updateSites(S2, newC2, S2->C2, 1); // neighbours
	// Update combined site for Site and neighbours
	Spointer S3, S4;
	S3 = moveIt(S1, -1); S4 = moveIt(S2, 1);
	updateCombinedSites(stt);
	updateCombinedSites(S1); updateCombinedSites(S2);
	updateCombinedSites(S3); updateCombinedSites(S4);
	// add ring counts
	m_pah->m_rings7_Embedded++;
	//printSites(stt);
	OpenBabel::OBMol newmol = passPAH();
	newmol = optimisePAH(newmol);
	passbackPAH(newmol);
}

// ************************************************************
// ID35- R7 growth on embedded-obstructed R5
// ************************************************************
void PAHProcess::proc_GR7_FEACR5(Spointer& stt, Cpointer C_1, Cpointer C_2) {
	proc_GR7_R5R6AC(stt, C_1, C_2);
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
	OpenBabel::OBMol mol = passPAH();
	mol = mol = optimisePAH(mol);
	passbackPAH(mol);
	proc_L6_BY6(stt, C_1, C_2);
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
	OpenBabel::OBMol mol = passPAH();
	mol = mol = optimisePAH(mol);
	passbackPAH(mol);
	proc_L6_BY6(stt, C_1, C_2);
	m_pah->m_rings--; m_pah->m_rings7_Embedded++;
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