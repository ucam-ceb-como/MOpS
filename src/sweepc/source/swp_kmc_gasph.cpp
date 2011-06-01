/*!
  * \Author     Zakwan Zainuddin (zz260)
  * \file       swp_gasph.cpp
  *
  * \brief        Implementation of the KMCGasph class declared in swp_kmc_gasph.h
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Stores the gas phase profiles (time, pressure, temperature, species concentrations)

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

#include <iostream>
#include <string>
#include <map>
#include "swp_kmc_gasph.h"
#include "swp_kmc_typedef.h"
#include "swp_kmc_jump_process.h"
#include "swp_params.h"

#include "rng.h"
#include "choose_index.hpp"
#include "csv_io.h"
#include "string_functions.h"
#include "linear_interpolator.hpp"
#include "choose_index.hpp"

using namespace std;
using namespace Sweep;
using namespace Sweep::KMC_ARS;

// Constructors
// Default
KMCGasph::KMCGasph() {
}
//! Copy Constructor
KMCGasph::KMCGasph(KMCGasph& m) {
    m_profile_val = m.m_profile_val;
    m_profile_labels = m.m_profile_labels;
    setInterpolators(m_profile_val, m_profile_labels, m_profile_labels[0]);
    m_gpoint.m_data = m.m_gpoint.m_data;
    jplist = m.jplist;
}
//! From external file containing profiles
KMCGasph::KMCGasph(const std::string& filename) {
    cout<<"kMC Model initiated..\n Profile from "<<filename<<".\n\n";
    m_csvfilename = filename;
}
//! Destructor
KMCGasph::~KMCGasph() {
    m_gasProfile.clear();
}

//! Write Processes
//! Set csv input file name
void KMCGasph::setCSVname(const std::string& filename) {
    m_csvfilename = filename;
}
//! Load processes from process list
void KMCGasph::loadProcesses(std::vector<JumpProcess*> (*jp)(const KMCGasPoint&)) {//how##
    jplist = jp(this->m_gpoint);
}
//! Loads profile of specific profile names from an external file
void KMCGasph::loadProfileCSV() {
    // open csv file containing profile
    int op = m_csvfile.Open(m_csvfilename, false);
    // check if opening successful
    if(op!=0) {
        cout<<"ERROR: Failed to open "<<m_csvfilename<<'\n';
		std::ostringstream msg;
        msg << "Not able to open gas profile values from "
            << m_csvfilename
            << " (Sweep::KMC_ARS::KMCGasph::LinearInterpolator)";
        return;
    }else {
        cout<<"File "<< m_csvfilename << " opened successfully. Preparing to load values..\n";
    }
    // read profile labels from csv file
    //cout<<"\nReading profiles...\n";
    m_csvfile.Read(m_profile_labels);
    // remove whitespaces from each name and create empty vectors for each.
    for(int i=0; i!= (int)m_profile_labels.size(); i++) {
        m_profile_labels[i] = Strings::removeWhiteSpace(m_profile_labels[i]);
        m_profile_val.push_back(rvector());
    }
    // initialising gaspoints
    getIt();
    // create a temporary vector to store each line in csv file
    rvector temp;
    m_csvfile.Read(temp);
    do {
        for(int i=0; i!=(int)temp.size(); i++) { // read each line and puts in vectors of column values
            std::string now = m_profile_labels[i];
            real y;
            y=temp[i];
            m_profile_val[i].push_back(y);
        }
        m_csvfile.Read(temp); // temp is cleared each time Read is called and reads the next line
    }
    while (!temp.empty()); // Read returns empty vector when end of line is reached
	if((int)m_profile_val[0].size() > (int)m_profile_val[1].size())
		m_profile_val[0].pop_back(); // delete last empty line in .csv file
    // create profile interpolators for each 
    setInterpolators(m_profile_val, m_profile_labels, m_profile_labels[m_gpoint.Time]);
    m_csvfile.Close();
    // initialise m_gaspoint with values at t=0
    interpolateProfiles(0, true, 1);
    //cout<<"\nFinished reading profiles.\n";
} 

//! Creates profile interpolators
void KMCGasph::setInterpolators(std::vector<rvector>& values, // vector of values
                                std::vector<std::string>& labels, // vector of profile labels
                                std::string& Xp) {// profile describing position of data (independent profile, X)
    // find the index for the independent profile X
    int Xp_index = Strings::findinlist(Xp, labels);
    // obtain vector of values for X
    rvector Xp_val = values[Xp_index];
    // create interpolators against X for each profile
    for(int i=0; i!=(int) labels.size(); i++) {
        Utils::LinearInterpolator<real,real> interp(values[Xp_index], values[i]);
        // store interpolator in m_gasProfile vector of interpolators
        m_gasProfile.push_back(interp);
        // displays label of profile successfully loaded
       // cout<<labels[i]<<" profile loaded..\n";
    }
}
// KMC Algorithm and Read Processes
////! Calculates each jump rates and stores in a vector
//rvector KMCGasph::calculateRates(std::vector<JumpProcess*>& p) const {
//    // temporary vector to store rates of each JumpProcess
//    rvector temp; 
//    for(int i=0; i!= (int)p.size(); i++) {
//        real a = p[i]->getRate(); // gets rate
//        temp.push_back(a); // store in temporary vector
//    }
//    return temp;
//}
//! Calculates total jump rates
//real KMCGasph::getTotalRate(const rvector& jrate) {
//    // calculates sum of all elements of jrate
//    real sum=0;
//    for(int i=0; i!=(int) jrate.size(); i++) {
//        sum += jrate[i];
//    }
//    return sum;
//}
//! Choosing a reaction to be taken place, returns vector index of process
int KMCGasph::chooseReaction(const rvector& jrate, real (*rand_u01)()) const {
    // chooses index from a vector of weights (real number in this case) randomly
    int ind = chooseIndex<real>(jrate, rand_u01);
    return ind;
}
//! Interpolates all profiles and calculates gas conc. at time tnow, multiplied by a 
//! a reduction factor r_factor.
void KMCGasph::interpolateProfiles(const real& tnow, bool interpolate_gas, real r_factor) {
    //if(label == "") return 1; <-- moved to getVal
    // temporary storage for profile labels
    std::vector<std::string> vec = m_profile_labels;
    // finds indices for temperature, time and pressure (to calculate concentration)
    int Tindex = m_profile_it[1];
    int Timeindex = m_profile_it[0];
    int Pindex = m_profile_it[m_profile_it.size()-2];
    //// checks if the indices above were successfully found
    //if(Tindex<0) {
    //    cout<<"ERROR: Temperature label \""<<Tlabel<<"\" cannot be found...";
    //    return;
    //}
    //if(Timeindex<0) {
    //    cout<<"ERROR: Time label \""<<Timelabel<<"\" cannot be found...";
    //    return;
    //}
    //if(Pindex<0) {
    //    cout<<"ERROR: Pressure label \""<<Plabel<<"\" cannot be found...";
    //    return;
    //}
    // P and T values
    real P_num = m_gasProfile[m_profile_it[13]].interpolate(tnow);
    real T_num = m_gasProfile[Tindex].interpolate(tnow);
    // for conversion factor calculation
    real Rc = 8.3142;
    real A = 1/*.01325*//(Rc*10);
    intvector& p = m_profile_it;
    if(interpolate_gas) {
        for(int i=0; i!=(int)p.size()-1; i++) {
            int a = p[i];
            if(a != -1) {
            // get label index of data value
            //int index = Strings::findinlist(label, vec); <-- replaced with i
            real y; // conversion factor for data value
            // checks if data value is not a value for T, P or t, i.e. mol fraction of gas specified as 'label'
            if(a!= Tindex && a!= Pindex && a!= Timeindex) {
                // calculates conversion factor according to ideal gas law (n/V) = (P/RT)
                y = P_num*A*r_factor/T_num;
            }else y=1; // otherwise if data value is value for T, P or t, conversion not required
            // update m_interpolated with converted value
            real gc = m_gasProfile[a].interpolate(tnow);
            m_gpoint.m_data[i]= gc * y;
            }
        }
    }else {
        m_gpoint.m_data[m_gpoint.T] = m_gasProfile[Tindex].interpolate(tnow);
        m_gpoint.m_data[m_gpoint.P] = m_gasProfile[Pindex].interpolate(tnow);
    }
}
//! Gets t_max
real KMCGasph::gett_stop() {
    int index = Strings::findinlist(Timelabel, m_profile_labels);
    // Get the last value for time
    real temp = m_profile_val[index].back();
    return temp;
}
/*
void changeValues(intvector& vc, int& a, int& v, const string& s, int m) {
    a = v;
    vc[m] = v;
    cout<<"Profile "<<s<<" on column "<<(a+1)<<'\n';
}*/

void KMCGasph::getIt() {
    initIt();
    intvector& p = m_profile_it;
    int i;
    std::vector<std::string> vec = m_profile_labels;
    for(i=0; i<(int) vec.size(); i++) {
        if(vec[i] == Timelabel) p[0] = i; //time
        else if(vec[i] == Tlabel) p[1] = i; //T
        else if(vec[i] == "H2") p[2] = i; //H2
        else if(vec[i] == "H") p[3] = i; //H
        else if(vec[i] == "O2") p[4] = i; //O2
        else if(vec[i] == "OH") p[5] = i; //OH
        else if(vec[i] == "C2H2") p[6] = i; //C2H2
        else if(vec[i] == "C2H6") p[7] = i; //C2H6
        else if(vec[i] == "A1") p[8] = i; //C6H6
        else if(vec[i] == "H2O") p[9] = i; //H2O
        else if(vec[i] == "CH4") p[10] = i; //CH4
        else if(vec[i] == "CO") p[11] = i; //CO
        else if(vec[i] == "CO2") p[12] = i; //CO2
        else if(vec[i] == Plabel) p[13] = i; //P
    }
    int s = (int) vec.size();
    p[m_gpoint.total-1] = s;
}

void KMCGasph::initIt() {
    m_profile_it.clear();
    for(int i=0; i!=m_gpoint.total; i++) m_profile_it.push_back(-1);
}