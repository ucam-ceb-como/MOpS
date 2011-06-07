/*!
  * \Author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_simulator.cpp
  *
  * \brief      Implementation for swp_kmc_simulator.h
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Simulator for the kMC Model.

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

#include "swp_cell.h"
#include "mt19937.h"
#include "csv_io.h"
#include "rng.h"
#include <string>
#include <iostream>
#include <time.h>
#include <math.h>
#include <cstdlib>
#include "string_functions.h"
#include "swp_kmc_simulator.h"


using namespace std;
using namespace Sweep;
using namespace Sweep::KMC_ARS;
using namespace Strings;

// Default names for CSV outputs if not specified
std::string default_timer_csv = "KMC_Model/PAH_loop_timer.csv";
std::string default_rxncount_csv = "KMC_Model/PAH_reaction_count.csv";
std::string default_pahlist_csv = "KMC_Model/PAH_CH_site_list.csv";

// Vector of all site types
std::vector<kmcSiteType> allSiteType = vectSiteType();

//! Default Constructor
KMCSimulator::KMCSimulator() {
    m_t = 0;
	//input format is kept unchanged. use C2H4_profiles_0p080.csv which is the same as gasphase.inp
	std::string p = "0p080";//added by dongping 12.04
	// Gas profile (csv)
    std::string profileInputName = "C2H4_profiles_"+p+".csv";
	// Name of dot file (PAH structure)
    std::string DOToutputName = "DOT files/KMC_"+p+"_";
	// Total actual times
    std::string CSVtimerName = "Results/PAH_loop_timer_"+p+"_cpp.csv";
	// Reactions counts
    std::string CSVreactioncountName = "Results/Reaction_Count_"+p+"_cpp.csv";
	// PAH details (N(C), N(H), N(sites))
    std::string CSVpahlistName = "Results/PAH_CH_Site_List_"+p+"_cpp.csv";
	// Rates of reactions for last run
    std::string CSVrates = "Results/PAH_rates_cpp.csv";
	std::string CSVtimestep="zztimestep.csv";
	

    setCSVinputName(profileInputName);
    setDOToutputName(DOToutputName);
    setCSVtimerName(CSVtimerName);
    setCSVreactioncountName(CSVreactioncountName);
    setCSVpahlistName(CSVpahlistName);
    setCSVratesName(CSVrates);
	setCSVtimestep(CSVtimestep);
	m_simGas.setCSVname(m_csv_in);
    m_simGas.loadProfileCSV();
	m_simGas.loadProcesses(JumpProcessList::obtainJumpProcess);
}

//! Copy Constructor
KMCSimulator::KMCSimulator(KMCSimulator& s) {
    m_t = s.m_t;
    KMCGasph m_simGas(s.m_simGas);
    //m_simPAH = new PAHStructure(*(s.m_simPAH));
    m_simPAHp = PAHProcess(*m_simPAH);
}

//! Default Destructor
KMCSimulator::~KMCSimulator() {
    //delete m_simPAH;
}

//! Set PAH to be simulated
void KMCSimulator::targetPAH(PAHStructure& pah) {
    m_simPAH = &pah;
    m_simPAHp = PAHProcess(*m_simPAH);
}
//! KMC algorithm to calculate timestep
real KMCSimulator::timeStep(real totalRate, real (*rand_u01)()) const {
    // get an exponentially distributed random number
    real logrand = logrnd(rand_u01);
    // time step = -N(e) / R(total)
    real tau = -logrand/totalRate;
    return tau;
}
//! Run simulation for a set number of times until maxTime is reached.
//! If maxTime specified as 0, last time value from csv file is taken.
void KMCSimulator::runSimulation(const int& total_runs,
									const real& startTime, 
									const real& maxTime, 
									const int& steps,
									const bool CHcol,
									const bool save_rates,
									const bool save_dots,
									const std::string& CHoutput,
									const bool save_CH,
									int (*rand_int)(int,int),
									real (*rand_u01)()) {
	assert(total_runs>0);
	// Initialising output csv files
    initCSVIO();
    /*m_simGas.setCSVname(m_csv_in);
    m_simGas.loadProfileCSV();
    m_simGas.loadProcesses(JumpProcessList::obtainJumpProcess);*/
	//targetPAH(pah);
	KMCGasph gasph_copy(m_simGas);
    // get t_max for simulation, if t_max not specified, set to final value in time profile
    real t_max;
	//const int run=1;
    if(maxTime == 0) t_max = m_simGas.gett_stop();
    else t_max = maxTime;
    //bool save_rates = false;
    m_simPAHp.m_rates_save = save_rates;
	// set minimum number of steps
    int min_steps = steps;
    // waiting time
    real t_step_max = t_max/min_steps;
	//initialise structure
	m_simPAH = new PAHStructure();
	m_simPAHp.setPAH(*m_simPAH);
	m_simPAHp.initialise(PYRENE);
    //--------for C-H counts per interval----
	CSV_data CHdata(*this);
    int maxInterval = min_steps;
	if(save_CH) {
    CHdata.setName(CHoutput);
    CHdata.initData(total_runs, maxInterval, t_max, m_simPAHp.getCHCount(), m_simGas);
	}
    //---------------------------------------
	for(int run=1; run<= total_runs; run++) {
		unsigned long sd = 103000+run;
		Sweep::init_genrand(sd);
		//initialise structure
		m_simPAHp.initialise(PYRENE);
        // Start timer
        clock_t timerStart = clock();
        // set current time to startTime
        m_t = startTime;
        // name for dot files
        std::string dotfile; 
        // interval number count
        int interval_num = -1; // -1 - starting
        // time of next iteration
        real t_next = 0;
        // loop counts
        int count = 0;
        // Store jump rates into a vector
        rvector jump_rates_vector((int)m_simGas.jplist.size());
        // find label for T
        int T_label = m_simGas.m_gpoint.T;
        // set reaction counts to 0
        initReactionCount();
        if(save_CH) CHdata.m_intervalcount = 0;
        bool savedotperinterval = false;
		vector<double> time_step_vector;//##
		// for cloned PAH
		bool copied = false;
		PAHStructure* m_copyPAH;
		real starttime = 0;

        while (m_t < t_max) {
            // simulator will not run if csv filename for gas profiles not specified
            if (m_csv_in.length() == 0) {
                cout << "ERROR: No CSV input file specified, stopping simulator.\n";
                return;
            }

			if(m_t > (t_max/2) && !copied) {
				m_copyPAH = m_simPAHp.clonePAH();
				starttime = m_t;
				copied = true;
				if(save_dots) {
					dotfile = m_dot_out+ Strings::cstr(run) + "_" + Strings::cstr(count)+ ".dot";
					m_simPAHp.saveDOT(dotfile);//to obtain the PAH structure
					PAHProcess p(*m_copyPAH);
					p.saveDOT(std::string("DOT files/copiedstructure_start.dot"));
				}
			}
            // Interpolate profile values at time m_t
            m_simGas.interpolateProfiles(m_t, true, 1);

            // Calculate rates of each jump process and store in vector & store total jump rate
            real totalJumpRate = JumpProcessList::calculateRates(m_simGas.m_gpoint, m_simPAHp, m_t, m_simGas.jplist, jump_rates_vector);
            rvector temp = jump_rates_vector;
            
            // Calculate time step, update time
            real t_step;
            if(!save_rates) t_step = timeStep(totalJumpRate, rand_u01);
            else t_step = t_step_max;
			time_step_vector.push_back(t_step);
			//writetimestepCSV(time_step_vector);
            t_next = m_t+t_step;

            if(t_next < t_max && t_step < t_step_max) {
                //saveDOTperXsec(t_step_max, sd, m_t, t_max, gasph_copy, interval_num);
                // Choose reaction according to rates
                int chosen_proc = m_simGas.chooseReaction(jump_rates_vector, rand_u01);
                // Get site type and structure change process from the jump process class
                JumpProcess* jp_perf = m_simGas.jplist[chosen_proc];
                kmcSiteType jp_site = jp_perf->getSiteType();
                int jp_id = jp_perf->getID();
                

                //--------for C-H counts per interval----
                if(save_CH) CHdata.addData(m_simPAHp.getCHCount(), t_next, run, m_simPAHp, savedotperinterval);
                //---------------------------------------
                // Update data structure
                if(!save_rates) {
                    //cout<<m_simGas.jplist[chosen_proc]->getName()<<'\t'<<m_simPAHp.getCHCount().first << '\n';
                    bool process_success = m_simPAHp.performProcess(jp_site, jp_id, rand_int);
					/*if(!m_simPAHp.checkCoordinates()) {
						cout<<"Coordinates of structure did not pass test. Aborting..\n";
						std::ostringstream msg;
						msg << "ERROR: Structure did not pass PAHProcess::checkCoordinates."
							<< " (Sweep::KMC_ARS::KMCSimulator::runSimulation)";
						throw std::runtime_error(msg.str());
						assert(false);
					}
					if(!m_simPAHp.checkSiteContinuity()) {
						cout<<"Structure has non-continuous sites. Aborting..\n";
						std::ostringstream msg;
						msg << "ERROR: Structure did not pass PAHProcess::checkSiteContinuity."
							<< " (Sweep::KMC_ARS::KMCSimulator::runSimulation)";
						throw std::runtime_error(msg.str());
						assert(false);
					}*/
                    if(process_success) {
                    m_rxn_count[chosen_proc]++;
                    }
                }
            }else {
                t_next = m_t+t_step_max;
                //else m_t = t_max;
            }
            // For Rates
            if(save_rates) writeRatesCSV(run, m_t, jump_rates_vector);
            m_t = t_next;
            
            //if(run>1){
                    //saveDOTperXLoops(1, count, run);
            //}else saveDOTperXLoops(10, count, run);
            count++;
            // Save 1500th run into DOT file
            /*if(count==1500) {
            dotfile = m_dot_out+ Strings::cstr(run) + "_" + Strings::cstr(count)+ ".dot";
            m_simPAH->saveDOT(dotfile);
            }*/
        }
        //--------for C-H counts per interval----
		if(save_CH) {
			while(CHdata.m_intervalcount < maxInterval) {
				CHdata.addData(m_simPAHp.getCHCount(), t_max, run, m_simPAHp, savedotperinterval);
			}
		}
        //---------------------------------------
        // calculate time elapsed
        clock_t timerStop = clock();
        double timerElapsed = double(timerStop-timerStart)/CLOCKS_PER_SEC;

        // save structure into dot file
        //dotfile = m_dot_out+ Strings::cstr((int)run) + "_" + Strings::cstr(count)+ ".dot";
        if(save_dots) {
			dotfile = m_dot_out+ Strings::cstr(run) + "_" + Strings::cstr(count)+ ".dot";
			m_simPAHp.saveDOT(dotfile);//to obtain the PAH structure
		}

        // for gif animation:
        //saveDOTperXsec(t_step_max, sd, t_max, t_max, gasph_copy, interval_num);

        // show time elapsed on display
        cout<<"\nSimulator ran for run "<<run<<" looping "<<count<<" number of times...\n";
        cout<<"Ran for "<<timerElapsed<<" seconds.\n\n";
        if(m_simPAHp.getCHCount().first < 12 && save_CH) { // omit benzene results
            cout<<"C count out of range. Redoing run "<<run<<"\n\n";
            CHdata.delData(run);
			run--;
        } else {
        // save csv files for pah information
        writeTimerCSV(count, timerElapsed);
        writeRxnCountCSV(m_rxn_count);
        writeCHSiteCountCSV();
		//writetimestep(time_step_vector);
		}
		delete m_simPAH;
		updatePAH(m_copyPAH, starttime, (t_max - starttime), steps/2, rand_int, rand_u01, 1,1);
		if(save_dots) {
			dotfile = "DOT files/copiedstructure_end.dot";
			m_simPAHp.saveDOT(dotfile);//to obtain the PAH structure
		}
	}
    if(save_CH) CHdata.writeCSV(CHcol, false);
    // close csv files
    m_timer_csv.Close();
    m_rxn_csv.Close();
    m_pah_csv.Close();
    m_rates_csv.Close();
	delete m_simPAH;
}

//! Update structure of PAH after time dt
real KMCSimulator::updatePAH(PAHStructure* pah, 
							const real tstart, 
							const real dt,  
							const int waitingSteps,  
							int (*rand_int)(int,int), 
							real (*rand_u01)(), 
							real r_factor,
							int PAH_ID) {
   
	vector<int> m_rxn_count(18,0);
	real m_t = tstart;
    real t_max = m_t + dt;
    targetPAH(*pah);
	/*if(m_simPAHp.checkCoordinates())
		cout<<"Coordinates of structure OK. Commencing updatePAH..\n";
	else {
		cout<<"Coordinates of structure did not pass test. Aborting..\n";
		std::ostringstream msg;
		msg << "ERROR: Trying updatePAH on a structure which did not pass PAHProcess::checkCoordinates."
			<< " (Sweep::KMC_ARS::KMCSimulator::updatePAH)";
		throw std::runtime_error(msg.str());
		assert(false);
	}
	if(m_simPAHp.checkSiteContinuity())
		cout<<"Site continuity of structure OK. Commencing updatePAH..\n";
	else {
		cout<<"Structure has uncontinuous site. Aborting..\n";
		std::ostringstream msg;
		msg << "ERROR: Trying updatePAH on a structure which did not pass PAHProcess::checkSiteContinuity."
			<< " (Sweep::KMC_ARS::KMCSimulator::updatePAH)";
		throw std::runtime_error(msg.str());
		assert(false);
	}*/
    rvector jump_rates_vector((int)m_simGas.jplist.size());
    real t_next = m_t;
    real t_step_max = dt/waitingSteps;
    real oldtnext;
	int loopcount=0;

    while (m_t < t_max) {
		//this->m_simPAHp.printStruct();//++++
        m_simGas.interpolateProfiles(m_t, true, r_factor);
		loopcount++;

        // Calculate rates of each jump process and store in vector & store total jump rate
        real totalJumpRate = JumpProcessList::calculateRates(m_simGas.m_gpoint, m_simPAHp, m_t, m_simGas.jplist, jump_rates_vector);
        rvector temp = jump_rates_vector;
        // Calculate time step, update time
        real t_step = timeStep(totalJumpRate, rand_u01);
        t_next = m_t+t_step;
        if(t_next < t_max) {
			//if (pah->numofC()>5000&&pah->numofC()<6000)//||pah->havebridgeC()
			//if (PAH_ID==224835||PAH_ID==1011604||PAH_ID==1111604)
            // saveDOTperLoop(100000*tstart,loopcount,PAH_ID);
            // Choose reaction according to rates
            int chosen_proc = m_simGas.chooseReaction(jump_rates_vector, rand_u01);
            // Get site type and structure change process from the jump process class
            JumpProcess* jp_perf = m_simGas.jplist[chosen_proc];
            kmcSiteType jp_site = jp_perf->getSiteType();
            int jp_id = jp_perf->getID();
            //cout<<m_simGas.jplist[chosen_proc]->getName()<<'\n';//++++

            // Update data structure
            bool process_success = m_simPAHp.performProcess(jp_site, jp_id, rand_int);
			/*if(m_simPAH->m_parent->ID() % 100000 == 609) {
			if(!m_simPAHp.checkCoordinates()) {
				cout<<"ERROR: Invalid coordinates produced after performing process "
					<<jp_perf->getName()<<" (ID" <<jp_id<<")\n";
			}
			}*/
			if(process_success) {
                m_rxn_count[chosen_proc]++;
            }
        }else {
            oldtnext = t_next;
            t_next = m_t+t_step_max;
        }
		m_t = t_next;
        //saveDOTperXLoops(1, count, run);
       //count++;
    }
    return m_t;
}

/*void KMCSimulator::testSimulation(PAHStructure& pah, const unsigned long seed, int totalruns) {
    for(int run=1; run<=totalruns; run++) {
        clock_t timerStart = clock();
        int total_interval = 10;
        real maxtime = 0.0222029;
        real step = maxtime/total_interval;
        real tnow=0;
        // loop counts
        int count = 0;
        // Start timer
        initReactionCount();
        m_simPAHp.initialise(PYRENE);
        while(tnow<maxtime) {
            //tnow = i*step;
            tnow = updatePAH(pah, tnow, step, m_simGas, 10, count);
        }
        //updatePAH(pah, tnow, maxtime, m_simGas);
        m_simPAHp.saveDOT(std::string("DOT files/Test_update_PAH.dot"));
        // calculate time elapsed
        clock_t timerStop = clock();
        double timerElapsed = double(timerStop-timerStart)/CLOCKS_PER_SEC;
        writeTimerCSV(count, timerElapsed);
        writeRxnCountCSV(m_rxn_count);
        writeCHSiteCountCSV();
        cout<<"\nSimulator ran for run "<<run<<" looping "<<count<<" number of times...\n";
        cout<<"Ran for "<<timerElapsed<<" seconds.\n\n";
    }
    
    m_timer_csv.Close();
    m_rxn_csv.Close();
    m_pah_csv.Close();
}*/

//! Set csv filename for gas profiles
void KMCSimulator::setCSVinputName(const std::string& filename) {
    m_csv_in = filename;
}

//! Set output DOT file name "filename"_runs_finalloopnum.dot
void KMCSimulator::setDOToutputName(const std::string& filename) {
    m_dot_out = filename;
}
//! Set output CSV file name to keep track of timer counts
void KMCSimulator::setCSVtimerName(const std::string& filename) {
    m_timer_name = filename;
}
//! Set output CSV file name to keep track of reaction counts
void KMCSimulator::setCSVreactioncountName(const std::string& filename) {
    m_rxncount_name = filename;
}
//! Set output CSV file name to keep track of CH and site counts
void KMCSimulator::setCSVpahlistName(const std::string &filename) {
    m_pahlist_name = filename;
}

//! Set output CSV file name to keep track of time step//##
void KMCSimulator::setCSVtimestep(const std::string &filename) {
    m_timestep_name = filename;
}
//! Set output CSV file name to keep track of CH and site counts
void KMCSimulator::setCSVratesName(const std::string &filename) {
    m_rates_name = filename;
}
//! Writes data for timeCount.csv
void KMCSimulator::writeTimerCSV(const int& loop, const double& elapsedTime) {
    vector<double> temp;
    // writes: | Total Loop | Total time elapsed |
    temp.push_back((double)loop);
    temp.push_back(elapsedTime);
    m_timer_csv.Write(temp);
}
void KMCSimulator::writetimestep(const std::vector<double>& timestep){
	m_timestep_csv.Write(timestep);
}
//! Writes data for reaction_count.csv
void KMCSimulator::writeRxnCountCSV(const std::vector<int>& rc) {
    // change int vector to float
    std::vector<float> temp;
    for(int i=0; i<(int) rc.size(); i++) {
        temp.push_back((float)rc[i]);
    }
    m_rxn_csv.Write(temp);
}
//! Writes data for CH_site_list.csv
void KMCSimulator::writeCHSiteCountCSV() {
    std::vector<float> temp;
    // get CH count
    intpair CH = m_simPAHp.getCHCount();
    temp.push_back((float)CH.first);
    temp.push_back((float)CH.second);
    // get counts for all site types
    for(int i=0; i<(int)allSiteType.size(); i++) {
        int scount = m_simPAHp.getSiteCount(allSiteType[i]);
        temp.push_back((float) scount);
    }
    m_pah_csv.Write(temp);
}
//! Writes data for rates count (csv)
void KMCSimulator::writeRatesCSV(int runNo, real& time, rvector& v_rates) {
    //int convC[] = {1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19};
    int convM[] = {1, 2,14,15, 8,10,11,12,13, 3, 7, 9, 5, 4, 6,22,24,16,21};
    int ID;
    if(runNo==1) {
    std::vector<real> temp(25,0);
    temp[0] = time;
    for(int i=0; i!=(int)v_rates.size(); i++) {
        ID = m_simGas.jplist[i]->getID();
        temp[convM[ID-1]] = v_rates[i];
    }
    m_rates_csv.Write(temp);
    }
}
//! Initialise CSV_IOs
void KMCSimulator::initCSVIO() {
    // Check if CSV file names specified, if not put them as default
    if(m_timer_name.length() == 0) {
        cout<<"WARNING: Output CSV name for time count is not specified. Defaulting to "<<default_timer_csv<<"\n";
        m_timer_name = default_timer_csv;
    }
    if(m_rxncount_name.length() == 0) {
        cout<<"WARNING: Output CSV name for reaction count is not specified. Defaulting to "<<default_rxncount_csv<<"\n";
        m_rxncount_name = default_rxncount_csv;
    }
    if(m_pahlist_name.length() == 0) {
        cout<<"WARNING: Output CSV name for CH and site counts is not specified. Defaulting to "<<default_pahlist_csv<<"\n";
        m_pahlist_name = default_pahlist_csv;
    }
    // Open csv file
    m_timer_csv.Open(m_timer_name, true);
    m_rxn_csv.Open(m_rxncount_name, true);
    m_pah_csv.Open(m_pahlist_name, true);
    m_rates_csv.Open(m_rates_name, true);
	 m_timestep_csv.Open(m_timestep_name, true);//##
    // Write column headings for CSV files
    writeCSVlabels();
}
//! Initialise reaction count
void KMCSimulator::initReactionCount() {
    m_rxn_count.clear();
    for(int i=0; i<(int)m_simGas.jplist.size(); i++) {
        m_rxn_count.push_back(0);
    }
}
//! Load gas profiles
void KMCSimulator::loadGasProfiles() {
	m_simGas.setCSVname(m_csv_in);
    m_simGas.loadProfileCSV();
    m_simGas.loadProcesses(JumpProcessList::obtainJumpProcess);
}
//! Write column headings for CSV files
void KMCSimulator::writeCSVlabels() {
    // Write headings for timer
    std::vector<string> timer_headings;
    timer_headings.push_back("Total Loops"); timer_headings.push_back("Time Elapsed");
    m_timer_csv.Write(timer_headings);
    // Write headings for reaction count
    std::vector<string> rxn_headings;
    for(int i=0; i<(int)m_simGas.jplist.size(); i++) {
        // gets name of each jump process and puts them in a row
        rxn_headings.push_back(m_simGas.jplist[i]->getName());
    }
    m_rxn_csv.Write(rxn_headings);
    // Write headings for CH and site list
    std::vector<string> pah_headings;
    // write headings for N_C and N_H
    pah_headings.push_back("N_C");
    pah_headings.push_back("N_H");
    // write headings for number of sites for all site types "N(sitetype)"
    for(int i=0; i<(int) allSiteType.size(); i++) {
        std::string header = "N(";
        header = header.append(kmcSiteName(allSiteType[i]));
        header = header.append(")");
        pah_headings.push_back(header);
    }
    m_pah_csv.Write(pah_headings);
}
//! Save the structure DOT file after every X loops
void KMCSimulator::saveDOTperXLoops(int X, int &loopcount, int& runcount) {
    int m = loopcount % X;
    if(m==0) {
        string filename = "KMC_DEBUG/Run_";
        filename.append(Strings::cstr(runcount));
        filename.append("_Loop_");
        filename.append(Strings::cstr(loopcount));
        filename.append(".dot");
        m_simPAHp.saveDOT(filename);
    }
}
//added by dongping 18/04
void KMCSimulator::saveDOTperLoop(int LOOPcount,int loopcount, int PAH_ID) {
        string filename = "KMC_DEBUG/ID_";
		filename.append(Strings::cstr(PAH_ID));
		filename.append("_Run_");
        filename.append(Strings::cstr(LOOPcount));
        filename.append("_Loop_");
        filename.append(Strings::cstr(loopcount));
        filename.append(".dot");
        m_simPAHp.saveDOT(filename);
}
//! Save the structure DOT file after every X simulation sec interval
void KMCSimulator::saveDOTperXsec(const real& X, const int& seed, const real &time, const real &time_max, KMCGasph& copyMod, int& intervalcount) {
    int interval = (int) ceil(time/X);
    std::string graphTitle;
    if(intervalcount == -1) {
        copyMod.interpolateProfiles(0, false, 1);
        std::string temp = Strings::cstr((int)ceil(copyMod.m_gpoint.m_data[copyMod.m_gpoint.T]))+"K";
        string filename = "KMC_DEBUG/";
        //filename.append(Strings::cstr(runcount));
        filename = filename+Strings::cstr(seed)+"-0.00000_s__"+temp+".dot";
        //graphTitle = "0.00000s";
        m_simPAHp.saveDOT(filename);
        intervalcount = 0;
    }
    while(interval > intervalcount || time == time_max) {
        real timenow = intervalcount * X;
        int sec = (int) floor(timenow);
        int dec = (int) (floor((timenow - (real) sec)*100000));
        std::string dec_str;
        if(dec<10) dec_str = "0000";
        else if(dec<100) dec_str = "000";
        else if(dec<1000) dec_str = "00";
        else if(dec<10000) dec_str = "0";
        else dec_str = "";
        dec_str = dec_str.append(Strings::cstr(dec));
        string filename = "KMC_DEBUG/";
        //filename.append(Strings::cstr(runcount));
        //filename.append("_");
        copyMod.interpolateProfiles(timenow, false, 1);
        std::string temp = Strings::cstr((int)ceil(copyMod.m_gpoint.m_data[copyMod.m_gpoint.T]))+"K";
        filename = filename+Strings::cstr(seed)+"-"+Strings::cstr(sec)+"."+dec_str+"_s__"+temp+".dot";
        //graphTitle = Strings::cstr(sec)+"."+dec_str+"s";
        m_simPAHp.saveDOT(filename);
        intervalcount++;
        if(time==time_max) break;
    }
}

// CSV data ----
CSV_data::CSV_data(KMCSimulator& st) {
    m_sim = &st;
}
CSV_data::~CSV_data() {}
void CSV_data::initData(int max_runs, int no_of_interv, real max_time, intpair N_CH_initial, KMCGasph& gasph) {
    cout<<"Initialising CH_data vector...\n";
    // vector of zeros for each run
    intvector zeros(no_of_interv+1, 0);
    m_time.clear();
    m_T.clear();
    m_dataC.clear();
    m_dataH.clear();
    // calculate length of an interval
    m_dt = max_time/no_of_interv;
    // reference to gaspoint
    KMCGasPoint& gp = gasph.m_gpoint;
    // initialising time values
    for(int i=0; i<=no_of_interv; i++){
        real timetemp = m_dt*i;
        m_time.push_back(timetemp);
        gasph.interpolateProfiles(timetemp, false, 1);
        real Ttemp = gp.m_data[gp.T];
        m_T.push_back(Ttemp);
    }
    // setting all values to zero for all runs
    for(int i=0; i<max_runs; i++) {
        m_dataC.push_back(zeros);
        m_dataH.push_back(zeros);
        // to include t = 0
        m_dataC[i][0] = N_CH_initial.first;
        m_dataH[i][0] = N_CH_initial.second;
    }
    m_intervalcount = 0;
    
    cout<<"CH_data initialised!!\n";
}
// Compares time and adds data if interval reached
void CSV_data::addData(intpair N_CH, real time, int runNo, PAHProcess& pp, bool savedot) {
    // Zakwan's code
    int& c = m_intervalcount;
    // how many intervals have passed
    int interv_now = (int) floor(time/m_dt);
    // check if new interval has reached
    if(interv_now > c) {
        // to find how many data points skipped if time step larger than interval
        int intervals_jumped = 0;
        real Temp;
        if((interv_now-c) > 1) {
            
            intervals_jumped = interv_now - c -1;
            // fill skipped data points with last data point before jump
            for(int i=1; i<=intervals_jumped; i++) {
                Temp = m_T[c+i];
                m_dataC[runNo-1][c+i] = m_dataC[runNo-1][c+i-1];
                m_dataH[runNo-1][c+i] = m_dataH[runNo-1][c+i-1];
                if(savedot) {
                std::string filename = "KMC_DEBUG/";
                real timenow = (c+i)*m_dt;
                int sec = (int) floor(timenow);
                int dec = (int) (floor((timenow - (real) sec)*100000));
                std::string dec_str;
                if(dec<10) dec_str = "0000";
                else if(dec<100) dec_str = "000";
                else if(dec<1000) dec_str = "00";
                else if(dec<10000) dec_str = "0";
                else dec_str = "";
                dec_str = dec_str.append(Strings::cstr(dec));
                filename = filename+Strings::cstr(sec)+"."+dec_str+"s__"+Strings::cstr(Temp)+"K.dot";
                pp.saveDOT(filename);
                }
            }
            // update current data point
            c = interv_now;
            

            Temp = m_T[c];
            m_dataC[runNo-1][c] = N_CH.first;
            m_dataH[runNo-1][c] = N_CH.second;
            if(savedot){
            std::string filename = "KMC_DEBUG/";
            int sec = (int) floor(time);
            int dec = (int) (floor((time - (real) sec)*100000));
            std::string dec_str;
            if(dec<10) dec_str = "0000";
            else if(dec<100) dec_str = "000";
            else if(dec<1000) dec_str = "00";
            else if(dec<10000) dec_str = "0";
            else dec_str = "";
            dec_str = dec_str.append(Strings::cstr(dec));
            filename = filename+Strings::cstr(sec)+"."+dec_str+"s__"+Strings::cstr(Temp)+"K.dot";
                pp.saveDOT(filename);
            }
        }else {
            // update current data point
            c = interv_now;
            
            Temp = m_T[c];
            m_dataC[runNo-1][c] = N_CH.first;
            m_dataH[runNo-1][c] = N_CH.second;
            if(savedot) {
            std::string filename = "KMC_DEBUG/";
            int sec = (int) floor(time);
            int dec = (int) (floor((time - (real) sec)*100000));
            std::string dec_str;
            if(dec<10) dec_str = "0000";
            else if(dec<100) dec_str = "000";
            else if(dec<1000) dec_str = "00";
            else if(dec<10000) dec_str = "0";
            else dec_str = "";
            dec_str = dec_str.append(Strings::cstr(dec));
            filename = filename+Strings::cstr(sec)+"."+dec_str+"s__"+Strings::cstr(Temp)+"K.dot";
                pp.saveDOT(filename);
            }
        }
    }else return;
    /* Abhijeet's
    real c = m_intervalcount*m_dt;
    // how many intervals have passed
    //int interv_now = (int) floor(time/m_dt);
    if(time >= c) {
        m_dataC[runNo-1][m_intervalcount] = N_CH.first;
        m_dataH[runNo-1][m_intervalcount] = N_CH.second;
        m_intervalcount++;
    }*/
}
// delete data of run
void CSV_data::delData(int runNo) {
    int s = (int) m_dataC[runNo-1].size();
    for(int i=1; i!=s; i++) { // exclude first data point
        m_dataC[runNo-1][i] = 0;
        m_dataH[runNo-1][i] = 0;
    }
}
// Set name of output csv file containing C-H values
void CSV_data::setName(const std::string filename){
    m_name = filename;
    //std::string num = Strings::cstr(runNo);
    //m_name = m_name.append(num);
    //m_name = m_name.append(".csv");
}
// write CSV file of data. Option is given to write PAH data in columns
// or rows (since prior to Excel 2007 only 256 columns can be viewed in an Excel spreadsheet)
void CSV_data::writeCSV(bool col, bool keep_data) {
    CSV_IO csvfile(m_name, true);
    std::vector<std::string> tempLine;
    // writes data in columns
    if(col) {
        for(int i=0; i<(int)m_time.size(); i++) {
            tempLine.clear();
            tempLine.push_back(Strings::cstr(m_time[i]));
            tempLine.push_back(Strings::cstr(m_T[i]));
            //m_time.erase(m_time.begin());
            for(int runNo=0; runNo<(int) m_dataC.size(); runNo++) {
                tempLine.push_back(Strings::cstr(m_dataC[runNo][i]));
                tempLine.push_back(Strings::cstr(m_dataH[runNo][i]));
                //m_dataC[runNo].erase(m_dataC[runNo].begin());
                //m_dataH[runNo].erase(m_dataH[runNo].begin());
            }
            csvfile.Write(tempLine);
        }
    } else{ // writes data in rows
        std::vector<std::string> tempLine2;
        csvfile.Write(m_time);
        csvfile.Write(m_T);
        for(int runNo=0; runNo<(int) m_dataC.size(); runNo++) {
            tempLine.clear();
            tempLine2.clear();
            // CSV_IO class is not able to write ints (maybe an addition should be made to the template?)
            // so convert ints into strings first
            for(int i=0; i<(int) m_dataC[runNo].size(); i++) {
                tempLine.push_back(Strings::cstr(m_dataC[runNo][i]));
                tempLine2.push_back(Strings::cstr(m_dataH[runNo][i]));
            }
            csvfile.Write(tempLine);
            csvfile.Write(tempLine2);
        }
    }
    // option to keep C-H data after completed writing
    if(!keep_data) {
        m_time.clear();
        m_dataC.clear();
        m_dataH.clear();
    }
}

