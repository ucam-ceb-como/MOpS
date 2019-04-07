/*!
  * \author     Zakwan Zainuddin (zz260)
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

#include "gpc_mech.h"
#include "gpc_mech_io.h"

#include "csv_io.h"

#include <string>
#include <iostream>
#include <time.h>
#include <math.h>
#include <cstdlib>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>

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
static std::vector<kmcSiteType> allSiteType = vectSiteType();

//! Default Constructor
KMCSimulator::KMCSimulator():
     m_gasprof(), m_mech(), m_gas(), m_simPAH(), m_t(), m_fromfile(false), m_kmcmech(),m_simPAHp()
{
}

//! Constructor from chemkin and gasphase files
KMCSimulator::KMCSimulator(const std::string gasphase, const std::string chemfile, const std::string thermfile)
{
    m_t=0;
    m_gasprof = new Sweep::GasProfile();
    LoadGasProfiles(gasphase, chemfile, thermfile);
    m_fromfile = true;
}
//! Constructor from a GasProfile object
KMCSimulator::KMCSimulator(Sweep::GasProfile& gprofile):
	m_gasprof(), m_mech(), m_gas(), m_simPAH(), m_t(0.0), m_fromfile(false), m_kmcmech(), m_simPAHp()
{
    std::cout << this << endl;
    m_gasprof = &gprofile;
    m_gas = new KMCGasPoint(gprofile, *gprofile[0].Gas.Species());
    m_mech = NULL;
}

//! Copy Constructor
KMCSimulator::KMCSimulator(KMCSimulator& s):
		m_gasprof(), m_mech(), m_gas(), m_simPAH(), m_t(s.m_t), m_fromfile(false),
		m_kmcmech(s.m_kmcmech),m_simPAHp()

{
    m_gasprof = s.m_gasprof;
    m_gas = new KMCGasPoint(*s.m_gas);
    //m_simPAH = new PAHStructure(*(s.m_simPAH));
    m_simPAHp = PAHProcess(*m_simPAH);
}

//! Default Destructor
KMCSimulator::~KMCSimulator() {
    //delete m_simPAH;
    if(m_fromfile) {
        delete m_gasprof;
        delete m_mech;
    }
    delete m_gas;
}

//! Set PAH to be simulated
void KMCSimulator::targetPAH(PAHStructure& pah) {
    m_simPAH = &pah;
    m_simPAHp = PAHProcess(*m_simPAH);
}

/*!
 * @param[in,out]    pah             PAH structure KMC-ARS jump process will be performed on.
 * @param[in]        tsart           The latest time the PAH was updated.
 * @param[in]        dt              The time over which the PAH is to be updated.
 * @param[in]        waitingSteps    Adjusts the size of dt. The larger it is the smaller dt is. Results were found to be insensitive to waitingSteps therefore it was hardcoded as 1 (tests by Dongping, dc516@cam.ac.uk).
 * @param[in]        rng             Random number generator.
 * @param[in]        r_factor        Model parameter: Growth factor, a multiplier that is applied to the growth rate of PAHs within primary particles when the number of PAHs exceeds a critical number of PAHs.
 * @param[in]        PAH_ID          "Unique" identification number attached to this PAH.
 */
double KMCSimulator::updatePAH(PAHStructure* pah, 
                            const double tstart, 
                            const double dt,  
                            const int waitingSteps,  
							const int maxloops,
                            rng_type &rng,
                            double r_factor,
                            int PAH_ID,
							bool calcrates,
							double ratefactor) {
	// wjm34: remove call to initReaction count to save time in updating.
    //initReactionCount();
    m_t = tstart;
    double t_max = m_t + dt;
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
    double t_next = m_t;
    double t_step_max = dt/waitingSteps;
    //double oldtnext;
    int loopcount=0;
	bool proceed = true;
	Spointer Sp1;
	int test;
    while (m_t < t_max && proceed) {
        //this->m_simPAHp.printStruct();// print out structure of this pah on the screen
        //m_simGas.interpolateProfiles(m_t, true, r_factor);
        loopcount++;

        // Calculate rates of each jump process
		if (calcrates){
			m_gas->Interpolate(m_t, r_factor);
			m_kmcmech.calculateRates(*m_gas, m_simPAHp, m_t);
		}

        // Calculate time step, update time
        typedef boost::exponential_distribution<double> exponential_distrib;
        exponential_distrib waitingTimeDistrib(m_kmcmech.TotalRate()*ratefactor);
        boost::variate_generator<rng_type &, exponential_distrib> waitingTimeGenerator(rng, waitingTimeDistrib);
        double t_step = waitingTimeGenerator();
        t_next = m_t+t_step;
        if(t_next < t_max && t_step < t_step_max) {

			//if (PAH_ID == 677){
			//	std::list<Site> tester = pah->GetSiteList();
			//	cout << "Check start " << PAH_ID << endl << t_next << endl;
			//	cout << "Rings = "<< pah->numofRings() << endl;;
			//	for (Sp1 = tester.begin(); Sp1 != tester.end(); ++Sp1){
			//		cout << (int)(Sp1->type) << " " << (int)(Sp1->comb) << endl;
			//	}
			//}

			//cout << "Check start " << PAH_ID << endl << t_next << endl;

            //if (pah->numofC()>5000&&pah->numofC()<6000)//||pah->havebridgeC()    //if (PAH_ID==224835)
            //if we want to check a PAH with specified ID or number of Carbon, 
            //its structure can be drawed by this function, used for tracking suspicious PAH.
            // saveDOTperLoop(100000*tstart,loopcount,PAH_ID);
            // Choose jump according to rates
            ChosenProcess jp_perf = m_kmcmech.chooseReaction(rng);
            //cout<<jp_perf.first->getName()<<'\n';//++++
            //std::ostringstream dotname;
            //dotname << "KMC_DEBUG/p_" << (t_next*1e8) << ".dot";
            //m_simPAHp.saveDOT(dotname.str());
            // Update data structure
            m_simPAHp.performProcess(*jp_perf.first, rng, PAH_ID);

			//Hard cut-off for PAHs. Cannot have less than one ring. Set number of carbons to 1 so that it will be invalidated
			//Set t_next to t_max so the updatePAH routine will be exited
			//if (pah->numofRings() < 4){
			//	proceed = false;
			//	//t_next = t_max;
			//}

			//if (PAH_ID == 677){
			//	std::list<Site> tester = pah->GetSiteList();
			//	cout << "Check end " << PAH_ID << endl << t_next << endl;
			//	cout << "Rings = " << pah->numofRings() << endl;;
			//	for (Sp1 = tester.begin(); Sp1 != tester.end(); ++Sp1){
			//		cout << (int)(Sp1->type) << " " << (int)(Sp1->comb) << endl;
			//	}
			//}

			//cout << "Check end " << PAH_ID << endl << t_next << endl;

            /*if(m_simPAH->m_parent->ID() % 100000 == 609) {
            if(!m_simPAHp.checkCoordinates()) {
                cout<<"ERROR: Invalid coordinates produced after performing process "
                    <<jp_perf->getName()<<" (ID" <<jp_id<<")\n";
            }
            }*/
            /*
             * wjm34: remove updating of jump counts to save time
            if(process_success) {
                m_rxn_count[jp_perf.second]++;
            }*/
        }else {
            //oldtnext = t_next;
            t_next = t_max;
        }
		if (loopcount == maxloops) proceed = false; //If maxloops is set to 0, this condition will never be true
        m_t = t_next;
    }
	return m_t;
}

//! Outputs rates into a csv file (assuming all site counts as 1)
void KMCSimulator::TestRates(const double tstart, const double tstop, const int intervals) {
    // set name of output file
    //setCSVratesName(filename);
    std::cout << "Saving Rates...\n";
    m_rates_csv.Open(m_rates_name, true);
    rvector rates(m_kmcmech.JPList().size(), 0);
    double dt = (tstop-tstart)/intervals;
    m_simPAHp.m_rates_save = true;
    // for each interval
    for(double t=tstart; t<= tstop; t+=dt) {
        // interpolate & calculate rates
        m_gas->Interpolate(t);
        m_kmcmech.calculateRates(*m_gas, m_simPAHp, t);
        rates = m_kmcmech.Rates();
        writeRatesCSV(t, rates);
    }
    m_simPAHp.m_rates_save = false;
    std::cout<<"Finished calculating rates for kMC mechanism. Results are saved in "
        <<m_rates_name<<"\n\n";
}

//! Obtains rates of PAH reactions with the current structure
rvector KMCSimulator::CurrentRates(PAHStructure* pah, double t) {
    m_simPAHp.setPAH(*pah);
    rvector rates(m_kmcmech.JPList().size(), 0);
    m_gas->Interpolate(t);
    m_kmcmech.calculateRates(*m_gas, m_simPAHp, t);
    rates = m_kmcmech.Rates();
    return rates;
}
//! Outputs gas concentrations into a csv file
void KMCSimulator::TestConc(const double& t_start, const double& t_stop, const int intervals, const std::string& filename) {
    CSV_IO csvio(filename, true);
    double dt = (t_stop-t_start)/intervals;
    std::vector<string> species(m_gas->m_total-2);
    rvector temp(m_gas->m_total-2); //exclude T & P
    for(size_t i=1; i<(temp.size()); i++) {
        species[i] = m_gas->SpNames()[i+1];
    }
    csvio.Write(species);
    for(double t=t_start; t<=t_stop; t+=dt) {
        temp[0] = t;
        m_gas->Interpolate(t);
        for(size_t i=1; i<(temp.size()-1); i++) {
            temp[i] = (*m_gas)[i+1];
        }
        csvio.Write(temp);
        
    }
}

/*void KMCSimulator::testSimulation(PAHStructure& pah, const unsigned long seed, int totalruns) {
    for(int run=1; run<=totalruns; run++) {
        clock_t timerStart = clock();
        int total_interval = 10;
        double maxtime = 0.0222029;
        double step = maxtime/total_interval;
        double tnow=0;
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
void KMCSimulator::writeRxnCountCSV() {
    // change int vector to float
    std::vector<std::string> temp;
    for(size_t i=0; i<m_rxn_count.size(); i++)
        temp.push_back(Strings::cstr(m_rxn_count[i]));
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
void KMCSimulator::writeRatesCSV(double& time, rvector& v_rates) {
    int convC[] = {1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21};
    //int convM[] = {1, 2,14,15, 8,10,11,12,13, 3, 7, 9, 5, 4, 6,22,24,16,21};
    int ID;
    //if(runNo==1) {
    int total_jp = 21;
    std::vector<double> temp(total_jp+1,0);
    temp[0] = time;
    for(int i=0; i!=(int)v_rates.size(); i++) {
        ID = m_kmcmech.JPList()[i]->getID();
        temp[convC[ID-1]] = v_rates[i];
    }
    m_rates_csv.Write(temp);
    //}
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
    m_rxn_count.assign(m_kmcmech.JPList().size(), 0);
}
//! Reads chemical mechanism / profile (if not obtained from Mops)
//! Similar function as Sweep::Flamesolver::LoadGasProfile
void KMCSimulator::LoadGasProfiles(const std::string gasphase, const std::string chemfile, const std::string thermfile) {
    m_mech = new Sprog::Mechanism();
    Sprog::IO::MechanismParser::ReadChemkin(chemfile, *m_mech, thermfile, 0);
    map<double,double> alpha_prof;

    // Clear the current gas-phase profile.
    m_gasprof->clear();

    // Open the file to read.
    ifstream fin(gasphase.c_str(), ios::in);
    if (!fin.good()) {
        throw runtime_error("Unable to open gas profile input "
                            "file (Sweep, KMCSimulator::LoadGasProfile).");
    }

    // Variables read from file.
    vector<string> subs;
    string delim = ",\t \r"; // Possible file delimiters (comma, tab and space).
    string line;

    // Get the first line (should be the header defining the columns).
    if (!getline(fin, line).eof()) {

        // Split the line to get the column headings.
        split(line, subs, delim);

        // Get important column indices (time, temperature and pressure).
        int tcol=-1, Tcol=-1, Pcol=-1;
        tcol = findinlist(string("Time"), subs);

        Tcol = findinlist(string("T"), subs);
        if(Tcol < 0)
            Tcol = findinlist(string("T[K]"), subs);

        Pcol = findinlist(string("P"), subs);

        // Columns to ignore, but which are useful to have in files for brush compatibility
        int Xcol = findinlist(string("X[cm]"), subs);
        int Dcol = findinlist(string("RHO[g/cm3]"), subs);
        int Vcol = findinlist(string("V[cm/s]"), subs);
        int Gcol = findinlist(string("GradT"), subs);
        int Acol = findinlist(string("Alpha"), subs);
        int Rcol = findinlist(string("wdotA4"), subs);


        // Check that the file contains required columns.
        if (tcol < 0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no Time "
                                "column (Sweep::KMCSimulator::LoadGasProfiles).");
        }
        if (Tcol < 0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no temperature "
                                "column (Sweep::KMCSimulator::LoadGasProfiles).");
        }
        if (Pcol < 0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no pressure "
                                "column (Sweep::KMCSimulator::LoadGasProfiles).");
        }

        // All other columns are chemical species.  Add them, and note
        // their columns.
        map<unsigned int,int> spcols;
        for (int i=0; (unsigned)i!=subs.size(); ++i) {
            if ((i!=tcol) && (i!=Tcol) && (i!=Pcol) && (i!=Acol) && (i!=Rcol) &&
                (i!=Xcol) && (i!=Dcol) && (i!=Vcol) && (i!=Gcol)) {
                // Try to find this species in the mechanism
                const int speciesMechIndex = m_mech->FindSpecies(subs[i]);

                if(speciesMechIndex < 0) {
                    std::ostringstream msg("Failed to find species ");
                    msg << subs[i] << " in mechanism (Sweep::KMCSimulator::LoadGasProfiles).";
                    throw std::runtime_error(msg.str());
                }
                // Found species
                spcols[i] = speciesMechIndex;
            }
        }

        // Now we can read the profile.
        while(!getline(fin, line).eof()) {
            // Set t=0 and create a new IdealGas object.
            double t = 0.0;
            double T = 0.0;
            double P = 0.0;
            double alpha = 0.0;
            double PAHRate = 0.0;
            GasPoint gpoint(m_mech->Species());

            // Split the line by columns.
            split(line, subs, delim);

            // Check that the mole fractions sum to 1
            double checkSum = 0.0;

            // Loop over all the elements in the line and save them
            // to the correct gas-phase variable.
            for (int i=0; (unsigned)i!=subs.size(); ++i) {
                if (i==tcol) {
                    // This is the Time column.
                    t = cdble(subs[i]);
                    gpoint.Time = t;
                } else if (i==Tcol) {
                    // This is the temperature column.
                    T = cdble(subs[i]);
                } else if (i==Pcol) {
                    // This is the pressure column.
                    P = cdble(subs[i]);
                } else if (i==Acol) {
                    alpha = cdble(subs[i]);
                } else if (i==Rcol) {
                    PAHRate = cdble(subs[i]);
                } else {
                    // This is a gas-phase species column.
                    map<unsigned int,int>::iterator isp = spcols.find(i);
                    if (isp != spcols.end()) {
                        const double frac = cdble(subs[i]);
                        assert(isp->second >= 0);
                        gpoint.Gas.RawData()[isp->second] = frac;
                        checkSum += frac;
                    }
                }
            }

            if((checkSum < 0.997) || checkSum > 1.003) {
                std::ostringstream msg;
                msg << "Mole fractions sum to " << checkSum
                    << ", but should sum to 1.000 (KMCSimulator::LoadGasProfiles)";
                throw std::runtime_error(msg.str());
            }

            // Set up the gas-phase by setting temperature, pressure and
            // normalising the mixture fractions.
            // TODO:  This will give the wrong component densities
            //        unless all species are specified!
            gpoint.Gas.SetTemperature(T);
            gpoint.Gas.SetPressure(P*1.0e5);
            gpoint.Gas.Normalise();
            gpoint.Gas.SetPAHFormationRate(PAHRate*1E6);

            // Add the profile point.
            alpha_prof[t] = alpha;
            m_gasprof->push_back(gpoint);
        }

        // Close the input file.
        fin.close();

        // Sort the profile by time.
        SortGasProfile(*m_gasprof);
        m_gas = new KMCGasPoint(*m_gasprof, m_mech->Species());
    } else {
        // There was no data in the file.
        fin.close();
        throw runtime_error("Input file contains no data "
                            "(Sweep::KMCSimulator::LoadGasProfiles).");
    }
}
//! Write column headings for CSV files
void KMCSimulator::writeCSVlabels() {
    // Write headings for timer
    std::vector<string> timer_headings;
    timer_headings.push_back("Total Loops"); timer_headings.push_back("Time Elapsed");
    m_timer_csv.Write(timer_headings);
    // Write headings for reaction count
    std::vector<string> rxn_headings;
    for(size_t i=0; i<m_kmcmech.JPList().size(); i++) {
        // gets name of each jump process and puts them in a row
        rxn_headings.push_back(m_kmcmech.JPList()[i]->getName());
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
    // Write headings for rate file
    std::vector<std::string> rates_header;
    rates_header.push_back("Time");
    for(size_t i=0; i<m_kmcmech.JPList().size()+1; i++) {
        std::ostringstream ID_name("ID");
        ID_name << (i+1);
        rates_header.push_back(ID_name.str());
    }
    m_rates_csv.Write(rates_header);
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
void KMCSimulator::saveDOTperXsec(const double& X, const int& seed, const double &time, const double &time_max, KMCMechanism& copyMod, int& intervalcount) {
    int interval = (int) ceil(time/X);
    std::string graphTitle;
    if(intervalcount == -1) {
        m_gas->Interpolate(0, 0);
        std::string temp = Strings::cstr((int)ceil((*m_gas)[m_gas->T]))+"K";
        string filename = "KMC_DEBUG/";
        //filename.append(Strings::cstr(runcount));
        filename = filename+Strings::cstr(seed)+"-0.00000_s__"+temp+".dot";
        //graphTitle = "0.00000s";
        m_simPAHp.saveDOT(filename);
        intervalcount = 0;
    }
    while(interval > intervalcount || time == time_max) {
        double timenow = intervalcount * X;
        int sec = (int) floor(timenow);
        int dec = (int) (floor((timenow - (double) sec)*100000));
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
        m_gas->Interpolate(timenow, 0);
        std::string temp = Strings::cstr((int)ceil((*m_gas)[m_gas->T]))+"K";
        filename = filename+Strings::cstr(seed)+"-"+Strings::cstr(sec)+"."+dec_str+"_s__"+temp+".dot";
        //graphTitle = Strings::cstr(sec)+"."+dec_str+"s";
        m_simPAHp.saveDOT(filename);
        intervalcount++;
        if(time==time_max) break;
    }
}

// CSV data ----
CSV_data::CSV_data(KMCSimulator& st):
		m_sim(NULL), m_name(), m_dataC(), m_dataH(), m_time(), m_T(),
		m_intervalcount(0), m_dt(0.0) {
    m_sim = &st;
}
CSV_data::~CSV_data() {}
void CSV_data::initData(int max_runs, int no_of_interv, double max_time, intpair N_CH_initial, KMCGasPoint& gp) {
    cout<<"Initialising CH_data vector...\n";
    // vector of zeros for each run
    intvector zeros(no_of_interv+1, 0);
    m_time.clear();
    m_T.clear();
    m_dataC.clear();
    m_dataH.clear();
    // calculate length of an interval
    m_dt = max_time/no_of_interv;
    // initialising time values
    for(int i=0; i<=no_of_interv; i++){
        double timetemp = m_dt*i;
        m_time.push_back(timetemp);
        gp.Interpolate(timetemp);
        double Ttemp = gp[gp.T];
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
void CSV_data::addData(intpair N_CH, double time, int runNo, PAHProcess& pp, bool savedot) {
    // Zakwan's code
    int& c = m_intervalcount;
    // how many intervals have passed
    int interv_now = (int) floor(time/m_dt);
    // check if new interval has reached
    if(interv_now > c) {
        // to find how many data points skipped if time step larger than interval
        int intervals_jumped = 0;
        double Temp;
        if((interv_now-c) > 1) {
            
            intervals_jumped = interv_now - c -1;
            // fill skipped data points with last data point before jump
            for(int i=1; i<=intervals_jumped; i++) {
                Temp = m_T[c+i];
                m_dataC[runNo-1][c+i] = m_dataC[runNo-1][c+i-1];
                m_dataH[runNo-1][c+i] = m_dataH[runNo-1][c+i-1];
                if(savedot) {
                std::string filename = "KMC_DEBUG/";
                double timenow = (c+i)*m_dt;
                int sec = (int) floor(timenow);
                int dec = (int) (floor((timenow - (double) sec)*100000));
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
            int dec = (int) (floor((time - (double) sec)*100000));
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
            int dec = (int) (floor((time - (double) sec)*100000));
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
    double c = m_intervalcount*m_dt;
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

// Test KMCGasPoint object. Linear interpolates values at 5 points from t = 0 to 0.005s
// and writes values of profile on console.
void KMCSimulator::TestGP() {
    cout<<endl<<"---(Sweep, KMC_ARS::KMCSimulator) Testing KMCGasPoint---"<<endl;
    for(double t=0; t<0.005; t+=0.001) {
        cout << "--At time "<<t <<"--"<<endl;
        m_gas->Interpolate(t);
        for(int i=0; i<m_gas->m_total; i++) {
            cout<<m_gas->SpNames()[i]<<"\t"<<(*m_gas)[i]<<endl;
        }
        cout<<endl;
    }
    cout<<"---(Sweep, KMC_ARS::KMCSimulator) Finished testing..."<<endl<<endl;
}
