/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_simulator.h
  *
  * \brief       Simulator for the kMC Model
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

#ifndef SWP_KMC_SIMULATOR_H
#define SWP_KMC_SIMULATOR_H

#include "swp_kmc_mech.h"
#include "swp_kmc_gaspoint.h"
#include "swp_kmc_pah_structure.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_typedef.h"
#include "swp_PAH_primary.h"
#include "swp_gas_profile.h"

#include "string_functions.h"
#include "csv_io.h"
#include <string>
#include <time.h>
#include <fstream>
#include <iostream>

//forward declaration
namespace Sprog {
    class Mechanism;
}

namespace Sweep{
    namespace KMC_ARS{
        class CSV_data;
        class KMCSimulator {
        public:
            friend class CSV_data;
            //! Constructor from a GasProfile object
            KMCSimulator(Sweep::GasProfile &gprofile);
            //! Constructor from chemkin and gasphase files
            KMCSimulator(const std::string gasphase, const std::string chemfile, const std::string thermfile);
            //! Copy Constructor
            KMCSimulator(KMCSimulator& s);
            //! Destructor
            virtual ~KMCSimulator();
            //! Initialise simulator from starting structure
            //void initialise(PAHStructure m_PAH);
            //! Initialise simulator
            //void initialise();
            // test gas profile
            void TestGP();
            //! Set PAH to be simulated
            void targetPAH(PAHStructure& pah);
            //! Set csv filename for gas profiles
            void setCSVinputName(const std::string& filename);
            //! Set output DOT file name "filename"_runs_finalLoopNum.dot
            void setDOToutputName(const std::string& filename);
            //! Set output CSV file name to keep track of timer counts
            void setCSVtimerName(const std::string& filename);
            //! Set output CSV file name to keep track of reaction counts
            void setCSVreactioncountName(const std::string& filename);
            //! Set output CSV file name to keep track of CH and site counts
            void setCSVpahlistName(const std::string& filename);
            //! Set output CSV file name to keep track of rates for one run
            void setCSVratesName(const std::string& filename);
            //! Writes data for timeCount.csv
            void writeTimerCSV(const int& loop, const double& elapsedTime);
            //! Writes data for reaction_count.csv
            void writeRxnCountCSV();
            //! Writes data for CH_site_list.csv
            void writeCHSiteCountCSV();
            //! Writes data for rates count (csv)
            void writeRatesCSV(double& time, rvector& v_rates);
            //! Initialise CSV_IOs
            void initCSVIO();
            //! Initialise reaction count
            void initReactionCount();
            //! Reads chemical mechanism / profile (if not obtained from Mops)
            void LoadGasProfiles(const std::string gasphase, const std::string chemfile, const std::string thermfile);
            //! Write column headings for CSV files
            void writeCSVlabels();
            //! Save the structure DOT file after every X loops
            void saveDOTperXLoops(int X, int& loopcount, int& runcount);
            //! Save the structure DOT file for particular PAH (ID)
            void saveDOTperLoop(int LOOPcount,int loopcount, int PAH_ID);
            //! Save the structure DOT file after every X simulation sec interval
            void saveDOTperXsec(const double& X, const int& seed, const double& time, const double &time_max, KMCMechanism& copyMod, int& intervalcount);
            //! Update structure of PAH after time dt.
            double updatePAH(PAHStructure* pah,         //! structure of pah.
                           const double tstart,       //! start time.
                           const double dt,           //! growth time.
                           const int waitingSteps,    //! waiting step used to calculate maximum time interval, currently use 1.
						   const int maxloops,        //! maximum number of loops to take. 0 means no limit
                           rng_type &rng,             //! random number generator.
                           double r_factor,           //! growth factor g, one important parameter used in this model.
                           int PAH_ID,                //! ID of this pah, used for debugging.
						   bool calcrates,
						   double ratefactor);        
            //! Outputs rates into a csv file (assuming all site counts as 1)
            void TestRates(const double tstart, const double tstop, const int intervals);
            //! Obtains rates of PAH reactions with the current structure
            rvector CurrentRates(PAHStructure* pah, double t);
            //! Outputs gas concentrations into a csv file
            void TestConc(const double& t_start, const double& t_stop, const int intervals, const std::string& filename);
            //void testSimulation(PAHStructure& pah, const unsigned long seed, int totalruns);
            void writetimestep(const std::vector<double>& timestep);
            void setCSVtimestep(const std::string &filename);
            std::string    m_timestep_name;
            //! returns the reaction counts
            std::vector<int> ReactionCounts() const; 
            //! CSV input filename
             std::string m_csv_in;
            //! DOT output filename
             std::string m_dot_out;
            //! CSV timer filename
             std::string m_timer_name;
            //! CSV reaction count filename
             std::string m_rxncount_name;
            //! CSV CH and site counts filename
             std::string m_pahlist_name;
            //! CSV rates counts filename (for one run)
             std::string m_rates_name;
            //! CSV io object for timer
             CSV_IO m_timer_csv;
            //! CSV io object for reaction count
             CSV_IO m_rxn_csv;
            //! CSV io object for CH and site counts
             CSV_IO m_pah_csv;
            //! CSV io object for rates counts (for one run)
             CSV_IO m_rates_csv;
            //! CSV io object for time step
             CSV_IO m_timestep_csv;
        private:
            //! Make default constructor private
            KMCSimulator();
            //! Pointer to gasprofile object
            Sweep::GasProfile* m_gasprof;
            //! Pointer to Sprog mechanism
            Sprog::Mechanism* m_mech;
            //! Gaspoint object
            KMCGasPoint* m_gas;
            //static int LOOPcount;
            //! The PAH data structure
            PAHStructure* m_simPAH;
            //! Current simulation time
            double m_t;
            //! Check if profile obtained from file rather than mops
            bool m_fromfile;
            //! KMC Mechanism
            KMCMechanism m_kmcmech;
            //! PAH process
            PAHProcess m_simPAHp;
            //! Reaction Counter
            std::vector<int> m_rxn_count;
        };

        // to write C-H data for PAHs after intervals
        class CSV_data {
        public:
            CSV_data(KMCSimulator& st);
            virtual ~CSV_data();
            // initialise data storage
            void initData(int max_runs, int no_of_interv, double max_time, intpair N_CH_initial, KMCGasPoint& gp);
            // Compares time and adds data if interval reached
            void addData(intpair N_CH, double time, int runNo, PAHProcess& pp, bool savedot);
            // delete data of run
            void delData(int runNo);
            // Set name of output csv file containing C-H values
            void setName(const std::string filename);
            // write CSV file of data. Option is given to write PAH data in columns
            // or rows (since prior to Excel 2007 only 256 columns can be viewed in a spreadsheet)
            void writeCSV(bool col, bool keep_data);

            // pointer to Simulator to extract data from
            KMCSimulator* m_sim;
            // name of output file
            std::string m_name;
            // stores C and H data
            std::vector<intvector> m_dataC;
            std::vector<intvector> m_dataH;
            // stores time
            rvector m_time;
            // stores temperature
            rvector m_T;
            // specifies the number of times data should be saved for whole simulation
            int m_intervalcount;
            // the length of an interval (max time / interval count)
            double m_dt;
        };
    }
}

#endif
