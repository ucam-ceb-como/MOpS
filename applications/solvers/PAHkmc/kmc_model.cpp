/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_model.cpp
  *
  * \brief        Example application of the KMC model
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    A main file with implementations of KMC Model

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
#include <sstream>
#include "sweep.h"

#include "swp_kmc_jump_process.h"
#include "swp_kmc_mech.h"
#include "swp_kmc_pah_structure.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_reaction.h"
#include "swp_kmc_simulator.h"
#include "swp_kmc_structure_comp.h"
#include "swp_kmc_typedef.h"
#include "swp_kmc_gaspoint.h"
#include "swp_params.h"

#include "linear_interpolator.hpp"
#include <map>
#include <vector>
#include "csv_io.h"
#include "string_functions.h"
#include <time.h>
#include "camxml.h"
#include <stdexcept>

using namespace std;
using namespace Sweep::KMC_ARS;
using namespace Strings;

int main(int argc, char *argv[])
{
    //Simulator->TestGP();
    std::string input_name = "kmc.inx";

    std::string gasphase_filename = "gasphase.inp";
    clock_t timerStart = clock();
    // read xml input file
    CamXML::Document doc;
    const CamXML::Element *root, *node;
    const CamXML::Attribute *attr;
    // output CH counts in columns
    bool save_CH = false;
    std::string CH_mode = "C";
    std::string dot_mode = "END";
    std::string rates_mode = "1";
    bool save_rates = false;
    bool save_dots = false;
    std::string CHOutput;
    int total_runs=0;
    double t_start=0;
    double t_end=0;
    int no_of_steps = 100;
    int no_of_rxnsteps = 500;
    StartingStructure startStruct = NONE;
    std::string startStruct_str = "";
    int R6_num = 0;
    int R5_num_Lone = 0;
	int R5_num_Embedded = 0;
	int R7_num_Lone = 0;
	int R7_num_Embedded = 0;
	int num_C = 0;
	int num_H = 0;
	std::list<cpair> intCarbons;
    std::string site_savemode = "END";

    KMCSimulator* Simulator;

    std::string rate_filename;

    if (doc.Load(input_name) == 0) {
        root = doc.Root();
        // get simulation parameters
        node = root->GetFirstChild("gasFile");
        gasphase_filename = node->Data();
        Simulator = new KMCSimulator (gasphase_filename,
        std::string("chem.inp"), std::string("therm.dat"));

        node = root->GetFirstChild("runs");
        total_runs = (int) Strings::cdble(node->Data());
        
        node = root->GetFirstChild("steps");
        no_of_steps = (int) Strings::cdble(node->Data());

        node = root->GetFirstChild("rxntimesteps");
        no_of_rxnsteps = (int) Strings::cdble(node->Data());

        node = root->GetFirstChild("mintime");
        t_start = Strings::cdble(node->Data());

        node = root->GetFirstChild("maxtime");
        t_end = Strings::cdble(node->Data());

        node = root->GetFirstChild("startStruct");
        if(node->Data() == "BENZENE") startStruct = BENZENE_C;
		else if(node->Data() == "TOLUENE") startStruct = TOLUENE_C;
        else if(node->Data() == "NAPHTHALENE") startStruct = NAPHTHALENE_C;
        else if(node->Data() == "PYRENE") startStruct = PYRENE_C;
		else if(node->Data() == "METHYLPYRENE") startStruct = METHYLPYRENE_C;
		else if(node->Data() == "BENZOPYRENE") startStruct = BENZOPYRENE_C;
        else if(node->Data() == "CORONENE") startStruct = CORONENE_C;
        else if(node->Data() == "TEST") startStruct = TEST_STRUCT;
        else {
            startStruct_str = node->Data();
            attr = node->GetAttribute("R6");
            R6_num = (int) Strings::cdble(attr->GetValue());
			attr = node->GetAttribute("R5_Lone");
			R5_num_Lone = (int)Strings::cdble(attr->GetValue());
			attr = node->GetAttribute("R5_Embedded");
			R5_num_Embedded = (int)Strings::cdble(attr->GetValue());
			attr = node->GetAttribute("R7_Lone");
			R7_num_Lone = (int)Strings::cdble(attr->GetValue());
			attr = node->GetAttribute("R7_Embedded");
			R7_num_Embedded = (int)Strings::cdble(attr->GetValue());
			attr = node->GetAttribute("R7_Embedded");
			num_C = (int)Strings::cdble(attr->GetValue());
			attr = node->GetAttribute("R7_Embedded");
			num_H = (int)Strings::cdble(attr->GetValue());
        }

        // get input file names
        node = root->GetFirstChild("gasFile");
        Simulator->setCSVinputName(node->Data());

        // get output file names
        node = root->GetFirstChild("dotFiles");
        Simulator->setDOToutputName(node->Data());
        attr = node->GetAttribute("mode");
        if(attr->GetValue() == "ON") save_dots = true;
        attr = node->GetAttribute("save");
        dot_mode = attr->GetValue();

        node = root->GetFirstChild("timer");
        Simulator->setCSVtimerName(node->Data());

        node = root->GetFirstChild("rxncounts");
        Simulator->setCSVreactioncountName(node->Data());

        node = root->GetFirstChild("pahlist");
        Simulator->setCSVpahlistName(node->Data());
        attr = node->GetAttribute("save");
        if(attr->GetValue() == "STEP") site_savemode = attr->GetValue();

        node = root->GetFirstChild("rates");
        rate_filename = node->Data();
        Simulator->setCSVratesName(rate_filename);
        attr = node->GetAttribute("mode");
        if(attr->GetValue() == "ON") save_rates = true;
        attr = node->GetAttribute("save");
        rates_mode = attr->GetValue();

        node = root->GetFirstChild("chcounts");
        CHOutput = node->Data();
        attr = node->GetAttribute("save");
        if(attr->GetValue() == "CH" || attr->GetValue() == "C" || attr->GetValue() == "C_H" || attr->GetValue() == "MW") {
            CH_mode = attr->GetValue();}
        attr = node->GetAttribute("mode");
        if(attr->GetValue() == "ON") save_CH = true;

        //KMC_Simulator.loadGasProfiles();
    }else {
        std::cout << "ERROR: Cannot open "<<input_name<<", aborting simulation...\n";
        return -1;
    }

    // Starting structure: Pyrene
    //KMC_Simulator.initialise(PYRENE);
    // Run simulation
    
    try{
        //Simulator->TestGP();
        
        if(save_rates && rates_mode == "1") {
            Simulator->TestRates(t_start, t_end, no_of_steps);
            Simulator->TestConc(t_start, t_end, no_of_steps, std::string("Results/gasConcentrations.csv"));
            return 0;
        }
        std::vector<PAHStructure*> pah(total_runs);
        PAHProcess pahp;
        /*for(size_t i=0; i<pah.size(); ++i) {
            pah[i] = new PAHStructure();
            pahp.setPAH(*pah[i]);
            pahp.initialise(startStruct);
        }*/
        int ID = 10000;
        Simulator->initCSVIO();

        // Create a random number generator object with a fixed seed,
        // if this code is moved inside a loop, the seed will need to
        // be made unique to the loop.
        Sweep::rng_type rng(123);

        double step_size=(t_end-t_start)/no_of_steps;

        std::vector< std::vector< pair<int,int> > > PAH_CH;
        std::vector<double> PAH_time;
        
        if(save_CH) {
            PAH_time.push_back(0);
            for(int j=0; j<no_of_steps; j++) {
                PAH_time.push_back(t_start+(j*step_size)+step_size);
            }
        }
            
        for(size_t i=0; i<pah.size(); i++) {
            pah[i] = new PAHStructure;
            pahp.setPAH(*pah[i]);
            std::cout << "Starting growth on PAH "<<i<<"...\n";
            if(startStruct != NONE)
                pahp.initialise(startStruct);
            else
				pahp.initialise(startStruct_str, R6_num, R5_num_Lone, R5_num_Embedded, R7_num_Lone, R7_num_Embedded, num_C, num_H, intCarbons);
            //pahp.setPAH(*pah[0]);
            std::cout << "Pointer to PAH:"<<pah[i]<<"\n";
            Sweep::rng_type rng2(1+i);
            std::vector< pair<int,int> > CH_counts;
            // initialising CH counts
            if(save_CH) {
                CH_counts.push_back(pahp.getCHCount());
            }
            // initialising PAH_rate file, write header files
                

            for(int j=0; j<no_of_steps; j++) {
                double t_now = t_start+(j*step_size);
                if(save_rates && rates_mode == "R") {
                    std::vector<double> rates = Simulator->CurrentRates(pah[i],t_now);
                    Simulator->writeRatesCSV(t_now,
                        rates);
                }
				Simulator->updatePAH(pah[i],
					t_now, step_size,
					1,//no_of_steps,
					0,
					rng,
					1,
					ID + i,
					true,
					1.0);
                if(save_CH) {
                    CH_counts.push_back(pahp.getCHCount());
                }
                if(save_dots && dot_mode == "STEP") {
                    std::ostringstream dotname;
                    dotname << "DOT files/PAH_" << i <<'_'<< j <<".dot";
                    pahp.saveDOT(dotname.str());
                }
                if(i == 0 && site_savemode == "STEP") {
                    Simulator->writeCHSiteCountCSV();
                }
            }
            PAH_CH.push_back(CH_counts);
            std::cout<<"done!\n";
            Simulator->writeRxnCountCSV();
            Simulator->writeCHSiteCountCSV();
            if(save_dots) {
                std::ostringstream dotname;
                dotname << "DOT files/PAH" << i <<".dot";
                pahp.saveDOT(dotname.str());
            }
            std::cout<<"Done simulation for PAH "<<ID+i<<std::endl<<endl;
            std::cout<<"PAH Structure: "<<pahp.SiteString(',')<<std::endl;
            std::cout<<"PAH Ring Counts ~ R6:"<<std::get<0>(pahp.getRingsCount())<<" -- R5:"
				<< std::get<1>(pahp.getRingsCount()) << std::endl;
            std::cout<<"PAH Atom Counts ~ C:"<<pahp.getCHCount().first<<" -- H:"
                <<pahp.getCHCount().second<<std::endl;
            std::cout<<"No of Edge C ~ "<<pahp.CarbonListSize()<<std::endl;
            std::cout<<"No of Bridges ~ "<<pahp.numberOfBridges()<<std::endl;
            std::cout<<"No of Sites ~ "<<pahp.SiteListSize()<<std::endl;

            //pahp.clearStructure();
            
        }
        
        if(save_CH) {
            CSV_IO file_CH(CHOutput, true);
            std::vector<std::string> header;
            header.push_back("Time");
            for(size_t j=0; j<PAH_CH.size(); j++) {
                if(CH_mode == "CH" || CH_mode == "C" || CH_mode == "MW" || CH_mode == "C_H") {
                    std::ostringstream headC("PAH-");
                    headC << (1+j);
                    if(CH_mode != "MW")
                        headC << "_C";
                        header.push_back(headC.str());
                    if(CH_mode == "CH") {
                        std::ostringstream headH("PAH-");
                        headH << (1+j);
                        headH << "_H";
                        header.push_back(headH.str());
                    }
                }
            }
            if(CH_mode == "C_H") {
                for(size_t j=0; j<PAH_CH.size(); j++) {
                    std::ostringstream headH("PAH-");
                    headH << (1+j);
                    headH << "_H";
                    header.push_back(headH.str());
                }
            }
            file_CH.Write(header);
            for(size_t t=0; t<PAH_time.size(); t++) {
                std::vector<std::string> CHline;
                std::ostringstream dat;
                // first put in the time
                dat << PAH_time[t];
                CHline.push_back(dat.str());
                for(size_t j=0; j<PAH_CH.size(); j++) {
                    dat.str("");
                    if(CH_mode == "C" || CH_mode == "C_H" || CH_mode == "CH")
                        dat << PAH_CH[j][t].first;
                    if(CH_mode == "CH")
                        dat << ',' << PAH_CH[j][t].second;
                    if(CH_mode == "MW")
                        dat << (PAH_CH[j][t].first*12 + PAH_CH[j][t].second);
                    if(dat.str() != "")
                        CHline.push_back(dat.str());
                }
                dat.str("");
                if(CH_mode == "C_H")
                    for(size_t j=0; j<PAH_CH.size(); j++) {
                        dat.str("");
                        dat << PAH_CH[j][t].second;
                        CHline.push_back(dat.str());
                    }
                file_CH.Write(CHline);
            }
            file_CH.Close();
        }
        // delete all
        std::cout << endl;
        for(size_t i=0; i<pah.size(); i++) {
            /*PAHProcess pahp(*pah[i]);
            std::cout << "PAHpointer:" << pah[i] << endl;
            std::cout << "\tm_cpositions: " << pah[i]->m_cpositions.size() << endl;
            std::cout << "\tm_carbonList: " << pahp.CarbonListSize() << endl;
            int sum_ = 0;
            for(std::list<Site>::const_iterator k=pahp.SiteList().begin(); k!=pahp.SiteList().end(); k++) {
                switch(k->type) {
                case FE:
                case R5:
                    sum_++;
                    break;
                case ZZ:
                case RFE:
                    sum_+=2; break;
                case AC:
                case RZZ:
                case RFER:
                    sum_+=3; break;
                case BY5:
                case RAC:
                case RZZR:
                    sum_+=4; break;
                case BY6:
                case RBY5:
                case RACR:
                    sum_+=5; break;
                default:
                    std::cout<<"ERROR: invalid SITE TYPE!\n";
                    break;
                }
            }
            std::cout<< "\tnumber of site members: " << sum_ <<endl;*/
            delete pah[i];
        }
        pah.clear();
        
    }
    catch(std::runtime_error &re) {
        std::cout << re.what();
        return -1;
    }
    //KMC_Simulator.m_simPAHp.initialise(PYRENE);
    //KMC_Simulator.testSimulation(*(KMC_Simulator.m_simPAH),1000UL, 3000);
    delete Simulator;
    // Stop timer
    clock_t timerStop = clock();
    
    // Calculate total time elapsed
    double timerElapsed = double(timerStop-timerStart)/CLOCKS_PER_SEC;

    std::cout<<"\n\n\nEnding Program..\n";
    std::cout<<"Ran for "<<timerElapsed<<" seconds.\n\n";

    //// Ends program in 5 seconds
    //timerStart = clock();
    //int secElapsed = 0;
    //int countdown = 5;
    //std::cout<<"Ending Program in "<<countdown<<" seconds..\n";
    //while(secElapsed<countdown) {
    //    timerStop = clock();
    //    if(secElapsed < (int)double(timerStop-timerStart)/CLOCKS_PER_SEC) {
    //        std::cout<<"Ending Program in "<<(countdown - ++secElapsed)<<" seconds..\n";
    //    }
    //}
    return 0;
}
