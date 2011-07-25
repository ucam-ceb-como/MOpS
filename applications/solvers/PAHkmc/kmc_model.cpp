/*!
  * \Author     Zakwan Zainuddin (zz260)
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

#include "mt19937.h"
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
    KMCSimulator Simulator(std::string("gasphase.inp"),
		std::string("chem.inp"), std::string("therm.dat"));
	Simulator.TestGP();
	std::string input_name = "kmc.inx";

    clock_t timerStart = clock();
	// read xml input file
	CamXML::Document doc;
	const CamXML::Element *root, *node;
	const CamXML::Attribute *attr;
	// output CH counts in columns
	bool save_CH = false;
	bool CH_col = false;
	bool save_rates = false;
	bool save_dots = false;
	std::string CHOutput;
	int total_runs=0;
	Sweep::real t_start=0;
	Sweep::real t_end=0;
	int no_of_steps = 100;
	int no_of_rxnsteps = 500;
	StartingStructure startStruct = NONE;

	if (doc.Load(input_name) == 0) {
		root = doc.Root();
		// get simulation parameters
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
		if(node->Data() == "PYRENE") startStruct = PYRENE;
		else if(node->Data() == "BENZENE") startStruct = BENZENE;

		// get input file names
		node = root->GetFirstChild("gasFile");
		//KMC_Simulator.setCSVinputName(node->Data());

		// get output file names
		node = root->GetFirstChild("dotFiles");
		//KMC_Simulator.setDOToutputName(node->Data());
		attr = node->GetAttribute("mode");
		if(attr->GetValue() == "ON") save_dots = true;

		node = root->GetFirstChild("timer");
		//KMC_Simulator.setCSVtimerName(node->Data());

		node = root->GetFirstChild("rxncounts");
		//KMC_Simulator.setCSVreactioncountName(node->Data());

		node = root->GetFirstChild("pahlist");
		//KMC_Simulator.setCSVpahlistName(node->Data());

		node = root->GetFirstChild("rates");
		//KMC_Simulator.setCSVratesName(node->Data());
		attr = node->GetAttribute("mode");
		if(attr->GetValue() == "ON") save_rates = true;

		node = root->GetFirstChild("chcounts");
		CHOutput = node->Data();
		attr = node->GetAttribute("col");
		if(attr->GetValue() == "true") CH_col = true;
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
        PAHStructure pah;
        PAHProcess pahp(pah);
        pahp.initialise(PYRENE);
        int ID = 10000;
        Simulator.updatePAH(&pah,
            t_start, (t_end-t_start),
            no_of_steps,
            Sweep::genrand_int, Sweep::genrand_real1,
            1,
            ID);
        pahp.saveDOT("DOT files/TestSimulator.dot");
	}
	catch(std::runtime_error &re) {
		std::cout << re.what();
		return -1;
	}
    //KMC_Simulator.m_simPAHp.initialise(PYRENE);
    //KMC_Simulator.testSimulation(*(KMC_Simulator.m_simPAH),1000UL, 3000);

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
	//	timerStop = clock();
	//	if(secElapsed < (int)double(timerStop-timerStart)/CLOCKS_PER_SEC) {
	//		std::cout<<"Ending Program in "<<(countdown - ++secElapsed)<<" seconds..\n";
	//	}
	//}
    return 0;
}