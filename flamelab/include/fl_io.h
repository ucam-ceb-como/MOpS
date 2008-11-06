/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This class manages the reading of input file flame.xml
  Licence:
    This file is part of "flameLab".

    flameLab is free software; you can redistribute it and/or
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
#ifndef FL_IO_H
#define FL_IO_H
#include<string>
#include<vector>
#include "fl_reactor.h"
#include "fl_initial.h"
#include "fl_solver_control.h"
#include "fl_params.h"
#include "fl_premix.h"
#include "console_io.h"
#include "camxml.h"
#include "gpc_mech.h"
#include "gpc_mixture.h"
using namespace std;
namespace FlameLab{

	class Premix;//forward declaration

	class FlameLabIO{
		vector<string> monitor;
		vector<unsigned int> consoleMask;
		Console_IO flameLabConsole;
		int monitorSwitch;
	public:
		enum{
			ON,
			OFF
		};
		FlameLabIO(){}
		~FlameLabIO(){}
		void readInput(const string &fileName, Reactor &reac, SolverControl &solver);
		void readGeometry(Reactor &reac, const CamXML::Element &node);
		void readOPConditions(Reactor &reac, const CamXML::Element &node);
		void readInlet(Reactor &reac, const CamXML::Element &node);
		void readNozzleConditions(Reactor &reac, InitialConditions &nozzle, const CamXML::Element &node);
		void readSolverControl(SolverControl &solver, const CamXML::Element &node);
		void readInitialGuess(Reactor &reac, const CamXML::Element &node);
		void readMonitor(const CamXML::Element &node);
		void prepareConsole(Sprog::Mechanism &mech, FlameLab::Premix &flame);
		void writeToConsole(Reactor &reac) const;
		void setMonitorSwtich(int n);
		const int& getMonitorSwitch() const;
		
	};
}

#endif