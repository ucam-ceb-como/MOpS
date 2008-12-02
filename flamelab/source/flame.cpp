/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This is the main file for flameLab. This file drives both
	premix flame code as well as counter flow diffusion flame code
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

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "fl_io.h"
#include "fl_premix.h"
#include "fl_counter_diffusion.h"
#include "fl_solver_control.h"
#include "fl_error_handler.h"
#include "fl_params.h"
#include "gpc.h"

using namespace std;
using namespace FlameLab;
using namespace Sprog;

int main(){
	string fChem("chem.inp");
	string fThermo("therm.dat");
	string fTrans("tran.dat");
	string fFlame("flame.xml");
	
	Reactor *flame = NULL;
	SolverManager *solver = NULL;
	SolverControl *solverControl = new SolverControl();
	static Mechanism mech;
	
	
	try{
		// read chemkin input files
		IO::MechanismParser::ReadChemkin(fChem,mech,fThermo,fTrans);
		
		// create a reator object just for reading purpose		
		flame = new Reactor();
		// create flamelab object
		FlameLab::FlameLabIO *fio = new FlameLabIO();
		// read flamelab input file
		fio->readInput(fFlame,*flame,*solverControl);
		
		// create the mixture object
		Thermo::Mixture mix(mech.Species());
		//cout << boolalpha << (*flame==flame->PremixFlame) << endl;
		// if premix reactor create a premix object and assign to the base class reactor
		if( flame->getReactorModel() == flame->PremixFlame || flame->getReactorModel() == flame->Plug ){
			solver = new Premix(mech);
			//solver->initSolver(mech,mix,*solverControl,*flame);
			if(solverControl->getSolMode() == solverControl->preProcess)
				fio->writeGrid(*flame);
			else				
				solver->solve(mech,*solverControl,*flame,*fio);

		}else if(flame->getReactorModel() == flame->CDflmae) {
			solver = new CounterDiffusion(mech);
			if(solverControl->getSolMode() == solverControl->preProcess)
				fio->writeGrid(*flame);
			else
				solver->solve(mech,*solverControl,*flame,*fio);

		}

		
				
	}catch(ErrorHandler rh){
		cout << rh.errorString <<endl;
	}

	delete flame;
	delete solverControl;

	cout << "Flame code finished successfully\n";	
	return 0;
}