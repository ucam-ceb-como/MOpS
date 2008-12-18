/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This is an abstract class which must be implemented by both premix
	and couterflow flame class. 
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

#ifndef FL_SOLVER_MANAGER_H
#define FL_SOLVER_MANAGER_H
#include "gpc.h"
#include "fl_params.h"
#include "fl_reactor.h"
#include "fl_solver_control.h"
//#include "fl_io.h"
#include "nvector/nvector_serial.h"
#include "cvode_impl.h" // For CVodeMem.
#include "cvode_band_impl.h" // For band mat.
#include "cvode_dense_impl.h"
//#include<cvode/cvode_spgmr.h>
//#include<cvode/cvode_bandpre.h>
// kinsol header files
//#include "kinsol_impl.h"
//#include "kinsol_band_impl.h"
//#include <kinsol/kinsol_spgmr.h>
#include <vector>
#include <map>
namespace FlameLab{
	class FlameLabIO;
	class SolverManager{

	protected:
		void *cvode_mem, *kin_mem, *bp_data; // memory space
		//real *ptrToSlnVector; // pointer to solution vector
		real tMax,timeStep, currentTime;
		int currentLoop;
		

		N_Vector solVect; // solution vector
		N_Vector derivative; // derivatives for residual evaluation
		N_Vector templ; // work space for Newton
		vector<real> variables;
		vector<real> mcVariables;//variables for multi cell

		typedef struct{
			int lower;
			int upper;
		}JacBounds;

		map<int,JacBounds> JacInfo;

	public:

		virtual void initSolver(SolverControl &sc,Reactor &reac)=0;
		virtual void solve(Sprog::Mechanism &mech, SolverControl &sc, Reactor &reac, FlameLabIO &io)=0;


	};
}
#endif
