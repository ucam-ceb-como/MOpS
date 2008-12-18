/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This class can be used to control the behavior of the solver
	based on the controls defined in the input file
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
#ifndef FL_SOLVER_CONTROL_H
#define FL_SOLVER_CONTROL_H
#include "fl_params.h"
#include "fl_reactor.h"
#include <map>

namespace FlameLab{
	class SolverControl{

		real aTol,rTol,resTol; 				// aobsolute, relative, adn residual tolerances
		real iniStep, maxStep, minStep;  // inital step and max allowed step size
		real maxTime, timeStep;				//maximum integration time
		int nSweep;								// number of sweeps in marching problem


	public:
		enum SolutionMode{
			steadyState,
			transient,
			preProcess
		};

		enum reportMode{
			Intermediate, //report after every time interval
			Final				//report only the final steady state result
		};

		enum method{
			DIRECT,
			KRYLOV
		};

		SolverControl(){
			setReportMode(Final);
		}
		~SolverControl(){}

		//set the absolute teolerance
		void setATol(real aTol); 

		//set the relative tolerenance
		void setRTol(real rTol);

		//set the residual teolerance
		void setResTol(real rsTol);

		//intial step size specification
		void setIniStep(real iniStep);

		//maximum allowed stepsize that the solver can take
		void setMaxStep(real maxStep);

		//minimum allowed step size
		void setMinStep(real minStep);

		//maximum integration time
		void setMaxTime(real maxTime);

		void setTimeStep(real timeStep);

		//set the solution mode
		void setSolMode(int mode);

		//set the report mode
		void setReportMode(int mode);

		//set the number of sweep
		void setNSweep(int n);

		
		void setOutInterval(std::map<real,real> interval, real dt);
		const int getIntervalCount() const;
		//return the output interval map for transient integration
		const map<map<real,real>,real>& getOutputInterval() const;
	
		//get the absolute tolerance
		real getATol()const;

		//get the relative tolerance
		real getRTol() const;

		//get the residual tolerance
		real getResTol() const;

		//return the initial step
		real getIniStep() const;

		//return the max step size
		real getMaxStep() const;

		//return the minimum allowes step size
		real getMinStep() const;

		//return the maximum integration time
		real getMaxTime() const;

		//return the time steps
		real getTimeStep() const;

		//return the solution mode (steady, transient)
		int getSolMode() const;

		//return the report mode (intermediate or final)
		int getReportMode() const;

		//return the number of sweep
		int getNSweep() const;

		//bool operator ==(SolutionMode sm) ;
		// set the solution method (either direct or krylov
		void setMethod(int n);

		// return the iteration method i.e direct or krylov
		int getMethod();



		//void initSolver(Reactor &reac, Sprog::Mechanism &mech, Sprog::Thermo::Mixture &mix);

	protected:
		int solMode,rptMode;//solution mode and report mode
		int method; // direct or krylov iteration method
		std::map<map<real, real>,real> outInterval;


	
	};
}

#endif
