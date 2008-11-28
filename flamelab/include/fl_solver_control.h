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

		real aTol,rTol; // aobsolute and relative tolerances
		real iniStep, maxStep, minStep; // inital step and max allowed step size
		real maxTime, timeStep;
		


	public:
		enum SolutionMode{
			steadyState,
			transient
		};

		SolverControl(){}
		~SolverControl(){}
		void setATol(real aTol);
		void setRTol(real rTol);
		void setIniStep(real iniStep);
		void setMaxStep(real maxStep);
		void setMinStep(real minStep);
		void setMaxTime(real maxTime);
		void setTimeStep(real timeStep);
		void setSolMode(int mode);
		
		void setOutInterval(std::map<real,real> interval, real dt);
		const int getIntervalCount() const;
		//return the output interval map for transient integration
		const map<map<real,real>,real>& getOutputInterval() const;
	

		real getATol()const;
		real getRTol() const;
		real getIniStep() const;
		real getMaxStep() const;
		real getMinStep() const;
		real getMaxTime() const;
		real getTimeStep() const;
		int getSolMode() const;
		//bool operator ==(SolutionMode sm) ;



		//void initSolver(Reactor &reac, Sprog::Mechanism &mech, Sprog::Thermo::Mixture &mix);

	protected:
		int solMode;
		std::map<map<real, real>,real> outInterval;


	
	};
}

#endif