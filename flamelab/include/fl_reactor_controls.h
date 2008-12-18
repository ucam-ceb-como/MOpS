/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This class controls the operation of the reactor
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

#ifndef FL_REACTOR_CONTROLS_H
#define FL_REACTOR_CONTROLS_H
#include "fl_error_handler.h"
namespace FlameLab{
	class ReactorControls{
	public:
		enum SpaceToTime{
			ON,
			OFF
		}DiffusionSwitch;

		enum ReactorRunModel{
			NDS,	// no diffusion steady
			NDT,	// no diffusion transient
			WDS,	// with diffusion steady
			WDT		// with diffusion transient			
		};

		ReactorControls(){};
		~ReactorControls(){};

		//defines the change over from marching problem to transient problem
		void setSpaceToTime(int onOrOff);

		//return the present status of integration
		int getSpaceToTime();

		//set the reactor run model
		void setReactorRunModel(int n);

		// return the reactor run model
		int getReactorRunModel();

		// set the diffusion switch
		void setDiffusion(int n);

		// get the diffusion switch
		bool getDiffusion();


	protected:
		int spaceToTime;
		int runModel;
		int DiffnSwtch;	
	};
}

#endif
