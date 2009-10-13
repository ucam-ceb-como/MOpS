/*
  Author(s):      Markus Sander (ms785)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2009 Markus Sander.

  File purpose:


  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
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

#ifndef SWEEP_TRAJECTORY_H
#define SWEEP_TRAJECTORY_H


#include <map>
#include <vector>
#include <fstream>
#include <string>

namespace Sweep
{
class Trajectory 
{
public:
    // Constructors.
    Trajectory(void); // Default constructor.

    // Destructors.
    virtual ~Trajectory(void); // Default destructor.

	void LoadPAHProfile(const std::string &file);

	struct trajectory_base {
		std::vector<double> time;
		std::vector<int> n_carbon_t;
    };
	std::vector<trajectory_base> alltrajectories;

    int maxID();

    double StartTime();


private:
	// the id of the PAH. This number is increased by 
	//int ID;
    int m_maxID;
    double m_starttime;
	
};
};

#endif
