/*
  Author(s):      Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Markus Sander.

  File purpose:
    Implementation of the PAHPrimary class declared in the
    swp_PAH_primary.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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

#include "swp_PAH_trajectory.h"
#include <stdexcept>
#include "csv_io.h"
#include <string>
#include "string_functions.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>


using namespace Sweep;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.
Trajectory::Trajectory() 
{
}


// Default destructor.
Trajectory::~Trajectory()
{
}

int Trajectory::maxID() const
{
    return m_maxID;
}

double Trajectory::StartTime() const
{
    return m_starttime;
}

void Trajectory::LoadPAHProfile(const std::string &file)
{
   	CSV_IO *csvinput = new CSV_IO(file,false);
	std::vector<std::string> values;	
	values.clear();
	csvinput->Read(values);
    m_maxID=values.size()/2-1;
	for (int j=0;j<m_maxID;j++)
	{
		trajectory_base newtraj;
		alltrajectories.push_back(newtraj);
	}
	for (int j=0;j<10000;j++)
		{	
			if (values.size()==0)
				break;
			for (int ID=0;ID<m_maxID;ID=ID+1)
			{
				string tempstring=values[0];
				double time=atof(tempstring.c_str());
                if (ID==0 && j==0)
                    m_starttime=time;
				alltrajectories.at(ID).time.push_back(time);
				tempstring=values[2*ID+2];
				int ncarb=atoi(tempstring.c_str());
				alltrajectories.at(ID).n_carbon_t.push_back(ncarb);
			}
			csvinput->Read(values);
		}
}