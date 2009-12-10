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

/*!
 * Read in the life stories for a collection of PAH molecules
 *
 *@param[in]    file    Name of file from which to read the PAH life stories
 *
 *@exception    std::runtime_error      No PAH stories found
 */
void Trajectory::LoadPAHProfile(const std::string &file)
{
   	CSV_IO *csvinput = new CSV_IO(file,false);
	std::vector<std::string> values;	
	values.clear();

    // Read first line
	csvinput->Read(values);

    // One column for number of H atoms in the PAH molecules
    // and one column for C atoms.  First column is the times
    m_maxID = (values.size() - 1) / 2;

    // Create an empty trajectory (life story) for each PAH
    // that is being read.
	for (int j=0;j<m_maxID;j++)
	{
		trajectory_base newtraj;
		alltrajectories.push_back(newtraj);
	}

	for (int j=0;j<10000;j++)
		{	
			if ((j == 0) && (values.size() == 0)) {
                std::string msg("No PAHs found in ");
                msg += file;
                msg += "(Trajectory::LoadPAHProfile)";
                throw std::runtime_error(msg);
            }

            // Now read the next step of each PAH story/trajectory
            for (int ID=0;ID<m_maxID;ID=ID+1)
			{
				double time=atof(values[0].c_str());

                if (ID==0 && j==0)
                    m_starttime=time;

                alltrajectories.at(ID).time.push_back(time);

                int ncarb=atoi(values[2*ID+2].c_str());
				alltrajectories.at(ID).n_carbon_t.push_back(ncarb);
			}

            // Read the next line
			csvinput->Read(values);
            if(values.empty()) {
                //no more lines
                break;
            }
		}
}