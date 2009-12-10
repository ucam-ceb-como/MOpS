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

#include "swp_PAH_database.h"
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
PAH_database::PAH_database() 
{
}


// Default destructor.
PAH_database::~PAH_database()
{
}

Trajectory const *PAH_database::GetTrajectory(double t) const
{
    Trajectory const *result=NULL;
    for (int i=0;i<m_num_traj;i++)
    {
        result=&(m_Trajectories.at(i));
        if (i==m_num_traj-1) 
            break;
        if (t<m_Trajectories.at(i+1).StartTime())
            break;
    }
    return result;
}

void PAH_database::LoadPAHProfiles()
{
    m_num_traj=5;
    for (int i=1;i<=m_num_traj;i++)
    {
        // Read one file for each position in the flame
        string fname = "PAH_data" + cstr(i) + ".csv";

        // Trajectory is a collection of life stories for many
        // different PAH molecules
        Trajectory newtraj;
        
        newtraj.LoadPAHProfile(fname.c_str());
        m_Trajectories.push_back(newtraj);
    }
}