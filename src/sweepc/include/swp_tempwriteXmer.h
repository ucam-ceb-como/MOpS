/*!
  * \author     Dongping Chen (dc516)
  * \file       swp_tempwriteXmer.h
  *
  * \brief      Implementation for swp_kmc_simulator.h
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    This file is used to dump Xmer information to a csv file temporarily.

  Licence:
    This file is part of "sweep".

    Sweep is free software; you can redistribute it and/or
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

#ifndef SWEEP_WRITEXMER_H
#define SWEEP_WRITEXMER_H

#include "csv_io.h"
#include <string>
#include "string_functions.h"

namespace Sweep 
{

void writeMononer(fvector &out)
{
	CSV_IO MONOMER("MONOMER.csv", true);

	std::vector<std::string> tempLine;
// writes data in columns

    for(size_t i=0; i!=out.size(); i++) {
        tempLine.clear();
        tempLine.push_back(Strings::cstr(out[i]));
        MONOMER.Write(tempLine);
		}
	MONOMER.Close();
}

void writeDimer(fvector &out) 
{
	CSV_IO DIMER("DIMER.csv", true);

	std::vector<std::string> tempLine;
    // writes data in columns

    for(size_t i=0; i!=out.size(); i++) {
        tempLine.clear();
        tempLine.push_back(Strings::cstr(out[i]));
        DIMER.Write(tempLine);
		}
	DIMER.Close();
}

void writeParimary(std::vector<fvector> &out) 
{
	CSV_IO parimary("parimary.csv", true);

	std::vector<std::string> tempLine;
	// writes data in columns

    for(size_t i=0; i!=out.size(); i++)    parimary.Write(out[i]);

	parimary.Close();
}

}// namespace Sweep
#endif
