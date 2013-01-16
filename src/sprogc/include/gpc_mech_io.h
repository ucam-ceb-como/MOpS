/*
  Author(s):      Matthew Celnik (msc37),
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This class reads/writes Sprog mechanism from file.  Currently only CHEMKIN formatted
    files are supported.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
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

#ifndef GPC_MECH_PARSER_H
#define GPC_MECH_PARSER_H

#include <string>
#include "gpc_mech.h"
#include "chemkinReader.h"

namespace Sprog
{
namespace IO
{
//! Read species properties and reaction data from CHEMKIN format files
class MechanismParser
{
public:

	//! Read in a mechanism and species properties including transport properties
	static void ReadChemkin
	(
        const std::string &filename,
        Sprog::Mechanism &mech,
        const std::string &thermofile,
		const int verbose = 0,
        const std::string &transFile = "NOT READ"
	);

	/*
	* Added by mm864 to take into account surf chem.  
	*
	*/
	static void ReadChemkin
	(
        const std::string &filename,
		const std::string &Surffilename,
        Sprog::Mechanism &mech,
        const std::string &thermofile,
		const std::string &Surfthermofile,
		const int verbose = 0,
        const std::string &transFile = "NOT READ"
	);
	
		
private:

    // Reads a CHEMKIN input file.
    static void ReadChemkin
    (
        ::IO::ChemkinReader& chemkinReader,
        Sprog::Mechanism &mech
    );


    //! method to read transport data
    static void ReadTransport
    (
        ::IO::ChemkinReader& chemkinReader,
        Sprog::Mechanism &mech
    );

    // Parse the units data from the REACTION line in a chemkin formatted mechanism file.
    static void parseCK_Units(
        const std::string &rxndef,          // String containing the REACTION statement.
        Sprog::Kinetics::ARRHENIUS &scale); // Scaling factors.
   
    

};
};
};

#endif
