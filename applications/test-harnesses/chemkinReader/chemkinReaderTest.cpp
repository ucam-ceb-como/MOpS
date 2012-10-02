/*!
 *\file chemkinReaderTest.cpp
 *\author Laurence R. McGlashan
 *
 *\brief Test harness for the chemkin reader.
 *
 *  Copyright (C) 2011 Laurence R. McGlashan.
 *

 Licence:
    This file is part of "mops".

    brush is free software; you can redistribute it and/or
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
    Prof Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
 *
 *
 */

#include "chemkinReader.h"
#include "stringFunctions.h"
#include "boost/algorithm/string/erase.hpp"
#include <boost/assert.hpp>

int main(int argc, char* argv[])
{

    std::cout << "from_string test:\n\n";
    BOOST_ASSERT(IO::from_string<double>("9.44477437E10")==9.44477437E10);
    std::cout << "9.44477437E10 = " << IO::from_string<double>("9.44477437E10") << std::endl;
    std::cout << "1.23398640E 01 = " << IO::from_string<double>("1.23398640E 01") << std::endl;
    std::string sillyFormat("1.23398640E 01");
    std::cout << "1.23398640E 01 = " << IO::from_string<double>(boost::erase_all_copy(sillyFormat," ")) << std::endl;
    BOOST_ASSERT(IO::from_string<double>(boost::erase_all_copy(sillyFormat," "))==1.23398640E01);
    std::cout << "1.23398640G01 = " << IO::from_string<double>("1.23398640G01") << std::endl;

    std::cout << "\nChemkin Reader Test\n";

    // Using arguments like this is rather ugly, but it saves writing full argument
    // handling code in a test program.
    const std::string chemfile(argv[1]);
    const std::string thermfile(argv[2]);
    const std::string transfile(4 == argc?argv[3]:"NOT READ");

    IO::ChemkinReader chemkinReader(chemfile,thermfile,transfile);
    
    chemkinReader.read();
    chemkinReader.check();

    std::cout << chemkinReader.elements()[0].getName() << std::endl;

    std::cout << chemkinReader.species()[0].name() << std::endl;
    std::cout << chemkinReader.species()[0].thermo().getPhase() << std::endl;
    std::cout << chemkinReader.species()[0].transport().getCollisionDiameter() << std::endl;

    chemkinReader.setSpecies()[0].transport().setCollisionDiameter(1.0);

    std::cout << chemkinReader.species()[0].transport().getCollisionDiameter() << std::endl;

    return 0;

}
