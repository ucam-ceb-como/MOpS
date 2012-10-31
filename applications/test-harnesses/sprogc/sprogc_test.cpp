/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics) test program.

  File purpose:
    A test program for the sprogc library.
*/

#include "gpc.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace Sprog;
using namespace std;


int main()
{
    Mechanism mechIN, mechOUT;
    // Load the test mechanism.
 
 int surf_switch = 0; // change this for including / excluding surface chemistry files
 
 std::string NOT_READ = "NOT READ";
 
 
if (surf_switch == 1){

	try {

        Sprog::IO::MechanismParser::ReadChemkin("chem.inp", "surfchem.inp", mechOUT, "therm.dat", "surftherm.dat" ,1, NOT_READ);

	    std::ofstream ofs("serializeTest");
	    boost::archive::text_oarchive oa(ofs);
	    oa << mechOUT;
        ofs.close();

        cout << "Output stream." << endl;

        mechOUT.WriteDiagnostics("logOUT");

        std::ifstream ifs("serializeTest");
        boost::archive::text_iarchive oi(ifs);
        oi >> mechIN;

        cout << "Read stream." << endl;

        mechIN.WriteDiagnostics("logIN");

    } catch (exception &e) {
        cout << e.what();
        return 1;
    }




} else{

	try {

        Sprog::IO::MechanismParser::ReadChemkin("chem.inp", mechOUT, "therm.dat",1, "tran.dat");

	    std::ofstream ofs("serializeTest");
	    boost::archive::text_oarchive oa(ofs);
	    oa << mechOUT;
        ofs.close();

        cout << "Output stream." << endl;

        mechOUT.WriteDiagnostics("logOUT");

        std::ifstream ifs("serializeTest");
        boost::archive::text_iarchive oi(ifs);
        oi >> mechIN;

        cout << "Read stream." << endl;

        mechIN.WriteDiagnostics("logIN");

    } catch (exception &e) {
        cout << e.what();
        return 1;
    }

	}
	
    return 0;

}
