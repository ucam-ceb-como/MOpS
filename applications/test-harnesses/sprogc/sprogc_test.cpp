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


using namespace Sprog;
using namespace std;


int main()
{
    Mechanism mechIN, mechOUT;
	int k;
    // Load the test mechanism.
    try {

        Sprog::IO::MechanismParser::ReadChemkin("chem.inp", mechOUT, "therm.dat",1.0,"tran.dat");

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

    return 0;

}
