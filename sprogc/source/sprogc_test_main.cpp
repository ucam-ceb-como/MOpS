/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics) test program.

  File purpose:
    A test program for the sprogc library.
*/

#include "gpc.h"
#include <iostream>
#include <stdexcept>

using namespace Sprog;
using namespace std;


void main()
{
    Mechanism mech;

    // Load the test mechanism.
    try {
        IO::MechanismParser::ReadChemkin("chem.inp", mech, "therm.dat");
    } catch (exception &e) {
        cout << e.what();
        return;
    }


}