/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    This is the main file for the mops solver.  Includes the driver
    program.
*/

#include "mops.h"
#include "sprog.h"
#include "console_io.h"
#include <vector>

using namespace Mops;
using namespace std;

int main(void)
{
    Reactor *batch = NULL;
    Mechanism mech;
    Solver solver;
    vector<TimeInterval> times;

    // Read the chemical mechanism.
    Sprog::IO::MechanismParser::ReadChemkin("chem.inp", mech, "therm.dat");

    // Read the settings file.
    batch = Settings_IO::LoadFromXML_V1("mops.inx", batch, times, solver, mech);

    // Solve reactor.
    solver.SolveReactor(*batch, times);

    // Post-process.
    solver.PostProcess(solver.OutputFile());

    // Clear up memory.
    delete batch;
}