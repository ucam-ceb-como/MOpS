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
    Mixture mix;
    Reactor *batch;
    Mechanism mech;
    Settings settings;
    vector<TimeInterval> times;

    // Read the chemical mechanism.
    Sprog::IO::Mechanism_IO::ReadChemkin("chem.inp", mech, "therm.dat");

    // Read the settings file.
    batch = NULL;
    batch = Settings_IO::LoadFromXML_V1("mops.inx", batch, times, settings, mech);

    /*
    // Create initial conditions.
    vector<real> initconds(mech.Species().size(), 0.0);
    initconds[mech.FindSpecies("C2H4")] = 0.148936;
    initconds[mech.FindSpecies("O2")] = 0.178723;
    initconds[mech.FindSpecies("N2")] = 0.663830;
    initconds[mech.FindSpecies("AR")] = 0.008511;
    real initT = 1650.0;  // Temperature (K).
    real initP = 1.01325; // Pressure(bar).

    // Calculate initial density = P/RT.
    real initD = initP * 1.0e5 / (Sprog::R * initT);

    // Set initial mixture.
    mix.SetSpecies(&mech.Species());
    mix.SetFracs(initconds);
    mix.SetTemperature(initT);
    mix.SetDensity(initD);

    // Fill reactor with mixture and initialise.
    batch->Fill(&mix);
    */

    batch->Initialise(times[0].StartTime());

    // Set up output.
    Console_IO cio;
    vector<string> header;
    header.push_back(settings.ConsoleVariable(0));
    header.push_back(settings.ConsoleVariable(1));
    header.push_back(settings.ConsoleVariable(2));
    header.push_back(settings.ConsoleVariable(3));
    header.push_back(settings.ConsoleVariable(4));
    header.push_back(settings.ConsoleVariable(5));
    cio.PrintDivider();
    cio.PrintRow(header);
    cio.PrintDivider();
    vector<real> out(6);
    cio.SetAutoHeader(header);
    cio.EnableAutoDividers();

    int iC2H2 = mech.FindSpecies("C2H4");
    int iH    = mech.FindSpecies("H");
    int iA4   = mech.FindSpecies("A4");

    // Initial output.
    out[0] = batch->Time();
    out[1] = batch->Mixture()->MolarConc(iC2H2) * 1.0e-6;
    out[2] = batch->Mixture()->MolarConc(iH) * 1.0e-6;
    out[3] = batch->Mixture()->MolarConc(iA4) * 1.0e-6;
    out[4] = batch->Mixture()->Density() * 1.0e-6;
    out[5] = batch->Mixture()->Temperature();
    cio.SetColumnFormat(5, Console_IO::Float);
    cio.PrintRow(out);

    // Solve reactor.
    int nsteps = 100;
    real tstop = 0.05, t2;
    real dt = tstop / (real)nsteps;
    for (int i=0; i<nsteps; i++) {
        t2 = batch->Time() + dt;
        batch->Solve(t2);
        out[0] = batch->Time();
        out[1] = batch->Mixture()->MolarConc(iC2H2) * 1.0e-6;
        out[2] = batch->Mixture()->MolarConc(iH) * 1.0e-6;
        out[3] = batch->Mixture()->MolarConc(iA4) * 1.0e-6;
        out[4] = batch->Mixture()->Density() * 1.0e-6;
        out[5] = batch->Mixture()->Temperature();
        cio.PrintRow(out);
    }
}