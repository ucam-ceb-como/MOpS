// main.cpp : Defines the entry point for the console application.
//

#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
#include <iostream>
#include <stdio.h>
#include <tchar.h>
#include <exception>
#include <time.h>
#include "sweep.h"
#include "swpabf.h"
#include "swpoutput.h"

using namespace std;
using namespace Sweep;
using namespace Sweep::ABF;

int main(int argc, _TCHAR* argv[])
{
    // Program timer.
    clock_t ct1, ct2;

    // Variables for the solver and the solver object.
    int i=0, n=100, r=0, nruns=10;
    Sweep::real t = 0.0, dt = 0.0002;
    Sweep::Solver solver;
    char outfile[] = "sweep3-coag.dat";

    // File and console output.
    Sweep::SweepOutput output;

    // Create test case flame.
    PremixFlame flame;
//    flame.ReadProfile("chem-debug.dat");
    flame.ReadProfile("chem-debug.dat");

    // Create a mechanism.
    Sweep::Mechanism mech;
    mech.SetSpeciesList(flame.GetSpeciesList());
    ABFMech::InitHACA(flame.GetSpeciesList());

    // Read mechanism from file.
    try {
        Sweep::XMLIO::ReadMechanism("sweep.xml", mech);
    } catch (exception &e) {
        // Failed to read mecanism.
        output.PrintConsoleMsg(e.what());
        return -1;
    }
    // Create coagulation process.
    mech.AddCoagulation();

    // Initialise ensemble to correct size.
    flame.Ensemble().Initialise(65536, DefaultParticle::NCACHE+mech.ComponentCount()+mech.ValueCount());

    // Run the simulation.
    for (r=1; r<=nruns; r++) {
        ct1 = clock();
        printf("Run number: %d\n", r);
        flame.Reset(6.5336e11);
        t = 0.0;

        output.Open(outfile, r);
        output.PrintToConsole(t, flame, mech);

        for (i=1; i<=n; i++) {
            // Solve time step.
            solver.Run(&t, (Sweep::real)i * dt, flame, mech);

            // Provide some output.
            output.Write(t, flame, mech);
            output.PrintToConsole(t, flame, mech);
        }

        output.WritePSL(t, flame, 1);
        output.Close();
        ct2 = clock();
        std::cout << (double)(ct2-ct1) / CLOCKS_PER_SEC << " seconds for run." << endl;
    }

    for (i=0; i<(int)solver.m_processcounter.size(); i++) {
        cout << solver.m_processcounter[i] << endl;
    }

    // Post-process output.
    output.PostProcess(outfile, 1, nruns, n, mech);
    output.PostProcessPSL(outfile, 1, nruns, 1, mech);

	return 0;
}
