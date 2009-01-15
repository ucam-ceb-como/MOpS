/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics) test program.

  File purpose:
    A test program for the sprogc library.
*/

#include "gpc.h"
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace Sprog;
using namespace std;


void main()
{
    Mechanism mech;
	int k;
    // Load the test mechanism.
    try {
        IO::MechanismParser::ReadChemkin("chem.inp", mech, "therm.dat","tran.dat");
		int k;
		SpeciesPtrVector spv = mech.Species();
		Sprog::Thermo::Mixture mix(spv);
		cout << " number of species " << mech.SpeciesCount() << endl;
		vector<double> massFrac(mech.SpeciesCount()),cp;
		Sprog::Thermo::IdealGas ig(spv);

		for(int i = 0; i!= mech.SpeciesCount(); i++)
			massFrac[i] = 0.0;

		double temp=1073.15;
		double p = 101325.0;
		ig.CalcCps(temp,cp);

		
		k = mech.FindSpecies("H2");				
		massFrac[k] = 0.0;		
		k = mech.FindSpecies("CO");
		massFrac[k] = 0.0;
		k = mech.FindSpecies("CH4");
		massFrac[k] = 0.0;
		k = mech.FindSpecies("CO2");
		massFrac[k] = 0.0;
		k = mech.FindSpecies("O2");
		massFrac[k] = 0.2329;		
		k = mech.FindSpecies("H2O");
		massFrac[k] = 0.0;
		k = mech.FindSpecies("N2");
		massFrac[k] = 0.7671;
		
		mix.SetMassFracs(massFrac);			
		mix.SetTemperature(temp);
		
		vector<double> Dkm(mech.SpeciesCount());
		Dkm = mix.getMixtureDiffusionCoeff(p);
		k = mech.FindSpecies("H2");
		cout << Dkm[k] << endl;
		k = mech.FindSpecies("O2");
		cout << Dkm[k] << endl;
		k = mech.FindSpecies("CH4");
		cout << Dkm[k] << endl;
		k = mech.FindSpecies("H2O");
		cout << Dkm[k] << endl;
		k = mech.FindSpecies("CO");
		cout << Dkm[k] << endl;		
		k = mech.FindSpecies("CO2");
		cout << Dkm[k] << endl;
		k = mech.FindSpecies("N2");
		cout << Dkm[k] << endl;

		cout << "mixture specific heat " << mix.getSpecificHeatCapacity(temp) << endl;


		

		//for(int i = 0; i!= mech.SpeciesCount(); i++){
		//	cout << spv[i]->Name() << " " << spv[i]->getSelfDiffusion(temp,p) << endl;
		//}

		//k = mech.FindSpecies("CO2");
		//cout << "H2 index " << k << endl;
		//
		//cout << spv[k]->getThermalConductivity(temp,p,cp[k]) << endl;
		
			

    } catch (exception &e) {
        cout << e.what();
        return;
    }


}