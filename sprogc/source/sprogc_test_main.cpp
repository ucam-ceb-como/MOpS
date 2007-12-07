/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics) test program.

  File purpose:
    A test program for the sprogc library.
*/

#include "gpc.h"
#include <iostream>

using namespace Sprog;

/*
int test_element(void);
int test_elements_in_mech(void);
int test_species(void);
void create_water_mech(Mechanism &mech);
*/

void main()
{
    Mechanism mech;

    try {
        IO::Mechanism_IO::ReadChemkin("chem.inp", mech, "therm.dat");
    } catch (exception &e) {
        cout << e.what();
        return;
    }       
}

/*
// Tests the basic functionality of the Elements class.
int test_element(void)
{
    Element *pDef, *pH, *pO, *pHcopy, Hassign;

    try {
        // Default constructor.
        pDef = new Element();
    } catch (...) {
        cout << "Element default constructor fail." << endl;
        return 1;
    }

    try {
        // Initialising constructor.
        pH = new Element("H", 1.0079);
        pO = new Element("O", 15.9994);
    } catch (...) {
        cout << "Element initialising constructor fail." << endl;
        return 1;
    }

    try {
        // Copy constructor.
        pHcopy = new Element(*pH);
    } catch (...) {
        cout << "Element copy constructor fail." << endl;
        return 1;
   }

    try {
        // Assignment.
        Hassign = *pH;
    } catch (...) {
        cout << "Element assignment fail." << endl;
        return 1;
    }

    try {
        // Equality (Element-Element).
        if (!(*pH == Hassign) || (*pH == *pO)) {
            cout << "Element-Element equality operator not working correctly.\n";
            return 1;
        }
    } catch (...) {
        cout << "Element-Element equality operator fail." << endl;
        return 1;
    }

    try {
        // Equality (Element-String).
        if (!(*pH == "H") || (*pH == "O")) {
            cout << "Element-String equality operator not working correctly.\n";
            return 1;
        }
    } catch (...) {
        cout << "Element-String equality operator fail." << endl;
        return 1;
    }

    try {
        // Inequality (Element-Element).
        if ((*pH != Hassign) || !(*pH != *pO)) {
            cout << "Element-Element inequality operator not working correctly.\n";
            return 1;
        }
    } catch (...) {
        cout << "Element-Element inequality operator fail." << endl;
        return 1;
    }

    try {
        // Inequality (Element-String).
        if ((*pH != "H") || !(*pH != "O")) {
            cout << "Element-String inequality operator not working correctly.\n";
            return 1;
        }
    } catch (...) {
        cout << "Element-String inequality operator fail." << endl;
        return 1;
    }

    // Get/Set name.
    try {
        // This block uses Name() and SetName().
        cout << "Setting name of element " << pH->Name() << " to Jam.\n";
        pH->SetName("Jam");
        cout << "What was H is now: " << pH->Name() << endl;
        pH->SetName("H");
    } catch (...) {
        cout << "Element Get/Set name fail.\n";
    }

    // Get/Set molecular weight.
    try {
        // This block uses MolWt() and SetMolWt().
        cout << "Mol. wt. of " << pH->Name() << " is: " << pH->MolWt() << endl;

        // Setting a normal value.
        cout << "Setting mol. wt. to 2.0 ... ";
        pH->SetMolWt(2.0);
        cout << "success.\n";

        // Setting zero: unphysical, should be an error.
        cout << "Setting mol. wt. to 0.0, expecting error ... ";
        try {
            pH->SetMolWt(0.0);
            cout << "error not given!\n";
            return 1;
        } catch (out_of_range &e) {
            // Correctly thrown out_of_range.
            cout << "error thrown correctly with message:\n  \"" << e.what() << "\"\n";
        } catch (...) {
            cout << "error not thrown correctly!\n";
            return 1;
        }
    } catch (...) {
        cout << "Element Get/Set molecular weight fail.\n";
    }

    // Clear up memory (tests destructor);
    delete pDef, pH, pO, pHcopy;

    return 0;
}

// Tests the functionality of the elements in a mechanism object.
int test_elements_in_mech(void)
{
    Mechanism mech;
    Element *pH, *pO, *pOins;

    try {
        // Adding empty element
        pH = mech.AddElement();
        pH->SetName("H");
        pH->SetMolWt(1.0079);
    } catch (...) {
        cout << "Adding empty element fail." << endl;
        return 1;
    }

    try {
        // Adding a copy of existing element.        
        pOins = new Element("O", 15.9994);
        pO = mech.AddElement(*pOins);

        // Sort out memory.
        delete pOins;
    } catch (...) {
        cout << "Adding element copy fail." << endl;
        return 1;
    }

    // Get/Set name.
    try {
        // Change name of H to O.  This should cause an error because
        // O is already in use.
        try {
            cout << "Attempt to set name of H to O, expecting error ... ";
            pH->SetName("O");
            cout << "error not given!\n";
            return -1;
        } catch (invalid_argument &e) {
            cout << "error thrown correctly with message:\n  \"" << e.what() << "\"\n";
        } catch (...) {
            cout << "error not thrown correctly!\n";
            return 1;
        }
    } catch (...) {
        cout << "Element Get/Set name fail.\n";
    }

    // Don't need to clear up memory, handled automatically by Mechanism object.
    return 0;
}


// Tests the basic functionality of the Species class.
int test_species(void)
{
    Species *CHO, *CHOcopy, *OH;

    // Create some elements.
    ElementPtrVector els;
    els.push_back(new Element("H", 1.0079));
    els.push_back(new Element("O", 15.9994));
    els.push_back(new Element("C", 12.011));
    Element *H=els[0], *O=els[1], *C=els[2];

    try {
        // Default constructor.
        CHO = new Species();
    } catch (...) {
        cout << "Species default constructor fail." << endl;
        return 1;
    }

    try {
        // Set up a species
        CHO->SetName("CHO");

        // Deliberately attempt to add elements before setting
        // elements vector.
        try {
            CHO->AddElement("C", 1); // Should throw exception.
            cout << "Did not catch error when adding element before setting elements vector!\n";
            return 1;
        } catch (exception &e) {
            // Everything went ok.
        }

        // Set the elements.
        CHO->SetElements(&els);

        // Deliberately attempt to add element with incorrect name.
        try {
            CHO->AddElement("Q", 1); // Should throw exception.
            cout << "Did not catch error when adding element with unknown name!\n";
            return 1;
        } catch (invalid_argument &e) {
            // Everything went ok.
        }

        // Deliberately attempt to add element with incorrect index.
        try {
            CHO->AddElement(300, 1); // Should throw exception.
            cout << "Did not catch error when adding element with index out of range!\n";
            return 1;
        } catch (out_of_range &e) {
            // Everything went ok.
        }

        // Add valid elements using all overloads of AddElement.
        CHO->AddElement("C", 1);       // Hopefully C!
        CHO->AddElement(ElComp(0, 1)); // Hopefully H!
        CHO->AddElement(1, 1);         // Hopefully O!
    } catch (...) {
        cout << "Species set up fail." << endl;
        return 1;
    }

    try {
        // Copy constructor.
        CHOcopy = new Species(*CHO);
    } catch (...) {
        cout << "Species copy constructor fail." << endl;
        return 1;
    }

    try {
        // Assignment.
        Species CHOassign = *CHO;
    } catch (...) {
        cout << "Species assignment fail." << endl;
        return 1;
    }

    // Create a second species for comparison.
    OH = new Species();
    OH->SetName("OH");
    OH->SetElements(&els);
    OH->AddElement("O", 1);
    OH->AddElement("H", 1);

    try {
        // Equality (Species-Species).
        if (!(*CHO == *CHOcopy) || (*CHO == *OH)) {
            cout << "Species-Species equality operator not working correctly.\n";
            return 1;
        }
    } catch (...) {
        cout << "Species-Species equality operator fail." << endl;
        return 1;
    }

    try {
        // Equality (Species-String).
        if (!(*CHO == "CHO") || (*CHO == "gobbledegook")) {
            cout << "Species-String equality operator not working correctly.\n";
            return 1;
        }
    } catch (...) {
        cout << "Species-String equality operator fail." << endl;
        return 1;
    }

    try {
        // Inequality (Species-Species).
        if ((*CHO != *CHOcopy) || !(*CHO != *OH)) {
            cout << "Species-Species inequality operator not working correctly.\n";
            return 1;
        }
    } catch (...) {
        cout << "Species-Species inequality operator fail." << endl;
        return 1;
    }

    try {
        // Inequality (Element-String).
        if ((*CHO != "CHO") || !(*CHO != "jibberish")) {
            cout << "Species-String inequality operator not working correctly.\n";
            return 1;
        }
    } catch (...) {
        cout << "Species-String inequality operator fail." << endl;
        return 1;
    }

    try {
        // Check mol. wt. was calculated correctly.
        real molwt = C->MolWt() + O->MolWt() + H->MolWt();
        cout << "Mol. wt. of CHO is:      " << CHO->MolWt() << "\n";
        cout << "Mol. wt of C + H + O is: " << molwt << "\n";
        if (CHO->MolWt() != molwt) {
            cout << "Mol. wt not calculated correctly!\n";
            return 1;
        }
    } catch (...) {
        cout << "Molecular weight fail.\n";
    }

    // Clear up memory (tests destructor);
    delete CHO, CHOcopy, OH;

    return 0;
}

THERMO_PARAMS buildThermo(real a1, real a2, real a3, real a4, real a5, real a6, real a7);

void create_water_mech(Mechanism &mech)
{
    try {
        // Define elemental hydrogen and oxygen.
        Element *H = mech.AddElement(Element("H", 1.0079));
        Element *O = mech.AddElement(Element("O", 15.9994));

        // Define hydrogen radical.
        Species *Hrad = mech.AddSpecies();
        Hrad->SetName("H");
        Hrad->AddElement("H", 1);

        // Define hydrogen molecule.
        Species *H2 = mech.AddSpecies();
        H2->SetName("H2");
        H2->AddElement("H", 2);
        H2->SetThermoStartTemperature(200.0);        
        H2->AddThermoParams(1000.0, buildThermo(2.34433112, 7.98052075e-3, -1.947815100e-5, 2.01572094e-8, -7.37611761e-12, -917.935173,  0.683010238));        
        H2->AddThermoParams(3500.0, buildThermo(3.3372792, -4.94024731e-5, 4.99456778e-7, -1.79566394e-10, 2.00255376e-14, -950.158922, -3.20502331));

        // Define oxygen radical.
        Species *Orad = mech.AddSpecies();
        Orad->SetName("O");
        Orad->AddElement("O", 1);

        // Define oxygen molecule.
        Species *O2 = mech.AddSpecies();
        O2->SetName("O2");
        O2->AddElement("O", 2);

        // Define OH radical.
        Species *OH = mech.AddSpecies();
        OH->SetName("OH");
        OH->AddElement("O", 1);
        OH->AddElement("H", 1);

        // Reaction: H + O2 = O + OH.
        Reaction &r1 = *mech.AddReaction();
        r1.SetName("H + O2 = O + OH");
        r1.AddReactant("H", (unsigned int)1);
        r1.AddReactant("O2", (unsigned int)1);
        r1.AddProduct("O", (unsigned int)1);
        r1.AddProduct("OH", (unsigned int)1);
        r1.SetArrhenius(ARRHENIUS(8.3e13, 0.0, 14413.0));
        r1.SetReversible(true);

    } catch (exception &e) {
        cout << e.what();
    }
}

THERMO_PARAMS buildThermo(real a1, real a2, real a3, real a4, real a5, real a6, real a7)
{
    THERMO_PARAMS thermo;
    thermo.Count = 7;
    thermo.Params[0] = a1;
    thermo.Params[1] = a2;
    thermo.Params[2] = a3;
    thermo.Params[3] = a4;
    thermo.Params[4] = a5;
    thermo.Params[5] = a6;
    thermo.Params[6] = a7;
    return thermo;
}
*/