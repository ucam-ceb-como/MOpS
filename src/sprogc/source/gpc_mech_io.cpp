/*
  Author(s):      Matthew Celnik (msc37),
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the MechanismParser class declared in the
    gpc_mech_io.h header file.

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

//#include "comostrings.h"
#include "gpc_mech_io.h"
#include "gpc_mech.h"

#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <cstdio>

using namespace Sprog::IO;
using namespace std;
using namespace Sprog::Kinetics;


/*!
 * @param[in]   filename        Path to main CHEMKIN mechanism file
 * @param[out]  mech            Mechanism into which to place newly read data
 * @param[in]   thermofile      Thermodynamic database file
 * @param[in]   transfile       Species transport properties
 * @param[in]   verbose         Level of informational output (higher is more, default is 0)
 */
void MechanismParser::ReadChemkin
(
    const std::string &filename,
    Sprog::Mechanism &mech,
    const std::string &thermofile,
    const int verbose,
    const std::string &transFile
)
{
    // member function added by vniod to enable the reading of transport data.
    // This function utilizes the ReadChemkin function which was provided earlier.

    ::IO::ChemkinReader chemkinReader(filename,thermofile,transFile);
    chemkinReader.read();
    if (verbose >= 1) chemkinReader.check();

    ReadChemkin(chemkinReader, mech);
    if (transFile != "NOT READ")
    {
        ReadTransport(chemkinReader, mech);
    }
    if (verbose >= 1) mech.WriteDiagnostics("mech.log");

}


void MechanismParser::ReadChemkin
(
    ::IO::ChemkinReader& chemkinReader,
    Sprog::Mechanism &mech
)
{

    size_t numElements = chemkinReader.elements().size();
    size_t numSpecies = chemkinReader.species().size();
    size_t numReactions = chemkinReader.reactions().size();

    // Clear current mechanism.
    mech.Clear();
    // CHEMKIN files are read in CGS units.
    mech.SetUnits(CGS);

    // Read Elements.
    for (size_t i=0; i<numElements; ++i)
    {
        Element *el = mech.AddElement();
        el->SetName(chemkinReader.elements()[i].getName());
        el->SetMolWtFromLibrary();
        delete el;
    }

    // Read Species.
    for (size_t i=0; i<numSpecies; ++i)
    {
        Species *sp = mech.AddSpecies();
        sp->SetName(chemkinReader.species()[i].name());

        std::map<std::string, int> elementMap = chemkinReader.species()[i].thermo().getElements();
        std::map<std::string, int>::const_iterator iter;
        for (iter = elementMap.begin(); iter != elementMap.end(); ++iter)
        {
            sp->AddElement(iter->first,iter->second);
        }

        sp->SetThermoStartTemperature(chemkinReader.species()[i].thermo().getTLow());
        Sprog::Thermo::THERMO_PARAMS lp, up;
        lp.Count = 7;
        up.Count = 7;
        for (size_t j=0; j<lp.Count; ++j)
        {
            lp.Params[j] = chemkinReader.species()[i].thermo().getLowerTemperatureCoefficients()[j];
            up.Params[j] = chemkinReader.species()[i].thermo().getUpperTemperatureCoefficients()[j];
        }
        sp->AddThermoParams(chemkinReader.species()[i].thermo().getTCommon(), lp);
        sp->AddThermoParams(chemkinReader.species()[i].thermo().getTHigh(), up);
        delete sp;
    }

    // Read Reactions.
    for (size_t i=0; i<numReactions; ++i)
    {
        Kinetics::Reaction *rxn = new Kinetics::Reaction;
        rxn->SetMechanism(mech);

        std::multimap<std::string, double> reactants = chemkinReader.reactions()[i].getReactants();
        std::multimap<std::string, double> products = chemkinReader.reactions()[i].getProducts();
        std::multimap<std::string, double>::const_iterator iter;
        for (iter = reactants.begin(); iter != reactants.end(); ++iter)
        {
            rxn->AddReactant(iter->first, iter->second);
        }
        for (iter = products.begin(); iter != products.end(); ++iter)
        {
            rxn->AddProduct(iter->first, iter->second);
        }

        rxn->SetArrhenius
        (
            ARRHENIUS
            (
                chemkinReader.reactions()[i].getArrhenius().A,
                chemkinReader.reactions()[i].getArrhenius().n,
                chemkinReader.reactions()[i].getArrhenius().E*4.184e7  // cal/mol to erg/mol
            )
        );
        rxn->SetReversible(chemkinReader.reactions()[i].isReversible());
        rxn->SetUseThirdBody(chemkinReader.reactions()[i].hasThirdBody());

        if
        (
            chemkinReader.reactions()[i].hasThirdBody()
         || chemkinReader.reactions()[i].isPressureDependent()
        )
        {

            if (chemkinReader.reactions()[i].getFallOffBody() != "")
            {
                cout << "Fall off body is " << chemkinReader.reactions()[i].getFallOffBody() << endl;
                rxn->SetFallOffThirdBody(mech.FindSpecies(chemkinReader.reactions()[i].getFallOffBody()));
            }

            multimap<string, double> thirdBodies
                = chemkinReader.reactions()[i].getThirdBodies();
            multimap<string, double>::const_iterator iter;
            for (iter = thirdBodies.begin(); iter != thirdBodies.end(); ++iter)
            {
                cout << "Add third body " << iter->first << " " << iter->second << endl;
                rxn->AddThirdBody(iter->first,iter->second);
            }

            if (chemkinReader.reactions()[i].hasLOW())
            {
                cout << "Set LOW." << endl;
                rxn->SetLowPressureLimit
                (
                    ARRHENIUS
                    (
                        chemkinReader.reactions()[i].getLOW()[0],
                        chemkinReader.reactions()[i].getLOW()[1],
                        chemkinReader.reactions()[i].getLOW()[2]*4.184e7  // cal/mol to erg/mol
                    )
                );
            }

            size_t sizeSRI = chemkinReader.reactions()[i].getSRI().size();
            if (sizeSRI == 3 && chemkinReader.reactions()[i].hasSRI())
            {
                real paramSRI[3];
                for (size_t j=0; j<3; ++j)
                {
                    paramSRI[j] = chemkinReader.reactions()[i].getSRI()[j];
                }
                rxn->SetFallOffParams(SRI, paramSRI);
            }

            if (chemkinReader.reactions()[i].hasTROE())
            {
                size_t sizeTROE = chemkinReader.reactions()[i].getTROE().size();
                if (sizeTROE == 3)
                {
                    real paramTROE[3];
                    for (size_t j=0; j<3; ++j)
                    {
                        paramTROE[j] = chemkinReader.reactions()[i].getTROE()[j];
                    }
                    rxn->SetFallOffParams(Troe3, paramTROE);
                }
                else if (sizeTROE == 4)
                {
                    real paramTROE[4];
                    for (size_t j=0; j<4; ++j)
                    {
                        paramTROE[j] = chemkinReader.reactions()[i].getTROE()[j];
                    }
                    rxn->SetFallOffParams(Troe4, paramTROE);
                }
                else
                {
                    throw std::range_error("Size of TROE is not 3 or 4.");
                }
            }

        }

        //rxn->SetRevLTCoeffs(LTCOEFFS(cdble(params[0]), cdble(params[1])));
        //rxn->SetLTCoeffs(LTCOEFFS(cdble(params[0]), cdble(params[1])));

        mech.AddReaction(rxn);
        delete rxn;
    }

    mech.BuildStoichXRef();

    mech.SetUnits(SI);

}


void MechanismParser::ReadTransport
(
    ::IO::ChemkinReader& chemkinReader,
    Sprog::Mechanism &mech
)
{

    for (size_t i=0; i!=mech.SpeciesCount(); ++i)
    {
        if (chemkinReader.species()[i].name() == mech.GetSpecies(i)->Name())
        {
            mech.GetSpecies(i)->setTransportData(chemkinReader.species()[i].transport());
            mech.GetSpecies(i)->getTransportData().setReducedDipoleMoment();
        }
        else
        {
            throw runtime_error
            (
                "Species : "
                + mech.GetSpecies(i)->Name()
                + " not found in transport data \n"
            );
        }
    }

}



// CHEMKIN PARSING ROUTINES.


// Parses the units from the REACTION line of a CHEMKIN formatted file.
void Sprog::IO::MechanismParser::parseCK_Units(const std::string &rxndef,
                                               Sprog::Kinetics::ARRHENIUS &scale)
{

    throw std::runtime_error("This doesn't work!");

    // Split the string into parts.
    vector<string> parts;
    //split(rxndef, parts, " ");

    // Default scaling.
    scale.A = 1.0;
    scale.n = 1.0;
    scale.E = 4.184e7; // Convert cal/mol -> erg/mol.

    vector<string>::iterator i;
    for (i=parts.begin(); i!=parts.end(); ++i) {
        if ((*i).compare("CAL/MOLE")==0) {
            // Default scaling.
        } else if ((*i).compare("KCAL/MOLE")==0) {
            scale.E = 4184.0e7;
        } else if ((*i).compare("JOULES/MOLE")==0) {
            scale.E = 1.0e7;
        } else if ((*i).compare("KJOULES/MOLE")==0) {
            scale.E = 1.0e10;
        } else if ((*i).compare("KJOU/MOLE")==0) {
            scale.E = 1.0e10;
        } else if ((*i).compare("KJOU/MOL")==0) {
            scale.E = 1.0e10;
        } else if ((*i).compare("KELVINS")==0) {
            scale.E = R;
        } else if ((*i).compare("EVOLTS")==0) {
            scale.E = 1.60217646e-12;
        } else if ((*i).compare("MOLES")==0) {
            // Default scaling.
        } else if ((*i).compare("MOLECULES")==0) {
            scale.A = 1.0 / NA;
        }
    }

}
