/*
  Author(s):      Martin Martin (mm864)
  Project:        sprog (gas-phase and surface chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2012 Martin Martin.

  File purpose:
    Implementation of the Stoichiometry class declared in the
    gpc_Delta_stoich.h header file.

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

#include "gpc_delta_stoich.h"
#include "gpc_reaction.h"
#include "string_functions.h"
#include <string>
#include <cmath>
#include <stdexcept>
using namespace Strings;
using namespace Sprog::Kinetics; 
using namespace std;


// CONSTRUCTORS AND DESTRUCTORS.
// Default constructor.
DeltaStoich::DeltaStoich(void)
{ 
  m_spName = "";
  total_stoich_sp = 0.0;
  reac_stoich_sp = 0.0;
  prod_stoich_sp = 0.0;
}


// Copy constructor.
DeltaStoich::DeltaStoich(const Sprog::Kinetics::DeltaStoich &delta_stoich)
{
  m_spName = delta_stoich.m_spName;
  total_stoich_sp = delta_stoich.total_stoich_sp; 
  reac_stoich_sp = delta_stoich.reac_stoich_sp; 
  prod_stoich_sp = delta_stoich.prod_stoich_sp;     
}

DeltaStoich::DeltaStoich(std::string &name) // Initialising constructor.
{
  m_spName = name;
  total_stoich_sp = 0.0;
  reac_stoich_sp = 0.0;
  prod_stoich_sp = 0.0;
}

DeltaStoich::DeltaStoich(std::istream &in) // Stream in initilisation
{
Deserialize(in);
}

// Destructor.
DeltaStoich::~DeltaStoich(void)
{
   // There is nothing special to destruct here.
}

// OPERATOR OVERLOADING.

// Assignment operator.
DeltaStoich &DeltaStoich::operator=(const Sprog::Kinetics::DeltaStoich &delta_stoich)
{
    // Check for self-assignment!
    if (this != &delta_stoich) {
        m_spName = delta_stoich.m_spName;
        total_stoich_sp = delta_stoich.total_stoich_sp;
        reac_stoich_sp = delta_stoich.reac_stoich_sp; 
	prod_stoich_sp = delta_stoich.prod_stoich_sp;
    }

    return *this;
}

// Comparison operator:  Compares a DeltaStoich object to a string.  Returns
// true if the DeltaStoich name is the same as the string.
bool DeltaStoich::operator==(const std::string &name) const
{
    return m_spName.compare(name)==0;
}


// Sets the species name.
void Sprog::Kinetics::DeltaStoich::SetSpeciesName(const std::string &name)
{
  m_spName = name;
}



// Sets/change the total stoichiometry change of this species. 
void Sprog::Kinetics::DeltaStoich::IncrementTotalStoich(const double &value) // value can be + or -
{
total_stoich_sp += value;
}
    
// Sets/change the reactant stoichiometry change of this species. 
void Sprog::Kinetics::DeltaStoich::IncrementReacStoich(const double &value) // value can be + or -
{
reac_stoich_sp += value;
}    

// Sets/change the product stoichiometry change of this species. 
void Sprog::Kinetics::DeltaStoich::IncrementProdStoich(const double &value) // value can be + or -
{
prod_stoich_sp += value;
}


// Sets the parent reaction.
void Sprog::Kinetics::DeltaStoich::SetReaction(Sprog::Kinetics::Reaction &rct)
{
    m_react = &rct;
}


// Prints a diagnostic output file containing all the stoichiometry change for each phase.  This is used to debug.
void Sprog::Kinetics::DeltaStoich::WriteDiagnostics(std::ostream &out) const
{
  string data = "";
    double val = 0.0;

 if (out.good()) {
        // Species Name.
        data = m_spName + " ";
        out.write(data.c_str(), data.length());

       // Total stoichiometry
        val = total_stoich_sp;
        data = cstr(val) + " ";
        out.write(data.c_str(), data.length());

	// Reactant stoichiometry
        val = reac_stoich_sp;
        data = cstr(val) + " ";
        out.write(data.c_str(), data.length());

	// Product stoichiometry
        val = prod_stoich_sp;
        data = cstr(val) + " ";
        out.write(data.c_str(), data.length());


    }

}

// Writes the phase to a binary data stream.
void Sprog::Kinetics::DeltaStoich::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialisation version number to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write the length of the species name to the stream.
        unsigned int n = m_spName.length();
        out.write((char*)&n, sizeof(n));

        // Write the species name to the stream.
        out.write(m_spName.c_str(), n);

        // Write the total stoich to the stream. 
        double totl = (double)total_stoich_sp;
        out.write((char*)&totl, sizeof(totl));

	// Write the reactant stoich to the stream. 
        double rea = (double)reac_stoich_sp;
        out.write((char*)&rea, sizeof(rea));

	// Write the prod stoich to the stream. 
        double pro = (double)prod_stoich_sp;
        out.write((char*)&pro, sizeof(pro));


    } else {
        throw invalid_argument("Output stream not ready (Sprog, DeltaStoich::Serialize).");
    }
}

void Sprog::Kinetics::DeltaStoich::Deserialize(std::istream &in)
{
    // Clear the phase of all current data.
    m_spName  = "";
    total_stoich_sp = 0.0;
    reac_stoich_sp = 0.0;
    prod_stoich_sp = 0.0;
    m_react  = NULL; 

   

    if (in.good()) {
        // Read the serialized phase version number.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0; // Need for reading name length.
        char *name = NULL;
        double totl= 0.0;
	double rea = 0.0;
	double pro = 0.0;
        switch (version) {
            case 0:
                // Read the length of the delta stoich name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the delta stoich name.
                name = new char[n];
                in.read(name, n);
                m_spName.assign(name, n);
                delete [] name;

		// Read the total stoich.
                in.read(reinterpret_cast<char*>(&totl), sizeof(totl));
				total_stoich_sp = (double)totl;

		// Read the react stoich.
                in.read(reinterpret_cast<char*>(&rea), sizeof(rea));
				reac_stoich_sp = (double)rea;

		 // Read the prod stoich.
                in.read(reinterpret_cast<char*>(&pro), sizeof(pro));
				prod_stoich_sp = (double)pro;
                break;
            default:
                throw runtime_error("DeltaStoich serialized version number "
                                    "is unsupported (Sprog, DeltaStoich::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, DeltaStoich::Deserialize).");
    }
}

