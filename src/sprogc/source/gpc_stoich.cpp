/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Stoichiometry class declared in the
    gpc_stoich.h header file.

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

#include "gpc_stoich.h"

using namespace Sprog;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Stoichiometry::Stoichiometry(void)
{
    m_species = -1;
    //m_phaseName = ""; // Added by mm864
}

// Copy constructor.
Stoichiometry::Stoichiometry(const Sprog::Stoichiometry &s)
{
    m_species = s.m_species;
    m_stoich = s.m_stoich;
    // m_phaseName = s.m_phaseName; // Added by mm864
}

// Initialising constructor.
Stoichiometry::Stoichiometry(unsigned int isp, const double &mu)
{
    m_species = isp;
    m_stoich = mu;
    //m_phaseName = ""; 
}

// OPERATOR OVERLOADING.

// Assignment operator.
Stoichiometry &Stoichiometry::operator=(const Sprog::Stoichiometry &s)
{
    // Check for self-assignment!
    if (this != &s) {
        m_species = s.m_species;
        m_stoich = s.m_stoich;
        //m_phaseName = s.m_phaseName; 
    }

    return *this;
}


// SPECIES DATA.

// Returns the index of the associated species.
int Stoichiometry::Index() const
{
    return m_species;
}

// Sets the associated species.
void Stoichiometry::SetSpecies(const unsigned int &sp)
{
    m_species = sp;
}

/*
// Sets the associated speciePhase
void Stoichiometry::SetSpeciesPhase(const Sprog::Species &species)
{
  m_phaseName = species.PhaseName();
}
*/


// STOICHIOMETRY VALUE.

// Returns the stoichiometry value.
const double &Stoichiometry::Mu() const
{
    return m_stoich;
}

// Sets the stoichiometry value.
void Stoichiometry::SetMu(const double &mu)
{
    m_stoich = mu;
}
// Add the stoichiometry
void Stoichiometry::IncMu(const double &mu)
{
    m_stoich += mu;
}
