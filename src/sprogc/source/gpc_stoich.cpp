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
template<class T>
Stoichiometry<T>::Stoichiometry(void)
{
    m_species = -1;
}

// Copy constructor.
template<class T>
Stoichiometry<T>::Stoichiometry(const Sprog::Stoichiometry<T> &s)
{
    m_species = s.m_species;
    m_stoich = s.m_stoich;
}

// Initialising constructor.
template<class T>
Stoichiometry<T>::Stoichiometry(unsigned int isp, const typename Stoichiometry<T>::stoich_val &mu)
{
    m_species = isp;
    m_stoich = mu;
}

// Destructor.
template<class T>
Stoichiometry<T>::~Stoichiometry()
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
template<class T>
Stoichiometry<T> &Stoichiometry<T>::operator=(const Sprog::Stoichiometry<T> &s)
{
    // Check for self-assignment!
    if (this != &s) {
        m_species = s.m_species;
        m_stoich = s.m_stoich;
    }

    return *this;
}


// SPECIES DATA.

// Returns the index of the associated species.
template<class T>
int Stoichiometry<T>::Index() const
{
    return m_species;
}

// Sets the associated species.
template<class T>
void Stoichiometry<T>::SetSpecies(const unsigned int &sp)
{
    m_species = sp;
}


// STOICHIOMETRY VALUE.

// Returns the stoichiometry value.
template<class T>
const T &Stoichiometry<T>::Mu() const
{
    return m_stoich;
}

// Sets the stoichiometry value.
template<class T>
void Stoichiometry<T>::SetMu(const T &mu)
{
    m_stoich = mu;
}

template<class T>
void Stoichiometry<T>::IncMu(const T &mu)
{
    m_stoich += mu;
}


template class Stoichiometry<int>;
template class Stoichiometry<real>;
