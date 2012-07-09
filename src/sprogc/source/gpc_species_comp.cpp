/*
  Author(s):      Martin Martin (mm864)
  Project:        sprog (gas-phase and surface chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2012 Martin Martin.

  File purpose:
    Inline function definitions for Species class.

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


#include "gpc_species_comp.h"
#include "gpc_species.h"

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SpComp::SpComp(void) 
{
    m_index = -1; // Invalid index for uninitialised SpComps.
    m_count = 0;
}

// Copy constructor.
SpComp::SpComp(const SpComp &sp)
{
    m_index = -1;
    *this   = sp; // Use operator=.
}

// Initialising constructor.
SpComp::SpComp(unsigned int i, unsigned int n)
{
    m_index = i;
    m_count = n;
}

// Destructor.
SpComp::~SpComp(void)
{
    // Nothing to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
SpComp &SpComp::operator=(const Sprog::SpComp &sp)
{
    if (this!=&sp) {
        m_index = sp.m_index;
        m_count = sp.m_count;
    }
    return *this;
}

// Increment operator:  Adds n to the species count.
SpComp &SpComp::operator+=(const unsigned int n)
{
    m_count += n;
    return *this;
}

// Decrement operator:  Removes n from the species count to a zero minimum.
SpComp &SpComp::operator-=(const unsigned int n)
{
    if (n > m_count) {
        m_count = 0;
    } else {
        m_count -= n;
    }
    return *this;
}

// Addition operator:  Returns a new SpComp object whose species count
// is this object's count + n.
const SpComp SpComp::operator+(const unsigned int n)
{
    SpComp sp(*this);
    return sp += n;
}

// Subtraction operator:  Returns a new SpComp object whose species count
// is this object's count - n to a minimum of zero.
const SpComp SpComp::operator-(const unsigned int n)
{
    SpComp sp(*this);
    return sp -= n;
}

// Comparison operator:  Returns true if both SpComp objects point to
// the same species.
bool SpComp::operator==(const SpComp &sp) const
{
    return (m_index == sp.m_index) && (m_index>=0);
}

// Inequality operator:  Returns false if both SpComp objects point to
// the same species.
bool SpComp::operator!=(const SpComp &sp) const
{
    return !(*this==sp);
}


// SPECIES INDEX.

// Returns a index to the species referred to by this SpComp object.
int SpComp::Index(void) const
{
    return m_index;
}

// Sets the index to the species referred to by this SpComp object.
void Sprog::SpComp::SetIndex(unsigned int sp)
{
    m_index = sp;
}


// SPECIES COUNT.

// Returns the species count for this composition.
unsigned int Sprog::SpComp::Count() const
{
    return m_count;
}

// Sets the species count for this composition.
void Sprog::SpComp::SetCount(unsigned int n)
{
    m_count = n;
}
