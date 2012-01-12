/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Mixture class declared in the mops_mixture.h
    header file.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#include "mops_mixture.h"

using namespace Mops;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
/*Mixture::Mixture(void)
{
}*/

// Default constructor (public, requires species list).
Mixture::Mixture(const Sweep::ParticleModel &model)
: Sweep::Cell(model)//,gp(NULL)
{
}

// Stream-reading constructor.
Mixture::Mixture(std::istream &in, const Sweep::ParticleModel &model)
: Sweep::Cell(in, model)//,gp(NULL)
{
}

// Default destructor.
Mixture::~Mixture(void)
{
}

//void Mixture::SetGasphaseProfile(Sweep::GasProfile* gasphase) {
//    gp=gasphase;
//}
//
//Sweep::GasProfile* Mixture::Gasphase(void) const{
//    return gp;
//}
// OPERATORS.

// Assignment operator.
Mixture &Mixture::operator =(const Mixture &rhs)
{
    if (this != &rhs) {
        // Invoke operator of base class.
        Sweep::Cell::operator=(rhs);
    }
    return *this;
}


// READ/WRITE/COPY.

// Creates a clone of the mixture.
Mixture *const Mixture::Clone() const
{
    return new Mixture(*this);
}
