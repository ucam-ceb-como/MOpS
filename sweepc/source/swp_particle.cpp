/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Particle class declared in the
    swp_particle.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#include "swp_particle.h"
#include "swp_particle_model.h"
#include "swp_submodel.h"
#include <cmath>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Particle::Particle(void)
: m_ensemble(NULL)
{
}

// Initialising constructor.
Particle::Particle(real time, const Sweep::ParticleModel &model)
: SubParticle(time, model), m_ensemble(NULL)
{
}

// Initialising constructor (from Primary particle).
Particle::Particle(Sweep::Primary &pri)
: SubParticle(pri), m_ensemble(NULL)
{
}

// Copy constructor.
Particle::Particle(const Sweep::Particle &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
Particle::Particle(std::istream &in, const Sweep::ParticleModel &model)
: SubParticle(in, model)
{
}

// Default destructor.
Particle::~Particle()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
Particle &Particle::operator=(const Sweep::Particle &rhs)
{
    if (this != &rhs) {
        SubParticle::operator=(rhs);
        m_ensemble = rhs.m_ensemble;
    }
    return *this;
}

// Addition-assignment operator.  This implements coagulation.
Particle &Particle::operator+=(const Sweep::Particle &rhs)
{
    SubParticle::operator+=(rhs);
    return *this;
}


// Addition operator.  This also implements coagulation.
const Particle Particle::operator +(const Sweep::Particle &rhs) const
{
    return Particle(*this) += rhs;
}


// PARENT ENSEMBLE.

// Returns the parent ensemble.
const Sweep::Ensemble *const Particle::Ensemble(void) const {return m_ensemble;}

// Sets the parent ensemble.
void Particle::SetEnsemble(Sweep::Ensemble &ens) {m_ensemble = &ens;}


// READ/WRITE/COPY.

// Creates a clone of the particle.
Particle *const Particle::Clone() const
{
    return new Particle(*this);
}
