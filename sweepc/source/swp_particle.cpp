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
#include "string_functions.h"
#include <cmath>
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Particle::Particle(void)
: m_Position(0.0)
, m_PositionTime(0.0)
{
}

// Initialising constructor.
Particle::Particle(real time, const Sweep::ParticleModel &model)
: SubParticle(time, model)
, m_Position(0.0)
, m_PositionTime(0.0)
{
}

// Initialising constructor (from Primary particle).
Particle::Particle(Sweep::Primary &pri)
: SubParticle(pri)
, m_Position(0.0)
, m_PositionTime(0.0)
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

/*!
 * @param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 */
Particle* Particle::createFromXMLNode(const CamXML::Element& xml, const Sweep::ParticleModel& model,
                                      int (*rand_int)(int, int))
{
    // Read initial particle composition.
    vector<CamXML::Element*> subitems; 
    xml.GetChildren("component", subitems);
    fvector components(model.ComponentCount(), 0);
    
    for (vector<CamXML::Element*>::iterator j=subitems.begin(); j!=subitems.end(); ++j) {
        // Get component ID.
        string str = (*j)->GetAttributeValue("id");
        int id = model.ComponentIndex(str);

        if (id >= 0) {
            // Get component value (XML uses dx to match format for inception).
            str = (*j)->GetAttributeValue("dx");
            components[id] = Strings::cdble(str);
        } else {
            // Unknown component in mechanism.
            throw std::runtime_error(str + ": Component not found in mechanism \
                                       (Sweep, Particle::createFromXMLNode).");
        }
    }

    // Read initial tracker variable values.
    xml.GetChildren("track", subitems);
    fvector trackers(model.TrackerCount(), 0);
    
    for (vector<CamXML::Element*>::iterator j=subitems.begin(); j!=subitems.end(); j++) {
        // Get tracker ID.
        string str = (*j)->GetAttributeValue("id");
        int id = model.GetTrackerIndex(str);

        if (id >= 0) {
            // Get tracker value (XML uses dx to match format for inception).
            str = (*j)->GetAttributeValue("dx");
            trackers[id] = Strings::cdble(str);
        } else {
            // Unknown tracker variable in mechanism.
            throw std::runtime_error(str + ": Tracker variable not found in mechanism. \
                                       (Sweep, Particle::createFromXMLNode).");
        }
    }
    
    
    //Particle* const pNew = new Particle(0, model);
    Particle* const pNew = model.CreateParticle(0.0, rand_int);
    
    // Initialise the new particle.
    pNew->Primary()->SetComposition(components);
    pNew->Primary()->SetValues(trackers);
    pNew->UpdateCache();
    
    return pNew;
}



// Compound assignment (coagulation).
// OPERATOR OVERLOADING.

// Assignment operator.
Particle &Particle::operator=(const Sweep::Particle &rhs)
{
    if (this != &rhs) {
        SubParticle::operator=(rhs);
        m_Position = rhs.m_Position;
        m_PositionTime = rhs.m_PositionTime;
    }
    return *this;
}

/*!
 * Both position and the associated time must be updated together.  This is
 * because it makes no sense to specify a position without knowing when it
 * applies.
 *
 *@param[in]    x       New position of particle
 *@param[in]    t       Time at which new position is correct
 */
void Particle::setPositionAndTime(const real x, const real t) {
    m_Position = x;
    m_PositionTime = t;
}

// READ/WRITE/COPY.

// Creates a clone of the particle.
Particle *const Particle::Clone() const
{
    return new Particle(*this);
}

/*!
 * Perform checks on the internal data structure.  This is mainly for
 * testing and checking purposes; it should not be called from performance
 * critical sections.
 * 
 * @return	true iff internal structures pass tests
 */
bool Particle::IsValid() const {
	return SubParticle::IsValid();
}
