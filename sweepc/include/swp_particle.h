/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Particle class represents the top node in the sub-particle tree.  This
    object is placed in the Ensemble and exhibits an interface which allows
    bulk particle properties to be summed and stored in the ensemble binary tree.

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

#ifndef SWEEP_PARTICLE_H
#define SWEEP_PARTICLE_H

#include "swp_params.h"
#include "swp_particle_model.h"
#include "swp_subparticle.h"
#include "camxml.h"
#include <vector>
#include <list>
#include <iostream>

namespace Sweep
{
class Ensemble;

/*!
 * \brief Particle which can move around, coagulate and have internal structure
 * 
 * The interface that this class should offer is not entirely clear.  This class
 * may eventually become a template parameter so that simulations can be compiled
 * for different particle types.  Alternatively an interface may be defined and then
 * subclassed.  This will require some thought. 
 */

class Particle : public SubParticle
{
public:
	// Constructors.
    Particle(                             // Initialising constructor.
        real time,                        // Create time.
        const Sweep::ParticleModel &model // Defining particle model.
        );
    Particle(Sweep::Primary &pri);        // Initialising constructor (from primary).
    Particle(const Particle &copy);       // Copy constructor.
    Particle(                             // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Model to which this particle subscribes.
        );

	// Destructor.
    ~Particle(void);
    
    //! Create a new particle using the model according to the xml data
    static Particle* createFromXMLNode(const CamXML::Element& xml, const Sweep::ParticleModel& model);

    // Operators.
    Particle &operator=(const Particle &rhs);
    Particle &operator+=(const Particle &rhs);
    const Particle operator+(const Particle &rhs) const;

    // POSITION DATA
    
    //! Get spatial position
    real getPosition() const {return m_Position;}
    
    //! Get time at which spatial position was valid
    real getPositionTime() const {return m_PositionTime;}
    
    //! Set spatial position of particle and time at which it applied
    void setPositionAndTime(const real x, const real t);

    // READ/WRITE/COPY.

    //! Clone the particle.
    Particle *const Clone() const;
    
    //! Internal consistency check
    bool IsValid() const;

    // Why is there no serialisation?

private:

    //! Spatial position of particle
    real m_Position;

    //! Time at which position was valid
    real m_PositionTime;

    // Can't create a particle without knowledge of the components
    // and the tracker variables.
    Particle(void);
};

typedef std::vector<Particle*> PartPtrVector;
typedef std::list<Particle*> PartPtrList;
}

#endif
