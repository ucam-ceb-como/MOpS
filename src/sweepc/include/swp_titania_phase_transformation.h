/*!
 * @file    swp_titania_phase_transformation.h
 * @author  Casper Lindberg
 * @brief   Titanium Dioxide crystal phase transformation
 *
 *   Author(s):      Casper Lindberg
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2017 Casper Lindberg
 *
 *   File purpose:
 *     
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public License
 *      as published by the Free Software Foundation; either version 2
 *      of the License, or (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this program; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *      02111-1307, USA.
 *
 *   Contact:
 *      Prof Markus Kraft
 *      Dept of Chemical Engineering
 *      University of Cambridge
 *      New Museums Site
 *      Pembroke Street
 *      Cambridge
 *      CB2 3RA, UK
 *
 *      Email:       mk306@cam.ac.uk
 *      Website:     http://como.cheng.cam.ac.uk
*/
#ifndef SWP_TITANIA_PHASE_TRANSFORMATION_H
#define SWP_TITANIA_PHASE_TRANSFORMATION_H

#include "swp_params.h"
#include "swp_particle_process.h"
#include "swp_property_indices.h"
#include "swp_process_type.h"
#include "gpc_rate_params.h"

namespace Sweep {

// Forward declare Mechanism class.
class Mechanism;

namespace Processes {

class TitaniaPhaseTransformation: public ParticleProcess {

public:
	// Mechanism constructor
    TitaniaPhaseTransformation(const Sweep::Mechanism &mech);

	// Copy constructor
	TitaniaPhaseTransformation(const TitaniaPhaseTransformation &copy);

	// Stream-reading constructor
	TitaniaPhaseTransformation(std::istream &in,const Sweep::Mechanism &mech);

	// Assignment operator
	TitaniaPhaseTransformation &operator =(const TitaniaPhaseTransformation &rhs);

	// Creates a copy of the particle process
	TitaniaPhaseTransformation *const Clone(void) const;

	// PARTICLE PROPERTY ID.

    // Returns the ID number of the particle property to which
    // the rate of this process is proportional.
    Sweep::PropID PropertyID(void) const;

    //! ID number of the particle property to which the rate of this process is proportional.
    void SetPropertyID(
        Sweep::PropID pid
        );

	// RATE CONSTANT AND PARAMETERS.

    //! Returns the Arrhenius parameter.
    Sprog::Kinetics::ARRHENIUS &Arrhenius();
    const Sprog::Kinetics::ARRHENIUS &Arrhenius() const;

    //! Sets the Arrhenius parameters.
    void SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr);

	// Returns the number of rate terms for this process.
	unsigned int TermCount(void) const;

	//Returns majorant rate
	double MajorantRate(double t, const Cell &sys, const Particle &sp) const;

	//Passes the system rate to an iterator
	double RateTerms(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom,
                             fvector::iterator &iterm) const;
	//Returns the rate for the whole system
	double Rate(
        double t,
        const Cell &sys,
        const Geometry::LocalGeometry1d &local_geom
        ) const;

	//Calculates the single-particle rate
	double Rate(
        double t,
        const Cell &sys,
        const Particle &sp
        ) const;
	
	//Performs the interparticle process system-wide
	int Perform(double t, Sweep::Cell &sys,
		const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const;

	//Performs process on a given particle n times
	int Perform(double t, Cell &sys, Particle &sp, rng_type &rng,
                unsigned int n) const;

	// Returns the process type.  Used to identify different
    // processes and for serialisation.
    ProcessType ID(void) const;

	//Writes the object to a binary stream.
	void Serialize(std::ostream &out) const;

	//Reads object from a binary stream
	void Deserialize(
		std::istream &in,
        const Sweep::Mechanism &mech);

protected:
	//!Majorant parameter. 
    const static double m_majfactor;

	// Default constructor
	TitaniaPhaseTransformation();

	//! Arrhenius rate parameters.
    Sprog::Kinetics::ARRHENIUS m_arr;

	//! Particle property to which the rate of the process is proportional.
    Sweep::PropID m_pid;

};

}

}

#endif /* SWP_TITANIA_PHASE_TRANSFORMATION_H */
