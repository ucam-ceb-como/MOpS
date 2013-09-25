/*!
 * \file   swp_constant_inception.h
 * \author Robert I A Patterson
 *
 * \brief  Class for inception at a constant rate
 *
 *  Copyright (C) 2010 Robert I A Patterson.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#ifndef SWEEP_CONSTANT_INCEPTION_H
#define SWEEP_CONSTANT_INCEPTION_H

#include "swp_inception.h"

#include <utility>

namespace Geometry
{
    // Forward declaration
    class LocalGeometry1d;
}

namespace Sweep
{
// Forward declare the Mechanism class.
class Mechanism;

namespace Processes
{

/*!
 *  Incept particles at a constant rate, independent of conditions.  This
 *  is useful for testing.  The particle size does not have to be constant,
 *  currently log normal size distributions are supported, with constant
 *  size being handled as the zero variance case.
 */
class ConstantInception : public Inception
{
public: 
    //! Create as part of a mechanism
	ConstantInception(const Sweep::Mechanism &mech, const double rate,
                      const std::vector<double>& locations,
                      const std::vector<double>& scales);

	//! Create a copy
	ConstantInception(const ConstantInception &copy);

	//! Stream-reading constructor.
	ConstantInception(std::istream &in, const Sweep::Mechanism &mech);

    // Destructors.
    ~ConstantInception();

    //! Assign
    ConstantInception &operator=(const ConstantInception &rhs);

    //! Turn fixed position inception on or off
    void useFixedPosition(const bool b) {mUseFixedPosition = b;}

    //! Specify location for fixed position inception
    void setFixedPosition(const double r) {mFixedPosition = r;}

    // PERFORMING THE PROCESS.

    //! Perform a coagulation with particles chosen according to the additive kernel
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng
        ) const;

	// TOTAL RATE CALCULATIONS.

    //! Rate of the process for the given system.
    double Rate(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom) const;

	// RATE TERM CALCULATIONS.

    //! Number of rate terms for this process.
    unsigned int TermCount() const;

    // Calculates the rate terms given an iterator to a double vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all terms.
    double RateTerms(
        double t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        const Geometry::LocalGeometry1d &local_geom,                  // position information
        fvector::iterator &iterm // Iterator to the first term.
        ) const;


    // READ/WRITE/COPY.

    //! Polymorphic copy
    ConstantInception *const Clone() const;

    //! Identify process type for serialisation
    ProcessType ID() const;

    //! Write the object to a binary stream.
    void Serialize(std::ostream &out) const;

    //! Read the object from a binary stream.
    void Deserialize(std::istream &in, const Sweep::Mechanism &mech);

protected:
    //! Cannot have an inception without a mechanism
    ConstantInception(void);

private:
    //! Indicate the all particles are to be incepted at exactly the same fixed position
    bool mUseFixedPosition;

    //! Locate for fixed position inception see mUseFixedPosition
    double mFixedPosition;

    //! mean (1st) and stddev (2nd) of logarithm of the components
    std::vector<std::pair<double, double> > mComponentDistributions;
};
}
}

#endif
