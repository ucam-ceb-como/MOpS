/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of an inception process.  Inceptions are modelled as the coagulation
    of 2 gas-phase species.  The rate is given by the collision kernel of these species,
    therefore one must provide their masses and diameters.  Kernel parameters are
    calculated internally.

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

#ifndef SWEEP_PAH_INCEPTION_H
#define SWEEP_PAH_INCEPTION_H

#include "swp_inception.h"

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
 *@brief Inception by dimerisation of PAH molecules
 *
 * The rate is calculated from a special column of chemical data for use with the PAH-PP model.
 */
class PAHInception : public Inception
{
public: 
    // Constructors.

    //! Initialising constructor.
    PAHInception(const Sweep::Mechanism &mech,
                 const EnvironmentInterface::PropertyIndex pah_index);

    PAHInception(const PAHInception &copy);        // Copy constructor.
    PAHInception(                               // Stream-reading constructor.
        std::istream &in,                    //  - Input stream.
        const Sweep::Mechanism &mech         //  - Parent mechanism.
        );

    // Destructors.
    ~PAHInception(void);

    // Operators.
    PAHInception &operator=(const PAHInception &rhs);

    // PERFORMING THE PROCESS.

    //! Perform a coagulation with particles chosen according to the additive kernel
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const;

    //! Perform an event to transfer mass from gasphase to particle pahse. This function is only used for PAH-PP model
    int AddInceptedPAH(
        int i,
        double t,
        Cell &sys,
        rng_type &rng) const;

	// TOTAL RATE CALCULATIONS.

    // Returns rate of the process for the given system.
    double Rate(
        double t,          // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d &local_geom
        ) const;

	// RATE TERM CALCULATIONS.

    // Returns the number of rate terms for this process.
    unsigned int TermCount(void) const;

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

    // Creates a copy of the inception.
    PAHInception *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    ProcessType ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Rate parameters.

    // Default constructor is protected to prevent an inception being
    // defined without knowledge of the parent mechanism.
    PAHInception(void);

private:
    //! Index for PAH formation rate in gas phase data
    EnvironmentInterface::PropertyIndex mPAHFormationIndex;

};
}
}

#endif
