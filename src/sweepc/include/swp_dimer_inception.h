/*
  Author(s):      Robert Patterson and Markus Sander (ms785)
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
    Prof Markus Kraft
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

#ifndef SWEEP_DIMER_INCEPTION_H
#define SWEEP_DIMER_INCEPTION_H

#include "swp_inception.h"

namespace Sweep
{
// Forward declare the Mechanism class.
class Mechanism;
// Forward declare the interface to the gas phase mixture
class EnvironmentInterface;

namespace Processes
{

class DimerInception : public Inception
{
public:
    // Constructors.
    DimerInception(const Sweep::Mechanism &mech); // Initialising constructor.
    DimerInception(const DimerInception &copy);        // Copy constructor.
    DimerInception(                               // Stream-reading constructor.
        std::istream &in,                    //  - Input stream.
        const Sweep::Mechanism &mech         //  - Parent mechanism.
        );

    // Destructors.
    ~DimerInception(void);

    // Operators.
    DimerInception &operator=(const DimerInception &rhs);


	// TOTAL RATE CALCULATIONS.

    // Returns rate of the process for the given system.
    double Rate(
        double t,          // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d& local_geom // Information regarding surrounding cells and boundaries
        ) const;

    // PERFORMING THE PROCESS.

    //! Performs the process on the given system.
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const;


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


    // Sets the coagulation kernel constants given incepting species
    // masses and diameters.
    void SetInceptingSpecies(
        double m1, // Mass of species 1.
        double m2, // Mass of species 2.
        double d1, // Collision diameter of species 1.
        double d2  // Collision diameter of species 2.
        );

    //! Inception only according to free mol kernel with these sizes
    void SetInceptingSpeciesFreeMol(double m1, double m2, double d1, double d2);

    // READ/WRITE/COPY.

    // Creates a copy of the inception.
    DimerInception *const Clone(void) const;

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
    // Default constructor is protected to prevent an inception being
    // defined without knowledge of the parent mechanism.
    DimerInception(void);

    //! A faster rate calculation routine for Inception events only.
    double Rate(
        const EnvironmentInterface &gas,
        double sqrtT,           // Square-root of temperature.
        double MFP,             // Gas mean free path.
        double vol              // Particle ensemble sample volume.
        ) const;

    // Calculates the gas-phase chemistry contribution to the rate
    // expression.  This is overloaded as Avogadro's number must be
    // included in the terms for inception processes.
    double chemRatePart(
        const EnvironmentInterface &gas
        ) const;

private:
    // Rate parameters.
    double m_kfm;          // Free-molecular kernel parameter.
    double m_ksf1, m_ksf2; // Slip-flow kernel parameters.

    // Free-molecular enhancement factor.
    const double m_efm;

};
}
}

#endif
