/** 
 * File:   swp_constcoag.h
 * Author: riap2
 * Project:        sweep (population balance solver)
 * Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2011 Robert I A Patterson.

  File purpose:
    Declaration of the coagulation process for the constant kernel.

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

#ifndef SWP_SYMMETRICFRAG_H
#define	SWP_SYMMETRICFRAG_H

#include "swp_fragmentation.h"

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Processes
{
class SymmetricFragmentation : public Fragmentation
{
public:
    // Constructors.
    SymmetricFragmentation(const Sweep::Mechanism &mech); // Default constructor.

    SymmetricFragmentation(                     // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );
    
    //* Virtual destructor
    virtual ~SymmetricFragmentation() {};
    
    virtual SymmetricFragmentation* const Clone() const {return new SymmetricFragmentation(*this);};

    //* Returns the process type for identification during serialisation
    virtual ProcessType ID(void) const {return Symmetric_Fragmentation_ID;};

    // TOTAL RATE CALCULATION.

    // Returns the rate of the process for the given system.
    virtual double Rate(double t,          // Time.
                      const Cell &sys, // System for which to calculate rate.
                      const Geometry::LocalGeometry1d &local_geom //position information
                      ) const;


    // RATE TERM CALCULATION.

    // Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a double vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all rate terms.
    virtual double RateTerms(
        double t,                  // Time.
        const Cell &sys,       // Indicates true kernel (not majorant).
        const Geometry::LocalGeometry1d &local_geom, //position information
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

    //! Perform a coagulation with particles chosen according to the additive kernel
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const;
    
private:
    // Coagulation rate types.  These define how the rate is 
    // calculated and how the particles are chosen.
    static const unsigned int TYPE_COUNT = 1;
    enum TermType {
        ConstantTimesConstant,
    };

    //* Calculate kernel between two particles
    virtual double FragKernel(const Particle &sp, const Cell& sys) const;

    //* Calculate majorant kernel between two particles
    virtual double MajorantKernel(const Particle &sp, const Cell& sys, const MajorantType maj) const;

    //* Arbitrary factor to give some headroom for LPDA
    static const double s_MajorantFactor;
};

} // namespace Processes

} // namespace Sweep

#endif	/* SWP_CONSTCOAG_H */

