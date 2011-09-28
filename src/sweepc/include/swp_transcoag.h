/* 
 * File:   swp_transcoag.h
 * Author: riap2
 * Project:        sweep (population balance solver)
 * Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of the coagulation process.  This coagulation process uses a transition
    kernel (REF).  It calculates the rates for the free molecular and slip-flow regimes.

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

#ifndef SWP_TRANSCOAG_H
#define	SWP_TRANSCOAG_H

#include "swp_coagulation.h"

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

// Forward declare class used for sums in the binary tree
class TreeWeightedCache;


namespace Transport
{
    // Forward declaration of unused argument type
    struct TransportOutflow;
}

namespace Processes
{

class TransitionCoagulation : public Coagulation
{
public:
    // Constructors.
    TransitionCoagulation(const Sweep::Mechanism &mech); // Default constructor.

    TransitionCoagulation(                     // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );
    
    //* Virtual destructor
    virtual ~TransitionCoagulation() {};
    
    virtual TransitionCoagulation* const Clone() const;

    //* Returns the process type for identification during serialisation
    virtual ProcessType ID(void) const {return Transition_Coagulation_ID;};


    // TOTAL RATE CALCULATION.

    // Returns the rate of the process for the given system.
    virtual real Rate(real t,          // Time.
                      const Cell &sys, // System for which to calculate rate.
                      const Geometry::LocalGeometry1d &local_geom
                      ) const;


    // RATE TERM CALCULATION.

    // Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all rate terms.
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        const Geometry::LocalGeometry1d &local_geom, // Spatial position information
        fvector::iterator &iterm // Iterator to the first term.
        ) const;


    // PERFORMING THE PROCESS.

    //! Perform a coagulation with particles chosen according to the transition regime kernel
    virtual int Perform(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const;

protected:
    //! Transition coagulation kernel between two particles
    virtual real CoagKernel(
        const Particle &sp1, // First particle.
        const Particle &sp2, // Second particle.
        const Cell &sys
        ) const;

    //! Majorant coagulation kernel between two particles
    virtual real MajorantKernel(
        const Particle &sp1, // First particle.
        const Particle &sp2, // Second particle.
        const Cell &sys,
        const MajorantType maj) const;

private:
        // Coagulation rate types.  These define how the rate is 
    // calculated and how the particles are chosen.
    static const unsigned int TYPE_COUNT = 6;
    enum TermType {
        SlipFlow1,
        SlipFlow2,
        SlipFlow3,
        SlipFlow4,
        FreeMol1,
        FreeMol2
    };
    

    // Free-molecular enhancement factor.  Currently hardcoded
    // for soot particles (m_efm = 2.2).
    static const real m_efm;

        // More efficient rate routine for coagulation only.  
    // All parameters required to calculate rate passed 
    // as arguments.
    real Rate(
        const TreeWeightedCache &data, // Particle model data.
        real n,     // Number of particles.
        real sqrtT, // Square root of the temperature.
        real T_mu,  // T / viscosity of air.
        real MFP,   // Gas mean-free path.
        real vol    // System sample volume.
        ) const;

    // More efficient rate routine for coagulation only.  
    // All parameters required to calculate rate terms
    // passed as arguments.
    real RateTerms(
        const TreeWeightedCache &data, // Particle model data.
        real n,     // Number of particles.
        real sqrtT, // Square root of the temperature
        real T_mu,  // T / viscosity of air.
        real MFP,   // Gas mean-free path.
        real vol,   // System sample volume.
        fvector::iterator &iterm // Iterator to first coagulation term.
        ) const;

        // COAGULATION KERNEL ROUTINES.

    // Returns the free-molecular coagulation kernel value for the 
    // two given particles.  Can return either the majorant or
    // true kernel.
    real FreeMolKernel(
        const Particle &sp1, // First particle.
        const Particle &sp2, // Second particle.
        real T,              // Temperature.
        real P,              // Pressure.
        const bool maj       // true=majorant kernel, false=true kernel.
        ) const;

    // Returns the slip-flow coagulation kernel value for the 
    // two given particles.  Can return either the majorant or
    // true kernel.
    real SlipFlowKernel(
        const Particle &sp1, // First particle.
        const Particle &sp2, // Second particle.
        real T,              // Temperature.
        real P,              // Pressure.
        const bool maj       // true=majorant kernel, false=true kernel.
        ) const;

}; //class TransitionCoagulation

} // namespace Processes

} // namespace Sweep



#endif	/*SWP_TRANSCOAG_H */

