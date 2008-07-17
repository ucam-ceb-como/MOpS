/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
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

#ifndef SWEEP_COAGULATION_H
#define SWEEP_COAGULATION_H

#include "swp_params.h"
#include "swp_process.h"
#include "swp_process_type.h"
#include "swp_particle_cache.h"
#include <vector>
#include <iostream>

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Processes
{
class Coagulation;
typedef std::vector<Coagulation*> CoagPtrVector;

class Coagulation : public Process
{
public:
    // Constructors.
    Coagulation(const Sweep::Mechanism &mech); // Default constructor.
    Coagulation(const Coagulation &copy);      // Copy-constructor.
    Coagulation(                     // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );        

    // Destructor.
    ~Coagulation(void);

    // Operators.
    Coagulation &operator=(const Coagulation &rhs);

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
    
    // Different types of transition coagulation 
    // kernel majorant types.
    enum MajorantType {
        None,    // Indicates true kernel (not majorant).
        FreeMol, // Free-molecular majorant.
        SlipFlow // Slip-flow majorant.
    };


    // TOTAL RATE CALCULATION.

    // Returns the rate of the process for the given system.
    real Rate(
        real t,         // Time.
        const Cell &sys // System for which to calculate rate.
        ) const;

/*
    // Calculates the process rate using the given 
    // chemical conditions, rather than those conditions in the
    // given system.
    real Rate(
        real t,         // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys // System for which to calculate rate.
        ) const;
*/

    // RATE TERM CALCULATION.

    // Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all rate terms.
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

/*
    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  The given chemical conditions are used instead of those
    // in the given system object.  Returns the sum of all rate terms.
    virtual real RateTerms(
        real t,                   // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys,          // System for which to calculate rate terms.
        fvector::iterator &iterm  // Iterator to the first term.
        ) const;
*/

    // PERFORMING THE PROCESS.

    // Performs the process on the given system. Must return 0
    // on success, otherwise negative.
    int Perform(
        real t,                // Time.
        Cell &sys,             // System to update.
        unsigned int iterm = 0 // The process term responsible for this event.
        ) const;


    // COAGULATION KERNEL ROUTINES.

    // Returns the transition coagulation kernel value for the 
    // two given particles.
    real CoagKernel(
        const Particle &sp1, // First particle.
        const Particle &sp2, // Second particle.
        real T,              // Temperature.
        real P,              // Pressure.
        MajorantType maj     // Type of majorant kernel to return (or true kernel).
        ) const;

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

    // READ/WRITE/COPY.

    // Creates a copy of the coagulation process.
    virtual Coagulation *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    virtual ProcessType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Free-molecular enhancement factor.  Currently hardcoded
    // for soot particles (m_efm = 2.2).
    static const real m_efm;

    // Default constructor is protected to prevent coagulations being
    // defined without knowledge of the parent mechanism.
    Coagulation(void);

    // More efficient rate routine for coagulation only.  
    // All parameters required to calculate rate passed 
    // as arguments.
    real Rate(
        const ParticleCache &data, // Particle model data.
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
        const ParticleCache &data, // Particle model data.
        real n,     // Number of particles.
        real sqrtT, // Square root of the temperature.
        real T_mu,  // T / viscosity of air.
        real MFP,   // Gas mean-free path.
        real vol,   // System sample volume.
        fvector::iterator &iterm // Iterator to first coagulation term.
        ) const;
};
};
};

#endif
