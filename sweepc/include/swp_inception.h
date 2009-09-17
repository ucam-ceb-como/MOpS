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

#ifndef SWEEP_INCEPTION_H
#define SWEEP_INCEPTION_H

#include "swp_params.h"
#include "swp_process.h"
#include "swp_process_type.h"
#include "swp_particle.h"
#include "swp_cell.h"
#include <vector>

namespace Sweep
{
// Forward declare the Mechanism class.
class Mechanism;

namespace Transport
{
    // Forward declaration of unused argument type
    struct TransportOutflow;
}

namespace Processes
{
class Inception;
typedef std::vector<Inception*> IcnPtrVector;

class Inception : public Process
{
public: 
    // Constructors.
    Inception(const Sweep::Mechanism &mech); // Initialising constructor.
    Inception(const Inception &copy);        // Copy constructor.
    Inception(                               // Stream-reading constructor.
        std::istream &in,                    //  - Input stream.
        const Sweep::Mechanism &mech         //  - Parent mechanism.
        );

    // Destructors.
    ~Inception(void);

    // Operators.
    Inception &operator=(const Inception &rhs);


    // RATE CONSTANT.

    // Returns the rate constant.
    real A(void) const;

    // Sets the rate constant.
    void SetA(real a);


    // INCEPTION KERNEL.

    // Sets the coagulation kernel constants given incepting species
    // masses and diameters.
    void SetInceptingSpecies(
        real m1, // Mass of species 1.
        real m2, // Mass of species 2.
        real d1, // Collision diameter of species 1.
        real d2  // Collision diameter of species 2.
        );


    // PROPERTIES OF INCEPTED PARTICLES.

    // Returns the composition vector of the new particle.
    const fvector &ParticleComp(void) const;

    // Returns the amount of the ith component of the new particle.
    real ParticleComp(unsigned int i) const;

    // Sets the particle composition vector.
    void SetParticleComp(const fvector &comp);

    // Sets the amount of the ith component in the new particle.
    void SetParticleComp(unsigned int i, real comp);

    // Returns the tracker variable vector of the new particle.
    const fvector &ParticleTrackers(void) const;

    // Returns the value of the ith tracker variable of the
    // new particle.
    real ParticleTrackers(unsigned int i) const;

    // Sets the new particle tracker variable vector.
    void SetParticleTrackers(const fvector &track);

    // Sets the value of the ith tracker variable in the
    // new particle.
    void SetParticleTracker(unsigned int i, real track);


	// TOTAL RATE CALCULATIONS.

    // Returns rate of the process for the given system.
    real Rate(
        real t,         // Time.
        const Cell &sys // System for which to calculate rate.
        ) const;

/*
	// Calculates the process rate using the given 
    // chemical conditions, rather than those conditions in the
    // given system.
    real Rate(
        real t,                   // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys           // System for which to calculate rate.
        ) const;
*/

    // Calculates the rate of multiple inceptions given a
    // vector of inceptions and an iterator to a vector of
    // reals for output.
    static real CalcRates(
        real t,                   // Time.
//        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys,          // System for which to calculate rates.
        const IcnPtrVector &icns, // Vector of inception processes.
        fvector &rates,           // Output rates vector.
        unsigned int start = 0    // Vector position to start at in vector rates.
        );


	// RATE TERM CALCULATIONS.

    // Returns the number of rate terms for this process.
    unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all terms.
    real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

/*
    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  The given chemical conditions are used instead of those
    // in the given system object.
    real RateTerms(
        real t,                   // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys,          // System for which to calculate rate terms.
        fvector::iterator &iterm  // Iterator to the first term.
        ) const;
*/

	// PERFORMING THE PROCESS.

    // Performs the process on the given system.  The responsible rate term is given
    // by index.  Returns 0 on success, otherwise negative.
    int Perform(
        real t,                // Time.
        Cell &sys,             // System to update.
        unsigned int iterm = 0,// The process term responsible for this event.
        Transport::TransportOutflow *out = 0 // Not used for this process
        ) const;


    // READ/WRITE/COPY.

    // Creates a copy of the inception.
    Inception *const Clone(void) const;

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
    real m_a;            // Rate constant.
    real m_kfm;          // Free-molecular kernel parameter.
    real m_ksf1, m_ksf2; // Slip-flow kernel parameters.

    // Properties of newly incepted particles.
    fvector m_newcomp; // Composition.
    fvector m_newvals; // Tracker variable values.

    // Free-molecular enhancement factor.  Currently hardcoded
    // for soot particles (m_efm = 2.2).
    static const real m_efm;

    // Default constructor is protected to prevent an inception being
    // defined without knowledge of the parent mechanism.
    Inception(void);

    // A faster rate calculation routine for Inception events only.  Requires all the
    // parameters that would otherwise be calculated by the routine to be passed as
    // arguments.
    real Rate(
        const fvector &fracs, // Gas-phase species mole fractions.
        real density,         // Gas-phase density.
        real sqrtT,           // Square-root of temperature.
        real T_mu,            // T / viscosity of air.
        real MFP,             // Gas mean free path.
        real vol              // Particle ensemble sample volume.
        ) const;

    // Calculates the gas-phase chemistry contribution to the rate
    // expression.  This is overloaded as Avogadro's number must be
    // included in the terms for inception processes.
    real chemRatePart(
        const fvector &fracs, // Species mole fractions in gas phase.
        real density          // Gas phase molar density.
        ) const;
};
};
};

#endif
