/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Defines a condensation reaction.  Condensation reactions are treated
    differently from surface reactions as they are modelled as free molecular
    collisions.  The assumptions made in this model are:

    1.  The colliding species is much smaller than the recipient particle.

    Before the condensation process can be used it must be provided the mass and
    diameter of the condensing species in order to calculate the rate terms.  The
    SetCondensingSpecies() function is used for this.
*/

#ifndef SWEEP_CONDENSATION_H
#define SWEEP_CONDENSATION_H

#include "swp_particleprocess.h"
#include "swp_processtype.h"
#include <iostream>

namespace Sweep
{
class Condensation : public ParticleProcess
{
public: 
    // Constructors.
    Condensation(void); // Default constructor.
    Condensation(const Condensation &copy); // Copy constructor.
    Condensation(std::istream &in); // Stream-reading constructor.

    // Destructor.
    ~Condensation(void);

    // Operators.
    Condensation &operator=(const Condensation &rhs);


    // RATE CONSTANT AND PARAMETERS.

    // Returns the fixed rate constant.
    real A() const;

    // Sets the fixed rate constant.
    void SetA(real a);

    // Sets the coagulation kernel parameters given the mass and
    // collision diameter of the condensing species.
    void SetCondensingSpecies(
        real m, // Mass of the condensing species.
        real d  // Diameter of the condensing species.
        );


    // TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

    // Returns rate of the process for the given system.
    real Rate(
        real t,         // Time.
        const Cell &sys // System for which to calculate rate.
        ) const;

	// Calculates the process rate using the given 
    // chemical conditions, rather than those conditions in the
    // given system.
    real Rate(
        real t,         // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys // System for which to calculate rate.
        ) const;


	// SINGLE PARTICLE RATE CALCULATIONS.

    // Returns the rate of the process for the given particle in
    // the system. Process must be linear in particle number.
    real Rate(
        real t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const;

	// Returns rate of the process for the given particle using the
    // given chemical conditions rather than those conditions in the
    // the given system.
    real Rate(
        real t,            // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys,   // System to which the particle belongs.
        const Particle &sp // Particle for which to calculate rate.
        ) const;

    // Returns majorant rate of the process for the given system.
    real MajorantRate(
        real t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const;

	// Calculates the majorant process rate using the given 
    // chemical conditions, rather than those conditions in the
    // given system.
    real MajorantRate(
        real t,            // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys,   // System to which the particle belongs.
        const Particle &sp // Particle for which to calculate rate.
        ) const;


	// RATE TERM CALCULATIONS.
    //   These routines return the individual rate terms for a 
    //   process, which may have multiple terms (e.g. condensation).

    // Returns the number of rate terms for this process.
    unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.
    real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

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


    // PERFORMING THE PROCESS.

    // Performs the process on the given system.  The responsible rate term is given
    // by index.  Returns 0 on success, otherwise negative.
    int Perform(
        real t,            // Time.
        Cell &sys,         // System to update.
        unsigned int iterm // The process term responsible for this event.
        ) const;

    // Performs the process on a given particle in the system.  Particle
    // is given by index.  The process is performed n times.
    int Perform(
        real t,        // Current time (s).
        Cell &sys,     // System to which the particle belongs.
        Particle &sp,  // Particle for which to perform process.
        unsigned int n // Number of times to perform the process.
        ) const;


    // READ/WRITE/COPY.
    
    // Creates a copy of the particle process.
    virtual Condensation *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    virtual ProcessType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(std::istream &in);

protected:
    // Number of terms in the condensation rate expression.
    static const unsigned int TERM_COUNT;

    // Condensation majorant parameter.  The true rate
    // is multiplied by this parameter to get the majorised rate.
    static const real m_majfactor;

    // Free-molecular enhancement factor.  Currently hardcoded
    // for soot particles (m_efm = 2.2).
    static const real m_efm;

    real m_a; // Rate constant.
    real m_kfm1, m_kfm2, m_kfm3; // Free-mol term parameters.
};
};

#endif
