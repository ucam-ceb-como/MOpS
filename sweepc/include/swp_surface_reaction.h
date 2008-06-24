/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a surface reaction process.  Surface reactions are single
    particle processes, hence are deferred by default.
*/

#ifndef SWEEP_SURFRXN_H
#define SWEEP_SURFRXN_H

#include "swp_params.h"
#include "swp_particle_process.h"
#include "swp_submodel_type.h"

namespace Sweep
{
// Forward declare mechanism class.
class Mechanism;

namespace Processes
{
class SurfaceReaction : public ParticleProcess
{
public:
    // Constructors.
    SurfaceReaction(const Sweep::Mechanism &mech); // Default constructor.
    SurfaceReaction(const SurfaceReaction &copy);  // Copy constructor.
    SurfaceReaction(                 // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism
        );

    // Destructor.
    virtual ~SurfaceReaction(void);

    // Operators.
    SurfaceReaction &operator=(const SurfaceReaction &rhs);


    // ARRHENIUS COEFFICIENTS.

    // Returns the Arrhenius parameter.
    Sprog::Kinetics::ARRHENIUS &Arrhenius();
    const Sprog::Kinetics::ARRHENIUS &Arrhenius() const;

    // Sets the Arrhenius parameters.
    void SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr);


    // PARTICLE PROPERTY ID.

    // Returns the ID number of the particle property to which
    // the rate of this process is proportional.
    unsigned int PropertyID(void) const;

    // Returns the ID number of the particle number for which the
    // PropertyID is valid.  The mechanism should check that this
    // model is enabled.
    SubModels::SubModelType ModelID(void) const;

    // Sets the ID number of the particle property to which
    // the rate of this process is proportional.
    void SetPropertyID(
        unsigned int i,   // ID number of particle property.
        SubModels::SubModelType modelid // The model for which this ID is valid. 
          = SubModels::BasicModel_ID //  - Default model is basic particle properties.
        );


    // TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

    // Returns rate of the process for the given system.
    virtual real Rate(
        real t,         // Time.
        const Cell &sys // System for which to calculate rate.
        ) const;

	// Calculates the process rate using the given 
    // chemical conditions, rather than those conditions in the
    // given system.
    virtual real Rate(
        real t,         // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys // System for which to calculate rate.
        ) const;


	// SINGLE PARTICLE RATE CALCULATIONS.

    // Returns the rate of the process for the given particle in
    // the system. Process must be linear in particle number.
    virtual real Rate(
        real t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const;

	// Returns rate of the process for the given particle using the
    // given chemical conditions rather than those conditions in the
    // the given system.
    virtual real Rate(
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
    virtual unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  The given chemical conditions are used instead of those
    // in the given system object.
    virtual real RateTerms(
        real t,                   // Time.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Cell &sys,          // System for which to calculate rate terms.
        fvector::iterator &iterm  // Iterator to the first term.
        ) const;


    // PERFORMING THE PROCESS.

    // Performs the process on the given system.  The responsible rate term is given
    // by index.  Returns 0 on success, otherwise negative.
    virtual int Perform(
        real t,                // Time.
        Cell &sys,             // System to update.
        unsigned int iterm = 0 // The process term responsible for this event.
        ) const;

    // Performs the process on a given particle in the system.  Particle
    // is given by index.  The process is performed n times.
    virtual int Perform(
        real t,        // Current time (s).
        Cell &sys,     // System to which the particle belongs.
        Particle &sp,  // Particle for which to perform process.
        unsigned int n // Number of times to perform the process.
        ) const;


    // READ/WRITE/COPY.
    
    // Creates a copy of the particle process.
    virtual SurfaceReaction *const Clone(void) const;

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
    // Surface reaction majorant parameter.  The true rate
    // is multiplied by this parameter to get the majorised rate.
    const static real m_majfactor;

    // Arrhenius rate parameters.
    Sprog::Kinetics::ARRHENIUS m_arr;

    // Particle property to which the rate of the process is
    // proportional.
    unsigned int m_pid;

    // Particle model for which the above particle property ID
    // is valid.
    SubModels::SubModelType m_modelid;

    // Default constructor is protected to prevent reactions being
    // defined without knowledge of the parent mechanism.
    SurfaceReaction(void);
};
};
};

#endif
