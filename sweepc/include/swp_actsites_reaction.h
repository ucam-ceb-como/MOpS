/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a surface reaction process which includes a term
    for active sites.  Active site densities are calculated using an
    ActiveSitesModel, which must be set before the reaction is used.
*/

#ifndef SWEEP_ACTSITE_RXN_H
#define SWEEP_ACTSITE_RXN_H

#include "swp_params.h"
#include "swp_surface_reaction.h"
#include "swp_actsites_model.h"

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Processes
{
class ActSiteReaction : public SurfaceReaction
{
public:
    // Constructors.
    ActSiteReaction(const Sweep::Mechanism &mech); // Default constructor.
    ActSiteReaction(const ActSiteReaction &copy);  // Copy constructor.
    ActSiteReaction(                 // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    // Destructor.
    virtual ~ActSiteReaction(void);

    // Operators.
    ActSiteReaction &operator=(const ActSiteReaction &rhs);


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


    // ACTIVE SITES MODEL.

    // Sets the active sites model.
    void SetModel(ActSites::ActSitesModel &model);

    // Returns the active sites model.
    ActSites::ActSitesModel *const Model(void) const;


    // READ/WRITE/COPY.
    
    // Creates a copy of the particle process.
    virtual ActSiteReaction *const Clone(void) const;

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
    // The active sites model instance.
    ActSites::ActSitesModel *m_asmodel;

    // Default constructor is protected to prevent reactions being
    // defined without knowledge of the parent mechanism.
    ActSiteReaction(void);
};
};
};

#endif
