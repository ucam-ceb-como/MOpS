/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The ActiveSitesModel is a base abstract class which defines the interface
    for active site models.  These models are used by the ActiveSitesReaction
    which is a specialization of the SurfaceReaction class.  The concept
    is that particle surfaces have active "sites" with which gas-phase species
    may react.  The ActiveSitesModel uses gas-phase and particle properties
    to calculate the concentration of these active sites.
*/

#ifndef SWEEP_ACTIVESITES_MODEL_H
#define SWEEP_ACTIVESITES_MODEL_H

#include "swp_params.h"
#include "swp_ensemble.h"
#include "swp_particle.h"
#include "swp_actsites_type.h"
#include "sprog.h"

namespace Sweep
{
namespace ActSites
{
class ActSitesModel
{
public:
    // Calculates the active site density for the given gas-phase
    // and particle ensemble.
    virtual real SiteDensity(
        real t,                     // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Ensemble &particles   // Particle ensemble.
        ) const = 0;

    // Calculates the active site density for the given gas-phase
    // and particle.
    virtual real SiteDensity(
        real t,                     // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Particle &part        // Particle.
        ) const = 0;

    // Returns the active-sites model type.
    virtual ActSitesType ID(void) const = 0;
protected:
    // ActiveSitesModels are singletons.
    ActSitesModel(void); // Default constructor.
    ActSitesModel(const ActSitesModel &copy); // Copy constructor.
    virtual ~ActSitesModel(void); // Destructor.
    ActSitesModel &operator=(const ActSitesModel &rhs); // Assignment.
};
};
};

#endif
