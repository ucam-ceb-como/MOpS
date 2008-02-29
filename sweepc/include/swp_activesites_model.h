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
#include "swp_modeltype.h"
#include "sprog.h"

namespace Sweep
{
class ActiveSitesModel
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
    virtual ModelType ID(void) const = 0;
protected:
    // ActiveSitesModels are singletons.
    ActiveSitesModel(void); // Default constructor.
    ActiveSitesModel(const ActiveSitesModel &copy); // Copy constructor.
    virtual ~ActiveSitesModel(void); // Destructor.
    ActiveSitesModel &operator=(const ActiveSitesModel &rhs); // Assignment.
};
};

#endif
