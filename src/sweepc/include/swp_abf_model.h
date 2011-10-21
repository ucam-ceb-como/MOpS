/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This class defines the ABF Hydrogen Abstraction - C2H2 Addtion (HACA) model
    for soot particles as discussed by Appel et al. (2000).  This model includes
    routines for calculating the fraction of radical sites on soot particle surfaces
    using the steady-state assumption and a routine for for returning the number of
    active sites using the alpha correlation given by Appel et al. (2000).  These
    routines are only valid for the surface growth model described in that paper.

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

#ifndef SWEEP_ABF_MODEL_H
#define SWEEP_ABF_MODEL_H

#include "swp_params.h"
#include "swp_actsites_model.h"
#include "swp_particle_model.h"
#include "swp_ensemble.h"
#include "swp_particle.h"
#include "swp_mechanism.h"

#include <map>

// forward declaration
namespace Sprog {
namespace Thermo {
    class IdealGas;
}
}

namespace Sweep
{
namespace ActSites
{
class ABFModel : public ActSitesModel
{
public:
    // Enumeration of the allowable forms of the alpha function.
    enum AlphaForm {
        AlphaProfile,     // Use profile of (time,alpha) coordinates.
        AlphaConst,       // Use a constant value of alpha, default = 1.0.
        AlphaCorrelation, // Use ABF alpha correlation.
    };


    // MODEL INITIALISATION.

    // Initialises the ABF model by saving the indices of
    // species required to calculate the steady-state.  Returns
    // <0 if model failed to initialise.
    int Initialise(const ParticleModel &model);

    // Adds the ABF surface reactions to a mechanism.
    static void AddToMech(Mechanism &mech);


    // SITE DENSITY CALCULATION.

    // Calculates the active site density for the given gas-phase
    // and particle ensemble.
    real SiteDensity(
        real t,                     // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Ensemble &particles   // Particle ensemble.
        ) const;

    // Calculates the active site density for the given gas-phase
    // and particle.
    real SiteDensity(
        real t,                     // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Particle &part        // Particle.
        ) const;


    // ALPHA CORRELATION.

    // Loads a time profile for alpha into the model.  Also sets 
    // the model to use this profile for calculations.
    void SetAlphaProfile(const std::map<real,real> &alphas);

    // Tells the model to use the alpha profile for calculations.
    void UseAlphaProfile(void);

    // Tells the model to use a constant value for alpha, and
    // sets its value.
    void SetAlphaConstant(real alpha);

    // Tells the model to use a constant value for alpha.
    void UseAlphaConstant(void);

    // Tells the model to use the ABF correlation for alpha.
    void UseAlphaCorrelation(void);


    // SINGLETON.

    // Returns the one and only instance of the ABF model.
    static ABFModel &Instance(void);

    // Returns the active-sites model type.
    ActSites::ActSitesType ID(void) const;

private:
    // Surface concentration of sites (active or otherwise).
    const static real m_sitedens;

    // Indices of required species in mechanism.
    int A4, C2H2, O2, OH, CO, H, H2, H2O;
    
    // Index of carbon in components.
    int iC;

    // Form of alpha to use.
    AlphaForm m_aform;

    // Constant alpha value.
    real m_aconst;

    // Profile of Alpha values.  Coordinates are (time,alpha).  If this
    // profile is empty then alpha is calculated using the ABF correlation.
    std::map<real,real> m_alpha_prof;

    // Singleton implementation.
    ABFModel(void);
    ABFModel(const ABFModel &copy);
    ~ABFModel(void);
    ABFModel &operator=(const ABFModel &rhs);

    // Returns the fraction of surface sites which are radicals.
    real radicalSiteFraction(const Sprog::Thermo::IdealGas &gas) const;

    // Returns alpha for a particle ensemble.
    real alpha(
        real t,                     // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Ensemble &particles   // Particle ensemble.
        ) const;

    // Returns alpha for a single particle.
    real alpha(
        real t,                     // Current time (s).
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        const Particle &part        // Particle.
        ) const;

    // Returns alpha using correlation.
    static real alpha(
        real T, // Temperature (K).
        real M1 // Reduced 1st C atom moment.
        );

    // Returns alpha linearly interpolated from the profile.
    real alpha(real t) const;
};
};
};

#endif
