/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ABF_Model class declared in the
    swp_abf_model.h header file.

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

#include "swp_abf_model.h"
#include "swp_params.h"
#include "swp_component.h"
#include "swp_particle_model.h"
#include "swp_actsites_type.h"
#include "swp_particle.h"
#include "swp_ensemble.h"
#include "swp_mechanism.h"
#include "sprog.h"
#include <map>

using namespace Sweep;
using namespace Sweep::ActSites;
using namespace std;

const real ABFModel::m_sitedens = 2.3e19; // sites/m2.

// SINGLETON IMPLEMENTATION.

// Default constructor.
ABFModel::ABFModel()
: A4(-1), C2H2(-1), O2(-1), OH(-1), CO(-1),
  H(-1), H2(-1), H2O(-1), iC(-1), m_aform(AlphaConst),
  m_aconst(1.0)
{
}

// Copy constructor.
ABFModel::ABFModel(const ABFModel &copy)
{
    *this = copy;
}

// Default destructor.
ABFModel::~ABFModel()
{
}

// Assignment operator.
ABFModel &ABFModel::operator =(const ABFModel &rhs)
{
    m_alpha_prof.clear();

    if (this != &rhs) {
        A4   = rhs.A4;
        C2H2 = rhs.C2H2;
        O2   = rhs.O2;
        OH   = rhs.OH;
        CO   = rhs.CO;
        H    = rhs.H;
        H2   = rhs.H2;
        H2O  = rhs.H2O;
        iC   = rhs.iC;
        m_aform  = rhs.m_aform;
        m_aconst = rhs.m_aconst;
        for (map<real,real>::const_iterator i = rhs.m_alpha_prof.begin();
             i!=rhs.m_alpha_prof.end(); ++i) {
            m_alpha_prof[i->first] = i->second;
        }
    }
    return *this;
}

// Returns the one and only instance of the ABF model.
ABFModel &ABFModel::Instance()
{
    static ABFModel inst;
    return inst;
}

// Returns the active-sites model type.
ActSitesType ABFModel::ID(void) const {return ABFSites_ID;}


// MODEL INITIALISATION.

// Initialises the ABF model by saving the indices of
// species required to calculate the steady-state.  Returns
// <0 if model failed to initialise.
int ABFModel::Initialise(const ParticleModel &model)
{
    // Invalidate all indices.
    A4   = -1;
    C2H2 = -1;
    O2   = -1;
    OH   = -1;
    CO   = -1;
    H    = -1;
    H2   = -1;
    H2O  = -1;
    iC   = -1;

    // Locate required species in mechanism.
    Sprog::SpeciesPtrVector::const_iterator i;
    unsigned int j = 0;
    for (i=model.Species()->begin(); i!=model.Species()->end(); ++i, ++j) {
        if ((*i)->Name().compare("A4")==0) {
            A4 = j;
        } else if  ((*i)->Name().compare("C2H2")==0) {
            C2H2 = j;
        } else if  ((*i)->Name().compare("O2")==0) {
            O2 = j;
        } else if  ((*i)->Name().compare("OH")==0) {
            OH = j;
        } else if  ((*i)->Name().compare("CO")==0) {
            CO = j;
        } else if  ((*i)->Name().compare("H")==0) {
            H = j;
        } else if  ((*i)->Name().compare("H2")==0) {
            H2 = j;
        } else if  ((*i)->Name().compare("H2O")==0) {
            H2O = j;
        }
    }

    // Locate carbon component in mechanism.
    CompPtrVector::const_iterator k;
    for (k=model.Components().begin(), j=0; k!=model.Components().end(); ++k, ++j) {
        if ((*k)->Name().compare("C") == 0) {
            iC = j;
            break;
        }
    }

    // Check that indices are valid.
    if ((A4>=0) && (C2H2>=0) && (O2>=0) && (OH>=0) &&
        (CO>=0) && (H>=0) && (H2>=0) && (H2O>=0) && (iC>=0)) {
        return 0;
    } else {
        return -1;
    }
};

// Adds the ABF surface reactions to a mechanism.
void ABFModel::AddToMech(Mechanism &mech)
{
}


// SITE DENSITY CALCULATION.

// Calculates the active site density for the given gas-phase
// and particle ensemble.
real ABFModel::SiteDensity(real t, const Sprog::Thermo::IdealGas &gas,
                           const Ensemble &particles) const
{
    return m_sitedens * alpha(t, gas, particles) * radicalSiteFraction(gas);
}

// Calculates the active site density for the given gas-phase
// and particle.
real ABFModel::SiteDensity(real t, const Sprog::Thermo::IdealGas &gas,
                           const Particle &part) const
{
    return m_sitedens * alpha(t, gas, part) * radicalSiteFraction(gas);
}

// Returns the fraction of surface sites which are radicals.
real ABFModel::radicalSiteFraction(const Sprog::Thermo::IdealGas &gas) const
{
    real r1f, r1b, r2f, r2b, r3f, r4f, r5f, rdenom;
    real T  = gas.Temperature();
    real RT = RCAL * T;

    // Calculate the forward and back reaction rates.
    r1f = 4.2e+07 * exp(-13.0/RT)				 * gas.MolarConc(H);
    r1b = 3.9e+06 * exp(-11.0/RT)				 * gas.MolarConc(H2);
    r2f = 1.0e+04 * exp(-1.43/RT) * pow(T,0.734) * gas.MolarConc(OH);
    r2b = 3.68e+2 * exp(-17.1/RT) * pow(T,1.139) * gas.MolarConc(H2O);
    r3f = 2.0e+07								 * gas.MolarConc(H);
    r4f = 8.0e+01 * exp( -3.8/RT) * pow(T,1.56)  * gas.MolarConc(C2H2);
    r5f = 2.1e+06 * exp( -7.47/RT)				 * gas.MolarConc(O2);
    rdenom = r1b+r2b+r3f+r4f+r5f;

    if (rdenom > 0.0) {
        real f = (r1f+r2f) / rdenom;
        return f / (f + 1.0);
    } else {
        return 0.0;
    }
};

// Returns alpha for a particle ensemble.
real ABFModel::alpha(real t, const Sprog::Thermo::IdealGas &gas,
                     const Ensemble &particles) const
{
    switch (m_aform) {
        case AlphaCorrelation:
            return alpha(gas.Temperature(),
                         particles.GetSums().Property(ParticleCache::iNumCarbon) /
                         particles.Count());
        case AlphaProfile:
            return alpha(t);
        case AlphaConst:
        default:
            return m_aconst;
    }
}

// Returns alpha for a single particle.
real ABFModel::alpha(real t, const Sprog::Thermo::IdealGas &gas,
                     const Particle &part) const
{
    switch (m_aform) {
        case AlphaCorrelation:
            return alpha(gas.Temperature(), part.Composition(iC));
        case AlphaProfile:
            return alpha(t);
        case AlphaConst:
        default:
            return m_aconst;
    }
}

// Returns alpha using correlation.
real ABFModel::alpha(Sweep::real T, Sweep::real M1)
{
    if (M1 > 0.0) {
        real a = 12.65 - (5.63e-3 * T);
        real b = -1.38 + (6.8e-4 * T);
        return max(0.0, tanh((a / log10(M1)) + b));
    } else {
        return 0.0;
    }
}

// Returns alpha linearly interpolated from the profile.
real ABFModel::alpha(real t) const
{
    // Get the time point after the required time.
    map<real,real>::const_iterator j = m_alpha_prof.upper_bound(t);

    if (j == m_alpha_prof.begin()) {
        // This time is before the beginning of the profile.  Return
        // the first time point.
        return j->second;
    } else {
        // Get the time point before the required time.
        map<real,real>::const_iterator i = j; --i;

        // Get alpha at both points..
        real a1 = i->second;
        real a2 = j->second;

        // Calculate time interval between points i and j.
        real dt_pro = j->first - i->first;

        // Calculate time interval between point i and current time.
        real dt = t - i->first;

        // Now use linear interpolation to calculate alpha.
        return a1 + ((a2 - a1) * dt / dt_pro);
    }
}


// ALPHA CORRELATION.

// Loads a time profile for alpha into the model.  Also sets
// the model to use this profile for calculations.
void ABFModel::SetAlphaProfile(const std::map<real,real> &alphas)
{
    m_alpha_prof.clear();
    for (map<real,real>::const_iterator i = alphas.begin();
         i!=alphas.end(); ++i) {
        m_alpha_prof[i->first] = i->second;
    }
    m_aform = AlphaProfile;
}

// Tells the model to use the alpha profile for calculations.
void ABFModel::UseAlphaProfile(void)
{
    m_aform = AlphaProfile;
}

// Tells the model to use a constant value for alpha, and
// sets its value.
void ABFModel::SetAlphaConstant(real alpha)
{
    m_aconst = alpha;
    m_aform  = AlphaConst;
}

// Tells the model to use a constant value for alpha.
void ABFModel::UseAlphaConstant(void)
{
    m_aform = AlphaConst;
}

// Tells the model to use the ABF correlation for alpha.
void ABFModel::UseAlphaCorrelation(void)
{
    m_aform = AlphaCorrelation;
}
